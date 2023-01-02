#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#include <ctime>
#include <cassert>
#include "./lrc_module.hpp"
#include "./maths_module.hpp"
#include "./md_lj_module.hpp"
#include "./averages_module.hpp"
#include "./config_io_module.hpp"

/*  Takes in a configuration of atoms (positions, velocities)
    Cubic periodic boundary conditions
    Conducts molecular dynamics using BAOAB algorithm of BJ Leimkuhler and C Matthews
    Appl. Math. Res. eXpress 2013, 34–56 (2013); J. Chem. Phys. 138, 174102 (2013)
    Uses no special neighbour lists

    Reads several variables and options from standard input using JSON format
    Leave input empty "{}" to accept supplied defaults

    Positions r are divided by box length after reading in and we assume mass=1 throughout
    However, input configuration, output configuration, most calculations, and all results 
    are given in simulation units defined by the model
    For example, for Lennard-Jones, sigma = 1, epsilon = 1

    Despite the program name, there is nothing here specific to Lennard-Jones
    The model is defined in md_lj_module*/

#define nblock        10
#define nstep         1000
#define temperature   1.0
#define r_cut         2.5
#define dt            0.005
#define gamma         1.0

std::vector<VariableType> calc_variables(PotentialType tot, double** r, int n, double box ){
/*
     --- some notes from Allen-Tildesley ---
     Variables of interest, of class VariableType, containing three attributes:
        .val: the instantaneous value
        .nam: used for headings
        .method: indicating averaging method
      If not set below, .method adopts its default value of avg
      The .nam and some other attributes need only be defined once, at the start of the program,
      but for clarity and readability we assign all the values together below */


    // Preliminary calculations (n,v,f,total are taken from the calling program)
    double vol = pow(box,3);                        //  Volume
    double rho = n / vol;                           //  Density
    double kin = 0.5*                               //  Kinetic energy
    double fsq = pow(force (tot, n, box, r_cut, r),2);   //  Total squared force from md_lj_module
    // std::cout << " ---- force:   " << fsq << " ---- \n";

    // Initial energy and overlap check
    potential (tot,n, box, r_cut, r);

    // Internal energy (cut-and-shifted) per atom
    // Total KE plus total cut-and-shifted PE divided by N
    VariableType e_s;
    e_s.nam = "E/N cut&shifted";
    e_s.val = (kin + tot.pot)/n;

    // Internal energy (full, including LRC) per atom
    // LRC plus total KE plus total cut (but not shifted) PE divided by N
    VariableType e_f;
    e_f.nam = "E/N full";
    e_f.val = potential_lrc(rho,r_cut) + ((kin+tot.cut)/n);

    // Pressure (cut-and-shifted)
    // Ideal gas contribution plus total virial divided by V
    VariableType p_c;
    p_c.nam = "P cut&shifted";
    p_c.val = rho*temperature + tot.vir/vol;

    // Pressure (full, including LRC)
    // LRC plus ideal gas contribution plus total virial divided by V
    VariableType p_f;
    p_f.nam = "P full";
    p_f.val = pressure_lrc(rho,r_cut) + rho*temperature + tot.vir/vol;

    // Kinetic temperature
    // Momentum is not conserved, hence 3N degrees of freedom
    VariableType t_k;
    t_k.nam = "T kinetic";
    t_k.val = val = 2.0*kin/(3*n);

    // Configurational temperature
    // Total squared force divided by total Laplacian
    VariableType t_c;
    t_c.nam = "T config";
    t_c.val = fsq/tot.lap;

    // Heat capacity (cut-and-shifted)
    // Total energy divided by temperature and sqrt(N) to make result intensive;
    VariableType c_s;
    c_s.nam = "Cv/N cut&shifted";
    c_s.val = (kin+tot.cut)/(temperature*sqrt(n));
    c_s.method = msd;
    c_s.instant = false;

    // Heat capacity (full)
    // Total energy divided by temperature and sqrt(N) to make result intensive; LRC does not contribute
    VariableType c_f;
    c_f.nam = "Cv/N full";
    c_f.val = (kin+tot.pot)/(temperature*sqrt(n));
    c_f.method = msd;
    c_f.instant = false;

    // list the VariableType objects
    std::vector<VariableType> variables;

    variables.push_back(e_s);
    variables.push_back(e_f);
    variables.push_back(p_s);
    variables.push_back(p_f);
    variables.push_back(t_c);
    variables.push_back(t_f);
    variables.push_back(c_s);
    variables.push_back(c_f);

    return variables;

}

void a_propagator(double time, int n, double box, double** r, double** vel){
/*    A: drift step propagator.
    t is the time over which to propagate (typically dt/2). */

    scalar2DArrayMultip(n,3,t/box,vel);    // (t/box)*velocity
    sum2DArrays(n,3,vel,r);                // Positions in box=1 units
    rint2D(n,3,r);                         // Periodic boundaries
}

void b_propagator(potentialType tot, double time, int n, double box, double** vel, double** r, double** f){
/*  B: kick step propagator.
    t is the time over which to propagate (typically dt/2).
    v is accessed from the calling program. */
    f = force (tot, n, box, r_cut, r);
    scalar2DArrayMultip(n,3,t,f);    // t * force
    matMultip(n,3,vel,f,vel);    // t * force
    sum2DArrays(n,3,vel,vel);
}

void o_propagator ( double t, double** vel, int n ){
/* O: friction and random contributions propagator.

    t is the time over which to propagate (typically dt).
    v, n, temperature, and gamma are accessed from the calling program. */

    double** rnd;
    rand2DArray(n,3,rnd);
    double x = gamma*t
    double c;
    if (x > 0.0001)
        c = 1-exp(-2*x);
    else 
        c = -2/3*pow(x,4)+4/3*pow(x,3)-2.0*pow(x,2)+2.0*x;
    
    matMultip(n,3,c,c);
    scalar2DArrayMultip(n,3,exp(-x),vel);
    matMultip(n,3,sqrt(temperature),c);
    scalar2DArrayMultip(n,3,c,rnd);
    sum2DArrays(n,3,rnd,vel);
}



int main(){

    // initial time for calculating the processing time
    std::clock_t ti = std::clock();

    // Preliminary calculations (n,r,total are taken from the calling program)
    const char* file = "cnf.inp";

    // Read in initial configuration
    std::ifstream input(file);

    int n;
    double box;

    input >> n;
    input >> box;
    input.close();

    int na   = 3;
    int ndim = 3;
    double** r      = allocate2DArray(n,3);
    double** e      = allocate2DArray(n,4);
    double** vel    = allocate2DArray(n,3);
    double** angvel = allocate2DArray(n,3);
    bool quaternion = false;
    bool with_v     = true;

    read_cnf_mols(file, quaternion, with_v, r, e, vel, angvel);

    scalar2DArrayDivision(n,3,box,r);    // Convert positions to box units
    rint2D(n,3, r);                      // Periodic boundaries

    // Initial forces, potential, etc plus overlap check
    PotentialType total;

    potential (total,n, box, r_cut, r);

    assert (!total.ovr); 
    std::cout << "No overlap in initial configuration! \n";
        
    std::cout << '\n';
    std::cout << "bd_nvt_lj \n";
    std::cout << "Brownian dynamics, constant-NVT ensemble \n";
    std::cout << "Particle mass=1 throughout \n";
    std::cout << '\n';

    introduction()

    //  Write out parameters
    std::cout << '\n';
    printf("%16s %42d   \n", "Number of blocks",            nblock);
    printf("%25s %33d   \n", "Number of steps per block",   nstep);
    printf("%25s %33d   \n", "Potential cutoff distance",   r_cut);
    printf("%20s %36.6f \n", "Time Step",                   dt);
    printf("%20s %36.6f \n", "Friction coeffcient",         gamma);
    printf("%20s %37.6f \n", "Specified temperature",       temperature);
    printf("%7s  %50.6f \n", "Ideal diffusion coeff",       temperature/gamma);
    std::cout << '\n'; 

    std::cout << '\n';
    printf("%10s  %39d   \n", "Number of particles",n);
    printf("%10s  %48.6f \n", "Box length", box);
    printf("%7s   %50.6f \n", "Density", n/pow(box,3));
    std::cout << '\n'; 

    BlockVar blk_var;
    
    // Initialize arrays for averaging and write column headings
    run_begin (calc_variables(total, r, n, box), blk_var, ti);

    for (int blk{0}; blk < nblock;++blk){                                // Loop over blocks

	  blk_begin(n_avg,blk_var); // from md_lj_module

        for (int stp{0}; stp<nstep;++stp){                               // Loop over steps    

        b_propagator ( dt/2 ) // B kick half-step
        a_propagator ( dt/2 ) // A drift half-step
        o_propagator ( dt )   // O random velocities and friction step
        a_propagator ( dt/2 ) // A drift half-step

          }
    }
}
