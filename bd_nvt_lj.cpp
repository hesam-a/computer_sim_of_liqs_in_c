#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <ctime>
#include <cassert>
#include <iomanip>
#include <time.h>
#include "./lrc_module.hpp"
#include "./maths_module.hpp"
#include "./md_lj_module.hpp"
#include "./averages_module.hpp"
#include "./config_io_module.hpp"

// Brownian dynamics, NVT ensemble.

#define nblock        10
#define nstep         1000
#define temperature   1.0
#define r_cut         2.5
#define dt            0.005
#define gamma         1.0

std::vector<VariableType> calc_variables(PotentialType tot, double** r, int n, double box, double** vel, double** f ){
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
    double vol      = pow(box,3);                        //  Volume
    double rho      = n / vol;                           //  Density
    double** f_sq   = allocate2DArray(n,3) ;
    double** vel_sq = allocate2DArray(n,3);
    elementWise2DProduct(n,3,vel,vel,vel_sq);
    double vel_sum  = elementSum2D(n,3, vel_sq);
    double kin      = 0.5 * vel_sum;                     //  Kinetic energy    

    force (tot, n, box, r_cut, r, f);              
    elementWise2DProduct(n,3,f,f,f_sq);                  //  Total squared force from md_lj_module
    double fsq = elementSum2D(n,3, f_sq);

    // Internal energy (cut-and-shifted) per atom
    // Total KE plus total cut-and-shifted PE divided by N
    VariableType e_s;
    e_s.nam = "E/N cut&shifted";
    e_s.val = (kin + tot.pot)/n;

    // Internal energy (full, including LRC) per atom
    // LRC plus total KE plus total cut (but not shifted) PE divided by N
    VariableType e_f;
    e_f.nam = "E/N full";
    e_f.val = potential_lrc(rho,r_cut) + (kin+tot.cut)/n;

    // Pressure (cut-and-shifted)
    // Ideal gas contribution plus total virial divided by V
    VariableType p_s;
    p_s.nam = "P cut&shifted";
    p_s.val = rho*temperature + tot.vir/vol;

    // Pressure (full, including LRC)
    // LRC plus ideal gas contribution plus total virial divided by V
    VariableType p_f;
    p_f.nam = "P full";
    p_f.val = pressure_lrc(rho,r_cut) + rho*temperature + tot.vir/vol;

    // Kinetic temperature
    // Momentum is not conserved, hence 3N degrees of freedom
    VariableType t_k;
    t_k.nam = "T kinetic";
    t_k.val = 2.0*kin/(3*n);

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
    variables.push_back(p_s);
    variables.push_back(e_f);
    variables.push_back(p_f);
    variables.push_back(t_k);
    variables.push_back(t_c);
    variables.push_back(c_s);
    variables.push_back(c_f);

    free2DArray(n,f_sq);
    free2DArray(n,vel_sq);

    return variables;

}

void deletePointer(std::vector<VariableType> vars, BlockVar blk_var){
//  A function for deleting the VariableType objects called by calc_variables

    // delete blockVariable's instances
    delete [] blk_var.run_avg;
    delete [] blk_var.run_err;
    delete [] blk_var.blk_avg;
    delete [] blk_var.blk_msd;
    delete [] blk_var.values;
    delete [] blk_var.addd;
    delete [] blk_var.mask;
    delete [] blk_var.methodd;

}

void a_propagator(double t, int n, double box, double** r, double** vel){
/*    A: drift step propagator.
    t is the time over which to propagate (typically dt/2). */

    double** vel_p = allocate2DArray(n,3);
    scalar2DArrayMultip(n,3,t/box,vel,vel_p);    // (t/box)*velocity
    sum2DArrays(n,3,vel_p,r);                    // Positions in box=1 units
    rint2D(n,3,r);                               // Periodic boundaries

    free2DArray(n,vel_p);
}

void b_propagator(double t, int n, double box, double** r, double** vel, double** f){
/*  B: kick step propagator.
    t is the time over which to propagate (typically dt/2).
    v is accessed from the calling program. */

    double** f_p = allocate2DArray(n,3);
    scalar2DArrayMultip(n,3,t,f,f_p);         // t * force
    sum2DArrays(n,3,f_p,vel);

    free2DArray(n,f_p);
}

void o_propagator ( double t, int n, double** vel ){
/* O: friction and random contributions propagator.

    t is the time over which to propagate (typically dt).
    v, n, temperature, and gamma are accessed from the calling program. */

    double     x = gamma*t;
    double    c1 = 2.0;
    double    c2 = -2.0;
    double    c3 = 4.0/3.0;
    double    c4 = -2.0/3.0;
    double    c;
    double** rnd = allocate2DArray(n,3);

    if (x > 0.0001)
        c = 1-exp(-2*x);
    else 
        c = x * ( c1 + x * ( c2 + x * ( c3 + x * c4 ) ) );
//        c = -2/3*pow(x,4)+4/3*pow(x,3)-2.0*pow(x,2)+2.0*x;
    
    c = sqrt(c);
    rand2DArray(n,3,rnd);
    scalar2DArrayMultip(n,3,exp(-x),vel,vel);
    scalar2DArrayMultip(n,3,c*sqrt(temperature),rnd,rnd);
    sum2DArrays(n,3,rnd,vel);

    free2DArray(n,rnd);
}


int main(){

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
    The model is defined in md_lj_module. */

    // initial time for calculating the processing time
    std::clock_t ti = std::clock();

    // Preliminary calculations (n,r,total are taken from the calling program)
    const char* file = "cnf_vel.inp";

    // Read in initial configuration
    std::ifstream input(file);

    int n;
    double box;

    input >> n;
    input >> box;
    input.close();

    double** r      = allocate2DArray(n,3); // position
    double** f      = allocate2DArray(n,3); // force
    double** e      = allocate2DArray(n,4); // just satisfy the "read_cnf_mols" function
    double** vel    = allocate2DArray(n,3); // velocity
    double** angvel = allocate2DArray(n,3); // just satisfy the "read_cnf_mols" function
    bool quaternion = false;
    bool with_v     = true;

    read_cnf_mols(file, quaternion, with_v, r, e, vel, angvel);

    scalar2DArrayDivision(n,3,box,r);    // Convert positions to box units
    rint2D(n,3, r);                      // Periodic boundaries

    
    std::cout << '\n';
    std::cout << "bd_nvt_lj \n";
    std::cout << "Brownian dynamics, constant-NVT ensemble \n";
    std::cout << "Particle mass=1 throughout \n";
    std::cout << '\n';

    introduction();

    //  Write out parameters
    std::cout << '\n';
    printf("%16s %42d   \n", "Number of blocks",            nblock);
    printf("%25s %33d   \n", "Number of steps per block",   nstep);
    printf("%25s %33f   \n", "Potential cutoff distance",   r_cut);
    printf("%9s  %48.6f \n", "Time Step",                   dt);
    printf("%19s %39.6f \n", "Friction coeffcient",         gamma);
    printf("%20s %37.6f \n", "Specified temperature",       temperature);
    printf("%7s  %36.6f \n", "Ideal diffusion coeff",       temperature/gamma);
    std::cout << '\n'; 

    std::cout << '\n';
    printf("%10s  %38d   \n", "Number of particles",n);
    printf("%10s  %47.6f \n", "Box length", box);
    printf("%7s   %49.6f \n", "Density", n/pow(box,3));
    std::cout << '\n'; 

    BlockVar blk_var;

    // Initial forces, potential, etc plus overlap check
    PotentialType total;

    force (total, n, box, r_cut, r, f);

    assert (!total.ovr); 
    std::cout << "No overlap in initial configuration! \n";
    std::cout << '\n';

    int n_avg = calc_variables(total, r, n, box, vel, f).size();

    // Initialize arrays for averaging and write column headings
    run_begin (calc_variables(total, r, n, box, vel, f), blk_var, ti);


    for (int blk{0}; blk < nblock;++blk){                                // Loop over blocks

	blk_begin(n_avg,blk_var);

        for (int stp{0}; stp<nstep;++stp){                               // Loop over steps    
	    
 	    b_propagator ( dt/2, n, box, r, vel, f );    // B kick half-step
            a_propagator ( dt/2, n, box, r, vel );       // A drift half-step
            o_propagator ( dt,   n, vel );               // O random velocities and friction step
            a_propagator ( dt/2, n, box, r, vel );       // A drift half-step

            force (total, n, box, r_cut, r, f);

            assert (!total.ovr); 
            //std::cout << "No overlap in configuration! \n";

            b_propagator ( dt/2, n, box, r, vel, f ); // A drift half-step

            blk_add (calc_variables(total, r, n, box, vel, f),blk_var);

        }
	
        blk_end ( blk, n_avg, blk_var);                         // Output block averages
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(3) << std::to_string(blk+1);
        std::string sav_tag(ss.str());                         // Number configuration by block
        double** out_r = allocate2DArray(n,3);
        scalar2DArrayMultip(n,3,box,r,out_r);
        //write_cnf_mols ("cnf."+sav_tag, n, box, quaternion, with_v, out_r, e, vel, angvel );    // Save configuration
        free2DArray(n,out_r);
    }
    run_end (calc_variables(total, r, n, box, vel, f), blk_var, ti);
    force (total,n, box, r_cut, r, f);

    assert (!total.ovr); 
    std::cout << "No overlap in final configuration! \n";

    double** out_r = allocate2DArray(n,3);
    scalar2DArrayMultip(n,3,box,r,out_r);
    //write_cnf_mols ("cnf.out", n, box, quaternion, with_v, out_r, e, vel, angvel );    // Save configuration
    free2DArray(n,out_r);
    deletePointer(calc_variables(total, r, n, box, vel, f), blk_var);

    conclusion();
    free2DArray(n,r);
    free2DArray(n,f); 
    free2DArray(n,e);
    free2DArray(n,vel); 
    free2DArray(n,angvel);

}
