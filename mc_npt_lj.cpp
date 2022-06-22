#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#include <ctime>
#include <cassert>
#include "./lrc_module.hpp"
#include "./math_module.hpp"
#include "./mc_lj_module.hpp"
#include "./averages_module.hpp"
#include "./read_atoms.hpp"

#define nblock        10
#define nstep         100
#define temperature   1.0
#define pressure      0.69
#define r_cut         2.5
#define dr_max        0.15
#define db_max        0.025



std::vector<VariableType> calc_variables(PotentialType tot, double** r, int n, double box, double m_ratio,double v_ratio){
/*  ---- Some notes from Allen-Tildesley ----
  
    Calculates all variables of interest.

    They are collected and returned as a vector, for use in the main program.

     In this example we simulate using the cut (but not shifted) potential
     Accordingly, < p_c > should match the input pressure and the values
     of < p_c >, < e_c > and density should be consistent (for this potential)
     For comparison, long-range corrections are also applied to give
     estimates of < e_f > and < p_f > for the full (uncut) potential. */


    // Preliminary calculations (n,r,total are taken from the calling program)
    double vol = pow(box,3);                     //  Volume
    double rho = n / vol;                        //  Density
    double fsq = force_sq (n, box, r_cut, r );   // Total squared force

    // Initial energy and overlap check
    potential (tot,n, box, r_cut, r);

    // Move acceptance ratio
    VariableType* m_r;
    m_r = new VariableType();
    m_r->nam = "Move ratio";
    m_r->val = m_ratio;
    m_r->instant = false;

    // Volume acceptance ratio
    VariableType* v_r;
    v_r = new VariableType();
    v_r->nam = "Volume ratio";
    v_r->val = v_ratio;
    v_r->instant = false;

    // Density
    VariableType* density;
    density = new VariableType();
    density->nam = "Density";
    density->val = rho;

    // Internal energy per atom for simulated, cut, potential
    // Ideal gas contribution plus cut (but not shifted) PE divided by N
    VariableType* e_c;
    e_c = new VariableType();
    e_c->nam = "E/N cut";
    e_c->val = 1.5*temperature + tot.pot/n;

    // Internal energy per atom for full potential with LRC
    // LRC plus ideal gas contribution plus cut (but not shifted) PE divided by N
    VariableType* e_f;
    e_f = new VariableType();
    e_f->nam = "E/N full";
    e_f->val = potential_lrc(rho,r_cut) + 1.5*temperature + tot.pot/n;

    // Pressure for simulated, cut, potential
    // Delta correction plus ideal gas contribution plus total virial divided by V
    VariableType* p_c;
    p_c = new VariableType();
    p_c->nam = "P cut";
    p_c->val = pressure_delta(rho,r_cut) + rho*temperature + tot.vir/vol;

    // Pressure for full potential with LRC
    // LRC plus ideal gas contribution plus total virial divided by V
    VariableType* p_f;
    p_f = new VariableType();
    p_f->nam = "P full";
    p_f->val = pressure_lrc(rho,r_cut) + rho*temperature + tot.vir/vol;

    // Configurational temperature
    // Total squared force divided by total Laplacian
    VariableType* t_c;
    t_c = new VariableType();
    t_c->nam = "T config";
    t_c->val = fsq/tot.lap;

/*  Heat capacity (cut but not shifted)
    MSD of excess "enthalpy" divided by temperature and sqrt(N) to make result intensive
    NB this is not really the excess Cp/NkB, it simply omits the kinetic energy fluctuations
    i.e. we add the ideal gas part of Cv/NkB, 1.5, to get total Cp/NkB */
    double enp = tot.pot+pressure*vol;
    VariableType* c_c;
    c_c = new VariableType();
    c_c->nam     = "Cp/N cut";
    c_c->val     = enp/(temperature*sqrt(n));
    c_c->method  = msd;
    c_c->add     = 1.5;
    c_c->instant = false;

/*  Heat capacity (full)
    MSD of excess "enthalpy" divided by temperature and sqrt(N) to make result intensive
    NB this is not really the excess Cp/NkB, it simply omits the kinetic energy fluctuations
    i.e. we add the ideal gas part of Cv/NkB, 1.5, to get total Cp/NkB */
    double enpf = n*potential_lrc(rho,r_cut)+tot.pot+pressure*vol;
    VariableType* c_f;
    c_f = new VariableType();
    c_f->nam     = "Cp/N full";
    c_f->val     = enpf/(temperature*sqrt(n));
    c_f->method  = msd;
    c_f->add     = 1.5;
    c_f->instant = false;

    // Volume MSD
    VariableType* vol_msd;
    vol_msd = new VariableType();
    vol_msd->nam     = "Volume MSD";
    vol_msd->val     = vol;
    vol_msd->method  = msd;
    vol_msd->instant = false;


    // list the VariableType objects
    std::vector<VariableType> variables;

    variables.push_back(*m_r);
    variables.push_back(*v_r);
    variables.push_back(*density);
    variables.push_back(*e_c);
    variables.push_back(*p_c);
    variables.push_back(*e_f);
    variables.push_back(*p_f);
    variables.push_back(*t_c);
    variables.push_back(*c_c);
    variables.push_back(*c_f);
    variables.push_back(*vol_msd);

    return variables;

    delete [] m_r;
    delete [] v_r;
    delete [] density;
    delete [] e_c;
    delete [] p_c;
    delete [] e_f;
    delete [] p_f;
    delete [] t_c;
    delete [] c_c;
    delete [] c_f;
    delete [] vol_msd;
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

    double** r = allocate2DArray(n,3);
    r = read_cnf_atoms(file,r);

    r = scalar2DArrayDivision(n,3,box,r);      // Convert positions to box units
    r = rint2D(n,3,r);                             // Periodic boundaries


    BlockVar blk_var;

    std::cout << '\n';
    std::cout << "mc_npt_lj \n";
    std::cout << "Monte Carlo, constant-NPT ensemble \n";
    std::cout << "Simulation uses cut (but not shifted) potential \n";
    std::cout << '\n';

    introduction();
//  Write out parameters
    std::cout << '\n';
    printf("%16s %23d   \n", "Number of blocks",          nblock);
    printf("%25s %14d   \n", "Number of steps per block", nstep);
    printf("%20s %18.6f \n", "Specified temperature",     temperature);
    printf("%25s %14.6f \n", "Potential cutoff distance", r_cut);
    printf("%20s %19.6f \n", "Maximum displacement",      dr_max);
    std::cout << '\n';

    std::cout << '\n';
    printf("%18s  %19d   \n", "Number of particles",n);
    printf("%18s  %20.6f \n", "Box length", box);
    printf("%18s  %20.6f \n", "Density", n/pow(box,3));
    std::cout << '\n';


    // Initial energy and overlap check
    PotentialType total;

    potential (total,n, box, r_cut, r);

    assert (!total.ovr);
    std::cout << "No overlap in initial configuration! \n";

    // zero out the move ratio
    double m_ratio = 0.0;
    double v_ratio;

    int n_avg = calc_variables(total, r, n, box, m_ratio, v_ratio).size();

    run_begin (calc_variables(total, r, n, box, m_ratio, v_ratio), blk_var, ti);

    for (int blk{0}; blk < nblock;++blk){                                // Loop over blocks

        blk_begin(n_avg,blk_var);

        for (int stp{0}; stp<nstep;++stp){                               // Loop over steps

            double moves = 0.0;
            for (int atm{0};atm <n;++atm){                               // Loop over atoms

                double** rj  = allocate2DArray(n-1,3);
                double*  r_i = new double[3];
                double*  ri  = new double[3];
                rj = remove2DArray(n,atm,r);                             // Array of all the other atoms

                for (int j{0};j<3;++j)
                   r_i[j] = r[atm][j];

                PotentialType partial_old;
                potential_1 (partial_old, n-1, r_i, box, r_cut, rj );    // Old atom potential, virial etc

                ri = random_translate_vector(dr_max/box, r_i );          // Trial move to new position (in box=1 units)
                ri = rint1D(3,ri);                                       // Periodic boundary correction

                PotentialType partial_new;
                potential_1 (partial_new, n-1, ri, box, r_cut, rj );                       // New atom potential, virial etc

                    if (!partial_new.ovr){                                                 // Test for non-overlapping configuration
                        double delta;
                        delta = partial_new.pot - partial_old.pot;                         // Use cut (but not shifted) potential
                        delta = delta / temperature;

                        if (metropolis (delta)){                                           // Accept Metropolis test
                            total.pot = total.pot + partial_new.pot - partial_old.pot;     // Update total values
                            total.vir = total.vir + partial_new.vir - partial_old.vir;     // Update total values
                            total.lap = total.lap + partial_new.lap - partial_old.lap;     // Update total values
                            total.ovr = total.ovr + partial_new.ovr - partial_old.ovr;     // Update total values

                            update2DArray(n, 3,r_i,ri,r);                                  // Update position
                            moves = moves + 1;                                             // Increment move counter
                        }
                    }
                delete [] r_i;
                delete [] ri ;
                //free2DArray(n-1,rj);
            }

            m_ratio = moves / n;
    
            v_ratio   = 0.0;                                // Zero volume move counter
    	    double zeta      = (double) rand()/RAND_MAX;    // Uniform random number in range (0,1)
            zeta             = 2.0*zeta-1.0;                // Now in range (-1,+1)
            double box_scale = exp(zeta*db_max);            // Sampling log(box) and log(vol) uniformly
            double box_new   = box*box_scale;               // New box (in sigma units)
            double den_scale = 1.0 / pow(box_scale,3);      // Density scaling factor

            PotentialType total_new;    
            potential (total_new,n, box_new, r_cut, r);     // New total energy, virial etc
    
            if (!total_new.ovr){                                              // Test for non-overlapping configuration
    	    double delta;
                delta = total_new.pot - total.pot;                            // Use cut (but not shifted) potential
                delta = delta + pressure * ( pow(box_new,3) - pow(box,3));    // Add PV term
                delta = delta / temperature;                                  // Divide by temperature
                delta = delta + (n+1) * log(den_scale);                       // Factor (n+1) consistent with log(box) sampling
    
                if (metropolis(delta)){      // Accept Metropolis test
                    total   = total_new;   // Update total values
                    box     = box_new;     // Update box
                    v_ratio = 1.0;         // Set volume move counter
    
    	        }
    	    }

        blk_add (calc_variables(total, r, n, box, m_ratio, v_ratio),blk_var);
	}

        blk_end ( blk, n_avg, blk_var);
    }

    run_end (calc_variables(total, r, n, box, m_ratio, v_ratio), blk_var, ti);

    conclusion();

    free2DArray(n,r);
}
