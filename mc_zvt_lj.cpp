#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iomanip>
#include <cassert>
#include "./lrc_module.hpp"
#include "./maths_module.hpp"
#include "./mc_lj_module.hpp"
#include "./averages_module.hpp"
#include "./config_io_module.hpp"

#define nblock        10
#define nstep         1000
#define temperature   1.0
#define activity      0.079
#define prob_move     0.34
#define r_cut         2.5
#define dr_max        0.15


std::vector<VariableType> calc_variables(PotentialType tot, double** r, int n, double box, double m_ratio,double c_ratio,double d_ratio){
/*  ---- Some notes from Allen-Tildesley ----
  
    Calculates all variables of interest.

    They are collected and returned as a list, for use in the main program.

    In this example we simulate using the cut (but not shifted) potential
    The values of < p_c >, < e_c > and < density > should be consistent (for this potential)
    For comparison, long-range corrections are also applied to give
    estimates of < e_f > and < p_f > for the full (uncut) potential
    The value of the cut-and-shifted potential is not used, in this example. */

    // Preliminary calculations (n,r,total are taken from the calling program)
    double vol = pow(box,3);                     //  Volume
    double rho = n / vol;                        //  Density
    double fsq = force_sq (n, box, r_cut, r );   // Total squared force

    // Initial energy and overlap check
    potential (tot,n, box, r_cut, r);

    // Move acceptance ratio
    VariableType m_r;
    m_r.nam = "Move ratio";
    m_r.val = m_ratio;
    m_r.instant = false;

    // creation acceptance ratio
    VariableType c_r;
    c_r.nam = "Create ratio";
    c_r.val = c_ratio;
    c_r.instant = false;

    // destruction acceptance ratio
    VariableType d_r;
    d_r.nam = "Destroy ratio";
    d_r.val = d_ratio;
    d_r.instant = false;

    // Number
    VariableType number;
    number.nam = "N";
    number.val = static_cast<double>(n);

    // Density
    VariableType density;
    density.nam = "Density";
    density.val = rho;

    // Internal energy per atom for simulated, cut, potential
    // Ideal gas contribution plus cut (but not shifted) PE divided by N
    VariableType e_c;
    e_c.nam = "E/N cut";
    e_c.val = 1.5*temperature + tot.pot/n;

    // Internal energy per atom for full potential with LRC
    // LRC plus ideal gas contribution plus cut (but not shifted) PE divided by N
    VariableType e_f;
    e_f.nam = "E/N full";
    e_f.val = potential_lrc(rho,r_cut) + 1.5*temperature + tot.pot/n;

    // Pressure for simulated, cut, potential
    // Delta correction plus ideal gas contribution plus total virial divided by V
    VariableType p_c;
    p_c.nam = "P cut";
    p_c.val = pressure_delta(rho,r_cut) + rho*temperature + tot.vir/vol;

    // Pressure for full potential with LRC
    // LRC plus ideal gas contribution plus total virial divided by V
    VariableType p_f;
    p_f.nam = "P full";
    p_f.val = pressure_lrc(rho,r_cut) + rho*temperature + tot.vir/vol;

    // Configurational temperature
    // Total squared force divided by total Laplacian
    VariableType t_c;
    t_c.nam = "T config";
    t_c.val = fsq/tot.lap;

    // Number MSD
    VariableType n_msd;
    n_msd.nam     = "N MSD";
    n_msd.val     = static_cast<double>(n);
    n_msd.method  = msd;
    n_msd.instant = false;

    // list the VariableType objects
    std::vector<VariableType> variables;

    variables.push_back(m_r);
    variables.push_back(c_r);
    variables.push_back(d_r);
    variables.push_back(number);
    variables.push_back(density);
    variables.push_back(e_c);
    variables.push_back(p_c);
    variables.push_back(e_f);
    variables.push_back(p_f);
    variables.push_back(t_c);
    variables.push_back(n_msd);

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


/*  Takes in a configuration of atoms (positions)
    Cubic periodic boundary conditions
    Conducts grand canonical Monte Carlo at the given temperature and activity
    Uses no special neighbour lists
   
    Reads several variables and options from standard input using JSON format
    Leave input empty "{}" to accept supplied defaults
   
    Positions r are divided by box length after reading in
    However, input configuration, output configuration, most calculations, and all results
    are given in simulation units defined by the model
    For example, for Lennard-Jones, sigma = 1, epsilon = 1
   
    Note that long-range corrections are not included in the acceptance/rejection
    of creation and destruction moves, but are added to the simulation averages
   
    Despite the program name, there is nothing here specific to Lennard-Jones
    The model is defined in mc_lj_module. */


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

    scalar2DArrayDivision(n,3,box,r);      // Convert positions to box units
    rint2D(n,3,r);                             // Periodic boundaries

    BlockVar blk_var;

    std::cout << '\n';
    std::cout << "mc_zpt_lj \n";
    std::cout << "Monte Carlo, constant-zPT ensemble \n";
    std::cout << "Simulation uses cut (but not shifted) potential \n";
    std::cout << '\n';

    double prob_create = (1-prob_move)/2; // So that create and destroy have equal probabilities
    introduction();

//  Write out parameters
    std::cout << '\n';
    printf("%16s %43d   \n", "Number of blocks",              nblock);
    printf("%25s %34d   \n", "Number of steps per block",     nstep);
    printf("%20s %38.6f \n", "Specified temperature",         temperature);
    printf("%18s %41.6f \n", "Specified activity",            activity); 
    printf("%19s %40.6f \n", "Probability of move",           prob_move);  
    printf("%25s %30.6f \n", "Probability of create/destroy", prob_create);
    printf("%20s %39.6f \n", "Maximum displacement",          dr_max);
    std::cout << '\n';

    printf("%10s  %39d   \n", "Number of particles",n);
    printf("%10s  %48.6f \n", "Box length", box);
    printf("%7s   %50.6f \n", "Density", n/pow(box,3));
    std::cout << '\n';


    // Initial energy and overlap check
    PotentialType total;

    potential (total,n, box, r_cut, r);

    assert (!total.ovr);
    std::cout << "No overlap in initial configuration! \n";
    std::cout << '\n';

    // zero out the move ratio
    double m_ratio = 0.0;
    double c_ratio = 0.0;
    double d_ratio = 0.0;

    int n_avg = calc_variables(total, r, n, box, m_ratio, c_ratio, d_ratio).size();
                                                                         
    run_begin (calc_variables( total, r, n, box, m_ratio, c_ratio, d_ratio), blk_var, ti);

    for (int blk{0}; blk < nblock;++blk){                                // Loop over blocks

    blk_begin(n_avg,blk_var);

        for (int stp{0}; stp<nstep;++stp){                               // Loop over steps
    
            double  m_try = 0., m_acc = 0.;
            double  c_try = 0., c_acc = 0.;
            double  d_try = 0., d_acc = 0.;
            int     ntry = n;              // Each step consists of ntry tries (during which N might vary)

            for (int itry{0};itry <n;++itry){
    
                double zeta = (double)rand()/RAND_MAX; // Uniform random number in range (0,1)
    
                if (zeta < prob_move){ 

                    m_try = m_try + 1;
                    int i = rand() % n + 0;         // Choose moving particle at random
		    int nn = n-1;

                    double** rj  = allocate2DArray(nn,3);
                    double*  r_i = new double[3];
                    double*  ri  = new double[3];
    
                    for (int j{0};j<3;++j)
                       r_i[j] = r[i][j];

                    remove2DArray(n,i,r,rj);                               // Array of all the other atoms
    	  	    PotentialType partial_old;
                    potential_1 (partial_old, nn, r_i, box, r_cut, rj );    // Old atom potential, virial etc
                    assert (!partial_old.ovr);
    
                    random_translate_vector(dr_max/box, r_i, ri);          // Trial move to new position (in box=1 units)
                    rint1D(3,ri);                                       // Periodic boundary correction

                    PotentialType partial_new;
                    potential_1 (partial_new, nn, ri, box, r_cut, rj );     // New atom potential, virial etc
    
                    if (!partial_new.ovr){                                                 // Test for non-overlapping configuration
                        double delta;
                        delta = partial_new.pot - partial_old.pot;                         // Use cut (but not shifted) potential
                        delta = delta / temperature;
    
                        if (metropolis (delta)){                                           // Accept Metropolis test
                            total.pot = total.pot + partial_new.pot - partial_old.pot;     // Update total pot 
                            total.vir = total.vir + partial_new.vir - partial_old.vir;     // Update total vir
                            total.lap = total.lap + partial_new.lap - partial_old.lap;     // Update total lap
                            total.ovr = total.ovr + partial_new.ovr - partial_old.ovr;     // Update total ovr 
                            update2DArrayZVT(n,3,r_i,ri,r);                                  // Update position		
                            m_acc = m_acc + 1;                                             // Increment move counter
			}
		    }
		delete [] ri;
		delete [] r_i;
		free2DArray(nn,rj);
		}

	        else if (zeta < prob_move + prob_create){      // Try create

		    int nn = n+1;
                    double*  ri  = new double[3];
                    c_try = c_try + 1;
                    rand1DArray(3,ri);                                         // Three uniform random numbers in range (0,1)
	            scalar1DArraySubtract(3, 0.5, ri);                         // Now in range (-0.5,+0.5) for box=1 units
                    PotentialType partial_new;
                    potential_1 (partial_new, n, ri, box, r_cut, r);           // New atom potential, virial etc
       
                    if (!partial_new.ovr){                                     //Test for non-overlapping configuratio
			double delta;
                        delta = partial_new.pot / temperature;                 // Use cut (not shifted) potential
                        delta = delta - log ( activity * pow(box,3)/(n+1));    // Activity term for creation
                        if (metropolis (delta)){                               // Accept Metropolis test
	            	    r = createArray(n,3, ri, r);                           // Add new particle to r array
                            n = n + 1;                                       // New value of N

                            total.pot = total.pot + partial_new.pot;         // Update total pot 
                            total.vir = total.vir + partial_new.vir;         // Update total vir
                            total.lap = total.lap + partial_new.lap;         // Update total lap
                            total.ovr = total.ovr + partial_new.ovr;         // Update total ovr
                            c_acc = c_acc + 1;	                         // Increment creation move counter
			}
		    }
		delete [] ri;
	        }

	        else{                                                       // Try destroy

		    int nn = n-1;
                    double*  r_i = new double[3];
                    double** rj  = allocate2DArray(nn,3);
                    int atm = rand() % n + 0;         // Choose moving particle at random
                    d_try = d_try + 1;

                    for (int j{0};j<3;++j)
                       r_i[j] = r[atm][j];

                    remove2DArray(n,atm,r,rj);                               // Array of all the other atoms

                    PotentialType partial_old;
                    potential_1 (partial_old, nn, r_i, box, r_cut, rj );    // Old atom potential, virial etc
                    assert (!partial_old.ovr);

	            double delta;
                    delta = -partial_old.pot / temperature;                  // Use cut (not shifted) potential
                    delta = delta - log ( n / (activity * pow(box,3)) );     // Activity term for creation
                    if (metropolis ( delta )){                               // Accept Metropolis test

	                r = annihilateArray(n,3, r_i, r);                        // Delete particle from r array
                        n = n - 1;                                           // New value of N
                        total.pot = total.pot - partial_old.pot;             // Update total pot 
                        total.vir = total.vir - partial_old.vir;             // Update total vir
                        total.lap = total.lap - partial_old.lap;             // Update total lap
                        total.ovr = total.ovr - partial_old.ovr;             // Update total ovr
                        d_acc = d_acc + 1;                                    // Increment destruction move counter
	            }
		    delete [] r_i;
		    free2DArray(nn,rj);
	        }
	    }
	    if (m_try>0) m_ratio = m_acc/m_try; else m_ratio =0.0;
	    if (c_try>0) c_ratio = c_acc/m_try; else c_ratio =0.0;
	    if (d_try>0) d_ratio = d_acc/m_try; else d_ratio =0.0;
            blk_add (calc_variables(total, r, n, box, m_ratio, c_ratio, d_ratio),blk_var);
        }

        blk_end ( blk, n_avg, blk_var);
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(3) << std::to_string(blk+1);
        std::string sav_tag(ss.str());
        double** out_r = allocate2DArray(n,3);
        scalar2DArrayMultip(n,3,box,r,out_r);
        write_cnf_atoms ("cnf."+sav_tag, n, box,out_r );
        free2DArray(n,out_r);

    }

    run_end (calc_variables(total, r, n, box, m_ratio, c_ratio, d_ratio), blk_var, ti);

    potential (total,n, box, r_cut, r);
    assert (!total.ovr);
    std::cout << "No overlap in final configuration! \n";

    double** out_r = allocate2DArray(n,3);
    scalar2DArrayMultip(n,3,box,r,out_r);
    write_cnf_atoms ("cnf.out", n, box,out_r );
    free2DArray(n,out_r);
    deletePointer(calc_variables(total, r, n, box, m_ratio, c_ratio, d_ratio), blk_var);

    conclusion();

    free2DArray(n,r);

}
