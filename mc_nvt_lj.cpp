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
#include "./config_io_module.hpp"


#define nblock        10
#define nstep         1000
#define temperature   1.0
#define r_cut         2.5
#define dr_max        0.15


std::vector<VariableType> calc_variables(PotentialType tot, double** r, int n, double box, double m_ratio){

    // Preliminary calculations (n,r,total are taken from the calling program)
    double vol = pow(box,3);                     //  Volume
    double rho = n / vol;                        //  Density
    double fsq = force_sq (n, box, r_cut, r );   // Total squared force
    // std::cout << " ---- force:   " << fsq << " ---- \n";

    // Initial energy and overlap check
    potential (tot,n, box, r_cut, r);

    // Move acceptance ratio
    VariableType m_r;
    m_r.nam = "Move ratio";
    m_r.val = m_ratio;
    m_r.instant = false;

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

    // Heat capacity (full)
    // MSD potential energy divided by temperature and sqrt(N) to make result intensive; LRC does not contribute
    // We add ideal gas contribution, 1.5, afterwards
    VariableType c_f;
    c_f.nam = "Cv/N full";
    c_f.val = tot.pot/(temperature*sqrt(n));
    c_f.method = msd;
    c_f.add = 1.5;
    c_f.instant = false;

    // list the VariableType objects
    std::vector<VariableType> variables;

    variables.push_back(m_r);
    variables.push_back(e_c);
    variables.push_back(p_c);
    variables.push_back(e_f);
    variables.push_back(p_f);
    variables.push_back(t_c);
    variables.push_back(c_f);


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

    std::cout << '\n';
    printf("%10s  %39d   \n", "Number of particles",n);
    printf("%10s  %48.6f \n", "Box length", box);
    printf("%7s   %50.6f \n", "Density", n/pow(box,3));
    std::cout << '\n'; 
    
    scalar2DArrayDivision(n,3,box,r);    // Convert positions to box units
    rint2D(n,3, r);                      // Periodic boundaries

    // Initial energy and overlap check
    PotentialType total;

    potential (total,n, box, r_cut, r);

    assert (!total.ovr); 
    std::cout << "No overlap in initial configuration! \n";

    BlockVar blk_var;

    std::cout << '\n';
    std::cout << "mc_nvt_lj \n";
    std::cout << "Monte Carlo, constant-NVT ensemble \n";
    std::cout << "Simulation uses cut (but not shifted) potential \n";
    std::cout << '\n';

    introduction();

//  Write out parameters
    std::cout << '\n';
    printf("%16s %43d   \n", "Number of blocks",          nblock);
    printf("%25s %34d   \n", "Number of steps per block", nstep);
    printf("%20s %38.6f \n", "Specified temperature",     temperature);
    printf("%25s %34.6f \n", "Potential cutoff distance", r_cut);
    printf("%20s %39.6f \n", "Maximum displacement",      dr_max);
    std::cout << '\n';

    // zero out the move ratio
    double m_ratio = 0.0;
    int n_avg = calc_variables(total, r, n, box, m_ratio).size();

    run_begin (calc_variables(total, r, n, box, m_ratio), blk_var, ti);

    for (int blk{0}; blk < nblock;++blk){                                // Loop over blocks

	blk_begin(n_avg,blk_var);

        for (int stp{0}; stp<nstep;++stp){                               // Loop over steps

            double moves = 0.0;
            for (int atm{0};atm <n;++atm){                               // Loop over atoms

		double** rj = allocate2DArray(n-1,3);
	        double*  r_i = new double[3];
	        double*  ri  = new double[3];

		remove2DArray(n,atm,r,rj);

                for (int j{0};j<3;++j){
                   r_i[j] = r[atm][j];
		}

		PotentialType partial_old;
	        potential_1 (partial_old, n-1, r_i, box, r_cut, rj );    // Old atom potential, virial etc

		random_translate_vector(dr_max/box, r_i, ri);          // Trial move to new position (in box=1 units)
		rint1D(3,ri);                                       // Periodic boundary correction

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
		free2DArray(n-1,rj);
	        delete [] r_i;
	        delete [] ri ;
	    }

        m_ratio = moves / n;

	blk_add (calc_variables(total, r, n, box, m_ratio),blk_var);
	}

    blk_end ( blk, n_avg, blk_var);
    }

    run_end (calc_variables(total, r, n, box, m_ratio), blk_var, ti);

    deletePointer(calc_variables(total, r, n, box, m_ratio), blk_var);

    conclusion();

    free2DArray(n,r);
}
