#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#include <ctime>
#include <iomanip>
#include <cassert>
#include "./lrc_module.hpp"
#include "./maths_module.hpp"
#include "./mc_poly_lj_module.hpp"
#include "./averages_module.hpp"
#include "./config_io_module.hpp"

#define nblock        10
#define nstep         10
#define temperature   1.0
#define r_cut         2.5
#define dr_max        0.05
#define de_max        0.05

std::vector<VariableType> calc_variables(PotentialType tot, int n, double box, double m_ratio){

    // Preliminary calculations (n,r,total are taken from the calling program)
    double vol = pow(box,3);                     //  Volume
    double rho = n / vol;                        //  Density

    // Move acceptance ratio
    VariableType m_r;
    m_r.nam = "Move ratio";
    m_r.val = m_ratio;
    m_r.instant = false;

    // Internal energy per molecule (shifted-force potential)
    // Ideal gas contribution (assuming nonlinear molecules) plus total PE divided by N
    VariableType e_sf;
    e_sf.nam = "E/N shifted force";
    e_sf.val = 3.0*temperature + tot.pot/n;

    // Pressure (shifted-force potential)
    // Ideal gas contribution plus total virial divided by V
    VariableType p_sf;
    p_sf.nam = "P shifted force";
    p_sf.val = rho*temperature + tot.vir/vol;

    // list the VariableType objects
    std::vector<VariableType> variables;

    variables.push_back(m_r);
    variables.push_back(e_sf);
    variables.push_back(p_sf);

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

    // Bond vectors in body-fixed frame
    // Isosceles triangle, 3 sites, with unit bond length and bond angle alpha
    // which we set to 75 degrees here
    double** db    = allocate2DArray(3,3);
    double** db2   = allocate2DArray(3,3);
    double alpha   = 75.0 * M_PI / 180.0;
    double alpha2  = alpha / 2.0;
    
    double dbb[3][3] = { -sin(alpha2) , 0.0, -cos(alpha2)/3.0,0.0, 0.0, 2.0*cos(alpha2)/3.0, sin(alpha2), 0.0,-cos(alpha2)/3.0};
    
    for (int i{0};i<3;++i){
        for (int j{0};j<3;++j){
            db[i][j] = dbb[i][j];
        }
    }
    
    elementWise2DProduct(3,3,db,db,db2);
    double* db_sum = new double[3];
    sum2DArrays(3,3,db2,db_sum,1);
    double max = max1DArray(3,db_sum);
    double diameter = 2.0 * sqrt(max); // Molecular diameter

    // Preliminary calculations (n,r,total are taken from the calling program)
    const char* file = "cnf_poly.inp";

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
    double** ea     = allocate2DArray(3,3);
    double** e      = allocate2DArray(n,4);
    double** vel    = allocate2DArray(n,3);
    double** angvel = allocate2DArray(n,3);
    double*** d     = allocate3DArray(n,na,3);
    bool quaternion = true;
    bool with_v     = false;

    read_cnf_mols(file, quaternion, with_v, r, e, vel, angvel);

    
    scalar2DArrayDivision(n,3,box,r);    // Convert positions to box units
    rint2D(n,3, r);                      // Periodic boundaries
    
    // Calculate all bond vectors
    for (int i{0};i<n;++i){   

        double* ei  = new double[4];
	for (int j{0}; j<4; ++j){
	    ei[j] = e[i][j];
	}
	q_to_a(ei,ea);                  // Rotation matrix for i
	for (int j{0}; j<3; ++j){
	    for (int k{0}; k<3; ++k){
	        for (int l{0}; l<3; ++l){
		    d[i][j][k] += db[j][l] * ea[l][k];      // NB: equivalent to ai_T*db, ai_T=transpose of ai
		}
	    }
	}
        delete [] ei;
    }

    // Initial energy and overlap check
    PotentialType total;

    potential (total, n, box, r_cut, diameter, r, d);

    assert (!total.ovr); 
    std::cout << "No overlap in initial configuration! \n";

    BlockVar blk_var;

    std::cout << '\n';
    std::cout << "mc_nvt_poly_lj \n";
    std::cout << "Monte Carlo, constant-NVT ensemble, polyatomic molecule \n";
    std::cout << "Simulation uses cut and not shifted potential \n";
    std::cout << '\n';

    introduction(na, db, diameter); 
    srand (time(NULL));

//  Write out parameters
    std::cout << '\n';
    printf("%16s %42d   \n", "Number of blocks",          nblock);
    printf("%25s %33d   \n", "Number of steps per block", nstep);
    printf("%20s %37.6f \n", "Specified temperature",     temperature);
    printf("%20s %36.6f \n", "Maximum r displacement",    dr_max);
    printf("%20s %36.6f \n", "Maximum e displacement",    de_max);
    printf("%10s %39d   \n", "Number of particles",n);
    printf("%10s %48.6f \n", "Box length", box);
    printf("%7s  %50.6f \n", "Density", n/pow(box,3));
    std::cout << '\n'; 

    // zero out the move ratio
    double m_ratio = 0.0;
    int n_avg = calc_variables(total,n,box,m_ratio).size();

    run_begin  (calc_variables(total,n,box,m_ratio), blk_var, ti);

    for (int blk{0}; blk < nblock;++blk){                                // Loop over blocks

	blk_begin(n_avg,blk_var);

        for (int stp{0}; stp<nstep;++stp){                               // Loop over steps

            double moves = 0.0;
            for (int atm{0};atm <n;++atm){                               // Loop over atoms

		double** rj  = allocate2DArray(n-1,ndim);
		double*** dj = allocate3DArray(n-1,ndim,ndim);
		double** d_i = allocate2DArray(ndim,ndim);
		double** di  = allocate2DArray(ndim,ndim);
	        double*  r_i = new double[ndim];
	        double*  ei  = new double[4];
	        double*  e_i = new double[4];
	        double*  ri  = new double[ndim];

		remove2DArray(n,atm,r,rj);
		remove3DArray(n, ndim, atm, d, dj);

                for (int j{0};j<ndim;++j)
                   r_i[j] = r[atm][j];

                for (int j{0};j<4;++j)
                   e_i[j] = e[atm][j];
		
                for (int j{0};j<ndim;++j){
                    for (int k{0};k<ndim;++k){
                        d_i[j][k] = d[atm][j][k];
		    }
		}

		PotentialType partial_old;
		potential_1 (partial_old, n-1, r_cut, diameter, r_i, d_i, box, rj, dj );  // Old atom potential, virial etc

		random_translate_vector(dr_max/box, r_i, ri);       // Trial move to new position (in box=1 units)
		rint1D(3,ri);                                       // Periodic boundary correction
		random_rotate_quaternion (de_max, e_i, ei);         //Trial rotation
                q_to_a ( ei, ea );                                  // Rotation matrix for i
		matMultip(ndim, ndim, db, ea, di);                  // NB: equivalent to ai_T*db, ai_T=transpose of ai

		PotentialType partial_new;
		potential_1 (partial_new, n-1, r_cut, diameter, ri, di, box, rj, dj );    // New atom potential, virial etc
		
                if (!partial_new.ovr){                                               // Test for non-overlapping configuration
                    double delta;
                    delta = partial_new.pot - partial_old.pot;                       // Use cut (but not shifted) potential
                    delta = delta / temperature;

                    if (metropolis (delta)){                                         // Accept Metropolis test

                        total.pot = total.pot + partial_new.pot - partial_old.pot;   // Update total values
                        total.vir = total.vir + partial_new.vir - partial_old.vir;   // Update total values
                        total.ovr = total.ovr + partial_new.ovr - partial_old.ovr;   // Update total values

    	                update2DArray(ndim,atm,ri,r);                                // Update position
    	                update2DArray(4,atm,ei,e);                                   // Update quaternion
			update3DArray(ndim,ndim,atm,di, d);                          // Update bond vectors
                        moves = moves + 1;                                           // Increment move counter
		        }
		    }
                    free2DArray(n-1,rj);
                    free3DArray(n-1,ndim,dj);
                    free2DArray(ndim,d_i);
                    free2DArray(ndim,di);
                    delete []  r_i;
                    delete []  ei;
                    delete []  e_i;
                    delete []  ri;
	    }
            m_ratio = moves / n;

	    blk_add (calc_variables(total,n,box,m_ratio),blk_var);
	}

        blk_end ( blk, n_avg, blk_var);
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(3) << std::to_string(blk+1);
        std::string sav_tag(ss.str());
        double** out_r = allocate2DArray(n,3);
        scalar2DArrayMultip(n,3,box,r,out_r);
	write_cnf_mols ("cnf."+sav_tag, n, box, quaternion, with_v, out_r, e, vel, angvel);
        free2DArray(n,out_r);
    }

    run_end (calc_variables(total,n,box,m_ratio), blk_var, ti);

    potential (total, n, box, r_cut, diameter, r, d);
    assert (!total.ovr);
    std::cout << "No overlap in final configuration! \n";

    double** out_r = allocate2DArray(n,3);
    scalar2DArrayMultip(n,3,box,r,out_r);
    write_cnf_mols ("cnf.out", n, box, quaternion, with_v, out_r, e, vel, angvel);
    free2DArray(n,out_r);
    deletePointer(calc_variables(total,n,box,m_ratio), blk_var);

    conclusion();
    delete [] db_sum;
    free2DArray(n,r);
    free2DArray(3,ea);
    free2DArray(3,db); 
    free2DArray(3,db2); 
    free2DArray(n,e);
    free2DArray(n,vel);
    free2DArray(n,angvel);
    free3DArray(n,na,d);
}
