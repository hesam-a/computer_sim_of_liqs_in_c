#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <cassert>
#include <numeric>
#include <stdlib.h>
#include <iterator>
#include "./math_module.hpp"

// ------  Routines for MC simulation, polyatomic molecule, LJ atoms  --------


// Cutoff distance and force-shift parameters (all private) chosen as per the reference:
// S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)

int    na       = 3;
double r_cut    = 2.612;    // in sigma=1 units, where r_cut = 1.2616 nm, sigma = 0.483 nm
double sr_cut   = 1.0/r_cut;
double sr_cut6  = pow(sr_cut,6);
double sr_cut12 = pow(sr_cut6,2);
double lambda1  = 4.0*(7.0*sr_cut6-13.0*sr_cut12);
double lambda2  = -24.0*(sr_cut6-2.0*sr_cut12)*sr_cut;

class PotentialType{

    public:
	double pot;  // = 0.0; 
	double vir;  // = 0.0;
	double lap;  // = 0.0;
	bool   ovr;  // = false;

};

void introduction(int na,double** db,double diameter){ 

/*  Prints out introductory statements at start of run. */

    std::cout << "Lennard-Jones potential" << '\n'; 
    std::cout << "Cut-and-force-shifted) " << '\n'; 
    std::cout << "Diameter, sigma = 1"     << '\n'; 
    std::cout << "Well depth, epsilon = 1" << '\n'; 

    printf("%20s %30d \n", "Number of atoms per molecule", na);
    std::cout << '\n'; 
    
    for (size_t i{0}; i<3;++i){
        printf("%s %2ld", "Body-fixed atom vector",i+1);
        for (size_t j{0}; j<3;++j){
            printf("%15.6f ", db[i][j]);
        }
        std::cout << '\n';
    }
    std::cout << '\n'; 
    printf("%s %40.6f \n", "Molecular diameter", diameter);
    printf("%s %53.6f \n", "r_cut", r_cut);
    printf("%s %39.6f \n", "Force-shift lambda1", lambda1);
    printf("%s %39.6f \n", "Force-shift lambda2", lambda2);
    std::cout << '\n'; 

}

void conclusion(){
/*  Prints out concluding statements at end of run. */

	std::cout << "Program ends \n";
	std::cout << "\n";
	std::cout << "\n";

}

void potential_1 (PotentialType &partial,int mm, double r_cut, double diameter, double* ri, double** ddi, double box, double** r, double*** dd ){

/*  partial.pot is the nonbonded cut (not shifted) potential energy of atom ri with a set of other atoms
    partial.vir is the corresponding virial of atom ri
    partial.ovr is a flag indicating overlap (potential too high) to avoid overflow
    If this is True, the values of partial.pot etc should not be used
    In general, r & d will be subsets of the complete set of simulation coordinates & bond vectors
    and none of their rows should be identical to ri, di */

/*  It is assumed that r has been divided by box
    Results are in LJ units where sigma = 1, epsilon = 1
    Note that this is the force-shifted LJ potential with a linear smoothing term
    S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002) */


    double rij_sq,sr2,sr6,sr12,pot,vir,lap,rmag,virab;
    bool   ovr;
    double sr2_ovr       = 1.77;      // Overlap threshold (pot > 100)
    double rm_cut_box    = ( r_cut + diameter ) / box; // Molecular cutoff in box=1 units
    double rm_cut_box_sq = pow(rm_cut_box,2);
    double r_cut_sq      = pow(r_cut,2);
    int ndim = 3;
    assert (rm_cut_box<0.5); //rm_cut/box too large

    partial.pot = 0.0;
    partial.vir = 0.0;
    partial.lap = 0.0;
    partial.ovr = false;

    for (int j{0};j<mm;++j){

	double* r_j  = new double[ndim]{};
        double* rij  = new double[ndim]{};
        double* rij2 = new double[ndim]{};

	for (int i{0};i<ndim;++i)
	   r_j[i] = r[j][i];

        subtract1DArrays(3,ri,r_j,rij);               // Separation vector
        rint1D(ndim,rij);                             // Periodic boundary conditions in box=1 units
        elementWise1DProduct(3,rij,rij,rij2);
        rij_sq = elementSum1D(3,rij2);                // Squared separation

        if (rij_sq < rm_cut_box_sq){                  // Check within cutoff
	    scalar1DArrayMultip(ndim,box,rij,rij);    // Now in sigma=1 units

            for (int a{0};a<na;++a){
                for (int b{0};b<na;++b){
		    double* rab  = new double[ndim]{};
		    double* rab2 = new double[ndim]{};
		    double* fab  = new double[ndim]{};

		    for (int i{0};i<ndim;++i)
			rab[i] = rij[i] + ddi[a][i] - dd[j][b][i]; // Atom-atom vector, sigma=1 units
                    elementWise1DProduct(3,rab,rab,rab2);
		    double rab_sq = elementSum1D(3,rab2);   // Squared atom-atom separation, sigma=1 units 

                    if (rab_sq < r_cut_sq){       // Test within potential cutoff 
                        sr2    = 1.0 / rab_sq;   // (sigma/rab)**2
                        ovr    = sr2 > sr2_ovr;  // Overlap if too close
                        if (ovr)
                            partial.ovr=true;

                        rmag  = sqrt(rab_sq);
                        sr6   = pow(sr2,3);
                        sr12  = pow(sr6,2);
                        pot   = 4.0*(sr12-sr6) + lambda1 + lambda2*rmag; // LJ atom-atom pair potential (force-shifted)
                        virab = 24.0*(2.0*sr12-sr6) - lambda2*rmag;      // LJ atom-atom pair virial
	                scalar1DArrayMultip(ndim,virab*sr2,rab,fab);     // LJ atom-atom pair force
			elementWise1DProduct(ndim,rij,fab,fab);
		        vir = elementSum1D(3,fab);
		        
                        partial.pot = partial.pot + pot; 
                        partial.vir = partial.vir + vir;
                        partial.ovr = partial.ovr + ovr;
                        partial.lap = partial.lap + lap;
		    }
		    delete [] rab;
		    delete [] rab2;
		    delete [] fab;
		}
	    }
	}
        delete [] rij2;
	delete [] rij;
	delete [] r_j;
    }
	
    // Include numerical factors
    partial.vir = partial.vir  / 3.0;         // divide virial by 3

}

void potential (PotentialType &total, int mm, double box,double r_cut, double diameter, double** r, double*** dd){
/*  Takes in box and r & d arrays, and calculates total potential etc.

    The results are returned as total, a PotentialType variable.
    
    Actual calculation performed by function potential_1 */

    int ndim = 3 ;

    total.pot = 0.0;
    total.vir = 0.0;
    total.lap = 0.0;
    total.ovr = false;

    double m = mm -1;
    for(int i{0};i<mm-1;++i){

        PotentialType part;

	double*   r_i    = new double[ndim];
        double**  dd_i   = allocate2DArray(ndim,ndim);
	double**  r_all  = allocate2DArray(m,ndim);
	double*** dd_all = allocate3DArray(m,ndim,ndim);

	for (int j{0};j<ndim;++j)
	   r_i[j] = r[i][j];

        int ii = i;
        for(int j=0;j<m;j++){
                ii += 1;
            for(int k{0}; k<ndim; ++k){
                r_all[j][k] = r[ii][k];
            }
        }

        for (int j{0};j<ndim;++j){
            for (int k{0};k<ndim;++k){
                dd_i[j][k] = dd[i][j][k];
            }
        }

        ii = i;
        for(int j=0;j<m;j++){
                ii += 1;
            for(int k{0}; k<ndim; ++k){
                for (int l{0};l<ndim;++l){
                    dd_all[j][k][l] = dd[ii][k][l];
                }
            }
        }

	potential_1 (part, m,r_cut, diameter,r_i, dd_i, box, r_all, dd_all );
        if (part.ovr){
            total.ovr = true;
            break;
	}
        total.pot = total.pot + part.pot;
        total.vir = total.vir + part.vir;
        total.ovr = total.ovr + part.ovr;
        total.lap = total.lap + part.lap;

	m -= 1;

        delete [] r_i;
        free3DArray(m+1,ndim,dd_all);
        free2DArray(m+1,r_all); 
        free2DArray(ndim,dd_i);

    }
}
