#include <iostream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <cassert>
#include <numeric>
#include <stdlib.h>
#include <iterator>
#include "./math_module.h"

class PotentialType{

    public:
	double pot; 
	double vir;
	double lap;
	bool   ovr;

};

void introduction(){

/*  Prints out introductory statements at start of run. */

	std::cout << "Lennard-Jones potential" << '\n'; 
    	std::cout << "Cut (but not shifted) "  << '\n'; 
    	std::cout << "Diameter, sigma = 1"     << '\n'; 
    	std::cout << "Well depth, epsilon = 1" << '\n'; 
}

void conclusion():{
/*  Prints out concluding statements at end of run. */

	std::cout << "Program ends";
}

void potential (PotentialType &total, int mm, double box, double r_cut, double** r ){
/*  Takes in box, cutoff range, and coordinate array, and calculates total potential etc.
    The results are returned as total, a PotentialType variable. */

    // Actual calculation performed by function potential_1

    int d = 3;

    total.pot = 0.0;
    total.vir = 0.0;
    total.lap = 0.0;
    total.ovr = false;

    double m = mm -1;
    for(int i{0};i<mm-1;++i){

        PotentialType part;
	double* r_i    = new double[3];
	double** r_all = allocate2DArray(m,d);

	for (int j{0};j<d;++j){
	   r_i[j] = r[i][j];
	}

        int ii = i;
        for(int j=0;j<m;j++){
                ii += 1;
            for(int k{0}; k<d; ++k){
                r_all[j][k] = r[ii][k];
            }
        }

	potential_1 (part, m, r_i, box, r_cut, r_all );
        if (part.ovr){
            total.ovr = true;
            break;
	}
        total.pot = total.pot + part.pot;
        total.vir = total.vir + part.vir;
        total.lap = total.lap + part.lap;

	m -= 1;
    }
}


void potential_1 (PotentialType &partial,int mm, double* ri, double box, double r_cut, double** r){
/*  Takes in coordinates of an atom and calculates its interactions.

    Values of box, cutoff range, and partner coordinate array are supplied.
    The results are returned as partial, a PotentialType variable. */


    double rij_sq,sr2,sr6,sr12,pot,vir,lap;
    bool ovr;
    double sr2_ovr      = 1.77;      // Overlap threshold (pot > 100)
    double r_cut_box    = r_cut / box;
    double r_cut_box_sq = pow(r_cut_box,2);
    double box_sq       = pow(box,2);

    int d = 3;

    partial.pot = 0.0;
    partial.vir = 0.0;
    partial.lap = 0.0;
    partial.ovr = false;

    for (int i{0};i<mm;++i){

        double* rij = new double[3];
	double* r_j = new double[3];

	for (int j{0};j<d;++j){
	   r_j[j] = r[i][j];
	}

        rij = subtract1DArrays(3,ri,r_j);   // Separation vector
        rij = rint(d,rij);                   // Periodic boundary conditions in box=1 units
        rij_sq = elementSum1D(3,elementWise1DProduct(3,rij,rij));       // Squared separation
    
        if (rij_sq < r_cut_box_sq){    // Check within cutoff
            rij_sq = rij_sq * box_sq;  // Now in sigma=1 units
            sr2    = 1.0 / rij_sq;     // (sigma/rij)**2
            ovr    = sr2 > sr2_ovr;    // Overlap if too close
    
            if (ovr){
                partial.ovr=true;
            }
    
            sr6  = pow(sr2,3);
            sr12 = pow(sr6,2);
            pot  = sr12 - sr6;                       // LJ pair potential (cut but not shifted)
            vir  = pot + sr12;                       // LJ pair virial
            lap  = ( 22.0*sr12 - 5.0*sr6 ) * sr2;    // LJ pair Laplacian
    
            partial.pot = partial.pot + pot; 
            partial.vir = partial.vir + vir;
            partial.lap = partial.lap + lap;
	}
    }
    partial.pot = partial.pot * 4.0;        // 4*epsilon
    partial.vir = partial.vir * 24.0 / 3.0; // 24*epsilon and divide virial by 3
    partial.lap = partial.lap * 24.0 * 2.0;
}

void force_sq (int mm,double box, double r_cut, double** r){
//   Calculates total squared force.

    double rij_sq,sr2,sr6,sr12,pot,vir,lap;
    double r_cut_box    = r_cut / box;
    double r_cut_box_sq = pow(r_cut_box,2);
    double box_sq       = pow(box,2);

    int d = 3;

    double** f     = allocate2DArray(n,d);
    double* fij    = new double[3];
    double* rij    = new double[3];
    double* rij_p  = new double[3];

    f = zeroMatrix(mm,d); // Initialize

    for(int i{0};i<mm;++i){
        for(int j{i+1};j<mm ;++j){
	    for (int k{0};k<3;++k){
		rij[k] = b[i][k] - b[j][k];                     // Separation vector
	    }

                rij    = rint(3,rij);                           // Periodic boundary conditions in box=1 units
                rij_p  = elementWise1DProduct(3,rij,rij);
                rij_sq = elementSum1D(3,rij_p);

                if (rij_sq < r_cut_box_sq){                         // Check within cutoff
                    rij_sq = rij_sq * box_sq;                       // Now in sigma=1 units 
                    rij    = scalar1DArrayMultip(3,box,rij);        // Now in sigma=1 units
                    sr2    = 1.0 / rij_sq;  
                    sr6    = pow(sr2,3);
                    sr12   = pow(sr6,2);
                    fij    = scalar1DArrayMultip(3, (2.0*sr12 - sr6),rij);
                    fij    = scalar1DArrayMultip(3,sr2,fij);
                    for (int k{0};k<3;++k){
                        f[i][k] = f[i][k] + fij[k];
                        f[j][k] = f[j][k] - fij[k];
                    }
                }
        }
    }
	    
    f = scalar2DArrayMultip(mm,d,24.0,f);
    std::cout << "f " << '\n';
    print2DArray(mm,3,f);

    return f; 
}
