#include <iostream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <cassert>
#include <numeric>
#include <stdlib.h>
#include <iterator>
#include "./maths_module.hpp"

// Force routine for MD simulation, Lennard-Jones atoms.

class PotentialType{

    public:
    double cut;  // = 0.0;
	double pot;  // = 0.0; 
	double vir;  // = 0.0;
	double lap;  // = 0.0;
	bool   ovr;  // = false;

};

void introduction(){

/*  Prints out introductory statements at start of run. */

	  std::cout << "Lennard-Jones potential" << '\n'; 
    std::cout << "Cut-and-shifted version for dynamics"  << '\n'; 
    std::cout << "Diameter, sigma = 1"     << '\n'; 
    std::cout << "Well depth, epsilon = 1" << '\n';

}

void conclusion(){
/*  Prints out concluding statements at end of run. */

	std::cout << "Program ends \n";
	std::cout << "\n";
	std::cout << "\n";

}


double** force (int mm,double box, double r_cut, double** r){
//   Takes in box, cutoff range, and coordinate array, and calculates forces and potentials etc.

    double rij_sq;
    double sr2_ovr = 1.77;
    double r_cut_box    = r_cut / box;
    double r_cut_box_sq = pow(r_cut_box,2);
    double box_sq       = pow(box,2);

    // Calculate potential at cutoff
    double sr2     = 1.0 / pow(r_cut,2); // in sigma=1 units
    double sr6     = pow(sr2,3);
    double sr12    = pow(sr6,2);
    double pot_cut = sr12 - sr6;         // Without numerical factor 4
    bool ovr;
    double cut;
    double vir;
    double pot;
    double lap;

    int d = 3;

    double** f     = allocate2DArray(mm,d);
    double* fij    = new double[3];
    double* rij    = new double[3];
    double* rijb   = new double[3];
    double* rij_p  = new double[3];

    // Initialize f
    zeroMatrix(mm,d,f);
    
    PotentialType total;

    total.cut = 0.0;
    total.pot = 0.0;
    total.vir = 0.0;
    total.lap = 0.0;
    total.ovr = false;

    for(int i{0};i<mm;++i){
        for(int j{i+1};j<mm ;++j){
	        for (int k{0};k<3;++k){
		        rij[k] = r[i][k] - r[j][k];            // Separation vector
	        }
            rint1D(3,rij);                             // Periodic boundary conditions in box=1 units

            elementWise1DProduct(3,rij,rij, rij_p);
            rij_sq = elementSum1D(3,rij_p);            // Squared separation
            if (rij_sq < r_cut_box_sq){                // Check within cutoff
                rij_sq = rij_sq * box_sq;              // Now in sigma=1 units 
                scalar1DArrayMultip(3,box,rij,rijb);   // Now in sigma=1 units
                sr2    = 1.0 / rij_sq;                 // (sigma/rij)**2
                ovr = sr2 > sr2_ovr;                   // Overlap if too close

                sr6    = pow(sr2,3);                   
                sr12   = pow(sr6,2);                   
                cut    = sr12 - sr6;                            // LJ pair potential (cut but not shifted)
                vir    = cut + sr12;                            // LJ pair virial
                pot    = cut - pot_cut;                         // LJ pair potential (cut and shifted)          
                lap    = (22.0 * sr12 - 0.5 * sr6) * sr2;       // LJ pair laplacian
                scalar1DArrayMultip(3,(vir * sr2),rijb,fij);    // LJ pair forces
              
                total.cut = total.cut + cut;
                total.pot = total.pot + pot;
                total.vir = total.vir + vir;
                total.lap = total.lap + lap;
                total.ovr = total.ovr + ovr;
              
                for (int k{0};k<3;++k){
                    f[i][k] = f[i][k] + fij[k];
                    f[j][k] = f[j][k] - fij[k];
                }
            }
        }
    }
	    
    scalar2DArrayMultip(mm,d,24.0,f,f);     // 24 *epsilon
    total.cut = total.cut * 4.0;            // 4  *epsilon
    total.pot = total.pot * 4.0;            // 4  *epsilon
    total.vir = total.vir * 24.0 / 3.0;     // 24 *epsilon and divided by virial by 3
    total.lap = total.lap * 24.0 / 2.0;     // 4  *epsilon and factor 2 for ij and ji

    free2DArray(mm,f);

    delete [] fij; 
    delete [] rij;  
    delete [] rijb;
    delete [] rij_p;

    return f; 
}


double hessian (int mm,double box, double r_cut, double** r, double** f){
    // Calculates Hessian function (for 1/N correction to config temp).

/*  This routine is only needed in a constant-energy ensemble
    It is assumed that positions are in units where box = 1
    but the result is given in units where sigma = 1 and epsilon = 1
    It is assumed that forces have already been calculated in array f */ 

    double rij_sq;
    double r_cut_box    = r_cut / box;
    double r_cut_box_sq = pow(r_cut_box,2);
    double box_sq       = pow(box,2);

    // Calculate potential at cutoff
    double ff, rf, sr2, sr6, sr8, sr10, v1, v2; 
    int d = 3;
    double hes = 0.0;

    double** f     = allocate2DArray(mm,d);
    double* fij    = new double[3];
    double* rij    = new double[3];
    double* rijb   = new double[3];
    double* rij_p  = new double[3];

    for(int i{0};i<mm;++i){
        for(int j{i+1};j<mm ;++j){
	          for (int k{0};k<3;++k){
		           rij[k] = r[i][k] - r[j][k];            // Separation vector
	          }
            rint1D(3,rij);                             // Periodic boundary conditions in box=1 units
            elementWise1DProduct(3,rij,rij, rij_p);
            rij_sq = elementSum1D(3,rij_p);            // Squared separation
          
            if (rij_sq < r_cut_box_sq){                // Check within cutoff
                rij_sq = rij_sq * box_sq;              // Now in sigma=1 units 
                scalar1DArrayMultip(3,box,rij,rijb);   // Now in sigma=1 units
                for (int k{0};k<3;++k){
		            fij[k] = f[i][k] - f[j][k];         // Difference in forces
                }
                ff   = dotProduct2D(3,fij,fij)
                rf   = dotProduct2D(3,rij,fij)
                sr2  = 1.0 / rij_sq;
                sr6  = pow(sr2,3);
                sr8  = sr6 * sr2;
                sr10 = sr8 * sr2;
                v1   = 24.0 * ( 1.0 - 2.0 * sr6 ) * sr8;
                v2   = 96.0 * ( 7.0 * sr6 - 2.0 ) * sr10;
                hes  = hes + v1 * ff + v2 * pow(rf,2);
            }
        }
    }
  
    return hes
}
