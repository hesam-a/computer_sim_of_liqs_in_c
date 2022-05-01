// This program is a C++ version of the "t_tensor.f90" program from the book "Computer simulations of liquids"
// by Allen and Tildesley
// It is used for calculation of electrostatic interactions between two linear molecules including dipole-dipole,
// dipole-quadrupole and quadrupole-quadrupole, using T tensors and Euler angles.

// While the book itself does not provide detailed discussion about the equations, "The theory of intermolecular 
// forces" by Anthony Stones is a good source for derivations and further discussion.

// some expalanations from Allen-Tildesley:

/* The dipole moment of molecule 1 is aligned along the axial vector e1
   The quadrupole tensor, quad1, is diagonal and traceless with
   quad1_xx = -0.5*quad1_mag, quad1_yy = -0.5*quad1_mag, and quad1_zz = quad1_mag in the molecule-fixed system.
   Similarly for molecule 2
   The vector r12 = r1-r2 points from 2 to 1.
   
   Forces are calculated by differentiating the T-tensor, giving the next higher rank T-tensor
   Torques are calculated from the angular dependence of dipole, quadrupole etc.
   potential V = mu_i g_i => torque tau_i = -epsilon_ijk mu_j g_k (=-cross_product)
   potential V = Q_ij G_ij => torque tau_l = -2 epsilon_lij Q_ik G_jk
   where ijkl are Cartesian indices and epsilon is the Levi-Civita symbol
   It is just necessary to identify the constants g_i, G_ij, in terms of the T tensor and the
   multipole on the other molecule. */
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <mkl.h>

#define nCols ???
#define nRows ???
#define k ???



double alpha = 1.0, beta = 0.0;
double *C;
C = (double*) malloc(nRows * nCols * sizeof(double));
for (i = 0; i < nRows*nCols; i++) C[i] = 0.0;

cblas_dger(CblasRowMajor, nRows, nCols, alpha, A, 1, B, 1, C, nCols);



// The tensor functions are defined here in order

int t2_tensor (r,r3){
/* Returns second-rank 3x3 interaction tensor.
    Supplied arguments should be the unit vector from 2 to 1 and
    the third power of the modulus of that vector.*/
   
}

int main(int agrc, char* argv[]){

    free(C);

    return 0;
}
