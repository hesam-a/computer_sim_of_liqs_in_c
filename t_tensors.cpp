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

// to allocate memory for a matrix [m x n]

double** allocateMatrix(int m,int n) {
    double *space;
    double **array;
    array = new double*[m];
    space = new double[m*n];
    for (int i=0; i<m; i++) {
        array[i] = space;
        space += n;
    }
    return array;
}

//to deallocate the memory:
void freeMatrix(double** array) {
    free(array[0]);
    //free array;
}

// some functions for matrix operations

void randMatrix(int m,int n, double** A) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] = (double) rand()/RAND_MAX;
        }
    }
}

void zeroMatrix(int m,int n, double** A) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] = 0.;
        }
    }
}

void identMatrix(int m,int n, double** A) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] = 1.;
        }
    }
}

void scalarMultip(int m,int n,double p, double** A) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] *= p;
        }
    }
}

void matSubtract(int m,int n, double** A, double** B) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] -= B[i][j];
        }
    }
}

void scalarDivision(int m,int n,double p, double** A) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] /= p;
        }
    }
}

void printMatrix(int m,int n, double** A){
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            printf("%15.10f  ",A[i][j]);
            //std::cout << A[i][j] << " ";
        }
        std::cout << '\n';
    }
}

// outer function

double* outer(int nCols, int nRows, double** A, double** B, double** C){

    double alpha = 1.0, beta = 0.0;

    cblas_dger(CblasRowMajor, nRows, nCols, alpha, A, 1, B, 1, C, nCols);

    return C;
}

// The t2 tensor function

double t2_tensor (r,r3){
/* Returns second-rank 3x3 interaction tensor.
    Supplied arguments should be the unit vector from 2 to 1 and
    the third power of the modulus of that vector.*/

    int nRows=3;
    int nCols=3;

    double** A  = allocateMatrix(nRows,nCols);
    double** t2 = allocateMatrix(nRows,nCols);

    outer(nCols,nRows, mat, mat, t2);
    scalarMultip(nRows,nCols,3,t2);

    identMatrix(nRows,nCols,A);

    matSubtract(nRows, nCols, t2, A);

    scalarDivision(nRows,nCols, r3, t2);

    return t2;
}

int main()
{
    int nRows=3 ;
    int nCols=3;
    double r3 = 4.;

    double** r=allocateMatrix(nRows,nCols);
    double** B=allocateMatrix(nRows,nCols);
    double** C=allocateMatrix(nRows,nCols);

    randMatrix(nRows,nCols,r);

    printMatrix(nRows,nCols,r);
    std::cout << '\n';
    B = t2_tensor(r,r3);
    printMatrix(nRows,nCols,B);
    std::cout << '\n';

    freeMatrix(r);
    freeMatrix(B);
    freeMatrix(C);

    return 0;
}
