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

// to allocate memory for a 2D Array
double** allocate2DArray(int m,int n) {
//    double** array = new double*[m];
//    double*  space = new double[m*n];
//
//    for (int i{0}; i<m; i++) {
//        array[i] = space;
//        space += n;
//    }
    double** array = new double*[m];

    for (int i{0}; i<m; i++) {
        array[i] = new double[n];
    }
    return array;
}

// to allocate memory for a 3D Array
double*** allocate3DArray(int p,int q, int r) {
    double*** array = new double**[p];

    for (int i{0}; i<p; i++) {
        array[i] = new double*[q];
        for (int j{0}; j<q; ++j){
            array[i][j] = new double[r];
        }
    }
    return array;
}

// deacllocate the memory of a 4D matrix
void free4DArray(int p,int q,int m,double**** array) {
    for (int i{0}; i<p ; i++){
        for (int j{0}; j<q; j++) {
            for (int k{0}; k<m; k++) {
                delete[] array[i][j][k];
            }
            delete[] array[i][j];
        }
        delete[] array[i];
    }
    delete[] array;
}

// deallocate the memory of a 2D matrix
void free2DArray(int m,double **array) {
    for (int i{0}; i<m; ++i)
        delete [] array[i];
        delete [] array;
}

// deacllocate the memory of a 3D matrix
void free3DArray(int p,int q,int r,double ***array) {
    for (int i{0}; i<p ; i++){
        for (int j{0}; j<q; j++) {
            delete[] array[i][j];
        }
        delete[] array[i];
    }
    delete[] array;
}

// randomly initialize a 1D array
void rand1DArray(int m, double *A) {
    for (int i=0;i<m;i++) {
        A[i] = (double) rand()/RAND_MAX;
    }
}


// randomly initialize a 2D array
void rand2DArray(int m,int n, double **A) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] = (double) rand()/RAND_MAX;
        }
    }
}

// randomly initialize a 3D array
void rand3DArray(int p,int q,int r, double ***A) {
    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<r;k++) {
                A[i][j][k] = (double) rand()/RAND_MAX;
            }
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
            if (i==j){
                A[i][j] = 1.;
            }
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

void print1DArray(int m, double* A){
    for (int i=0;i<m;i++) {
        printf("%15.10f  ",A[i]);
    }
    std::cout << '\n';
}

void print2DArray(int m,int n, double** A){
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            printf("%15.10f  ",A[i][j]);
            //std::cout << A[i][j] << " ";
        }
        std::cout << '\n';
    }
}

void print3DArray(int p, int q, int r, double*** A){
    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<r;k++) {
                printf("%15.10f  ",A[i][j][k]);
            }
            std::cout << '\n';
        }
        std::cout << '\n';
    }
}

// a function for 2D array outer product 
double** outer2D(int nCols, double* A, double* B){

    double** C = allocate2DArray(nCols,nCols);

    for (int i{0};i<nCols;++i){
        for (int j{0};j<nCols;++j){
                C[i][j] = A[i] * B[j] ;
        }
    }
    return C;
}
// a function for 3D array outer product
double*** outer3D(int nRows, int nCols, int n3rd, double* A, double* B, double* C){

    double*** D = allocate3DArray(nRows,nCols,n3rd);

    for (int i{0};i<nRows;++i){
        for (int j{0};j<nCols;++j){
            for (int k{0};k<n3rd;++k){
                D[i][j][k] = A[i] * B[j] * C[k];
            }
        }
    }
    return D;
}

// The t2 tensor function
double** t2_tensor (double* mat,double r3){
/* Returns second-rank 3x3 interaction tensor.
    Supplied arguments should be the unit vector from 2 to 1 and
    the third power of the modulus of that vector.*/

    int nCols=3;

    double** A  = allocate2DArray(nCols,nCols);
    double** t2 = allocate2DArray(nCols,nCols);
    
    t2 = outer2D(nCols, mat, mat);
    scalarMultip(nCols,nCols,3.,t2);

    identMatrix(nCols,nCols,A);

    matSubtract(nCols, nCols, t2, A);

    scalarDivision(nCols,nCols, r3, t2);

    return t2;
}


int main()
{
    int nRows=1;
    int nCols=3;
    int d3rd =3;
    double r3 = 4.;
//    double r4 = 4.;

    double* r = new double[nCols];
    double**  D=allocate2DArray(nCols,nCols);
//    double*** C=allocate3DArray(nCols,nCols,d3rd);

    std::cout << "1D Array r:" << '\n';
    rand1DArray(nCols,r);
    print1DArray(nCols,r);
    std::cout << '\n';

//    std::cout << "3D Array C:" << '\n';
//    C = outer3D(nCols,nCols,d3rd,r,r,r);
//    print3DArray(nCols,nCols,d3rd,C);

    D = t2_tensor(r,r3);
    std::cout << "2D Array D:" << '\n';
    print2DArray(nCols,nCols,D);
    std::cout << '\n';

    delete [] r;
    free2DArray(nCols,D);
//    free3DArray(nCols,nCols,d3rd,C);

    return 0;
}

