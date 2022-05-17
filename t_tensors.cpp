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

double**** allocate4DArray(int p,int q, int m, int n) {
    double**** array = new double***[p];

    for (int i{0}; i<p; i++) {
        array[i] = new double**[q];
        for (int j{0}; j<q; ++j){
            array[i][j] = new double*[m];
            for (int k{0}; k<m; ++k){
                array[i][j][k] = new double[n];
            }
        }
    }
    return array;
}

double***** allocate5DArray(int p,int q, int m, int n, int s) {
    double***** array = new double**** [p];

    for (int i{0}; i<p; i++) {
        array[i] = new double*** [q];
        for (int j{0}; j<q; ++j){
            array[i][j] = new double** [m];
            for (int k{0}; k<m; ++k){
                array[i][j][k] = new double* [n];
                for (int l{0}; l<n; ++l){
                    array[i][j][k][l] = new double [s];
                }
            }
        }
    }
    return array;
}

// deallocate the memory of a 2D Array
void free2DArray(int m,double **array) {
    for (int i{0}; i<m; ++i)
        delete [] array[i];
        delete [] array;
}

// deacllocate the memory of a 3D Array
void free3DArray(int p,int q,int r,double ***array) {
    for (int i{0}; i<p ; i++){
        for (int j{0}; j<q; j++) {
            delete[] array[i][j];
        }
        delete[] array[i];
    }
    delete[] array;
}

// deacllocate the memory of a 4D Array
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

// randomly initialize a 4D array
void rand4DArray(int p,int q,int m, int n, double ****A) {
    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<m;k++) {
                for (int l=0;l<n;l++) {
                    A[i][j][k][l] = (double) rand()/RAND_MAX;
                }
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

void scalar2DArrayMultip(int m,int n,double p, double** A) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] *= p;
        }
    }
}

void scalar3DArrayMultip(int m,int n, int r, double p, double*** A) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            for (int k=0;k<r;k++) {
                A[i][j][k] *= p;
            }
        }
    }
}

void scalar4DArrayMultip(int p,int q, int m, int n, double b, double**** A) {
    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<m;k++) {
                for (int l=0;l<n;l++) {
                    A[i][j][k][l] *= b;
                }
            }
        }
    }
}

void scalar2DArraySubtract(int m,int n, double** A, double** B) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] -= B[i][j];
        }
    }
}

void scalar2DArrayDivision(int m,int n,double p, double** A) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] /= p;
        }
    }
}

void scalar3DArrayDivision(int m,int n,int r, double p, double*** A) {
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            for (int k=0;k<r;k++) {
                A[i][j][k] /= p;
            }
        }
    }
}

void scalar4DArrayDivision(int p,int q, int m,int n, double b, double**** A) {
    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<m;k++) {
                for (int l=0;l<n;l++) {
                    A[i][j][k][l] /= b;
                }
            }
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

void print4DArray(int p, int q, int m, int n, double**** A){
    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<m;k++) {
                for (int l=0;l<n;l++) {
                    printf("%15.10f  ",A[i][j][k][l]);
                }
                std::cout << '\n';
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

// a function for 4D array outer product
double**** outer4D(int nRows, int nCols, int n3rd, int n4th, double* A, double* B, double* C, double* D){

    double**** E = allocate4DArray(nRows,nCols,n3rd,n4th);

    for (int i{0};i<nRows;++i){
        for (int j{0};j<nCols;++j){
            for (int k{0};k<n3rd;++k){
                for (int l{0};l<n4th;++l){
                    E[i][j][k][l] = A[i] * B[j] * C[k] * D[l];
                }
            }
        }
    }
    return E;
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

    scalar2DArraySubtract(nCols, nCols, t2, A);

    scalarDivision(nCols,nCols, r3, t2);

    return t2;
}

double*** t3_tensor(double *mat3, double r4){
/*  returns third-rank 3x3x3 interaction tensor (note positive sign).

    Supplied arguments should be the unit vector from 2 to 1 and
    the fourth power of the modulus of that vector. */

    int nCols=3;

    double*** t3 = allocate3DArray(nCols,nCols,nCols);

    t3 = outer3D(nCols,nCols,nCols,mat3,mat3,mat3);

    scalar3DArrayMultip(nCols,nCols,nCols,15.,t3);

    for (int i{0};i<3;++i){
            t3[i][i][i] = t3[i][i][i] - 9.0 * mat3[i];
        for (int j{0};j<3;++j){
            if (j == i){
                    continue;
            }
             t3[i][i][j] = t3[i][i][j] - 3.0 * mat3[j];
             t3[i][j][i] = t3[i][j][i] - 3.0 * mat3[j];
             t3[j][i][i] = t3[j][i][i] - 3.0 * mat3[j];
        }
    }

    scalar3DArrayDivision(nCols,nCols,nCols, r4,t3);

    return t3;

}

double**** t4_tensor(double* mat4, double r5){
    /*Returns fourth-rank 3x3x3x3 interaction tensor

    Supplied arguments should be the unit vector from 2 to 1 and
    the fifth power of the modulus of that vector. */

    int nCols=3;

    double**** t4 = allocate4DArray(nCols,nCols,nCols,nCols);
    double**   I  = allocate2DArray(nCols,nCols);

    t4 = outer4D(nCols,nCols,nCols,nCols,mat4,mat4,mat4,mat4);
    scalar4DArrayMultip(nCols,nCols,nCols,nCols,105.,t4);

    identMatrix(nCols,nCols,I);

    for (int i{0};i<nCols;++i){
        for (int j{0};j<nCols;++j){
            for (int k{0};k<nCols;++k){
                for (int l{0};l<nCols;++l){

                    t4[i][j][k][l] = t4[i][j][k][l] - 15.0 * (
                           mat4[i] * mat4[j] * I[k][l] + mat4[i] * mat4[k] * I[j][l]
                         + mat4[i] * mat4[l] * I[j][k] + mat4[j] * mat4[k] * I[i][l]
                         + mat4[j] * mat4[l] * I[i][k] + mat4[k] * mat4[l] * I[i][j])
                         + 3.0 * ( I[i][j] * I[k][l] + I[i][k] * I[j][l] + I[i][l] * I[j][k]);
                }
            }
        }
    }

    scalar4DArrayDivision(nCols,nCols,nCols,nCols,r5,t4);

    return t4;
}


int main()
{
    int nRows=1;
    int nCols=3;
    double r3 = 4.;
    double r4 = 4.;
    double r5 = 4.;

    double*   r = new double[nCols];
    double**  D = allocate2DArray(nCols,nCols);
    double*** T = allocate3DArray(nCols,nCols,nCols);
    double**** F = allocate4DArray(nCols,nCols,nCols,nCols);

    std::cout << "1D Array r:" << '\n';
    rand1DArray(nCols,r);
    print1DArray(nCols,r);
    std::cout << '\n';

    D = t2_tensor(r,r3);
    std::cout << "2D Array D:" << '\n';
    print2DArray(nCols,nCols,D);
    std::cout << '\n';

    T = t3_tensor(r,r4);
    std::cout << "3D Array T:" << '\n';
    print3DArray(nCols,nCols,nCols,T);

    std::cout << "4D Array F:" << '\n';
    F = t4_tensor(r,r5);
    print4DArray(nCols,nCols,nCols,nCols,F);

    std::cout << '\n';
    delete [] r;
    free2DArray(nCols,D);
    free3DArray(nCols,nCols,T);
    free4DArray(nCols,nCols,nCols,F);

}
