#define _USE_MATH_DEFINES
 
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <cassert>
#include <random>

// This is a math module for some functions used in the programs. Although the current version is completely working for the programs,
// some functions might need slight modifications.


// random device class instance, source of 'true' randomness for initializing random seed
std::random_device rd;

// Mersenne twister PRNG, initialized with seed from previous random device instance
std::mt19937 gen(rd());


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
        array[i] = new double[n]{};
    }

    return array;
}

double*** allocate3DArray(int p,int q, int r) {
    double*** array = new double** [p];

    for (int i{0}; i<p; i++) {
        array[i] = new double*[q];
        for (int j{0}; j<q; ++j){
            array[i][j] = new double[r]{};
        }
    }
    return array;
}

double**** allocate4DArray(int p,int q, int m, int n) {
    double**** array = new double*** [p];

    for (int i{0}; i<p; i++) {
        array[i] = new double**[q];
        for (int j{0}; j<q; ++j){
            array[i][j] = new double*[m];
            for (int k{0}; k<m; ++k){
                array[i][j][k] = new double[n]{};
	    }
        }
    }
    return array;
}


double***** allocate5DArray(int p,int q, int m, int n, int s) {
    //double* array = new double [p * q * m * n * s];

    double***** array = new double**** [p];

    for (int i{0}; i<p; i++) {
        array[i] = new double***[q];
        for (int j{0}; j<q; ++j){
            array[i][j] = new double**[m];
            for (int k{0}; k<m; ++k){
                array[i][j][k] = new double*[n];
                for (int l{0}; l<n; ++l){
                    array[i][j][k][l] = new double[s]{};
		}
	    }
        }
    }
    return array;
}

// deallocate the memory of a 2D matrix
void free2DArray(int m,double** array) {
    for (int i{0}; i<m; ++i)
        delete [] array[i];

    delete [] array;
}

// deacllocate the memory of a 3D matrix
void free3DArray(int p,int q, double*** array) {
    for (int i{0}; i<p ; i++){
        for (int j{0}; j<q; j++) {
            delete[] array[i][j];
        }
        delete[] array[i];
    }
    delete[] array;
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


// deacllocate the memory of a 5D matrix
void free5DArray(int p,int q,int m, int n, double***** array) {
    for (int i{0}; i<p ; i++){
        for (int j{0}; j<q; j++) {
            for (int k{0}; k<m; k++) {
                for (int l{0}; l<n; l++) {
		    delete [] array[i][j][k][l];
		}
                delete[] array[i][j][k];
	    }
            delete[] array[i][j];
        }
        delete[] array[i];
    }
    delete[] array;
}

void rand1DArray(int m, double* arr) {

    for (int i=0;i<m;i++) {
        arr[i] = (double) rand()/RAND_MAX;
    }
}

void rand2DArray(int m,int n, double** arr) {

// a function for generating a random 2D array with mean 0.0 and stdev 1.0

    // instance of class std::normal_distribution with specific mean and stddev
    std::normal_distribution<float> d(0, 1);

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
        //    arr[i][j] = (double) rand()/RAND_MAX;
            arr[i][j] = d(gen);
        }
    }
}


void rand3DArray(int p,int q,int r, double*** arr) {


    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<r;k++) {
                arr[i][j][k] = (double) rand()/RAND_MAX;
            }
        }
    }
}

void rand4DArray(int p,int q,int m, int n, double**** arr) {

    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<m;k++) {
                for (int l=0;l<n;l++) {
                    arr[i][j][k][l] = (double) rand()/RAND_MAX;
		}
            }
        }
    }
}

void rand5DArray(int p,int q,int m, int n,int s, double***** arr) {

    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<m;k++) {
                for (int l=0;l<n;l++) {
                    for (int r=0;r<s;r++) {
                        arr[i][j][k][l][r] = (double) rand()/RAND_MAX;
		    }
		}
            }
        }
    }
}

void zeroVec(int m, double* arr) {
    
    for (int i{0};i<m;i++)
	arr[i] = 0.0;
}

void zeroMatrix(int m,int n, double** arr) {

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            arr[i][j] = 0.;
        }
    }
}

void identMatrix(int m,int n, double** arr) {

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
	    if (i==j){
	        arr[i][j] = 1.;
	    }
        }
    }
}

void matMultip(int m, int n, double** A, double** B, double** C){

    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) {
            for (int k=0;k<3;k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}
    

void sum1DArrays(int m, double* A, double* B, double* sum) {

    for (int i=0;i<m;i++) {
            sum[i] = A[i] + B[i];
    }
}

void sum2DArrays(int m, int n,double** A, double* sum,int axis) {

    assert(axis==0 or axis==1);

    if (axis == 0){ // column-wise sum
	for (int i=0;i<m;i++) {
	    sum[i] = 0.0;
	    for (int j=0;j<n;j++) {
                sum[i] += A[i][j];
	    }
	}
    }
    else if (axis == 1){ // row-wise sum
	for (int i=0;i<m;i++) {
	    sum[i] = 0.0;
	    for (int j=0;j<n;j++) {
                sum[j] += A[j][i];
	    }
	}
    }
}

double max1DArray(int m, double* arr){
 
    double max=0.;	
    for (int i=0;i<m;i++) {
	if (arr[i] > max){
	    max = arr[i];
	}
    }
    return max;
}

void scalar1DArrayMultip(int m,double p, double* A, double* arr) {

    for (int i=0;i<m;i++) {
            arr[i] = A[i] * p;
    }
}

void scalar2DArrayMultip(int m,int n,double p, double** A, double** product) {

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            product[i][j] = A[i][j] * p;
        }
    }
}

void scalar3DArrayMultip(int m,int n, int r, double p, double*** A, double*** product) {

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            for (int k=0;k<r;k++) {
                product[i][j][k] = A[i][j][k] * p;
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

void scalar5DArrayMultip(int p,int q, int m, int n, int s, double b, double***** A) {
    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<m;k++) {
                for (int l=0;l<n;l++) {
                    for (int r=0;r<s;r++) {
                        A[i][j][k][l][r] *= b;
		    }
		}
	    }
        }
    }
}

void rint1D(int m,double* A){

    for (int i=0;i<m;i++)
        A[i] -= round(A[i]);
}

void rint2D(int m,int n,double** A){

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] -= round(A[i][j]);
	}
    }
}

void sqrt1DArray(int m, double*A, double* sq ){

    for (int i{0}; i<m;++i)
	sq[i] = sqrt(A[i]);
}


void subtract1DArrays(int m, double* A, double* B, double* C) {

    for (int i=0;i<m;i++)
        C[i] = A[i] - B[i];
}

void subtract2DArrays(int m,int n, double** A, double** B, double** C) {

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

void scalar1DArraySubtract(int m, double p, double* A) {

    for (int i=0;i<m;i++)
        A[i] -= p;
}

void scalar1DArrayDivision(int m,double p, double* A, double* product) {

    for (int i=0;i<m;i++)
            product[i] = A[i] / p;
}

void scalar2DArraySubtract(int m,int n, double p, double** A) {

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            A[i][j] -= p;
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

void scalar5DArrayDivision(int p,int q, int m,int n, int s, double b, double***** A) {

    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<m;k++) {
                for (int l=0;l<n;l++) {
                    for (int r=0;r<s;r++) {
                        A[i][j][k][l][r] /= b;
		    }
		}
            }
        }
    }
}

void print1DArray(int m, double* A){
    for (int i=0;i<m;i++) {
        printf("%14.10f  ",A[i]);
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

void print5DArray(int p, int q, int m, int n, int s, double***** A){
    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<m;k++) {
                for (int l=0;l<n;l++) {
                    for (int r=0;r<s;r++) {
                    printf("%15.10f  ",A[i][j][k][l][r]);
                    //    printf("%15.10f  ",*(A + i*q*m*n*s + j*m*n*s+ k*n*s + l*s +r));
		    }
		    std::cout << '\n';
		}
                std::cout << '\n';
            }
            std::cout << '\n';
        }
        std::cout << '\n';
    }
}

// outer function 
void outer2D(int nCols, double* A, double* B, double** C){

    for (int i{0};i<nCols;++i){
	for (int j{0};j<nCols;++j){
		C[i][j] = A[i] * B[j] ;
	}
    }
}

void outer3D(int nRows, int nCols, int n3rd, double* A, double* B, double* C, double*** D){
	
    for (int i{0};i<nRows;++i){
	for (int j{0};j<nCols;++j){
            for (int k{0};k<n3rd;++k){
		D[i][j][k] = A[i] * B[j] * C[k];
	    }
	}
    }
}


void outer4D(int nRows, int nCols, int n3rd, int n4th, double* A, double* B, double* C, double* D, double**** E){
	
    for (int i{0};i<nRows;++i){
	for (int j{0};j<nCols;++j){
            for (int k{0};k<n3rd;++k){
                for (int l{0};l<n4th;++l){
  		    E[i][j][k][l] = A[i] * B[j] * C[k] * D[l];
		}
	    }
	}
    }
}


void outer5D(int nRows, int nCols, int n3rd, int n4th, int n5th, double* A, double* B, double* C, double* D, double* E, double***** F){
	
    for (int i{0};i<nRows;++i){
	for (int j{0};j<nCols;++j){
            for (int k{0};k<n3rd;++k){
                for (int l{0};l<n4th;++l){
                    for (int r{0};r<n5th;++r){
  		        F[i][j][k][l][r] = A[i] * B[j] * C[k] * D[l] * E[r];
		    }
		}
	    }
	}
    }
}

void t2_tensor (double* mat,double r3, double** t2){
/* Returns second-rank 3x3 interaction tensor.
    Supplied arguments should be the unit vector from 2 to 1 and
    the third power of the modulus of that vector.*/

    int nCols  = 3;
    double** A = allocate2DArray(nCols,nCols);
 
    outer2D(nCols,mat,mat,t2);

    scalar2DArrayMultip(nCols,nCols,3.,t2,t2);

    identMatrix(nCols,nCols,A);
   
    subtract2DArrays(nCols,nCols, t2, A, t2);
    
    scalar2DArrayDivision(nCols,nCols, r3, t2);
     
    free2DArray(nCols,A);
}

void t3_tensor(double *mat3, double r4, double*** t3){
/*  returns third-rank 3x3x3 interaction tensor (note positive sign).
    Supplied arguments should be the unit vector from 2 to 1 and
    the fourth power of the modulus of that vector. */
	
    int nCols = 3;

    outer3D(nCols,nCols,nCols,mat3,mat3,mat3,t3);

    scalar3DArrayMultip(nCols,nCols,nCols,15.,t3,t3);

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

}

void t4_tensor(double* mat4, double r5, double**** t4){
    /*Returns fourth-rank 3x3x3x3 interaction tensor
    Supplied arguments should be the unit vector from 2 to 1 and
    the fifth power of the modulus of that vector. */	

    int nCols  = 3;
    double** A = allocate2DArray(nCols,nCols);

    outer4D(nCols,nCols,nCols,nCols,mat4,mat4,mat4,mat4,t4);
    scalar4DArrayMultip(nCols,nCols,nCols,nCols,105.,t4);

    identMatrix(nCols,nCols,A);

    for (int i{0};i<nCols;++i){
        for (int j{0};j<nCols;++j){
            for (int k{0};k<nCols;++k){
                for (int l{0};l<nCols;++l){

                    t4[i][j][k][l] = t4[i][j][k][l] - 15.0 * (
                           mat4[i] * mat4[j] * A[k][l] + mat4[i] * mat4[k] * A[j][l]
                         + mat4[i] * mat4[l] * A[j][k] + mat4[j] * mat4[k] * A[i][l]
                         + mat4[j] * mat4[l] * A[i][k] + mat4[k] * mat4[l] * A[i][j])
                         + 3.0 * ( A[i][j] * A[k][l] + A[i][k] * A[j][l] + A[i][l] * A[j][k]);
                }
            }
        }
    }

    scalar4DArrayDivision(nCols,nCols,nCols,nCols,r5,t4);

    free2DArray(nCols,A);
}


void t5_tensor(double* mat5, double r6, double***** t5){
/*  Returns fifth-rank 3x3x3x3x3 interaction tensor
    Supplied arguments should be the unit vector from 2 to 1 and
    the fifth power of the modulus of that vector. */	

    int nCols  = 3;
    double** A = allocate2DArray(nCols,nCols);

    outer5D(nCols,nCols,nCols,nCols,nCols,mat5,mat5,mat5,mat5,mat5,t5);
    scalar5DArrayMultip(nCols,nCols,nCols,nCols,nCols,945.,t5);

    identMatrix(nCols,nCols,A);
    for (int i{0};i<nCols;++i){
        for (int j{0};j<nCols;++j){
            for (int k{0};k<nCols;++k){
                for (int l{0};l<nCols;++l){
                    for (int m{0};m<nCols;++m){

			t5[i][j][k][l][m] = t5[i][j][k][l][m] - 105.0 * (
                          mat5[i] * mat5[j] * mat5[k] * A[l][m] + mat5[i] * mat5[j] * mat5[l] * A[k][m]  
                        + mat5[i] * mat5[j] * mat5[m] * A[k][l] + mat5[i] * mat5[k] * mat5[l] * A[j][m]  
                        + mat5[i] * mat5[k] * mat5[m] * A[j][l] + mat5[i] * mat5[l] * mat5[m] * A[j][k]  
                        + mat5[j] * mat5[k] * mat5[l] * A[i][m] + mat5[j] * mat5[k] * mat5[m] * A[i][l]  
                        + mat5[j] * mat5[l] * mat5[m] * A[i][k] + mat5[k] * mat5[l] * mat5[m] * A[i][j] )
                        + 15.0 * ( 
                          mat5[i] * ( A[j][k] * A[l][m] + A[j][l] * A[k][m] + A[j][m] * A[k][l] )
                        + mat5[j] * ( A[i][k] * A[l][m] + A[i][l] * A[k][m] + A[i][m] * A[k][l] )
                        + mat5[k] * ( A[i][j] * A[l][m] + A[i][l] * A[j][m] + A[i][m] * A[j][l] )
                        + mat5[l] * ( A[i][j] * A[k][m] + A[i][k] * A[j][m] + A[i][m] * A[j][k] )
                        + mat5[m] * ( A[i][j] * A[k][l] + A[i][k] * A[j][l] + A[i][l] * A[j][k] ) );
		    }
		}
	    }
	}
    }
    scalar5DArrayDivision(nCols,nCols,nCols,nCols,nCols,r6,t5);

    free2DArray(nCols,A);
}

void skew(double** vec, double* b){
/*  Returns contraction of supplied 3x3 matrix with Levi-Civita tensor.*/

    b[0] = vec[1][2] - vec[2][1];
    b[1] = vec[2][0] - vec[0][2];
    b[2] = vec[0][1] - vec[1][0];
}


double dotProduct1D(int n,double* a){

    double product{1.0};
    
    for (int i{0}; i<n; i++)
	product *= a[i];

    return product;
}

double dotProduct2D(int n,double* a, double* b){

    double product{0.0};
    
    for (int i{0}; i<n; i++)
	product = product + (a[i] * b[i]);

    return product;
}


void crossProduct(double* a, double* b, double* product){

    product[0] = a[1]*b[2] - a[2]*b[1];
    product[1] = a[2]*b[0] - a[0]*b[2];
    product[2] = a[0]*b[1] - a[1]*b[0];

}

double elementSum1D(int m, double* a){
    
    double product{0.};

    for (int i{0};i<m;++i)
        product += a[i];

    return product;
}

double elementSum2D(int m,int n, double** a){
    
    double product{0.};

    for (int i{0};i<m;++i){
        for (int j{0};j<n;++j){
         product += a[i][j];
	}
    }
    return product;
}

void elementWise1DProduct(int m, double* a, double* b, double* product){

    for (int i{0};i<m;++i)
         product[i] = a[i] * b[i];
}

void elementWise2DProduct(int m, int n, double** a, double** b, double** product){
    
    for (int i{0};i<m;++i){
        for (int j{0};j<n;++j){
            product[i][j] = a[i][j] * b[i][j];
        }
    }
}

double contract_i_i (int n,double* a, double* b){
/*  Returns a zero-rank contraction of a first-rank tensor
    with a first-rank tensor. */

    double c;
    c = dotProduct2D(n,a,b);
    
    return c;
}


void contract_ij_j (double** a, double* b, double* product){
/*  Returns a first-rank contraction offirst-rank tensor
    with a first-rank tensor. */

    for(int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    product[i] = product[i] + a[i][j] *  b[j];
	}
    }   
}	


double contract_ij_ij (double** a, double** b){
/*  c ! Returns a zero-rank contraction of a second-rank tensor
    with another second-rank tensor. */

    double product{0};
    double** dot = allocate2DArray(3,3);
    elementWise2DProduct(3,3,a,b,dot);
//    print2DArray(3,3,dot);
    for (int i{0};i<3;++i){
        for (int j{0};j<3;++j){
            product = product + dot[i][j];
        }
    }
    free2DArray(3,dot);
    return product;

}

void contract_ik_jk (double** a, double** b, double** c){
/*  Returns a second-rank contraction of a second-rank tensor
    with another second-rank tensor. */

    for (int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    for (int k{0};k<3;++k){
		    c[i][j] = c[i][j] + a[i][k] * b[j][k];
	    }
	}
    }
}


void contract_ijk_k (double*** a, double* b, double** c){
/*  Returns a second-rank contraction of a third-rank tensor
    and a first-rank tensor. */

    for (int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    for (int k{0};k<3;++k){
                c[i][j] =  c[i][j] + a[i][j][k] * b[k];
	    }
	}
    }
}

void contract_ijk_jk (double*** a, double** b, double* c){
/*  Returns a second-rank contraction of a fourth-rank tensor
    and a second-rank tensor */

    for (int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    for (int k{0};k<3;++k){
                c[i] =  c[i] + a[i][j][k] * b[j][k];
	    }
	}
    }
}

void contract_ijkl_kl(double**** a, double** b, double** c){
/*  Returns a second-rank contraction of a fourth-rank tensor
    and a second-rank tensor. */

    for (int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    for (int k{0};k<3;++k){
	        for (int l{0};l<3;++l){
                    c[i][j] =  c[i][j] + a[i][j][k][l] * b[k][l];
		}
	    }
	}
    }
}


void contract_ijklm_lm(double***** a, double** b, double*** c){
/*  Returns a third-rank contraction of a fifth-rank tensor
    and a second-rank tensor. */

    for (int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    for (int k{0};k<3;++k){
	        for (int l{0};l<3;++l){
	            for (int m{0};m<3;++m){
                        c[i][j][k] = c[i][j][k] + a[i][j][k][l][m] * b[l][m];
		    }
		}
	    }
	}
    }
}

void random_vector(double* rand_vec){
//  Returns a random unit vector as a numpy array of 3 elements. 

    double* zeta = new double[2];
    rand1DArray(2,zeta);                      // Two uniformly sampled random numbers in range (0,1)
    double c         = 2.0*zeta[0] - 1.0;     // Random cos(theta) uniformly sampled in range (-1,+1)
    double s,phi;                             // desclare sin, angle phi

    if (c >= 1.0){                            // Guard against very small chance of roundoff error
        s = 0.0;                              // Set sin(theta) to zero
    }
    else{
        s = sqrt(1.0-pow(c,2));               // Calculate sin(theta) from cos(theta), always positive
    }
    phi = zeta[1] * 2.0*M_PI;                 //  Random angle uniformly sampled in range (0,2*pi)

    rand_vec[0] = s*cos(phi);
    rand_vec[1] = s*sin(phi);
    rand_vec[2] = c;

    delete [] zeta;
 }


void remove2DArray(int m, int p, double** A, double** B){

    if (p <m-1){
        int ii = 0;
        for (int i{0};i<m;i++){
            if (i == p){
                i++;
            }
            for (int j{0};j<3;++j){
                B[ii][j] = A[i][j];
           }
           ii++;
        }
    }
    else if (p == m-1){
        for (int i{0};i<m-1;i++){
            for (int j{0};j<3;++j){
                B[i][j] = A[i][j];
	    }
	}
    }
}

void remove3DArray(int m, int n, int p, double*** A, double*** B){

    if (p <m-1){
        int ii = 0;
        for (int i{0};i<m;i++){
            if (i == p){
                i++;
            }
            for (int j{0};j<n;++j){
                for (int k{0};k<n;++k){
                B[ii][j][k] = A[i][j][k];
		}
           }
           ii++;
        }
    }
    else if (p == m-1){
        for (int i{0};i<m-1;i++){
            for (int j{0};j<3;++j){
                for (int k{0};k<3;++k){
                    B[i][j][k] = A[i][j][k];
		}
	    }
	}
    }
}

void random_translate_vector (double dr_max, double* old,double* ri ){
/*  Returns a vector translated by a random amount.
 
     A randomly chosen vector is added to the old one. */

    double* zeta = new double[3];
    rand1DArray(3,zeta);             // Three uniformly sampled random numbers in range (0,1)
    scalar1DArrayMultip(3,2.0,zeta,zeta);
    scalar1DArraySubtract(3,1.0,zeta);             // Now in range (-1,+1)
    scalar1DArrayMultip(3,dr_max,zeta,zeta);
    sum1DArrays( 3, old, zeta,ri);             // Move to new position

    delete [] zeta;
}

void quatmul (double* a, double* b, double* c){
//  Returns quaternion product of two supplied quaternions.

    c[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
    c[1] = a[1]*b[0] + a[0]*b[1] - a[3]*b[2] + a[2]*b[3];
    c[2] = a[2]*b[0] + a[3]*b[1] + a[0]*b[2] - a[1]*b[3];
    c[3] = a[3]*b[0] - a[2]*b[1] + a[1]*b[2] + a[0]*b[3];
}

void rotate_quaternion (double angle, double* axis,double* old, double* e){
/*  Returns a quaternion rotated by angle about axis relative to old quaternion.

    Note that the axis vector should be normalized and we test for this
    In general, the old quaternion need not be normalized, and the same goes for the result
    although in our applications we only ever use unit quaternions (to represent orientations) */
    double* rot  = new double[3];
    double* rot2 = new double[4];

    // Standard formula for rotation quaternion, using half angles
    scalar1DArrayMultip(3,sin(0.5*angle),axis,rot);
    rot2[0] = cos(0.5*angle);
    for(int i{0};i<3;++i)
	rot2[i+1] = rot[i];
    quatmul (rot2, old, e); // Apply rotation to old quaternion
    delete [] rot;
    delete [] rot2;
}

void random_quaternion(double* randq){
//  Returns a random unit quaternion as a numpy array of 4 elements.

    double* zeta = new double[2];
    double* zeta2= new double[2];
    double* beta = new double[2];
    double* beta2= new double[2];
    double norm1,norm2,f;

    while (true){                            // Loop until within unit disk
        rand1DArray(2,zeta);             
        scalar1DArrayMultip(2,2.0,zeta,zeta);
        scalar1DArraySubtract(2,1.0,zeta);   // Two uniform random numbers between -1 and 1
	elementWise1DProduct(2,zeta, zeta,zeta2);
	norm1 = elementSum1D(2, zeta2);      // Squared magnitude
        if (norm1 < 1.0)                     // Test for within unit disk
            break;
    }

    while (true){                            // Loop until within unit disk
        rand1DArray(2,beta);             
        scalar1DArrayMultip(2,2.0,zeta,zeta);
        scalar1DArraySubtract(2,1.0,zeta);   // Two uniform random numbers between -1 and 1
	elementWise1DProduct(2,zeta, zeta,zeta2);
	norm2 = elementSum1D(2, zeta2);      // Squared magnitude
        if (norm2 < 1.0)                     // Test for within unit disk
            break;
    }
    f = sqrt((1.0 - norm1)/norm2);
    randq[0] = zeta[0];
    randq[1] = zeta[1];
    randq[2] = beta[0]*f;
    randq[3] = beta[1]*f;
    delete [] zeta;
    delete [] zeta2;
    delete [] beta;
    delete [] beta2;
}

void random_rotate_quaternion (double angle_max, double* old, double* e){
/*  Returns a unit quaternion rotated by a maximum angle (in radians) relative to the old quaternion.

    Note that the reference quaternion should be normalized and we test for this */
    double randd = (double) rand()/RAND_MAX;
    double* axis = new double[3];

    random_vector(axis);                               // Choose random unit vector
    double angle = ( 2.0* randd - 1.0 ) * angle_max;          // Uniform random angle in desired range
    rotate_quaternion (angle, axis, old, e);           // General rotation function
    delete [] axis;
}

bool metropolis (double delta ){
//  Conduct Metropolis test, with safeguards.

    double exponent_guard = 75.0;

    if (delta > exponent_guard){                // Too high, reject without evaluating
        return false;
    }
    else if (delta < 0.0){                      // Downhill, accept without evaluating
        return true;
    }
    else{
        double zeta = (double) rand()/RAND_MAX; // Uniform random number in range (0,1)
        return exp(-delta) > zeta;         // Metropolis test
    }
}

void update2DArrayZVT(int m, int n, double* ri, double* rj, double** A){

    for (int i{0};i<m;++i){
        for(int j{0};j<n;++j){
            if (ri[j] == A[i][j]){
                A[i][j] = rj[j];
            }
        }
    }
}

void update2DArray(int m, int atom, double* ri, double** A){

    for (int i{0};i<m;++i)
	A[atom][i] = ri[i];
}

//void update3DArray(int m,int n, int p, double** ri, double** rj, double*** A){
void update3DArray(int m,int n, int atom, double** ri, double*** A){

    for (int i{0};i<m;++i){
	for (int j{0};j<n;++j){
	    A[atom][i][j] = ri[i][j];
	}
    }
}

double** createArray(int m,int n, double* ri, double** A){

    size_t newSize = m+1;
    double** newArr = allocate2DArray(newSize,n);

    for (int i{0};i<m;++i){
        for(int j{0};j<n;++j){
            newArr[i][j] = A[i][j];
        }
    }

    for(int k{0};k<n;++k)
        newArr[m][k] = ri[k];

    m = newSize;
    free2DArray(m-1,A);
    A = newArr;

    return A;
}

double** annihilateArray(int m,int n, double* ri, double** A){

    int newSize = m - 1;
    double** newArr = new double*[newSize];

    for (int i{0}; i<m-1; i++)
        newArr[i] = new double[n];

    bool last = false;
    int count =0;
    for (int i{0};i<m-1;++i){
        for(int j{0};j<n;++j){
            if (A[i][j] == ri[j]){
                count = i;
                last = true;
                break;
            }
            else{
                newArr[i][j] = A[i][j];
            }
        }
    }

    if (last == false){
    }
    else{
        for (int i{count+1};i<m;++i){
            for(int j{0};j<n;++j){
                newArr[i-1][j] = A[i][j];
            }
        }
    }

    m = newSize;
    free2DArray(m+1,A);
    A = newArr;

    return A;

}

void q_to_a (double* q, double** a){
/*  Returns a 3x3 rotation matrix calculated from supplied quaternion.

    The rows of the rotation matrix correspond to unit vectors of the molecule in the space-fixed frame
    The third row  a(3,:) is "the" axis of the molecule, for uniaxial molecules
    Use a to convert space-fixed to body-fixed axes thus: db = np.dot(a,ds)
    Use transpose of a to convert body-fixed to space-fixed axes thus: ds = np.dot(db,a)

    The supplied quaternion should be normalized and we check for this
    assert np.isclose(np.sum(q**2),1.0), 'quaternion normalization error {} {} {} {}'.format(*q) */

    // Write out row by row, for clarity
    a[0][0] =  pow(q[0],2)+ pow(q[1],2)-pow(q[2],2)-pow(q[3],2);           
    a[0][1] =  2*(q[1]*q[2]+q[0]*q[3]);          
    a[0][2] =  2*(q[1]*q[3]-q[0]*q[2]);             
    a[1][0] =  2*(q[1]*q[2]-q[0]*q[3]);
    a[1][1] =  pow(q[0],2)-pow(q[1],2)+pow(q[2],2)-pow(q[3],2);
    a[1][2] =  2*(q[2]*q[3]+q[0]*q[1]); 
    a[2][0] =  2*(q[1]*q[3]+q[0]*q[2]);
    a[2][1] =  2*(q[2]*q[3]-q[0]*q[1]);
    a[2][2] =  pow(q[0],2)-pow(q[1],2)-pow(q[2],2)+pow(q[3],2); 
}
