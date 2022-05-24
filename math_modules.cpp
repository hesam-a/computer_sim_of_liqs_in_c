#define _USE_MATH_DEFINES
 
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

// This is a math module for some functions used in the programs. Although the current version is completely working for the programs,
// some functions might need 

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

double*** allocate3DArray(int p,int q, int r) {
    double*** array = new double** [p];

    for (int i{0}; i<p; i++) {
        array[i] = new double*[q];
        for (int j{0}; j<q; ++j){
            array[i][j] = new double[r];
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
                array[i][j][k] = new double[n];
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
                    array[i][j][k][l] = new double[s];
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

double* rand1DArray(int m) {

    double* arr = new double[m];

    for (int i=0;i<m;i++) {
        arr[i] = (double) rand()/RAND_MAX;
    }
    return arr;
}

double** rand2DArray(int m,int n) {

    double** arr = allocate2DArray(m,n);

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            arr[i][j] = (double) rand()/RAND_MAX;
        }
    }
    return arr;
}

double*** rand3DArray(int p,int q,int r) {

    double*** arr = allocate3DArray(p,q,r);

    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<r;k++) {
                arr[i][j][k] = (double) rand()/RAND_MAX;
            }
        }
    }
    return arr;
}

double**** rand4DArray(int p,int q,int m, int n) {

    double**** arr = allocate4DArray(p,q,m,n);

    for (int i=0;i<p;i++) {
        for (int j=0;j<q;j++) {
            for (int k=0;k<m;k++) {
                for (int l=0;l<n;l++) {
                    arr[i][j][k][l] = (double) rand()/RAND_MAX;
		}
            }
        }
    }
    return arr;
}

double***** rand5DArray(int p,int q,int m, int n,int s) {

    double***** arr = allocate5DArray(p,q,m,n,s);

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
    return arr;
}

double** zeroMatrix(int m,int n) {

    double** arr = allocate2DArray(m,n);

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            arr[i][j] = 0.;
        }
    }
    return arr;
}

double** identMatrix(int m,int n) {

    double** arr = allocate2DArray(m,n);

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
	    if (i==j){
	        arr[i][j] = 1.;
	    }
        }
    }
    return arr;
}

double* sum1DArrays(int m, double* A, double* B) {

    double* sum = new double [m];

    for (int i=0;i<m;i++) {
            sum[i] = A[i] + B[i];
    }
    return sum;
}

double* scalar1DArrayMultip(int m,double p, double* A) {

    double* arr = new double [3];

    for (int i=0;i<m;i++) {
            arr[i] = A[i] * p;
    }
    return arr;
}

double** scalar2DArrayMultip(int m,int n,double p, double** A) {

    double** product = allocate2DArray(m,n);

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            product[i][j] = A[i][j] * p;
        }
    }
    return product;
}

double*** scalar3DArrayMultip(int m,int n, int r, double p, double*** A) {

    double*** product = allocate3DArray(m,n,r);

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            for (int k=0;k<r;k++) {
                product[i][j][k] = A[i][j][k] * p;
	    }
        }
    }
    return product;
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

double* subtract1DArrays(int m, double* A, double* B) {

    double* sub = new double [m];

    for (int i=0;i<m;i++) {
            sub[i] = A[i] - B[i];
    }
    return sub;
}

double** subtract2DArrays(int m,int n, double** A, double** B) {

    double** sub = allocate2DArray(m,n);

    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            sub[i][j] = A[i][j] - B[i][j];
        }
    }
    return sub;
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
        printf("%10.10f  ",A[i]);
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
double** outer2D(int nCols, double* A, double* B){

    double** C = allocate2DArray(nCols,nCols);

    for (int i{0};i<nCols;++i){
	for (int j{0};j<nCols;++j){
		C[i][j] = A[i] * B[j] ;
	}
    }
    return C;
}

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


double***** outer5D(int nRows, int nCols, int n3rd, int n4th, int n5th, double* A, double* B, double* C, double* D, double* E){

    double***** F = allocate5DArray(nRows,nCols,n3rd,n4th,n5th);
	
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

    return F;
}

double** t2_tensor (double* mat,double r3){
/* Returns second-rank 3x3 interaction tensor.
    Supplied arguments should be the unit vector from 2 to 1 and
    the third power of the modulus of that vector.*/

    int nCols=3;
   
    double** A  = allocate2DArray(nCols,nCols);
    double** t2 = allocate2DArray(nCols,nCols);
 
    t2 = scalar2DArrayMultip(nCols,nCols,3.,outer2D(nCols, mat, mat));

    A  = identMatrix(nCols,nCols);
   
    t2 = subtract2DArrays(nCols, nCols, t2, A);
    
    scalar2DArrayDivision(nCols,nCols, r3, t2);
     
    return t2;
}

double*** t3_tensor(double *mat3, double r4){
/*  returns third-rank 3x3x3 interaction tensor (note positive sign).
    Supplied arguments should be the unit vector from 2 to 1 and
    the fourth power of the modulus of that vector. */
	
    int nCols=3;

    double*** t3 = allocate3DArray(nCols,nCols,nCols);

    t3 = outer3D(nCols,nCols,nCols,mat3,mat3,mat3);

    t3 = scalar3DArrayMultip(nCols,nCols,nCols,15.,t3);

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
    double**   A  = allocate2DArray(nCols,nCols);

    t4 = outer4D(nCols,nCols,nCols,nCols,mat4,mat4,mat4,mat4);
    scalar4DArrayMultip(nCols,nCols,nCols,nCols,105.,t4);

    A = identMatrix(nCols,nCols);

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

    return t4;
}


double***** t5_tensor(double* mat5, double r6){
/*    Returns fifth-rank 3x3x3x3x3 interaction tensor
    Supplied arguments should be the unit vector from 2 to 1 and
    the fifth power of the modulus of that vector. */	

    int nCols=3;

    double***** t5 = allocate5DArray(nCols,nCols,nCols,nCols,nCols);
    double**    A  = allocate2DArray(nCols,nCols);

    t5 = outer5D(nCols,nCols,nCols,nCols,nCols,mat5,mat5,mat5,mat5,mat5);
    scalar5DArrayMultip(nCols,nCols,nCols,nCols,nCols,945.,t5);

    A = identMatrix(nCols,nCols);
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
    return t5;
}

double* skew(double** vec){
/*  Returns contraction of supplied 3x3 matrix with Levi-Civita tensor.*/
    double* b = new double[3];
    b[0] = vec[1][2] - vec[2][1];
    b[1] = vec[2][0] - vec[0][2];
    b[2] = vec[0][1] - vec[1][0];

    return b;
}

double dotProduct(int n,double* a, double* b){

    double product{0};
    
    for (int i{0}; i<n; i++){
	product = product + (a[i] * b[i]);
    }
    return product;
}


double* crossProduct(double* a, double* b){

    double* product = new double [3];

    product[0] = a[1]*b[2] - a[2]*b[1];
    product[1] = a[2]*b[0] - a[0]*b[2];
    product[2] = a[0]*b[1] - a[1]*b[0];

    return product;
}

double** elementWiseProduct(int m, int n, double** a, double** b){
    
    double** product = allocate2DArray(3,3);
    for (int i{0};i<m;++i){
        for (int j{0};j<n;++j){
            product[i][j] = a[i][j] * b[i][j];
        }
    }
    return product;
}

double contract_i_i (int n,double* a, double* b){
/*  Returns a zero-rank contraction of a first-rank tensor
    with a first-rank tensor. */

    double c;
    c = dotProduct(n,a,b);
    
    return c;
}


double* contract_ij_j (double** a, double* b){
/*  Returns a first-rank contraction offirst-rank tensor
    with a first-rank tensor. */
    double* product = new double[3];
    for(int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    product[i] = product[i] + a[i][j] *  b[j];
	}
    }
    return product;
}	


double contract_ij_ij (double** a, double** b){
/*  c ! Returns a zero-rank contraction of a second-rank tensor
    with another second-rank tensor. */

    double product{0};
    double** dot;
    dot = elementWiseProduct(3,3,a,b);
//    print2DArray(3,3,dot);
    for (int i{0};i<3;++i){
        for (int j{0};j<3;++j){
            product = product + dot[i][j];
        }
    }
    return product;
}

double** contract_ik_jk (double** a, double** b){
/*  Returns a second-rank contraction of a second-rank tensor
    with another second-rank tensor. */

    double** D = allocate2DArray(3,3);

    for (int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    for (int k{0};k<3;++k){
		    D[i][j] = D[i][j] + a[i][k] * b[j][k];
	    }
	}
    }
    return D;
}


double** contract_ijk_k (double*** a, double* b){
/*  Returns a second-rank contraction of a third-rank tensor
    and a first-rank tensor. */

    double** D = allocate2DArray(3,3);

    for (int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    for (int k{0};k<3;++k){
                D[i][j] =  D[i][j] + a[i][j][k] * b[k];
	    }
	}
    }
    return D;
}

double* contract_ijk_jk (double*** a, double** b){
/*  Returns a second-rank contraction of a fourth-rank tensor
    and a second-rank tensor */

    double* D = new double[3];

    for (int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    for (int k{0};k<3;++k){
                D[i] =  D[i] + a[i][j][k] * b[j][k];
	    }
	}
    }
    return D;
}

double** contract_ijkl_kl(double**** a, double** b){
/*  Returns a second-rank contraction of a fourth-rank tensor
    and a second-rank tensor. */

    double** D = allocate2DArray(3,3);

    for (int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    for (int k{0};k<3;++k){
	        for (int l{0};l<3;++l){
                    D[i][j] =  D[i][j] + a[i][j][k][l] * b[k][l];
		}
	    }
	}
    }
    return D;
}


double*** contract_ijklm_lm(double***** a, double** b){
/*  Returns a third-rank contraction of a fifth-rank tensor
    and a second-rank tensor. */

    double*** D = allocate3DArray(3,3,3);

    for (int i{0};i<3;++i){
	for (int j{0};j<3;++j){
	    for (int k{0};k<3;++k){
	        for (int l{0};l<3;++l){
	            for (int m{0};m<3;++m){
                        D[i][j][k] = D[i][j][k] + a[i][j][k][l][m] * b[l][m];
		    }
		}
	    }
	}
    }
    return D;
}

double* random_vector(){
//  Returns a random unit vector as a numpy array of 3 elements. 

    double* rand_vec = new double [3];
    double* zeta     = rand1DArray(2);        // Two uniformly sampled random numbers in range (0,1)
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

    return rand_vec; 
 }
