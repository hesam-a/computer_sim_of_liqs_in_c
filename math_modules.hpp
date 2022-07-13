#include <iostream>

double**    allocate2DArray(int m,int n);

double***   allocate3DArray(int p,int q, int r);

double****  allocate4DArray(int p,int q, int m, int n);

double***** allocate5DArray(int p,int q, int m, int n, int s);

// deallocate the memory of a 2D matrix
void free2DArray(int m,double** array);

// deacllocate the memory of a 3D matrix
void free3DArray(int p,int q, double*** array);

// deacllocate the memory of a 4D matrix
void free4DArray(int p,int q,int m,double**** array);

// deacllocate the memory of a 5D matrix
void free5DArray(int p,int q,int m, int n, double***** array);

void rand1DArray(int m, double* arr);

void rand2DArray(int m,int n, double** arr);

void rand3DArray(int p,int q,int r, double ***arr);

void rand4DArray(int p,int q,int m, int n, double ****arr);

void rand5DArray(int p,int q,int m, int n,int s, double *****arr);

void zeroMatrix (int m,int n,double** arr);

void identMatrix(int m,int n,double** arr);

void matMultip(int m, int n, double** A, double** B, double** C);

void scalar1DArrayMultip(int m,double p, double* A, double* arr );

void scalar2DArrayMultip(int m,int n,double p, double** A, double** product);

void scalar3DArrayMultip(int m,int n, int r, double p, double*** A, double*** product);

void scalar4DArrayMultip(int p,int q, int m, int n, double b, double**** A, double**** product);

void scalar5DArrayMultip(int p,int q, int m, int n, int s, double b, double***** A, double***** product);

void scalar2DArraySubtract(int m,int n, double p, double** A, double** product);

void scalar1DArraySubtract(int m, double p, double* A);

void subtract1DArrays(int m, double* A, double* B, double* C);

void subtract2DArrays(int m,int n, double** A, double** B, double** C);

void sum1DArrays(int m, double* A, double* B, double* sum);

void sum2DArrays(int m, int n,double** A, double* sum,int axis);

double max1DArray(int m, double* arr);

void scalar1DArrayDivision(int m,double p, double* A, double* product);

void sqrt1DArray(int m, double*A, double* sq);

void scalar2DArrayDivision(int m,int n,double p, double** A);

void scalar3DArrayDivision(int m,int n,int r, double p, double*** A);

void scalar4DArrayDivision(int p,int q, int m,int n, double b, double**** A) ;

void scalar5DArrayDivision(int p,int q, int m,int n, int s, double b, double***** A);

void print1DArray(int m, double* A);

void print2DArray(int m,int n, double** A);

void print3DArray(int p, int q, int r, double*** A);

void print4DArray(int p, int q, int m, int n, double**** A);

void print5DArray(int p, int q, int m, int n, int s, double***** A);

// outer functions

double**    outer2D(int nCols, double* A, double* B);

double***   outer3D(int nRows, int nCols, int n3rd, double* A, double* B, double* C);

double****  outer4D(int nRows, int nCols, int n3rd, int n4th, double* A, double* B, double* C, double* D);

double***** outer5D(int nRows, int nCols, int n3rd, int n4th, int n5th, double* A, double* B, double* C, double* D, double* E);

double** t2_tensor (double* mat,double r3);

double*** t3_tensor(double *mat3, double r4);

double**** t4_tensor(double* mat4, double r5);

double***** t5_tensor(double* mat5, double r6);

double* skew(double** vec);

double dotProduct1D(int n,double* a);

double dotProduct2D(int n,double* a, double* b);

void elementWise2DProduct(int m, int n, double** a, double** b,double** product);

void elementWise1DProduct(int m, double* a,double* b, double* product);

double elementSum1D(int m, double* a);

double elementSum2D(int m, int n, double** a);

double* crossProduct(double* a, double* b);

double contract_i_i (int n,double* a, double* b);

double* contract_ij_j (double** a, double* b);

double contract_ij_ij (double** a, double** b);

double** contract_ik_jk (double** a, double** b);

double** contract_ijk_k (double*** a, double* b);

double* contract_ijk_jk (double*** a, double** b);

double** contract_ijkl_kl(double**** a, double** b);

double*** contract_ijklm_lm(double***** a, double** b);

void random_vector(double* rand_vec);

void rint1D(int m,double* A);

void rint2D(int m,int n, double** A);

void remove2DArray(int m, int p, double** A, double** B);

void remove3DArray(int m, int n, int p, double*** A, double*** B);

void random_translate_vector (double dr_max, double* old, double* ri);

void quatmul (double* a, double* b, double* c);

void rotate_quaternion (double angle, double* axis,double* old, double* e);

void random_quaternion(double* randq);

void random_rotate_quaternion (double angle_max, double* old, double* e);

bool metropolis (double delta );

void update2DArray(int n, int atom, double* ri, double** A);

void update3DArray(int m,int n, int atom, double** ri, double*** A);

double** createArray(int m,int n, double* ri, double** A);

double** annihilateArray(int m,int n, double* ri, double** A);

void q_to_a (double* q, double** a);
