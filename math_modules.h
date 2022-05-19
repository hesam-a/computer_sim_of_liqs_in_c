#include <iostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

double** allocate2DArray(int m,int n);

double*** allocate3DArray(int p,int q, int r);

double**** allocate4DArray(int p,int q, int m, int n);

double***** allocate5DArray(int p,int q, int m, int n, int s);

// deallocate the memory of a 2D matrix
void free2DArray(int m,double** array);

// deacllocate the memory of a 3D matrix
void free3DArray(int p,int q, double*** array);

// deacllocate the memory of a 4D matrix
void free4DArray(int p,int q,int m,double**** array);

// deacllocate the memory of a 5D matrix
void free5DArray(int p,int q,int m, int n, double***** array);

void rand1DArray(int m, double *A);

void rand2DArray(int m,int n, double **A);

void rand3DArray(int p,int q,int r, double ***A);

void rand4DArray(int p,int q,int m, int n, double ****A);

void rand5DArray(int p,int q,int m, int n,int s, double *****A);

void zeroMatrix(int m,int n, double **A);

void identMatrix(int m,int n, double **A);

void scalar2DArrayMultip(int m,int n,double p, double** A);

void scalar3DArrayMultip(int m,int n, int r, double p, double*** A);

void scalar4DArrayMultip(int p,int q, int m, int n, double b, double**** A);

void scalar5DArrayMultip(int p,int q, int m, int n, int s, double b, double***** A);

void scalar2DArraySubtract(int m,int n, double** A, double** B);

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

double dotProduct(double* a, double* b);

double** elementWiseProduct(double** a, double** b);

double* crossProduct(double** a, double** b);

double contract_i_i (int n,double* a, double* b);

double* contract_ij_j (double* a, double** b);

double contract_ij_ij (double** a, double** b);


