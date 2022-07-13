#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <cassert>
#include <numeric>
#include <stdlib.h>
#include <iterator>
#include "./math_module.hpp"

// ------  Routines for MC simulation, polyatomic molecule, LJ atoms  --------

class PotentialType{

    public:
	double pot;  // = 0.0; 
	double vir;  // = 0.0;
	double lap;  // = 0.0;
	bool   ovr;  // = false;

};

void introduction(int na,double** db,double diameter);

void conclusion();

void potential_1 (PotentialType &partial,int mm, double r_cut, double diameter, double* ri, double** ddi, double box, double** r, double*** dd );

void potential (PotentialType &total, int mm,  double box,double r_cut, double diameter, double** r, double*** dd);

