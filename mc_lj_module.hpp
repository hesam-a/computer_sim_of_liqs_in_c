#include <iostream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "./maths_module.hpp"

class PotentialType{

    public:
	double pot; 
	double vir;
	double lap;
	bool   ovr;

};

void introduction();

void conclusion();

void potential_1 (PotentialType &partial,int mm, double* ri, double box, double r_cut, double** r);

void potential   (PotentialType &total, int mm, double box, double r_cut, double** r );

double force_sq (int mm,double box, double r_cut, double** r);
