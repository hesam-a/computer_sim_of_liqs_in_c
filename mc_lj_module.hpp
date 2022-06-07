#include <iostream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "./math_module.h"

class PotentialType{

    public:
	double pot; 
	double vir;
	double lap;
	bool   ovr;

};

void introduction();

void conclusion();

void potential (PotentialType &total, int mm, double box, double r_cut, double** r );

void potential_1 (PotentialType &partial,int mm, double* ri, double box, double r_cut, double** r);

void force_sq (int mm,double box, double r_cut, double** r);
