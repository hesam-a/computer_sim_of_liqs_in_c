#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES


//  Long-range and delta corrections for potential energy and pressure.

/*  The purpose of this module is simply to gather in one place the common
    functions for long-range and delta corrections for the Lennard-Jones potential.
    If a different potential is used in the simulation, a different file
    (with the same module name) containing different expressions should be substituted. */

double potential_lrc (double density,double r_cut ){
//  Calculates long-range correction for Lennard-Jones potential per atom.

    // density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    double sr3 = 1.0 / pow(r_cut,3);
    return M_PI * ( (8.0/9.0)  * pow(sr3,3)  - (8.0/3.0)  * sr3 ) * density;
}

double pressure_lrc (double density, double r_cut ){
//  Calculates long-range correction for Lennard-Jones pressure.

    // density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    double sr3 = 1.0 / pow(r_cut,3);
    return M_PI * ( (32.0/9.0) * pow(sr3,3)  - (16.0/3.0) * sr3 ) * pow(density,2);
}

double pressure_delta (double density, double r_cut ){
//  Calculates correction for Lennard-Jones pressure due to discontinuity in the potential at r_cut.

    // density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    double sr3 = 1.0 / pow(r_cut,3);
    return M_PI * (8.0/3.0) * (pow(sr3,3)  - sr3 ) * pow(density,2);
}
