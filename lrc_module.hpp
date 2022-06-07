#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES


//  Long-range and delta corrections for potential energy and pressure.

/*  The purpose of this module is simply to gather in one place the common
    functions for long-range and delta corrections for the Lennard-Jones potential.
    If a different potential is used in the simulation, a different file
    (with the same module name) containing different expressions should be substituted. */

double potential_lrc (double density,double r_cut );

double pressure_lrc (double density, double r_cut );

double pressure_delta (double density, double r_cut );
