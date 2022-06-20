#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <stdio.h>
#include <vector>
#include "./math_module.hpp"


double** read_cnf_atoms(const char* filename,double** coord){

    std::ifstream input(filename);

    int natom;
    double box_size;

    input >> natom;
    input >> box_size;

    for (int i{0}; i< natom; ++i){
        input >> coord[i][0] >> coord[i][1] >> coord[i][2];
    }
    input.close();

    return coord;
}
