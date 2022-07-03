#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <stdio.h>
#include <vector>
#include "boost/format.hpp"
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


void write_cnf_atoms (std::string filename, int nn, double box, double** coord){
//  Write out atomic configuration (only position).
    std::ofstream output(filename);

    output << boost::format("%5d \n") %nn;
    output << boost::format("%10.6f \n") %box;

    for (int i{0};i<nn;++i){
        for (int j{0};j<3;++j){
            output << boost::format(" %10.6f ") %coord[i][j];
        }
        output << '\n';
    }
}
