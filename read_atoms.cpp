#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <stdio.h>
#include <vector>
#include "./math_module.hpp"


std::vector<std::vector<double>> read_cnf_atoms(const char* filename,std::vector<std::vector<double>> coord){

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

  
//    std::cout << "Number of atoms: " << natom << '\n';
//    std::cout << "Input Cartesian coordinates:\n";
//    for(int i=0; i < natom; i++)
//      printf("%s %20.12f %20.12f %20.12f\n",  atomName[i].c_str(), x[i], y[i], z[i]);
