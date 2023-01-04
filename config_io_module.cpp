#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <stdio.h>
#include <vector>
#include <cassert>
#include <stdarg.h>
#include "boost/format.hpp"
#include "./maths_module.hpp"

#define line_begin "%15s" 
#define coord_dist "%15.8f"

double** read_cnf_atoms(const char* filename,double** coord){
//  Read in atomic configuration.

    std::ifstream input(filename);

    int natom;
    double box_size;

    input >> natom;       // Number of atoms
    input >> box_size;    // Simulation box length (assumed cubic)

    for (int i{0}; i< natom; ++i){
        input >> coord[i][0] >> coord[i][1] >> coord[i][2];
    }
    input.close();

    return coord;
}

void read_cnf_mols(const char* filename, bool quaternion, bool with_v, double** coord, double** orient, double** vel, double** angvel){
//  Read in molecular configuration.

    std::ifstream input(filename);

    int natom;
    double box_size;

    input >> natom;      // Number of atoms
    input >> box_size;   // Simulation box length (assumed cubic)

    if (quaternion && with_v){
        for (int i{0}; i< natom; ++i)    //The rest of the file should be uniformly formatted
            input >> coord[i][0] >> coord[i][1] >> coord[i][2] >> orient[i][0] >> orient[i][1] >> orient[i][2] >> orient[i][3] >> vel[i][0] >> vel[i][1] >> vel[i][2] >> angvel[i][0] >> angvel[i][1] >> angvel[i][2];
    }
    else if (quaternion && !with_v){
        for (int i{0}; i< natom; ++i)    //The rest of the file should be uniformly formatted
            input >> coord[i][0] >> coord[i][1] >> coord[i][2] >> orient[i][0] >> orient[i][1] >> orient[i][2] >> orient[i][3];
    }
    else if (!quaternion && !with_v){
        for (int i{0}; i< natom; ++i)    //The rest of the file should be uniformly formatted
            input >> coord[i][0] >> coord[i][1] >> coord[i][2] >> orient[i][0] >> orient[i][1] >> orient[i][2];
    }
    else if (!quaternion && with_v){
        for (int i{0}; i< natom; ++i)    //The rest of the file should be uniformly formatted
            input >> coord[i][0] >> coord[i][1] >> coord[i][2] >> vel[i][0] >> vel[i][1] >> vel[i][2];
    }
}


void write_cnf_atoms (std::string filename, int nn, double box, double** coord){
//  Write out atomic configuration.
    std::ofstream output(filename);

    output << boost::format(" %15d \n") %nn;      // Number of atoms
    output << boost::format(" %15.8f \n") %box;  // Simulation box length (assumed cubic)

    for (int i{0};i<nn;++i){
        for (int j{0};j<3;++j){
            output << boost::format(" %15.10f ") %coord[i][j];
	}
	output << '\n';
    }
}

void write_cnf_mols (std::string filename,int nn,double box,bool quaternion,bool with_v,double** coord,double** orient,double** vel,double** angvel){
//  Write out atomic configuration.

    std::ofstream output(filename);

    output << boost::format(" %15d \n") %nn;      // Number of atoms
    output << boost::format(" %15.8f \n") %box;  // Simulation box length (assumed cubic)

    if (quaternion && with_v){
        for (int i{0};i<nn;++i){
            for (int j{0};j<3;++j)
                output << boost::format(" %15.10f ") %coord[i][j];
            for (int j{0};j<4;++j)
                output << boost::format(" %15.10f ") %orient[i][j];
            for (int j{0};j<3;++j)
                output << boost::format(" %15.10f ") %vel[i][j];
            for (int j{0};j<3;++j)
                output << boost::format(" %15.10f ") %angvel[i][j];
            output << '\n';
        }
    }
    else if (quaternion && !with_v){ 
        for (int i{0};i<nn;++i){
            for (int j{0};j<3;++j)
                output << boost::format(" %15.10f ") %coord[i][j];
            for (int j{0};j<4;++j)
                output << boost::format(" %15.10f ") %orient[i][j];
            output << '\n';
        }
    }
    else if (!quaternion && !with_v){ 
        for (int i{0};i<nn;++i){
            for (int j{0};j<3;++j)
                output << boost::format(" %15.10f ") %coord[i][j];
            for (int j{0};j<3;++j)
                output << boost::format(" %15.10f ") %orient[i][j];
            output << '\n';
        }
    }
    else if (!quaternion && with_v){ 
        for (int i{0};i<nn;++i){
            for (int j{0};j<3;++j)
                output << boost::format(" %15.10f ") %coord[i][j];
            for (int j{0};j<3;++j)
                output << boost::format(" %15.10f ") %vel[i][j];
            output << '\n';
        }
    }
}
