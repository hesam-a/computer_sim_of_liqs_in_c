#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <stdio.h>
#include "./math_module.hpp"

//class Molecule {
//
//    public:
//
//        int natom;       // Number of atoms
//        double box_size;      // Simulation box length (assumed cubic)
//        double** coord; 
//
////        Molecule();
////        ~Molecule();
//};


double** read_cnf_atoms(const char* filename,double** coord);

void read_cnf_mols(const char* filename, bool quaternion, bool with_v, double** coord, double** orient, double** vel, double** angvel);

void write_cnf_atoms (std::string filename, int nn, double box, double** coord);

void write_cnf_mols (std::string filename, int nn, double box,bool quaternion, bool with_v, double** coord, double** orient, double** vel, double** angvel);
