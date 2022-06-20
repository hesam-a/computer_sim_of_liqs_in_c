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
