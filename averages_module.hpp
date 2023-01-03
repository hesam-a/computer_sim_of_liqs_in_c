#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <ctime>
#include <vector>
#include "./maths_module.hpp" 

// Options for averaging methods
#define avg  0
#define msd  1
#define cke  2

class VariableType{
//    Class encapsulating the essential information for simulation averages.

    public:

        const char* nam; 
        double      val;
        double      add = 0.0;
        int         method   = avg;
        bool        e_format = false;
        bool        instant  = true; 
};

class BlockVar{

    public:
            double blc_nrm = 0.0;
            double run_nrm = 0.0;

            double* run_avg = new double[12];
            double* run_err = new double[12];
            double* blk_avg = new double[12];
            double* blk_msd = new double[12];
            double* values  = new double[12];
            double* addd    = new double[12];
            bool*   mask    = new bool[12];
            int*    methodd = new int[12];

};

void time_stamp(bool end, std::clock_t ti);

std::vector<std::string> word_splitter(std::string str);

void run_begin (std::vector<VariableType> vars, BlockVar &blk_var,std::clock_t ti);

void blk_begin(int n_avg, BlockVar &blk_var);

void blk_add (std::vector<VariableType> vars, BlockVar &blk_var);

void blk_end (int blk, int n_avg, BlockVar &blk_var);

void run_end (std::vector<VariableType> vars, BlockVar &blk_var, std::clock_t ti);
