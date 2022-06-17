#include <iostream>
#include <cmath>
#include <list>
#include <stdlib.h>
#include <string>
#include <vector>
#include "./math_module.hpp" 

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

            double* run_avg = new double[7];
            double* run_err = new double[7];
            double* blk_avg = new double[7];
            double* blk_msd = new double[7];
            double* values  = new double[7];
            double* addd    = new double[7];
            bool*   mask    = new bool[7];
            int*    methodd = new int[7];

};

//void time_stamp(){
////    Function to print date, time, and cpu time information.
//
//
//    printf("{:45}{}".format("Date:",time.strftime("%Y/%m/%d")))
//    printf("{:47}{}".format("Time:",time.strftime("%H:%M:%S")))
//    printf("{:40}{:15.6f}".format("CPU time:",time.process_time()))
//}

std::vector<std::string> word_splitter(std::string str);

void run_begin (std::vector<VariableType> vars, BlockVar &blk_var);

void blk_begin(int n_avg, BlockVar &blk_var);

void blk_add (std::vector<VariableType> vars, BlockVar &blk_var);

void blk_end (int blk, int n_avg, BlockVar &blk_var);

void run_end (std::vector<VariableType> vars, BlockVar &blk_var);
