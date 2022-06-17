#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <cassert>
#include <algorithm>
#include <functional>
#include "./math_module.hpp" 
#include "./read_atoms.hpp"

// Calculation of run averages with output to stdout.

// Options for averaging methods
#define avg  0
#define msd  1
#define cke  2


// Internal variable
int      n_avg;

class VariableType{
//    Class encapsulating the essential information for simulation averages.

    public:
        const char* nam; 
        double      val;
        double      add=0.0;
        int         method   = avg;
        bool        e_format = false;
        bool        instant  = true; 
};

class BlockVar{

    public:
	    double blk_nrm = 0.0;
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

std::vector<std::string> word_splitter(std::string str){
//   A function for splitting the words. 

    size_t      pos          = 0;
    std::string delimit      = " ";
    std::string split;
    std::vector<std::string> tok;

    while((pos = str.find(delimit)) != std::string::npos){

        split = str.substr(0, pos);
        tok.push_back(split);
        str.erase(0, pos + delimit.length());
    }

    tok.push_back(str);

    return tok;
}


void run_begin (std::vector<VariableType> vars, BlockVar &blk_var){ 
//  Set up averaging variable based on supplied list of names & other attributes.

    int n_avg = vars.size();

    bool need_header = true;
    for (auto variable: vars){
        if (variable.instant){
	    if (need_header){
		std::cout << "Initial values"<< "\n";
                need_header = false;
	    }
            printf("%10s %29.6f \n",variable.nam,variable.val);
	}
    }

    // Store method options and add-constants in module variable
        for(std::size_t i = 0; i < n_avg; ++i){
	    blk_var.methodd[i] = vars[i].method;
            blk_var.addd[i]    = vars[i].add;
    }
    
    std::cout << "\n";

    // Zero averages and error accumulators at start of run
    blk_var.run_nrm = 0.0;

    for (int i{0}; i< n_avg;++i)
	    blk_var.run_avg[i] = 0.0;

    for (int i{0}; i<n_avg ;++i)
	    blk_var.run_err[i] = 0.0;

    std::cout << '\n';
    std::cout << "Run begins \n";
    std::cout << "=============================================================================================================================== \n";
    printf("%15s %15s %15s %15s %15s %15s %15s %15s \n","Block","Move","E/N","P","E/N","P", "T","Cv/N");
    printf("%15s %15s %15s %15s %15s %15s %15s %15s \n","  ","ratio","cut","cut","full","full","config","full");
    std::cout << "------------------------------------------------------------------------------------------------------------------------------- \n";

}

void blk_begin(int n_avg, BlockVar &blk_var){
//  Zero average variable at start of each block.

    blk_var.blk_nrm = 0.0; 
    for (int i{0}; i< n_avg;++i)
	blk_var.blk_avg[i] = 0.0;

    for (int i{0}; i< n_avg;++i)
	blk_var.blk_msd[i] = 0.0;

}

void blk_add (std::vector<VariableType> vars, BlockVar &blk_var){ 
//    Increment block-average variable.
 
    int n_avg = vars.size();
    double* values2= new double[n_avg]; 

    if(!(vars.size() == n_avg)){
        std::cout << "Mismatched variable arrays";
    }

    //for (std::size_t i=0; i<n_avg;++i){
    for (int i=0; i<n_avg;++i){
        blk_var.values[i] = vars[i].val;
        //std::cout  << vars[i].nam << ":  " << vars[i].val << "\n";
        //std::cout << "blk_var.values: " << blk_var.values[i] << "\n";
    }

    for(int i{0};i<n_avg;++i)
	values2[i] = 0.0;


    blk_var.blk_avg = sum1DArrays(n_avg,blk_var.blk_avg,blk_var.values);
    values2         = elementWise1DProduct(n_avg,blk_var.values,blk_var.values);
    blk_var.blk_msd = sum1DArrays(n_avg,blk_var.blk_msd,values2);
    for (std::size_t i=0; i<n_avg;++i)
        std::cout << "blk_msd:  " << blk_var.blk_msd[i] << "\n";
    blk_var.blk_nrm = blk_var.blk_nrm + 1.0;

    delete [] values2;
}

void blk_end (int blk, int n_avg, BlockVar &blk_var){
//  Write out averages at end of every block.

    double* constraint = new double[7];
    double* blk_avg2   = new double[7];

    for (int i{0}; i< n_avg;++i)
	    blk_avg2[i] = 0.0;

    for (int i{0}; i< n_avg;++i)
	    constraint[i] = 0.0;

    if (!(blk_var.blk_nrm>0.5))
        std::cout << "Block accumulation error \n";
    
    for(int i=0;i<n_avg;++i){
        blk_var.blk_avg[i] = blk_var.blk_avg[i] / blk_var.blk_nrm;
    }

    for(int i=0;i<n_avg;++i){
        blk_var.blk_msd[i] = blk_var.blk_msd[i] / blk_var.blk_nrm;
    }

    // Replace blk_avg by mean-squared deviations plus optional constant where required

    for (int i{0};i<n_avg;++i){ 
        if (blk_var.methodd[i] == msd)
            blk_var.mask[i] = true;
	
	else if(blk_var.methodd[i] == cke)
	    blk_var.mask[i] = true;
	
	else
	    blk_var.mask[i] = false;	
    }

    std::cout << '\n';
    blk_avg2     = elementWise1DProduct(n_avg,blk_var.blk_avg,blk_var.blk_avg);
    constraint   = subtract1DArrays(n_avg,blk_var.blk_msd, blk_avg2);
    blk_var.addd = sum1DArrays(n_avg,blk_var.addd, constraint);

    for(int i{0};i<n_avg;++i){
	if (blk_var.mask[i]){
            blk_var.blk_avg[i] = blk_var.addd[i];
	}
	else
	    continue;
    }

    //for (auto i : method){
    //    if (i == cke){
    //       cke_calc();             // Call special routine for Cv from KE fluctuations
    //    }
    //}

    blk_var.run_avg = sum1DArrays(n_avg, blk_var.run_avg,blk_var.blk_avg);
    blk_var.run_err = sum1DArrays(n_avg, blk_var.run_err,blk_avg2);
    blk_var.run_nrm = blk_var.run_nrm + 1.0;        // Increment run normalizer

    // Write out block averages
    printf("%15i %15f %15f %15f %15f %15f %15f %15f\n", blk+1,blk_var.blk_avg[0],blk_var.blk_avg[1],blk_var.blk_avg[2],blk_var.blk_avg[3],blk_var.blk_avg[4],blk_var.blk_avg[5],blk_var.blk_avg[6]);

}

void run_end (std::vector<VariableType> vars, BlockVar &blk_var){ 
//  Write out averages and error estimates at end of run.

    int n_avg = vars.size();

    double* run_avg2    = new double[n_avg];
    double* constraint2 = new double[n_avg];

    for (int i{0}; i< n_avg;++i)
	    run_avg2[i] = 0.0;

    for (int i{0}; i< n_avg;++i)
	    constraint2[i] = 0.0;

    if (!(blk_var.run_nrm>0.5))
        std::cout << "Run accumulation error \n";

    // NB, these are the crudest possible error estimates, based on the wholly unjustified
    // assumption that the blocks are statistically independent
    // For a discussion of errors, see Chapter 8 and the error_calc.py example

    for(int i=0;i<n_avg;++i)
        blk_var.run_avg[i] = blk_var.run_avg[i] / blk_var.run_nrm;

    for(int i=0;i<n_avg;++i)
        blk_var.run_err[i] = blk_var.run_err[i] / blk_var.run_nrm;
 

    run_avg2        = elementWise1DProduct(n_avg,blk_var.run_avg,blk_var.run_avg);
    blk_var.run_err = subtract1DArrays(n_avg,blk_var.run_err, run_avg2);             // Compute fluctuations of block averages

    // Normalize and get estimated errors guarding against roundoff
    constraint2 = scalar1DArrayDivision(n_avg, blk_var.run_nrm,blk_var.run_err);
    constraint2 = sqrt1DArray(n_avg,constraint2);
    for(int i{0};i<n_avg;++i){
        if (blk_var.run_err[i]>0.0){
            blk_var.run_err[i] = constraint2[i];
	}
	else{
	    blk_var.run_err[i] = 0.0;
	}
    }

    std::cout << "------------------------------------------------------------------------------------------------------------------------------- \n";
    printf("%15s %15f %15f %15f %15f %15f %15f %15f \n","Run averages",blk_var.run_avg[0],blk_var.run_avg[1],blk_var.run_avg[2],blk_var.run_avg[3],blk_var.run_avg[4],blk_var.run_avg[5],blk_var.run_avg[6]);
    printf("%15s %15f %15f %15f %15f %15f %15f %15f \n","Run errors",  blk_var.run_err[0],blk_var.run_err[1],blk_var.run_err[2],blk_var.run_err[3],blk_var.run_err[4],blk_var.run_err[5],blk_var.run_err[6]);
    std::cout << "=============================================================================================================================== \n";
    std::cout << '\n';
    std::cout << "Run ends \n";
    std::cout << '\n';
    //time_stamp()
    std::cout << '\n';

    bool need_header = true;
    for (auto variable: vars){
        if (variable.instant){
	    if (need_header){
		std::cout << "Final values"<< "\n";
                need_header = false;
	    }
            printf("%10s %29.6f \n",variable.nam,variable.val);
	}
    }
}

//void cke_calc(){
////  Special fluctuation formula for microcanonical heat capacity.
//
//    // Locate variable corresponding to kinetic temperature
//
//    bool found=false;
//    for head,subhead,temperature in zip(headings,subheads,blk_avg){
//        if ("T" in head+subhead and "kin" in head+subhead){
//            found=true;
//            break
//	}
//    }
//
//    if (assert found){
//	std::cout << "Could not find T kin variable \n";
//    }
//
//    // Apply special fluctuation formula for microcanonical ensemble heat capacity
//    // blk_avg[i] should contain mean-squared total KE, divided by N
//
//    blk_avg = np.where ( method==cke, 9.0 / ( 6.0 - 4.0 * blk_avg / temperature**2 ), blk_avg )
//}
