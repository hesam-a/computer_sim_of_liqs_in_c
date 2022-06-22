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
#include "boost/format.hpp"
#include "./math_module.hpp" 
#include "./read_atoms.hpp"

// Calculation of run averages with output to stdout.

// Options for averaging methods
#define avg  0
#define msd  1
#define cke  2

// Internal variable
#define colf_fmt  "%15.6f"         // Format for floats; we assume that 6 dp will be sufficient
#define head_fmt  "%15s"           // Format for heading strings

int     n_avg;
int     col_width = 15;            // Must be large enough to allow sensible format
int     nam_width = 2*col_width+1; // At most two column widths plus spacer


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

	    double* run_avg = new double[11];
	    double* run_err = new double[11]; 
	    double* blk_avg = new double[11];
	    double* blk_msd = new double[11];
	    double* values  = new double[11];
	    double* addd    = new double[11];
	    bool*   mask    = new bool[11];
	    int*    methodd = new int[11];
};


void time_stamp(bool end, std::clock_t ti){
//  Function to print date, time, and cpu time information.
    time_t timer ;
    struct std::tm *tim;

    std::time(&timer);
    tim = localtime(&timer);

    std::cout << '\n';
    printf("Date:   %27d/%d/%d \n",1900+tim->tm_year,tim->tm_mon,tim->tm_mday);
    printf("Time:   %26d:%d:%d \n",tim->tm_hour,tim->tm_min,tim->tm_sec);

    std::clock_t tf = std::clock();
    double time_elapsed = (tf-ti) / CLOCKS_PER_SEC;

    if (end)
       printf("CPU time: %30.6f \n", time_elapsed);

}


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


void run_begin (std::vector<VariableType> vars, BlockVar &blk_var, std::clock_t ti){ 
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
    
    //First column plus a column for each variable; allow one space between columns
    int line_width = col_width + n_avg * ( col_width + 1 );

/*  Store variable names in module variables
    Attempt to split name string at first space
    Build up format string for line of averages*/
    std::vector<std::string> headings; 
    headings.push_back("Block");
    std::vector<std::string> subheads;
    subheads.push_back(" ");
    std::vector<std::string> line_fmt;

    for (auto variable : vars){
	std::vector<std::string> parts  = word_splitter(variable.nam);
        headings.push_back(parts[0]);
        if (parts.size() > 1)
            subheads.push_back(parts[1]);
        else
            subheads.push_back(" ");
    }
    
    // Zero averages and error accumulators at start of run
    blk_var.run_nrm = 0.0;

    for (int i{0}; i< n_avg;++i)
	    blk_var.run_avg[i] = 0.0;

    for (int i{0}; i<n_avg ;++i)
	    blk_var.run_err[i] = 0.0;

    time_stamp(false, ti)


    //Write headings
    std::cout << '\n';
    std::cout << "Run begins \n";
    std::cout << std::string(line_width,'=') << '\n';
    for (int i{0}; i< headings.size();++i)
	std::cout << boost::format(head_fmt) %headings[i];
    std::cout << '\n';
    for (int i{0}; i<subheads.size();++i)
	std::cout << boost::format(head_fmt) %subheads[i];
    std::cout << '\n';
    std::cout << std::string(line_width,'-') << '\n';
    
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

    for (int i=0; i<n_avg;++i)
        blk_var.values[i] = vars[i].val;

    for(std::size_t i{0};i<n_avg;++i)
	values2[i] = 0.0;

    blk_var.blk_avg = sum1DArrays(n_avg,blk_var.blk_avg,blk_var.values);
    values2         = elementWise1DProduct(n_avg,blk_var.values,blk_var.values);
    blk_var.blk_msd = sum1DArrays(n_avg,blk_var.blk_msd,values2);

    blk_var.blk_nrm = blk_var.blk_nrm + 1.0;

    delete [] values2;
}

void blk_end (int blk, int n_avg, BlockVar &blk_var){
//  Write out averages at end of every block.

    double* constraint = new double[11];
    double* blk_avg2   = new double[11];

    for (int i{0}; i< n_avg;++i)
	    blk_avg2[i] = 0.0;

    for (int i{0}; i< n_avg;++i)
	    constraint[i] = 0.0;

    assert (blk_var.blk_nrm>0.5);
//    std::cout << "Block accumulation error \n";

    blk_var.blk_avg = scalar1DArrayDivision(n_avg, blk_var.blk_nrm,blk_var.blk_avg);
    blk_var.blk_msd = scalar1DArrayDivision(n_avg, blk_var.blk_nrm,blk_var.blk_msd);
    
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
    constraint   = sum1DArrays(n_avg,blk_var.addd, constraint);

    for(int i{0};i<n_avg;++i){
	if (blk_var.mask[i]){
            blk_var.blk_avg[i] = constraint[i];
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

    std::vector<double> block_variables;
    block_variables.push_back(blk+1);

    for (int i{0}; i< n_avg;++i)
	block_variables.push_back(blk_var.blk_avg[i]);

    // Write out block averages
    std::cout << boost::format("%15.0f") %block_variables[0];
    for (int i{1}; i< n_avg+1;++i)
	std::cout << boost::format(colf_fmt) %block_variables[i];
    std::cout << '\n';

    delete [] constraint;
    delete [] blk_avg2  ;
}

void run_end (std::vector<VariableType> vars, BlockVar &blk_var, std::clock_t ti){ 
//  Write out averages and error estimates at end of run.

    int n_avg = vars.size();
    int line_width = col_width + n_avg * ( col_width + 1 );

    double* run_avg2    = new double[n_avg];
    double* constraint2 = new double[n_avg];

    for (int i{0}; i< n_avg;++i)
	    run_avg2[i] = 0.0;

    for (int i{0}; i< n_avg;++i)
	    constraint2[i] = 0.0;

    assert (blk_var.run_nrm>0.5);
    //std::cout << "Run accumulation error \n";

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
    for (int i{0};i<n_avg;++i)
	    constraint2[i] = abs(constraint2[i]);

    constraint2 = sqrt1DArray(n_avg,constraint2);
    for(int i{0};i<n_avg;++i){
        if (blk_var.run_err[i]>0.0){
            blk_var.run_err[i] = constraint2[i];
	}
	else{
	    blk_var.run_err[i] = 0.0;
	}
    }

    std::cout << '\n';
    std::cout << std::string(line_width,'-') << '\n';
    std::cout << boost::format("%15s") %"Run averages";
    for (int i{0}; i< n_avg;++i)
	std::cout << boost::format(colf_fmt) %blk_var.run_avg[i];
    std::cout << '\n';
    std::cout << boost::format("%15s") %"Run errors";
    for (int i{0}; i<n_avg;++i)
	std::cout << boost::format(colf_fmt) %blk_var.run_err[i];
    std::cout << '\n';
    std::cout << std::string(line_width,'=') << '\n';
    std::cout << '\n';
    std::cout << "Run ends \n";
    std::cout << '\n';

    time_stamp(true, ti)
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

    delete [] run_avg2  ; 
    delete [] constraint2;
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
