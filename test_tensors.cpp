// This program is used to test the function for operations on tensors.

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "./math_module.hpp"



int main()
{

    int nRows=1;
    int nCols=3;
    double contract2;
    double r3 = 4.;
    double r4 = 4.;
    double r5 = 4.;
    double r6 = 4.;

    double*     r = new double[nCols];
    double*     a = new double[nCols];
    double*     contract = new double[nCols];
    double**    contract3 = allocate2DArray(nCols,nCols);
    double**    contract4 = allocate2DArray(nCols,nCols);
    double**    contract5 = allocate2DArray(nCols,nCols);
    double***    contract6 = allocate3DArray(nCols,nCols,nCols);
    double**    D = allocate2DArray(nCols,nCols);
    double**    b = allocate2DArray(nCols,nCols);
    double***   T = allocate3DArray(nCols,nCols,nCols);
    double****  F = allocate4DArray(nCols,nCols,nCols,nCols);
    double***** G = allocate5DArray(nCols,nCols,nCols,nCols,nCols);

    std::cout << "1D Array r:" << '\n';
    r = rand1DArray(nCols);
    print1DArray(nCols,r);
    std::cout << '\n';

    D = t2_tensor(r,r3);
    std::cout << "2D Array D:" << '\n';
    print2DArray(nCols,nCols,D);
    std::cout << '\n';

    T = t3_tensor(r,r4);
    std::cout << "3D Array T:" << '\n';
    print3DArray(nCols,nCols,nCols,T);

    std::cout << "4D Array F:" << '\n';
    F = t4_tensor(r,r5);
    print4DArray(nCols,nCols,nCols,nCols,F);

    std::cout << "5D Array G:" << '\n';
    G = t5_tensor(r,r6);
    print5DArray(nCols,nCols,nCols,nCols,nCols,G);
    std::cout << '\n';

    b = rand2DArray(nCols,nCols);
    std::cout << "randomized 2D Array b:" << '\n';
    print2DArray(nCols,nCols,b);
    std::cout << "skew:" << '\n';
    a = skew(b);
    print1DArray(nCols,a);
    std::cout << '\n';

    std::cout << "contract:" << '\n';
    contract = contract_ij_j(b,a);
    print1DArray(nCols,contract);

    std::cout << "contract ij ij:" << '\n';
    contract2 = contract_ij_ij(D,b);
    std::cout << contract2 << '\n';
    std::cout << '\n';

    std::cout << "contract ik jk:" << '\n';
    contract3 = contract_ik_jk(b,b);
    print2DArray(nCols,nCols,contract3);
    std::cout << '\n';

    std::cout << "contract ijk k:" << '\n';
    contract4 = contract_ijk_k(T,contract);
    print2DArray(nCols,nCols,contract4);
    std::cout << '\n';

    std::cout << "contract ijkl kl:" << '\n';
    contract5 = contract_ijkl_kl(F,contract4);
    print2DArray(nCols,nCols,contract5);
    std::cout << '\n';

    std::cout << "contract ijklm lm:" << '\n';
    contract6 = contract_ijklm_lm(G,contract5);
    print3DArray(nCols,nCols,nCols,contract6);
    std::cout << '\n';

    delete [] r;
    delete [] a;
    delete [] contract;
    free2DArray(nCols,D);
    free2DArray(nCols,b);
    free2DArray(nCols,contract3);
    free2DArray(nCols,contract4);
    free2DArray(nCols,contract5);
    free3DArray(nCols,nCols,T);
    free3DArray(nCols,nCols,contract6);
    free4DArray(nCols,nCols,nCols,F);
    free5DArray(nCols,nCols,nCols,nCols,G);

    return 0;

}
