
// This program is a C++ version of the "t_tensor.f90" program from the book "Computer simulations of liquids"
// by Allen and Tildesley
// It is used for calculation of electrostatic interactions between two linear molecules including dipole-dipole,
// dipole-quadrupole and quadrupole-quadrupole, using T tensors and Euler angles.

// While the book itself does not provide detailed discussion about the equations, "The theory of intermolecular 
// forces" by Anthony Stones is a good source for derivations and further discussion.

// some expalanations from Allen-Tildesley:

/* The dipole moment of molecule 1 is aligned along the axial vector e1
   The quadrupole tensor, quad1, is diagonal and traceless with
   quad1_xx = -0.5*quad1_mag, quad1_yy = -0.5*quad1_mag, and quad1_zz = quad1_mag in the molecule-fixed system.
   Similarly for molecule 2
   The vector r12 = r1-r2 points from 2 to 1.
   
   Forces are calculated by differentiating the T-tensor, giving the next higher rank T-tensor
   Torques are calculated from the angular dependence of dipole, quadrupole etc.
   potential V = mu_i g_i => torque tau_i = -epsilon_ijk mu_j g_k (=-cross_product)
   potential V = Q_ij G_ij => torque tau_l = -2 epsilon_lij Q_ik G_jk
   where ijkl are Cartesian indices and epsilon is the Levi-Civita symbol
   It is just necessary to identify the constants g_i, G_ij, in terms of the T tensor and the
   multipole on the other molecule. */

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "./math_module.h"


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
    double**    D = allocate2DArray(nCols,nCols);
    double**    b = allocate2DArray(nCols,nCols);
    double***   T = allocate3DArray(nCols,nCols,nCols);
    double****  F = allocate4DArray(nCols,nCols,nCols,nCols);
    double***** G = allocate5DArray(nCols,nCols,nCols,nCols,nCols);

    std::cout << "1D Array r:" << '\n';
    rand1DArray(nCols,r);
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

    std::cout << "5D Array F:" << '\n';
    G = t5_tensor(r,r6);
    print5DArray(nCols,nCols,nCols,nCols,nCols,G);
    std::cout << '\n';

    rand2DArray(nCols,nCols,b);
    std::cout << "randomized 2D Array b:" << '\n';
    print2DArray(nCols,nCols,b);
    std::cout << "skew:" << '\n';
    a = skew(b);
    print1DArray(nCols,a);
    std::cout << '\n';
   
    std::cout << "contract:" << '\n';
    contract = contract_ij_j(a,b);
    print1DArray(nCols,contract);
   
    std::cout << "contract ij ij:" << '\n';
    contract2 = contract_ij_ij(D,b);
    std::cout << contract2 << '\n';
    std::cout << '\n';
   

    delete [] r;
    delete [] a;
    delete [] contract;
    free2DArray(nCols,D);
    free2DArray(nCols,b);
    free3DArray(nCols,nCols,T);
    free4DArray(nCols,nCols,nCols,F);
    free5DArray(nCols,nCols,nCols,nCols,G);

    return 0;

}
