#include <iostream>
#include <cmath>
#include <stdio.h>
#include <string>
#include <cstring>
#include <iomanip>
#include <time.h>
#include <numeric>
#include <stdlib.h>
#include "./math_module.hpp"


// Default parameters
#define  nCols      3   // Number of columns
#define  d_min      0.5 // Minimum separation
#define  d_max      1.5 // Maximum separation
#define  mu1_mag    1.0 // Dipole moment of molecule 1
#define  mu2_mag    1.0 // Dipole moment of molecule 2
#define  quad1_mag  1.0 // Quadrupole moment of molecule 1
#define  quad2_mag  1.0 // Quadrupole moment of molecule 2

int main()
{

    double  r12_mag, c1, c2, c12, v12t, v12e, h;

    double*     g       = new double[nCols];
    double*     h1      = new double[nCols];
    double*     h2      = new double[nCols];
    double*     e1      = new double[nCols];
    double*     e2      = new double[nCols];
    double*     mu1     = new double[nCols];
    double*     mu2     = new double[nCols];
    double*     t1t     = new double[nCols];
    double*     t1tt1e  = new double[nCols];
    double*     t2tt2e  = new double[nCols];
    double*     t2t     = new double[nCols];
    double*     t1e     = new double[nCols];
    double*     t2e     = new double[nCols];
    double*     f12t    = new double[nCols];
    double*     f12e    = new double[nCols];
    double*     f12tf12e= new double[nCols];
    double*     r12     = new double[nCols];
    double*     r12_hat = new double[nCols];
    double*     cross   = new double[nCols];
    double**    gg      = allocate2DArray(nCols,nCols);
    double**    tt2     = allocate2DArray(nCols,nCols);
    double**    quad1   = allocate2DArray(nCols,nCols);
    double**    quad2   = allocate2DArray(nCols,nCols);
    double**    ident   = allocate2DArray(nCols,nCols);
    double***   tt3     = allocate3DArray(nCols,nCols,nCols);
    double***   ggg     = allocate3DArray(nCols,nCols,nCols);
    double****  tt4     = allocate4DArray(nCols,nCols,nCols,nCols);    
    double***** tt5     = allocate5DArray(nCols,nCols,nCols,nCols,nCols);

    std::cout << '\n';
    std::cout <<"T-tensor" << '\n';
    std::cout <<"Calculation of electrostatic interactions between linear molecules" << '\n';
    std::cout <<"using T-tensors and Euler angles" << '\n';
    std::cout << '\n';
    
    // Write out parameters
    printf( "Min separation d_min            %15.6f \n", d_min); 
    printf( "Max separation d_max            %15.6f \n", d_max); 
    printf( "Dipole moment of molecule       %15.6f \n", mu1_mag); 
    printf( "Dipole moment of molecule       %15.6f \n", mu2_mag); 
    printf( "Quadrupole moment of molecule   %15.6f \n", quad1_mag); 
    printf( "Quadrupole moment of molecule   %15.6f \n", quad2_mag); 
    std::cout << '\n';
    
    random_vector(e1);
    random_vector(e2);
 
    // Place atom 2 at origin and atom 1 in a random direction within desired distance range
    std::cout << " --- r12_hat --- ";
    random_vector(r12_hat);
    print1DArray(nCols,r12_hat);
    r12_mag = (double) rand()/RAND_MAX;
    // printf("%10.10f \n ",r12_mag);
    r12_mag = d_min + (d_max-d_min)*r12_mag;              // Magnitude of r12
    // std::cout <<r12_mag << '\n';
    scalar1DArrayMultip(nCols,r12_mag,r12_hat,r12);       // Within desired range of origin

    c1  = dotProduct2D(nCols,e1,r12_hat);                 // Cosine of angle between e1 and r12
//    printf("c1 %10.6f \n",c1);
    c2  = dotProduct2D(nCols,e2,r12_hat);                 // Cosine of angle between e2 and r12
//    printf("c2 %10.6f \n",c2);
    c12 = dotProduct2D(nCols,e1,e2);                      // Cosine of angle between e1 and e2
//    printf("c12 %10.6f \n",c12);
 
    printf("%s %25.6f %10.6f %10.6f\n", "Displacement   r12", r12[0],r12[1],r12[2]);
    printf("%s %25.6f %10.6f %10.6f\n", "Orientation     e1",  e1[0], e1[1], e1[2]);
    printf("%s %25.6f %10.6f %10.6f\n", "Orientation     e2",  e2[0], e2[1], e2[2]);
    std::cout << '\n';
    
    // Dipole vectors in space-fixed frame
    scalar1DArrayMultip(nCols,mu1_mag,e1,mu1);
    scalar1DArrayMultip(nCols,mu2_mag,e2,mu2);

    // Quadrupole tensors in space-fixed frame (traceless)
    outer2D(nCols,e1,e1,quad1);
    scalar2DArrayMultip(nCols,nCols,1.5,quad1,quad1);
//    print2DArray(nCols,nCols,quad1);
    identMatrix(nCols,nCols,ident);
    scalar2DArrayMultip(nCols,nCols,0.5,ident,ident);
    subtract2DArrays(nCols,nCols,quad1,ident,quad1);
//    print2DArray(nCols,nCols,quad1);
    scalar2DArrayMultip(nCols,nCols,quad1_mag,quad1,quad1);
//    std::cout << " --- quad1 --- \n";
//    print2DArray(nCols,nCols,quad1);

    outer2D(nCols,e2,e2,quad2);
    scalar2DArrayMultip(nCols,nCols,1.5,quad2,quad2);
//    print2DArray(nCols,nCols,quad2);
    identMatrix(nCols,nCols,ident);
    scalar2DArrayMultip(nCols,nCols,0.5,ident,ident);
    subtract2DArrays(nCols,nCols,quad2,ident,quad2);
    scalar2DArrayMultip(nCols,nCols,quad2_mag,quad2,quad2);
//    std::cout << " --- quad2 ---\n";
//    print2DArray(nCols,nCols,quad2);

    // The T tensors of each rank: T2, T3, T4, T5
    t2_tensor (r12_hat, pow(r12_mag,3),tt2); 
    t3_tensor (r12_hat, pow(r12_mag,4),tt3);
    t4_tensor (r12_hat, pow(r12_mag,5),tt4);
    t5_tensor (r12_hat, pow(r12_mag,6),tt5);


    // Heading
    printf("%66s %40s %40s \n", ".....Result from T tensor", ".....Result from Euler angles",".........Difference");
    std::cout << '\n';

    std::cout << "Dipole-dipole" << '\n';

    // Calculate the dipole-dipole energy
    contract_ij_j(tt2, mu2,g);                                         // Contract T2 with dipole 2
    v12t = -1 * contract_i_i (nCols, mu1, g );                              // Contract result with dipole 1
    v12e = (mu1_mag*mu2_mag/pow(r12_mag,3)) * ( c12 - 3.0 * c1 * c2 );      // Compare result from angles
    printf ("Energy %59.6f %40.6f %40.6f \n", v12t, v12e, v12t-v12e);       

    // Calculate the dipole-dipole force
    contract_ijk_k(tt3, mu2,gg);                                         // Contract T3 with dipole 2
    contract_ij_j (gg , mu1,f12t);                                        // Contract result with dipole 1
    scalar1DArrayMultip(nCols,-1,f12t,f12t);

    double* f12e1 = new double[nCols];
    scalar1DArrayMultip(nCols,(c12-5.0*c1*c2),r12_hat,f12e);
    scalar1DArrayMultip(nCols,c2,e1,f12e1);
    sum1DArrays(nCols,f12e1,f12e,f12e);
    scalar1DArrayMultip(nCols,c1,e2,f12e1);
    sum1DArrays(nCols,f12e,f12e1,f12e);
    scalar1DArrayMultip(nCols,(3.0*mu1_mag*mu2_mag/pow(r12_mag,4)),f12e,f12e);   // Compare result from angles
    subtract1DArrays(nCols, f12t, f12e,f12tf12e);
    printf("Force  %37.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f\n", f12t[0],f12t[1],f12t[2], f12e[0], f12e[1],f12e[2], f12tf12e[0],f12tf12e[1],f12tf12e[2]);

    delete [] f12e1;

    // Calculate the dipole-dipole torques
    zeroVec(3,g);
    contract_ij_j(tt2, mu2,g);               // Contract T2 with dipole 2
    scalar1DArrayMultip(nCols,-1,g,g);
    crossProduct(mu1, g,t1t);                 // Cross-product result with dipole 1
    scalar1DArrayMultip(nCols,-1,t1t,t1t);
    scalar1DArrayMultip(nCols,3.*c2,r12_hat,g);
    subtract1DArrays(nCols,e2,g,g);         // Compare result from angles
    crossProduct(e1, g,cross);
    scalar1DArrayMultip(nCols,-(mu1_mag*mu2_mag/pow(r12_mag,3)),cross,t1e);
    subtract1DArrays(nCols, t1t, t1e,t1tt1e);
    printf("Torque on 1 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t1t[0], t1t[1],t1t[2], t1e[0], t1e[1],t1e[2], t1tt1e[0], t1tt1e[1],t1tt1e[2]);
   
    zeroVec(3,g);
    contract_ij_j ( tt2, mu1,g );              // Contract T2 with dipole 1
    scalar1DArrayMultip(nCols,-1,g,g);
    crossProduct ( mu2, g,t2t );                // Cross-product result with dipole 2
    scalar1DArrayMultip(nCols,-1,t2t,t2t);

    scalar1DArrayMultip(nCols,3.*c1,r12_hat,g);
    subtract1DArrays(nCols,e1,g,g);            // Compare result from angles
    crossProduct(e2, g,cross);
    scalar1DArrayMultip(nCols,-(mu1_mag*mu2_mag/pow(r12_mag,3)),cross,t2e);
    subtract1DArrays(nCols, t1t, t1e,t2tt2e);
    printf("Torque on 2 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t2t[0], t2t[1],t2t[2], t2e[0], t2e[1],t2e[2], t2tt2e[0], t2tt2e[1],t2tt2e[2]);

    std::cout << '\n' ;
    std::cout << "Dipole-quadrupole" << '\n';

    // Calculate the dipole-quadrupole energy
    zeroVec(3,g);
    contract_ijk_jk(tt3, quad2,g);                                               // Contract T3 with quadrupole 2
    v12t = -(1.0/3.0) *contract_i_i (nCols,mu1, g);                                   // Contract result with dipole 1
    v12e = (1.5*mu1_mag*quad2_mag/pow(r12_mag,4)) * ( c1*(1.0-5.0*c2*c2) + 2.0*c2*c12 );
    printf ("Energy %59.6f %40.6f %40.2f \n", v12t, v12e, v12t-v12e);       

    // Calculate the dipole-quadrupole force
    zeroMatrix(3,3,gg);
    contract_ijkl_kl (tt4, quad2,gg);                            // Contract T4 with quadrupole 2
    zeroVec(3,f12t);
    contract_ij_j ( gg, mu1,f12t);
    scalar1DArrayMultip(nCols,-(1.0/3.0),f12t,f12t);             // Contract result with dipole 1
    scalar1DArrayMultip(nCols,(35.0*c1*pow(c2,2)-10.0*c2*c12-5.0*c1),r12_hat,f12e);
    scalar1DArrayMultip(nCols,(1.0 - 5.0*pow(c2,2)),e1,h1);
    scalar1DArrayMultip(nCols,(2.0*c12-10.0*c1*c2),e2,h2);
    sum1DArrays(nCols,h1,f12e,f12e);
    sum1DArrays(nCols,h2,f12e,f12e);                                                  // Compare result from angles
    scalar1DArrayMultip(nCols,-(1.5*mu1_mag*quad2_mag/pow(r12_mag,5)),f12e,f12e);     // Compare result from angles
//    print1DArray(nCols,f12e);
    subtract1DArrays(nCols, f12t, f12e,f12tf12e);
    printf("Force  %37.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n", f12t[0],f12t[1],f12t[2], f12e[0], f12e[1],f12e[2], f12tf12e[0],f12tf12e[1],f12tf12e[2]);

    // Calculate the dipole-quadrupole torques
    double* g1 = new double[nCols];
    zeroVec(3,g);
    contract_ijk_jk( tt3, quad2,g);
    scalar1DArrayMultip(nCols,-(1.0/3.0),g,g);                  // Contract T3 with quadrupole 2
    crossProduct (mu1, g,t1t);
    scalar1DArrayMultip(nCols,-1.,t1t,t1t);                   // Cross-product result with dipole 1
    scalar1DArrayMultip(nCols,(1.0-5.0*pow(c2,2)),r12_hat,g);
    scalar1DArrayMultip(nCols, 2.0*c2, e2,g1);
    sum1DArrays(nCols,g,g1,g);                                  // Compare result from angles
    crossProduct(e1, g,cross);
    scalar1DArrayMultip(nCols,-(1.5*mu1_mag*quad2_mag/pow(r12_mag,4)),cross,t1e);
    subtract1DArrays(nCols, t1t, t1e,t1tt1e);
    printf("Torque on 1 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t1t[0], t1t[1],t1t[2], t1e[0], t1e[1],t1e[2], t1tt1e[0], t1tt1e[1],t1tt1e[2]);

    double* skg = new double[nCols];
    double** gg1 = allocate2DArray(nCols,nCols);
    zeroMatrix(3,3,gg);
    contract_ijk_k(tt3, mu1,gg);
    scalar2DArrayMultip(nCols,nCols,-(1.0/3.0),gg,gg);               // Contract T3 with dipole 1
    contract_ik_jk (quad2,gg,gg1);                                    // Contract result with quadrupole 2
    skew(gg1,skg);
    scalar1DArrayMultip(nCols,-2.0,skg,t2t);                         // Contract with Levi-Civita symbol
    scalar1DArrayMultip(nCols,(c12-5.0*c1*c2),r12_hat,g);
    scalar1DArrayMultip(nCols,c2, e1,g1);
    sum1DArrays(nCols,g,g1,g);                                        // Compare result from angles
    crossProduct(e2,g,cross);
    scalar1DArrayMultip(nCols,-(3.0*mu1_mag*quad2_mag/pow(r12_mag,4)),cross,t2e);
    subtract1DArrays(nCols, t1t, t1e,t2tt2e);
    printf("Torque on 2 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t2t[0], t2t[1],t2t[2], t2e[0], t2e[1],t2e[2], t2tt2e[0], t2tt2e[1],t2tt2e[2]);
    delete [] g1;

    std::cout << '\n' ;
    std::cout << "Quadrupole-dipole" << '\n';

    // Calculate the quadrupole-dipole energy
    
    zeroVec(3,g);
    contract_ijk_jk(tt3, quad1,g);                                         // Contract T3 with quadrupole 1
    v12t = (1.0/3.0) *contract_i_i (nCols,g,mu2);                                     // Contract result with dipole 2
    v12e = -(1.5*quad1_mag*mu2_mag/pow(r12_mag,4)) * (c2*(1.0-5.0*c1*c1) + 2.0*c1*c12);
    printf ("Energy %59.6f %40.6f %40.2f \n", v12t, v12e, v12t-v12e); 

    // Calculate the quadrupole-dipole force
    zeroMatrix(3,3,gg);
    contract_ijkl_kl(tt4, quad1,gg);                             // Contract T4 with quadrupole 1
    zeroVec(3,g1);
    contract_ij_j(gg, mu2,g1);
    scalar1DArrayMultip(nCols,(1.0/3.0),g1,f12t);                // Contract result with dipole 2
    scalar1DArrayMultip(nCols,(35.0*c2*pow(c1,2)-10.0*c1*c12-5.0*c2),r12_hat,f12e);
    scalar1DArrayMultip(nCols,(1.0 - 5.0*pow(c1,2)),e2,h1);
    scalar1DArrayMultip(nCols,(2.0*c12-10.0*c1*c2),e1,h2);
    sum1DArrays(nCols,h1,f12e,f12e);
    sum1DArrays(nCols,h2,f12e,f12e);    
    scalar1DArrayMultip(nCols,(1.5*mu1_mag*quad2_mag/pow(r12_mag,5)),f12e,f12e);    // Compare result from angles
    subtract1DArrays(nCols, f12t, f12e,f12tf12e);
    printf("Force  %37.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n", f12t[0],f12t[1],f12t[2], f12e[0], f12e[1],f12e[2], f12tf12e[0],f12tf12e[1],f12tf12e[2]);

    // Calculate the quadrupole-dipole torques
    zeroMatrix(3,3,gg);
    contract_ijk_k(tt3, mu2,gg);
    scalar2DArrayMultip(nCols,nCols,(1.0/3.0),gg,gg);         // Contract T3 with dipole 2
    zeroMatrix(3,3,gg1);
    contract_ik_jk(quad1,gg,gg1);                     // Contract result with quadrupole 1
    skew(gg1,skg);
    scalar1DArrayMultip(nCols,-2.0,skg, t1t);                                   // Contract with Levi-Civita symbol
    scalar1DArrayMultip(nCols,(c12-5.0*c1*c2),r12_hat,g);
    scalar1DArrayMultip(nCols,c1,e2,g1);
    sum1DArrays(nCols,g,g1,g);                       // Compare result from angles
    crossProduct(e1,g,cross);
    scalar1DArrayMultip(nCols,3.0*quad1_mag*mu2_mag/pow(r12_mag,4),cross,t1e);
    subtract1DArrays(nCols, t1t, t1e,t1tt1e);
    printf("Torque on 1 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t1t[0], t1t[1],t1t[2], t1e[0], t1e[1],t1e[2], t1tt1e[0], t1tt1e[1],t1tt1e[2]);

    zeroVec(3,g);
    contract_ijk_jk(tt3,quad1,g);
    scalar1DArrayMultip(nCols,(1.0/3.0),g,g);             // Contract T3 with quadrupole 1
    crossProduct(mu2,g,cross);
    scalar1DArrayMultip(nCols,-1,cross,t2t);            // Cross product result with dipole 2
    scalar1DArrayMultip(nCols,(1.0-5.0*pow(c1,2)),r12_hat,g);
    scalar1DArrayMultip(nCols,2.0*c1,e1,g1);
    sum1DArrays(nCols,g1,g,g);                          // Compare result from angles
    crossProduct(e2,g,cross);
    scalar1DArrayMultip(nCols,(1.5*quad1_mag*mu2_mag/pow(r12_mag,4)),cross,t2e);
    subtract1DArrays(nCols, t2t, t2e,t2tt2e);
    printf("Torque on 2 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t2t[0], t2t[1],t2t[2], t2e[0], t2e[1],t2e[2], t2tt2e[0], t2tt2e[1],t2tt2e[2]);

    std::cout << '\n' ;
    std::cout << "Quadrupole-quadrupole" << '\n';

    // Calculate the quadrupole-quadrupole energy
    zeroMatrix(3,3,gg);
    contract_ijkl_kl(tt4,quad2,gg);                                // Contract T4 with quadrupole 2
    v12t = (1.0/9.0) * contract_ij_ij (quad1,gg);                  // Contract result with quadrupole 1
    v12e = (0.75*quad1_mag*quad2_mag/pow(r12_mag,5)) * (1.0 - 5.0*pow(c1,2) - 5.0*pow(c2,2) + 2.0*pow(c12,2) + 35.0*pow((c1*c2),2) - 20.0*c1*c2*c12 ); // Compare result from angles
    printf ("Energy %59.6f %40.6f %40.2f \n", v12t, v12e, v12t-v12e); 

    // Calculate the quadrupole-quadrupole force
    contract_ijklm_lm(tt5,quad2,ggg);                           // Contract T5 with quadrupole 2
    zeroVec(3,g);
    contract_ijk_jk(ggg,quad1,g);
    scalar1DArrayMultip(nCols,(1.0/9.0),g,f12t);               // Contract result with quadrupole 1
    scalar1DArrayMultip(nCols,(5.0 - 35.0*pow(c1,2) - 35.0*pow(c2,2) + 10.0*pow(c12,2) + 315.0*pow((c1*c2),2) - 140.0*c1*c2*c12),r12_hat,f12e);
    scalar1DArrayMultip(nCols,(10.0*c1 - 70.0*c1*pow(c2,2) + 20.0*c2*c12),e1,h1);
    scalar1DArrayMultip(nCols,(10.0*c2 - 70.0*c2*pow(c1,2) + 20.0*c1*c12),e2,h2);
    sum1DArrays(nCols, h1,f12e,f12e);
    sum1DArrays(nCols, h2,f12e,f12e);
    scalar1DArrayMultip(nCols,(0.75*quad1_mag*quad2_mag/pow(r12_mag,6)),f12e,f12e);      // Compare result from angles
    subtract1DArrays(nCols, f12t, f12e,f12tf12e);
    printf("Force  %37.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n", f12t[0],f12t[1],f12t[2], f12e[0], f12e[1],f12e[2], f12tf12e[0],f12tf12e[1],f12tf12e[2]);


    // Calculate the quadrupole-quadrupole torques
    zeroMatrix(3,3,gg);
    contract_ijkl_kl(tt4,quad2,gg);
    scalar2DArrayMultip(nCols,nCols,(1.0/9.0),gg,gg);          // Contract T4 with quadrupole 2
    zeroMatrix(3,3,gg1);
    contract_ik_jk(quad1, gg,gg1);                              // Contract result with quadrupole 1
    skew(gg1,skg);
    scalar1DArrayMultip(nCols,-2.0,skg,t1t);                   // Contract with Levi-Civita symbol
    scalar1DArrayMultip(nCols,(2.5*(c1*(7.0*pow(c2,2)-1.0)-2.0*c2*c12)), r12_hat,g);
    scalar1DArrayMultip(nCols,(5.0*c1*c2-c12), e2,g1);
    subtract1DArrays(nCols,g,g1,g);                            // Compare result from angles
    crossProduct(e1,g,cross);
    scalar1DArrayMultip(nCols,-(3.0*quad1_mag*quad2_mag/pow(r12_mag,5)),cross,t1e );
    subtract1DArrays(nCols, t1t, t1e,t1tt1e);
    printf("Torque on 1 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t1t[0], t1t[1],t1t[2], t1e[0], t1e[1],t1e[2], t1tt1e[0], t1tt1e[1],t1tt1e[2]);

    zeroMatrix(3,3,gg);
    contract_ijkl_kl(tt4,quad1,gg);
    scalar2DArrayMultip(nCols,nCols,(1.0/9.0),gg,gg);          // Contract T4 with quadrupole 1
    zeroMatrix(3,3,gg1);
    contract_ik_jk(quad2, gg,gg1);                              // Contract result with quadrupole 2
    skew(gg1,skg);
    scalar1DArrayMultip(nCols,-2.0,skg,t2t);                   // Contract with Levi-Civita symbol
    scalar1DArrayMultip(nCols,2.5*(c2*(7.0*pow(c1,2)-1.0)-2.0*c1*c12),r12_hat, g);
    scalar1DArrayMultip(nCols,(5.0*c1*c2-c12),e1,g1);
    subtract1DArrays(nCols,g,g1,g);                           // Compare result from angles
    crossProduct(e2,g,cross);
    scalar1DArrayMultip(nCols,-(3.0*quad1_mag*quad2_mag/pow(r12_mag,5)),cross,t2e);
    subtract1DArrays(nCols, t2t, t2e, t2tt2e);
    printf("Torque on 2 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t2t[0], t2t[1],t2t[2], t2e[0], t2e[1],t2e[2], t2tt2e[0], t2tt2e[1],t2tt2e[2]);

    delete [] g;
    delete [] h1;
    delete [] h2;
    delete [] e1;
    delete [] e2;
    delete [] mu1;
    delete [] mu2;
    delete [] t1t;
    delete [] t1tt1e;
    delete [] t2tt2e;
    delete [] t2t;
    delete [] t1e;
    delete [] t2e;
    delete [] f12t;
    delete [] f12e;
    delete [] f12tf12e;
    delete [] r12;
    delete [] r12_hat;
    delete [] cross;
    free2DArray(nCols,gg);   
    free2DArray(nCols,tt2);
    free2DArray(nCols,quad1);
    free2DArray(nCols,quad2);
    free2DArray(nCols,ident);
    free3DArray(nCols,nCols,tt3); 
    free3DArray(nCols,nCols,ggg);
    free4DArray(nCols,nCols,nCols,tt4); 
    free5DArray(nCols,nCols,nCols,nCols,tt5);

    return 0;

}
