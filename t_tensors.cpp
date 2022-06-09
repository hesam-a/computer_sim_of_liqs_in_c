
// This program is a C++ version of the "t_tensor.f90" program from the book "Computer simulations of liquids"
// by Allen and Tildesley
// It is used for calculation of electrostatic interactions between two linear molecules including dipole-dipole,
// dipole-quadrupole and quadrupole-quadrupole, using T tensors and Euler angles.

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
#include <cmath>
#include <stdio.h>
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
    
    e1 = random_vector();
//    std::cout << "e1:" << '\n';
//    print1DArray(nCols,e1);
    e2 = random_vector();
//    std::cout << "e2:" << '\n';
//    print1DArray(nCols,e2);
 
    // Place atom 2 at origin and atom 1 in a random direction within desired distance range
    r12_hat = random_vector();
    //print1DArray(nCols,r12_hat);
    r12_mag = (double) rand()/RAND_MAX;
    // printf("%10.10f \n ",r12_mag);
    r12_mag = d_min + (d_max-d_min)*r12_mag;                      // Magnitude of r12
    // std::cout <<r12_mag << '\n';
    r12     = scalar1DArrayMultip(nCols,r12_mag,r12_hat);   // Within desired range of origin

    c1  = dotProduct(nCols,e1,r12_hat);                 // Cosine of angle between e1 and r12
//    printf("c1 %10.6f \n",c1);
    c2  = dotProduct(nCols,e2,r12_hat);                 // Cosine of angle between e2 and r12
//    printf("c2 %10.6f \n",c2);
    c12 = dotProduct(nCols,e1,e2);                      // Cosine of angle between e1 and e2
//    printf("c12 %10.6f \n",c12);
 
    printf("%s %25.6f %10.6f %10.6f\n", "Displacement   r12", r12[0],r12[1],r12[2]);
    printf("%s %25.6f %10.6f %10.6f\n", "Orientation     e1",  e1[0], e1[1], e1[2]);
    printf("%s %25.6f %10.6f %10.6f\n", "Orientation     e2",  e2[0], e2[1], e2[2]);
    std::cout << '\n';
    
    // Dipole vectors in space-fixed frame
    mu1 = scalar1DArrayMultip(nCols,mu1_mag,e1);
//    std::cout << "mu1: " << '\n';
//    print1DArray(nCols,mu1);
    mu2 = scalar1DArrayMultip(nCols,mu2_mag,e2);
//    std::cout << "mu2: " << '\n';
//    print1DArray(nCols,mu2);

    // Quadrupole tensors in space-fixed frame (traceless)

    quad1 = scalar2DArrayMultip(nCols,nCols,1.5,outer2D(nCols,e1,e1));
//    print2DArray(nCols,nCols,quad1);
    ident = identMatrix(nCols,nCols);
    ident = scalar2DArrayMultip(nCols,nCols,0.5,ident);
    quad1 = subtract2DArrays(nCols,nCols,quad1,ident);
//    print2DArray(nCols,nCols,quad1);
    quad1 = scalar2DArrayMultip(nCols,nCols,quad1_mag,quad1);
//    print2DArray(nCols,nCols,quad1);

    quad2 = scalar2DArrayMultip(nCols,nCols,1.5,outer2D(nCols,e2,e2));
//    print2DArray(nCols,nCols,quad2);
    ident = identMatrix(nCols,nCols);
    ident = scalar2DArrayMultip(nCols,nCols,0.5,ident);
    quad2 = subtract2DArrays(nCols,nCols,quad2,ident);
    quad2 = scalar2DArrayMultip(nCols,nCols,quad2_mag,quad2);
//    print2DArray(nCols,nCols,quad2);

    // The T tensors of each rank: T2, T3, T4, T5
    tt2 = t2_tensor(r12_hat, pow(r12_mag,3)); 
//    std::cout << "tt2 " << '\n';
//    print2DArray(nCols,nCols,tt2);
    tt3 = t3_tensor (r12_hat, pow(r12_mag,4));
    tt4 = t4_tensor (r12_hat, pow(r12_mag,5));
    tt5 = t5_tensor (r12_hat, pow(r12_mag,6));


    // Heading
    printf("%66s %40s %40s \n", ".....Result from T tensor", ".....Result from Euler angles",".........Difference");
    std::cout << '\n';

    std::cout << "Dipole-dipole" << '\n';

    // Calculate the dipole-dipole energy
    g    = contract_ij_j(tt2, mu2);                                         // Contract T2 with dipole 2
//    print1DArray(nCols,g);
    v12t = -1 * contract_i_i (nCols, mu1, g );                              // Contract result with dipole 1
    v12e = (mu1_mag*mu2_mag/pow(r12_mag,3)) * ( c12 - 3.0 * c1 * c2 );      // Compare result from angles
    printf ("Energy %59.6f %40.6f %40.6f \n", v12t, v12e, v12t-v12e);       

    // Calculate the dipole-dipole force
    gg       = contract_ijk_k(tt3, mu2);                                        // Contract T3 with dipole 2
    f12t     = contract_ij_j (gg , mu1);                                        // Contract result with dipole 1
    f12t     = scalar1DArrayMultip(nCols,-1,f12t);

    f12e     = sum1DArrays(nCols,scalar1DArrayMultip(nCols,(c12-5.0*c1*c2),r12_hat),scalar1DArrayMultip(nCols,c2,e1)) ;
    f12e     = sum1DArrays(nCols,f12e,scalar1DArrayMultip(nCols,c1,e2));
    f12e     = scalar1DArrayMultip(nCols,(3.0*mu1_mag*mu2_mag/pow(r12_mag,4)),f12e);   // Compare result from angles
    f12tf12e = subtract1DArrays(nCols, f12t, f12e);
    printf("Force  %37.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f\n", f12t[0],f12t[1],f12t[2], f12e[0], f12e[1],f12e[2], f12tf12e[0],f12tf12e[1],f12tf12e[2]);

    // Calculate the dipole-dipole torques
    g      = contract_ij_j(tt2, mu2);               // Contract T2 with dipole 2
    g      = scalar1DArrayMultip(nCols,-1,g);
//    print1DArray(nCols,g);
    t1t    = crossProduct(mu1, g);                 // Cross-product result with dipole 1
    t1t    = scalar1DArrayMultip(nCols,-1,t1t);
//    print1DArray(nCols,t1t);
    g      = scalar1DArrayMultip(nCols,3.*c2,r12_hat);
    g      = subtract1DArrays(nCols,e2,g);         // Compare result from angles
    t1e    = scalar1DArrayMultip(nCols,-(mu1_mag*mu2_mag/pow(r12_mag,3)), crossProduct(e1, g));
    t1tt1e = subtract1DArrays(nCols, t1t, t1e);
    printf("Torque on 1 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t1t[0], t1t[1],t1t[2], t1e[0], t1e[1],t1e[2], t1tt1e[0], t1tt1e[1],t1tt1e[2]);
   
    g      = contract_ij_j ( tt2, mu1 );              // Contract T2 with dipole 1
    g      = scalar1DArrayMultip(nCols,-1,g);
    t2t    = crossProduct ( mu2, g );                // Cross-product result with dipole 2
    t2t    = scalar1DArrayMultip(nCols,-1,t2t);

    g      = scalar1DArrayMultip(nCols,3.*c1,r12_hat);
    g      = subtract1DArrays(nCols,e1,g);            // Compare result from angles
    t2e    = scalar1DArrayMultip(nCols,-(mu1_mag*mu2_mag/pow(r12_mag,3)), crossProduct(e2, g));
    t2tt2e = subtract1DArrays(nCols, t1t, t1e);
    printf("Torque on 2 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t2t[0], t2t[1],t2t[2], t2e[0], t2e[1],t2e[2], t2tt2e[0], t2tt2e[1],t2tt2e[2]);

    std::cout << '\n' ;
    std::cout << "Dipole-quadrupole" << '\n';

    // Calculate the dipole-quadrupole energy
    g    = contract_ijk_jk(tt3, quad2);                                               // Contract T3 with quadrupole 2
    v12t = -(1.0/3.0) *contract_i_i (nCols,mu1, g);                                   // Contract result with dipole 1
    v12e = (1.5*mu1_mag*quad2_mag/pow(r12_mag,4)) * ( c1*(1.0-5.0*c2*c2) + 2.0*c2*c12 );
    printf ("Energy %59.6f %40.6f %40.2f \n", v12t, v12e, v12t-v12e);       

    // Calculate the dipole-quadrupole force
    gg   = contract_ijkl_kl (tt4, quad2);                                              // Contract T4 with quadrupole 2
    f12t = scalar1DArrayMultip(nCols,-(1.0/3.0),contract_ij_j ( gg, mu1));             // Contract result with dipole 1
    f12e = scalar1DArrayMultip(nCols,(35.0*c1*pow(c2,2)-10.0*c2*c12-5.0*c1),r12_hat);
    h1   = scalar1DArrayMultip(nCols,(1.0 - 5.0*pow(c2,2)),e1);
    h2   = scalar1DArrayMultip(nCols,(2.0*c12-10.0*c1*c2),e2);
    f12e = sum1DArrays(nCols,h1,f12e);
    f12e = sum1DArrays(nCols,h2,f12e);                                            // Compare result from angles
    f12e = scalar1DArrayMultip(nCols,-(1.5*mu1_mag*quad2_mag/pow(r12_mag,5)),f12e);     // Compare result from angles
//    print1DArray(nCols,f12e);
    f12tf12e = subtract1DArrays(nCols, f12t, f12e);
    printf("Force  %37.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n", f12t[0],f12t[1],f12t[2], f12e[0], f12e[1],f12e[2], f12tf12e[0],f12tf12e[1],f12tf12e[2]);

    // Calculate the dipole-quadrupole torques
    g   = scalar1DArrayMultip(nCols,-(1.0/3.0),contract_ijk_jk( tt3, quad2));           // Contract T3 with quadrupole 2
    t1t = scalar1DArrayMultip(nCols,-1.,crossProduct ( mu1, g ));                       // Cross-product result with dipole 1
    g   = scalar1DArrayMultip(nCols,(1.0-5.0*pow(c2,2)),r12_hat);
    g   = sum1DArrays(nCols,g,scalar1DArrayMultip(nCols, 2.0*c2, e2));            // Compare result from angles
    t1e = scalar1DArrayMultip(nCols,-(1.5*mu1_mag*quad2_mag/pow(r12_mag,4)),crossProduct(e1, g ));
    t1tt1e = subtract1DArrays(nCols, t1t, t1e);
    printf("Torque on 1 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t1t[0], t1t[1],t1t[2], t1e[0], t1e[1],t1e[2], t1tt1e[0], t1tt1e[1],t1tt1e[2]);

    gg  = scalar2DArrayMultip(nCols,nCols,-(1.0/3.0),contract_ijk_k(tt3, mu1));               // Contract T3 with dipole 1
    gg  = contract_ik_jk ( quad2, gg );                                                 // Contract result with quadrupole 2
    t2t = scalar1DArrayMultip(nCols,-2.0,skew ( gg ));                                                            // Contract with Levi-Civita symbol
    g   = scalar1DArrayMultip(nCols,(c12-5.0*c1*c2),r12_hat);
    g   = sum1DArrays(nCols,g,scalar1DArrayMultip(nCols,c2, e1));                                        // Compare result from angles
    t2e = scalar1DArrayMultip(nCols,-(3.0*mu1_mag*quad2_mag/pow(r12_mag,4)), crossProduct(e2,g));
    t2tt2e = subtract1DArrays(nCols, t1t, t1e);
    printf("Torque on 2 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t2t[0], t2t[1],t2t[2], t2e[0], t2e[1],t2e[2], t2tt2e[0], t2tt2e[1],t2tt2e[2]);

    std::cout << '\n' ;
    std::cout << "Quadrupole-dipole" << '\n';

    // Calculate the quadrupole-dipole energy
    
    g    = contract_ijk_jk(tt3, quad1);                                               // Contract T3 with quadrupole 1
    v12t = (1.0/3.0) *contract_i_i (nCols,g,mu2);                                     // Contract result with dipole 2
    v12e = -(1.5*quad1_mag*mu2_mag/pow(r12_mag,4)) * (c2*(1.0-5.0*c1*c1) + 2.0*c1*c12);
    printf ("Energy %59.6f %40.6f %40.2f \n", v12t, v12e, v12t-v12e); 

    // Calculate the quadrupole-dipole force
    gg   = contract_ijkl_kl(tt4, quad1);                                               // Contract T4 with quadrupole 1
    f12t = scalar1DArrayMultip(nCols,(1.0/3.0),contract_ij_j(gg, mu2));                // Contract result with dipole 2
    f12e = scalar1DArrayMultip(nCols,(35.0*c2*pow(c1,2)-10.0*c1*c12-5.0*c2),r12_hat);
    h1   = scalar1DArrayMultip(nCols,(1.0 - 5.0*pow(c1,2)),e2);
    h2   = scalar1DArrayMultip(nCols,(2.0*c12-10.0*c1*c2),e1);
    f12e = sum1DArrays(nCols,h1,f12e);
    f12e = sum1DArrays(nCols,h2,f12e);    
    f12e = scalar1DArrayMultip(nCols,(1.5*mu1_mag*quad2_mag/pow(r12_mag,5)),f12e);    // Compare result from angles
    f12tf12e = subtract1DArrays(nCols, f12t, f12e);
    printf("Force  %37.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n", f12t[0],f12t[1],f12t[2], f12e[0], f12e[1],f12e[2], f12tf12e[0],f12tf12e[1],f12tf12e[2]);

    // Calculate the quadrupole-dipole torques
    gg  = scalar2DArrayMultip(nCols,nCols,(1.0/3.0),contract_ijk_k(tt3, mu2));         // Contract T3 with dipole 2
    gg  = contract_ik_jk(quad1,gg);                                                    // Contract result with quadrupole 1
    t1t = scalar1DArrayMultip(nCols,-2.0,skew (gg));                                   // Contract with Levi-Civita symbol
    g   = scalar1DArrayMultip(nCols,(c12-5.0*c1*c2),r12_hat);
    g   = sum1DArrays(nCols,g,scalar1DArrayMultip(nCols,c1,e2));                       // Compare result from angles
    t1e = scalar1DArrayMultip(nCols,3.0*quad1_mag*mu2_mag/pow(r12_mag,4),crossProduct(e1,g));
    t1tt1e = subtract1DArrays(nCols, t1t, t1e);
    printf("Torque on 1 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t1t[0], t1t[1],t1t[2], t1e[0], t1e[1],t1e[2], t1tt1e[0], t1tt1e[1],t1tt1e[2]);

    g   = scalar1DArrayMultip(nCols,(1.0/3.0),contract_ijk_jk(tt3,quad1));             // Contract T3 with quadrupole 1
    t2t = scalar1DArrayMultip(nCols,-1,crossProduct(mu2,g));                           // Cross product result with dipole 2
    g   = scalar1DArrayMultip(nCols,(1.0-5.0*pow(c1,2)),r12_hat);
    g   = sum1DArrays(nCols,scalar1DArrayMultip(nCols,2.0*c1,e1),g);                                               // Compare result from angles
    t2e = scalar1DArrayMultip(nCols,(1.5*quad1_mag*mu2_mag/pow(r12_mag,4)),crossProduct(e2,g));
    t2tt2e = subtract1DArrays(nCols, t2t, t2e);
    printf("Torque on 2 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t2t[0], t2t[1],t2t[2], t2e[0], t2e[1],t2e[2], t2tt2e[0], t2tt2e[1],t2tt2e[2]);

    std::cout << '\n' ;
    std::cout << "Quadrupole-quadrupole" << '\n';

    // Calculate the quadrupole-quadrupole energy
    gg   = contract_ijkl_kl(tt4,quad2);                                                    // Contract T4 with quadrupole 2
    v12t = (1.0/9.0) * contract_ij_ij (quad1,gg);                                          // Contract result with quadrupole 1
    v12e = (0.75*quad1_mag*quad2_mag/pow(r12_mag,5)) * (1.0 - 5.0*pow(c1,2) - 5.0*pow(c2,2) + 2.0*pow(c12,2) + 35.0*pow((c1*c2),2) - 20.0*c1*c2*c12 ); // Compare result from angles
    printf ("Energy %59.6f %40.6f %40.2f \n", v12t, v12e, v12t-v12e); 

    // Calculate the quadrupole-quadrupole force
    ggg  = contract_ijklm_lm(tt5,quad2);                                                   // Contract T5 with quadrupole 2
    f12t = scalar1DArrayMultip(nCols,(1.0/9.0), contract_ijk_jk(ggg,quad1));               // Contract result with quadrupole 1
    f12e = scalar1DArrayMultip(nCols,(5.0 - 35.0*pow(c1,2) - 35.0*pow(c2,2) + 10.0*pow(c12,2) + 315.0*pow((c1*c2),2) - 140.0*c1*c2*c12),r12_hat);
    h1   = scalar1DArrayMultip(nCols,(10.0*c1 - 70.0*c1*pow(c2,2) + 20.0*c2*c12),e1);
    h2   = scalar1DArrayMultip(nCols,(10.0*c2 - 70.0*c2*pow(c1,2) + 20.0*c1*c12), e2);
    f12e = sum1DArrays(nCols, h1,f12e);
    f12e = sum1DArrays(nCols, h2,f12e);
    f12e = scalar1DArrayMultip(nCols,(0.75*quad1_mag*quad2_mag/pow(r12_mag,6)),f12e);      // Compare result from angles
    f12tf12e = subtract1DArrays(nCols, f12t, f12e);
    printf("Force  %37.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n", f12t[0],f12t[1],f12t[2], f12e[0], f12e[1],f12e[2], f12tf12e[0],f12tf12e[1],f12tf12e[2]);


    // Calculate the quadrupole-quadrupole torques
    gg  = scalar2DArrayMultip(nCols,nCols,(1.0/9.0),contract_ijkl_kl(tt4,quad2));          // Contract T4 with quadrupole 2
    gg  = contract_ik_jk(quad1, gg);                                                       // Contract result with quadrupole 1
    t1t = scalar1DArrayMultip(nCols,-2.0,skew(gg));                                        // Contract with Levi-Civita symbol
    g   = scalar1DArrayMultip(nCols,(2.5*(c1*(7.0*pow(c2,2)-1.0)-2.0*c2*c12)), r12_hat);
    g   = subtract1DArrays(nCols,g,scalar1DArrayMultip(nCols,(5.0*c1*c2-c12), e2));        // Compare result from angles
    t1e = scalar1DArrayMultip(nCols,-(3.0*quad1_mag*quad2_mag/pow(r12_mag,5)), crossProduct(e1,g));
    t1tt1e = subtract1DArrays(nCols, t1t, t1e);
    printf("Torque on 1 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t1t[0], t1t[1],t1t[2], t1e[0], t1e[1],t1e[2], t1tt1e[0], t1tt1e[1],t1tt1e[2]);

    gg  = scalar2DArrayMultip(nCols,nCols,(1.0/9.0),contract_ijkl_kl(tt4,quad1));          // Contract T4 with quadrupole 1
    gg  = contract_ik_jk(quad2, gg);                                                       // Contract result with quadrupole 2
    t2t = scalar1DArrayMultip(nCols,-2.0,skew(gg));                                        // Contract with Levi-Civita symbol
    g   = scalar1DArrayMultip(nCols,2.5*(c2*(7.0*pow(c1,2)-1.0)-2.0*c1*c12),r12_hat);
    g   = subtract1DArrays(nCols,g,scalar1DArrayMultip(nCols,(5.0*c1*c2-c12),e1));            // Compare result from angles
    t2e = scalar1DArrayMultip(nCols,-(3.0*quad1_mag*quad2_mag/pow(r12_mag,5)),crossProduct(e2,g));
    t2tt2e = subtract1DArrays(nCols, t2t, t2e);
    printf("Torque on 2 %32.6f %10.6f %10.6f %18.6f %10.6f %10.6f %18.6f %10.6f %10.6f \n",  t2t[0], t2t[1],t2t[2], t2e[0], t2e[1],t2e[2], t2tt2e[0], t2tt2e[1],t2tt2e[2]);


    delete [] g,h1,h2, e1,e2, mu1,mu2,t1t,t2t,t1e,t2e,f12t,f12e,r12,r12_hat,f12tf12e;
    delete [] t1tt1e;
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
