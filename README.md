This repository includes some C++ version (written by H.A) of some programs from the book "Computer simulation of liquids"
by Allen and Tildesley. Professor Allen is aware of this repository and has stargazed it. 

H.A is writing the codes ONLY for the purpose of learning the physics of the topics and practicing c/c++. The codes are not computationally efficient (neither vectorization nor a vectorization-enabled package is used), since H.A is interested in learning the linear algebra of the topics (especially electrostatic interactions where higher-order tensors are involved) beside the physics.

So far it is written for the t_tensor program (chapter one) that is used to calculate the electrostatic interactions between two linear atoms, and NVT, NPT and zVT (Metropolis) Monte Carlo programs for atoms, NVT (Metropolis) Monte Carlo using quaternion for molecules (chapter four) and NVT brownian dynamics (chapter twelve), all for Leonard-Jones particles.

There are two ways provided for the compilation of the programs. First, via Makefile which compiles all libraries (generates *.o files) and compiles the main programs one by one (comment should be removed for the intended program). Second, via compileall.sh to compile them all at once without generating the *.o files.

In case the boost library is not available in the user's c header files library, it can be installed with this command:

sudo apt-get install libboost-all-dev

that needs root privilege.

"Any feedback is welcome!"
