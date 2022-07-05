This repository includes some C++ version (written by H.A) of some programs from the book "Computer simulation of liquids"
by Allen and Tildesley.

So far it is written for the t_tensor program (chapter one) that is used to calculate the electrostatic interactions between two linear atoms, NVT, NPT and zVT (Metropolis) Monte Carlo programs (chapter four) for Leonard-Jones particles.

The cmake file for compilation of the codes will be provided soon.

Currently, the programs, Monte Carlo NVT, NPT and zPT can be compiled with this command: 

g++ -g mc_nvt_lj.cpp math_module.cpp mc_lj_module.cpp lrc_module.cpp averages_module.cpp config_io_module.cpp -o mc_nvt

where "mc_nvt_lj.cpp" should be replaced with "mc_npt_lj.cpp" or "mc_zvt_lj.cpp", respectively for NPT or zVT codes.

In case the boost library is not available in the user's c header files library, it can be installed with this command:

sudo apt-get install libboost-all-dev

that needs root privilege.
