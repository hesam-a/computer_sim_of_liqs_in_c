This repository includes some C++ version (written by H.A) of some programs from the book "Computer simulation of liquids"
by Allen and Tildesley.

So far it is written for the t_tensor program (chapter one) that is used to calculate the electrostatic interactions between two linear atoms, NVT and NPT (Metropolis) Monte Carlo programs (chapter four) for Leonard-Jones particles.

The cmake file for compilation of the codes will be provided soon.

Currently, the mc_nvt_lj.cpp program can be compiled with this command:

g++ -g mc_nvt_lj.cpp math_module.cpp mc_lj_module.cpp lrc_module.cpp averages_module.cpp config_io_module.cpp -o mc_nvt

and run with command:

./mc_nvt

Again, the mc_npt_lj.cpp program can be compiled with this command:

g++ -g mc_npt_lj.cpp math_module.cpp mc_lj_module.cpp lrc_module.cpp averages_module.cpp config_io_module.cpp -o mc_npt

and run with command:

./mc_npt

In case the boost library is not available in the user's c header library, it can be installed with this command:

sudo apt-get install libboost-all-dev

that needs root privilege.
