#! /bin/bash

g++ -std=c++14 t_tensor.cpp    maths_module.cpp -o t_tensor
g++ -std=c++14 mc_nvt_lj.cpp   maths_module.cpp config_io_module.cpp mc_lj_module.cpp averages_module.cpp lrc_module.cpp -o mc_nvt
g++ -std=c++14 mc_npt_lj.cpp   maths_module.cpp config_io_module.cpp mc_lj_module.cpp averages_module.cpp lrc_module.cpp -o mc_npt
g++ -std=c++14 mc_zvt_lj.cpp   maths_module.cpp config_io_module.cpp mc_lj_module.cpp averages_module.cpp lrc_module.cpp -o mc_zvt
g++ -std=c++14 mc_nvt_poly.cpp maths_module.cpp config_io_module.cpp mc_poly_lj_module.cpp averages_module.cpp           -o mc_poly
g++ -std=c++14 bd_nvt_lj.cpp   maths_module.cpp config_io_module.cpp md_lj_module.cpp averages_module.cpp  rc_module.cpp -o bd_nvt
