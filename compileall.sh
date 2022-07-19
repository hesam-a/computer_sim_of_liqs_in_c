#! /bin/bash

g++ t_tensor.cpp    math_module.cpp -o t_tensor
g++ mc_nvt_lj.cpp   math_module.cpp config_io_module.cpp mc_lj_module.cpp averages_module.cpp lrc_module.cpp -o mc_nvt
g++ mc_npt_lj.cpp   math_module.cpp config_io_module.cpp mc_lj_module.cpp averages_module.cpp lrc_module.cpp -o mc_npt
g++ mc_zvt_lj.cpp   math_module.cpp config_io_module.cpp mc_lj_module.cpp averages_module.cpp lrc_module.cpp -o mc_zvt
g++ mc_nvt_poly.cpp math_module.cpp config_io_module.cpp mc_poly_lj_module.cpp averages_module.cpp           -o mc_poly

