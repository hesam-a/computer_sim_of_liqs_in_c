#==================
#   Compiler
#==================

CC = g++

#==================
#   Libraries
#==================

ITINCLUDE3  =   #-I/... 
ITINCLUDE   = ${ITINCLUDE3}
ITLIB       = -L/usr/lib/ # directory where the required libraries are

#----------------------------
#   Compilation options
#----------------------------

INCLUDE_PATH = -I./ ${ITINCLUDE}
LDFLAGS      = ${ITLIB}
WARN         = -Wno-deprecated
#CPPFLAGS     = -ffast-math -funroll-loops ${WARN} ${INCLUDE_PATH}
DEFINES      = -Wall

.cpp.o: ; ${CC} ${DEFINES} ${INCLUDE_PATH} ${WARN} -g -c $*.cpp


#===========================
#   Files to be compiled
#===========================

TENSOR  = math_module.o
MC_LJ   = math_module.o config_io_module.o mc_lj_module.o   averages_module.o lrc_module.o
MC_POLY = math_module.o config_io_module.o mc_poly_lj_module.o averages_module.o 

#t_tensor:	${TENSOR}
#		${CC} -g t_tensor.cpp   ${TENSOR}   -o tensor ${LDFLAGS} ${LIBRARY}     

#mc_nvt_lj:	${MC_LJ}
#		${CC} -g mc_nvt_lj.cpp   ${MC_LJ}   -o mc_nvt ${LDFLAGS} ${LIBRARY}     

#mc_npt_lj:	${MC_LJ}
#		${CC} -g mc_npt_lj.cpp   ${MC_LJ}   -o mc_npt ${LDFLAGS} ${LIBRARY}     

#mc_zvt_lj:	${MC_LJ}
#		${CC} -g mc_zvt_lj.cpp   ${MC_LJ}   -o mc_zvt ${LDFLAGS} ${LIBRARY}     

mc_nvt_poly:	${MC_POLY}
		${CC} -g mc_nvt_poly.cpp ${MC_POLY} -o mc_poly ${LDFLAGS} ${LIBRARY}     


clean:
	/bin/rm -rf *.o *.exe *~ *++
