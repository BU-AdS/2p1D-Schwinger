#!/bin/bash

# Simple test script to demonstrate how to use the 2D U(1) code
# The size of the lattice (L) is hardcoded in the main.cpp file
# to make writing new code simpler. Please please edit and remake
# if you wish to vary L.

rm -rf {gauge,data}
mkdir -p {gauge,data/{data,plaq,creutz,polyakov,rect,top,pion,vacuum}}

# configure preamble
#---------------------------------------------------------------
LX=48
LY=48
LZ=3
# Construct the correct executable
cp main_template.cpp main.cpp
cp Makefile_template Makefile
sed -i '.bak' -e s/__LX__/${LX}/g main.cpp Makefile
sed -i '.bak' -e s/__LY__/${LY}/g main.cpp Makefile
sed -i '.bak' -e s/__LZ__/${LZ}/g main.cpp Makefile

rm *.bak
make
#---------------------------------------------------------------

# The value of the coupling in the U(1) 2D theory
BETA=$1

# The value of the coupling in the extra dim
BETAZ=$2

# The total number of HMC iterations to perform.
HMC_ITER=1000
# The number of HMC iterations for thermalisation.
HMC_THERM=25
# The number of HMC iterations to skip bewteen measurements.
HMC_SKIP=5
# Dump the gauge field every HMC_CHKPT iterations after thermalisation.
HMC_CHKPT=5000
# If non-zero, read in the HMC_CHKPT_START gauge field. 
HMC_CHKPT_START=0
# Number of HMC steps in the integration 
HMC_NSTEP=40
# HMC trajectory time
HMC_TAU=1.0

# Number of APE smearing hits to perform when measuring topology
APE_ITER=100
# The alpha value in the APE smearing
APE_ALPHA=0.5

# The RNG seed
RNG_SEED=1234

# DYNAMIC (1) or QUENCHED (0)
DYN_QUENCH=1

# Lock the Z gauge to unit (1) or allow z dynamics (0)
ZLOCKED=1

# Dynamic fermion parameters 
# Fermion mass
MASS=0.00
# Maximum CG iterations
MAX_CG_ITER=1000
# CG tolerance
CG_EPS=1e-16

# Eigensolver parameters
# Tolerance on the residual
TOL=1e-8
# Maximum ARPACK iterations
ARPACK_MAXITER=100000

#polyACC
USE_ACC=0
AMAX=11
AMIN=1.0
N_POLY=100

# Measuremets: 1 = measure, 0 = no measure
# Polyakov loops
MEAS_PL=1
# Wilson loops and Creutz ratios
MEAS_WL=1
# Pion Correlation function
MEAS_PC=1
# Vacuum trace
MEAS_VT=0

#./2p1D-Wilson-LX48-LY48-LZ3 5.0 1 1000 25 5 5000 0 40 1.0 5 0.5 1234 1 1 0.00 1000 1e-16 1e-8 100000 0 11 1.0 100 1 1 1 0

command="./2p1D-Wilson-LX$LX-LY$LY-LZ$LZ $BETA $BETAZ $HMC_ITER $HMC_THERM $HMC_SKIP 
	      $HMC_CHKPT $HMC_CHKPT_START $HMC_NSTEP $HMC_TAU $APE_ITER $APE_ALPHA 
	      $RNG_SEED $DYN_QUENCH $ZLOCKED $MASS $MAX_CG_ITER $CG_EPS $TOL 
	      $ARPACK_MAXITER $USE_ACC $AMAX $AMIN $N_POLY $MEAS_PL $MEAS_WL $MEAS_PC 
	      $MEAS_VT"

echo $command

$command
