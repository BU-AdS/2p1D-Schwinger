#!/bin/bash

# Simple test script to demonstrate how to use the 2D U(1) code
# The size of the lattice (L) is hardcoded in the main.cpp file
# to make writing new code simpler. Please please edit and remake
# if you wish to vary L.

# The value of the coupling in the U(1) 2D theory
BETA=3.0

# The total number of HMC iterations to perform.
HMC_ITER=10000
# The number of HMC iterations for thermalisation.
HMC_THERM=1
# The number of HMC iterations to skip bewteen measurements.
HMC_SKIP=1
# Dump the gauge field every HMC_CHKPT iterations after thermalisation.
HMC_CHKPT=1000
# If non-zero, read in the HMC_CHKPT_START gauge field. 
HMC_CHKPT_START=0
# Number of HMC steps in the integration 
HMC_NSTEP=100
# Ficticious time step size
HMC_DT=0.01

# Number of APE smearing hits to perfrom when measuring topology
APE_ITER=2
# The 'alpha' value in the APE smearing
APE_ALPHA=0.5

# The RNG seed
RNG_SEED=1234

# DYNAMIC (1) or QUENCHED (0)
DYN_QUENCH=1

# Fermion mass
MASS=0.1

# Maximum CG iterations
MAX_CG_ITER=1000

# CG tolerance
CG_EPS=1e-6

# ARPACK parameters
# Krylov space size
NKV=32
# Requested converged eigenpairs
NEV=8
# Tolerance on the residual
TOL=1e-12
# Maximum ARPACK iterations
ARPACK_MAXITER=1000000

# OMP therads
THREADS=1

./u1 $BETA $HMC_ITER $HMC_THERM $HMC_SKIP $HMC_CHKPT $HMC_CHKPT_START $HMC_NSTEP \
     $HMC_DT $APE_ITER $APE_ALPHA $RNG_SEED $DYN_QUENCH $MASS $MAX_CG_ITER $CG_EPS \
     $NKV $NEV $TOL $ARPACK_MAXITER $THREADS
