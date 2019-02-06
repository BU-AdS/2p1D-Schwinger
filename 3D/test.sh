#!/bin/bash

# Simple test script to demonstrate how to use the (2+1)D U(1) code
# The size of the lattice (L,Lz) is hardcoded in the main.cpp file
# to make writing new code simpler. Please please edit and remake
# if you wish to vary L or Lz

# The value of the coupling in the x-, y-dimensions
BETA=4.0

# The value of the coupling in the z-dimensions.
# By construction, beta(Lz-1) = 0 to implement
# open boundary conditions.
BETAZ=4.0

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
HMC_NSTEP=200
# Ficticious time step size
HMC_DT=0.005

# Number of APE smearing hits to perfrom when measuring topology
APE_ITER=2
# The 'alpha' value in the APE smearing
APE_ALPHA=0.5

# The RNG seed
RNG_SEED=1234

# DYNAMIC (1) or QUENCHED (0)
DYN_QUENCH=0

# Lock the Z gauge to unit (1) or allow z dynamics (0)
ZLOCKED=1

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
NEV=16
# Tolerance on the residual
TOL=1e-12
# Maximum ARPACK iterations
ARPACK_MAXITER=1000000

./u1 $BETA $BETAZ $HMC_ITER $HMC_THERM $HMC_SKIP $HMC_CHKPT $HMC_CHKPT_START $HMC_NSTEP \
     $HMC_DT $APE_ITER $APE_ALPHA $RNG_SEED $DYN_QUENCH $ZLOCKED $MASS $MAX_CG_ITER $CG_EPS \
     $NKV $NEV $TOL $ARPACK_MAXITER
