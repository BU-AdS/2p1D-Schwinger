#!/bin/bash

# Simple test script to demonstrate how to use the (2+1)D U(1) code
# The size of the lattice (L,Lz) is hardcoded in the main.cpp file
# to make writing new code simpler. Please please edit and remake
# if you wish to vary L or Lz

# The value of the coupling in the x-, y-dimensions
BETA=1.0

# The value of the coupling in the z-dimensions.
# By construction, beta(Lz-1) = 0 to implement
# open boundary conditions.
BETAZ=1.0

# The total number of HMC iterations to perform.
HMC_ITER=6000

# The number of HMC iterations for thermalisation.
HMC_THERM=100

# The number of HMC iterations to skip bewteen measurements.
HMC_SKIP=100

# Dump the gauge field every HMC_CHKPT iterations after thermalisation.
HMC_CHKPT=100000000

# If non-zero, read in the HMC_CHKPT_START gauge field. 
HMC_CHKPT_START=0

# Number of APE smearing hits to perfrom when measuring topology
APE_ITER=2

# The 'alpha' value in the APE smearing
APE_ALPHA=0.5

# The RNG seed
RNG_SEED=1234

./u1 $BETA $BETAZ $HMC_ITER $HMC_THERM $HMC_SKIP $HMC_CHKPT $HMC_CHKPT_START $APE_ITER $APE_ALPHA $RNG_SEED 
