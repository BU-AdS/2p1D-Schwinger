#!/bin/bash

# Simple test script to demonstrate how to use the 2D U(1) code
# The size of the lattice (L) is hardcoded in the main.cpp file
# to make writing new code simpler. Please please edit and remake
# if you wish to vary L.

# The value of the coupling in the U(1) 2D theory
BETA=5.0

# The total number of HMC iterations to perform.
HMC_ITER=110

# The number of HMC iterations for thermalisation.
HMC_THERM=50

# The number of HMC iterations to skip bewteen measurements.
HMC_SKIP=10

# Dump the gauge field every HMC_CHKPT iterations after thermalisation.
HMC_CHKPT=10

# If non-zero, read in the HMC_CHKPT_START gauge field. 
HMC_CHKPT_START=0

# Number of APE smearing hits to perfrom when measuring topology
APE_ITER=0

# The 'alpha' value in the APE smearing
APE_ALPHA=0.5

# The RNG seed
RNG_SEED=1234

./u1 $BETA $HMC_ITER $HMC_THERM $HMC_SKIP $HMC_CHKPT $HMC_CHKPT_START $APE_ITER $APE_ALPHA $RNG_SEED 
