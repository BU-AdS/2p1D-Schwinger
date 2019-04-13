#!/bin/bash

#define the data set
BETA=4  #bulk beta
LX=48   #X extent
LY=48   #Y extent
T=25    #Folded timeslices
B=1     #JK block size
HMC_STEP=30

PATH_TO_2p1D=/Users/deanhowarth/2p1D/2p1D-Schwinger
DATA=/Users/deanhowarth/2p1D/analysis/dynamic/3D/beta${BETA}p0_betaZ1p0/LX${LX}_LY${LY}_LZ3/data/pion

#jack-knife calculator
#ensure error file is up to date:
(cd ${PATH_TO_2p1D}/utils/jack_knife; make)
cp ${PATH_TO_2p1D}/utils/jack_knife/jk_error .
MEFF_FILE="m_eff.dat"

#remove previous chrial extrapolation data
rm chiralExtrap.dat

for M in -0.04 -0.02 0.00 0.02 0.04 0.06 0.08 0.10; do
    
    FILE=pion_LX${LX}_LY${LY}_B${BETA}.000000_M${M}0000_tau1.000000_nHMCstep${HMC_STEP}
    
    echo "Looking for ${DATA}/${FILE}.dat"
    if [ -f ${DATA}/${FILE}.dat ]; then
	
	#query the file to get the number pf data points
	LINES="$(grep -c ^ ${DATA}/${FILE}.dat)"
	echo "Found ${LINES} data points."
	./jk_error ${DATA}/${FILE}.dat ${LINES} ${T} ${B}
	
	cp plot.p tmp.p
	sed -i '.bak' -e s/__T__/${LY}/g tmp.p
	sed -i '.bak' -e s/__MASS__/${M}/g tmp.p
	sed -i '.bak' -e s/__BETA__/${BETA}.0/g tmp.p
	sed -i '.bak' -e s/__TITLE__/${FILE}/g tmp.p
	sed -i '.bak' -e s/__MFILE__/${MEFF_FILE}/g tmp.p
	gnuplot tmp.p

	mv ${MEFF_FILE} ${FILE}.dat

    else
	echo "${FILE}.dat not found in ${DATA}."
	echo "Listing files in ${DATA}:"
	echo "$(ls -lht ${DATA})"
    fi
done

cp chiralExtrap.p tmp.p
sed -i '.bak' -e s/__BETA__/${BETA}/g tmp.p
gnuplot tmp.p
open chiralExtrap.eps

rm *.bak
rm tmp.p
