#!/bin/bash

# configuration preamble to construct the correct executable
#---------------------------------------------------------------
PATH_TO_SOURCE=/Users/deanhowarth/2p1D/2p1D-Schwinger/wilson/3D
cp ${PATH_TO_SOURCE}/*template* .

cp main_template.cpp main.cpp
cp Makefile_template Makefile

LX=$1
LY=$2
LZ=$3

sed -i '.bak' -e s/__LX__/${LX}/g main.cpp Makefile
sed -i '.bak' -e s/__LY__/${LY}/g main.cpp Makefile
sed -i '.bak' -e s/__LZ__/${LZ}/g main.cpp Makefile

rm *.bak
rm *template*
make

mkdir -p {gauge,logs,data/{data,plaq,creutz,polyakov,rect,top,pion,vacuum}}
#---------------------------------------------------------------

BETA=$4
BETAZ=$5

HMC_STEP=32

for M in -0.04 -0.02 0.00 0.02 0.04 0.06 0.08 0.10 ; do
    
    ./launcher.sh ${BETA} ${BETAZ} ${LX} ${LY} ${LZ} ${HMC_STEP} ${M} >& logs/log_${LX}_${LY}_${LZ}_${BETA}_${BETAZ}_m${M}.log &
    
done
