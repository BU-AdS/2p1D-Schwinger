#!/bin/bash

# configuration preamble to construct the correct executable
#---------------------------------------------------------------
PATH_TO_SOURCE=/Users/deanhowarth/2p1D/2p1D-Schwinger/wilson/2D
cp ${PATH_TO_SOURCE}/*template* .

cp main_template.cpp main.cpp
cp Makefile_template Makefile

LX=$1
LY=$2

sed -i '.bak' -e s/__LX__/${LX}/g main.cpp Makefile
sed -i '.bak' -e s/__LY__/${LY}/g main.cpp Makefile

rm *.bak
rm *template*
make

mkdir -p {gauge,logs,data/{data,plaq,creutz,polyakov,rect,top,pion,vacuum}}
#---------------------------------------------------------------

BETA=$3
BETA_0=${BETA}

HMC_STEP_0=21

while [ ${BETA} -le ${BETA_0} ]; do
    
    HMC_STEP_E=$((${BETA}))
    HMC_STEP=$(($HMC_STEP_0 + $HMC_STEP_E))
    echo $HMC_STEP

    ./launcher.sh ${BETA}.0 ${LX} ${LY} ${HMC_STEP} >& logs/log_${LX}_${LY}_${BETA}p0.log &
    
    let BETA=BETA+1
    
done
