#!/bin/bash

L=$1
LZ=$2
BETA=$3
BETAZ=$4

for L in ${L}; do

    DIR=beta${BETA}p0_betaZ${BETAZ}p0/LX${L}_LY${L}_LZ${LZ}
    mkdir -p ${DIR}
    (cd ${DIR}; cp ../../looper.sh .; cp ../../launcher.sh .; ./looper.sh $L $L $LZ $BETA $BETAZ)
    
done
