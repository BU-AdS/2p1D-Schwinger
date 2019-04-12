#!/bin/bash

BETA=3
BETAZ=1

mkdir logs

while [ ${BETA} -le 3 ]; do
    
    #./launcher.sh ${BETA}.0 ${BETAZ} >& logs/log_${BETA}_${BETAZ}.log &
    ./launcher.sh ${BETA}.0 ${BETAZ} 
    
    let BETA=BETA+1
    
done
