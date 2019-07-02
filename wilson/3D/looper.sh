#!/bin/bash

BETA=5
BETAZ=1

mkdir logs

while [ ${BETA} -le 5 ]; do
    
    ./launcher.sh ${BETA}.0 ${BETAZ} 
    
    let BETA=BETA+1
    
done
