#!/bin/bash

BETA=4

mkdir logs

while [ ${BETA} -le 4 ]; do
    
    ./launcher.sh ${BETA}.0 >& logs/log_${BETA}.log &
    
    let BETA=BETA+1
    
done
