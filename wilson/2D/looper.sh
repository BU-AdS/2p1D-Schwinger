#!/bin/bash

BETA=2

mkdir logs

while [ ${BETA} -le 2 ]; do
    
    ./launcher.sh ${BETA}.0 >& logs/log_${BETA}.log &
    
    let BETA=BETA+1
    
done
