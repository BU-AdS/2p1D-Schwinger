#!/bin/bash

BETA=1
BETAZ=$1

mkdir logs

while [ ${BETA} -le 8 ]; do
    
    ./test.sh ${BETA}.0 ${BETAZ} >& logs/log_${BETA}_${BETAZ}.log &
    
    let BETA=BETA+1
    
done
