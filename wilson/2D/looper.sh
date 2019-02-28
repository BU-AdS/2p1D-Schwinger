#!/bin/bash

BETA=1

mkdir -p {gauge,logs,data/{data,plaq,creutz,polyakov,rect,top}}

while [ ${BETA} -le 8 ]; do
    
    ./test.sh ${BETA}.0 >& logs/log_${BETA}.log &
    
    let BETA=BETA+1
    
done
