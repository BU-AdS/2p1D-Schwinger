#!/bin/bash

BETA=1
BETA_Z=1.0

mkdir -p beta_${BETA_Z}/{gauge,data/{data,plaq,creutz,polyakov,rect,top}}

while [ ${BETA} -le 8 ]; do
    
    ./test.sh ${BETA} ${BETA_Z} >& beta_${BETA_Z}/log_${BETA}_${BETA_Z}.log &
    
    let BETA=BETA+1
    
done
