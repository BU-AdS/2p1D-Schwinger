#!/bin/bash

L=$1
BETA=$2

for L in ${L}
do
    
    mkdir LX${L}_LY${L}
    (cd LX${L}_LY${L}; cp ../looper.sh .; cp ../launcher.sh .; ./looper.sh $L $L $BETA)
    
done
