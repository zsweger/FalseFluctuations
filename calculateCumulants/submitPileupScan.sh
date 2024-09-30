#!/bin/bash

# Define an array of pileup rates
PILEUP_RATES=("0.00001" "0.00001" "0.00002" "0.00005" "0.0001" "0.0002" "0.0005" "0.001" "0.002" "0.005" "0.01")

# Loop over each pileup rate
for PILEUP_RATE in "${PILEUP_RATES[@]}"
do
    echo "******* WORKING ON PILEUP RATE ${PILEUP_RATE} *******"
    rm file.list
    realpath ../generateTrees/outdir/*${PILEUP_RATE}pileupRate.root > file.list
    echo "Using filelist: "
    cat file.list

    ./run.csh
    wait
done
