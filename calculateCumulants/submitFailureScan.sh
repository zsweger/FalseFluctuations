#!/bin/bash

# Define an array of failure rates
FAILURE_RATES=("0.00001" "0.00002" "0.00005" "0.0001" "0.0002" "0.0005" "0.001" "0.002" "0.005" "0.01")

# Loop over each failure rate
for FAILURE_RATE in "${FAILURE_RATES[@]}"
do
    echo "******* WORKING ON FAILURE RATE ${FAILURE_RATE} *******"
    rm file.list
    realpath ../generateTrees/outdir/*${FAILURE_RATE}failureRate.root > file.list
    echo "Using filelist: "
    cat file.list

    ./run.csh
    wait
done
