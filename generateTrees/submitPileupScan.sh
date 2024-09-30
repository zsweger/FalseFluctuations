#!/bin/bash

# Define an array of pileup rates
PILEUP_RATES=("0.00001" "0.00001" "0.00002" "0.00005" "0.0001" "0.0002" "0.0005" "0.001" "0.002" "0.005" "0.01")

# Loop over each pileup rate
for PILEUP_RATE in "${PILEUP_RATES[@]}"
do
    # Submit the job using star-submit
    star-submit-template  -u ie -template Analysis_ToyModel1.xml -entities PILEUP_RATE=$PILEUP_RATE
done
