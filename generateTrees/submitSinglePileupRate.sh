#!/bin/bash

# Set the pileup rate
PILEUP_RATE="0.002"

# Submit the job using star-submit
star-submit-template  -u ie -template Analysis_ToyModel1.xml -entities PILEUP_RATE=$PILEUP_RATE
