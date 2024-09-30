#!/bin/bash

# Set the failure rate
FAILURE_RATE="0.01"

# Submit the job using star-submit
star-submit-template  -u ie -template Analysis_ToyModel2.xml -entities FAILURE_RATE=$FAILURE_RATE
