#!/bin/bash

for i in `seq 0 26`;do
    echo "Starting Index ${i}, of 26",
    ./run ${i}
done

