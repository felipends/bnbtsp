#!/bin/bash

# Run ./build/bnb for all instances in the instances folder and write result in a file
for file in instances/*.tsp; do
    echo "Running $file"
    (echo "INSTÂNCIA $file" && ./build/bnb $file 1) >> results.txt
done