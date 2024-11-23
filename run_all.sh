#!/bin/bash

# Lista de números associados aos arquivos
numbers=(1610 2020 3323 148 937 2085 2707 1272 6859 7013)

index=0
# Run ./build/bnb for all instances in the instances folder and write result in a file
echo "" > results.txt
for file in instances/*.tsp; do
    number=${numbers[$index]}
    echo "Running $file"
    (echo "INSTÂNCIA $file" && ./build/bnb $file 0 1 "$number") >> results.txt
    index=$((index+1))
done