#!/bin/bash

# Author: Riya Kayal
# Created: 03/03/2025

## replace with your orca batch submit script name and path
suborca='/hsnfs/software/orca/suborca'

INPUT_DIR="orca_opt"


    if [ ! -d "$INPUT_DIR" ]; then
        echo "Directory $INPUT_DIR not found."
        exit 1
    fi

    cd $INPUT_DIR
    for file in *.inp; do
        [ -e "$file" ] || continue  # skip if no files found
        
        echo "current File: $(basename "$file")"
	basename="${file%.inp}"
	
	if  grep -q "****ORCA TERMINATED NORMALLY****" "${basename}.mpi8.out"; then
           continue
        fi
        echo "---------------------------"
    
    ## replace suborca with your orca sbatch submit script name	
    $suborca -q ep29th  "${basename}"
    done

cd ../
