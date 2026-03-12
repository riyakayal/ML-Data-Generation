#!/bin/bash

# Author: Riya Kayal
# Created: 05/03/2025

## replace with your orca batch submit script name and path
suborca='/hsnfs/software/orca/suborca'

help()
{
   echo ""
   echo "Usage: $0 -a Basis_set"
   echo -e "\t-a OUTPUT_DIR"
   exit 1 # Exit script after printing help
}

while getopts "a:" opt
do
   case "$opt" in
      a ) Basis_set="$OPTARG" ;;
      ? ) help ;; # Print help in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$Basis_set" ] 
then
   echo "missing argument for output directory name";
   help
fi

# Begin script in case all parameters are correct
echo "$Basis_set"

## change basis set name as needed, e.g. "SVP", "TZVP", "QZVP" 
INPUT_DIR="orca_opt"
OUTPUT_DIR="mp2_$Basis_set"

let create_inputs=1 submit_jobs=1 

#----------------------------------------------------
# Create inputs
#----------------------------------------------------

if ((create_inputs));then
    if [ ! -d "$INPUT_DIR" ]; then
        echo "Directory $INPUT_DIR not found."
        exit 1
    fi
    
    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "Directory $OUTPUT_DIR not found."
        mkdir $OUTPUT_DIR
    fi

    cd $INPUT_DIR
    rm -f log
    for file in *.inp; do
        [ -e "$file" ] || continue  # skip if no files found
       	
        echo "current File: $(basename "$file")"
	basename="${file%.inp}"
	outfile="${basename}.mpi8.out"
        echo "---------------------------"
       	
	if  ! grep -q "****ORCA TERMINATED NORMALLY****" "$outfile"; then
           echo " $file did not run properly" >> log
           continue
        fi
	
	if  grep -q "The optimization did not converge but reached the maximum" "$outfile"; then
           echo " $file maximum optimization cycly reached" >> log
           continue
        fi

    ## create the input file	
    cp $file temp
    cp temp temp1
    #echo "$(head -n1 temp)"
    header="! MP2 def2-SVP TightSCF"

    #sed -i 's/B3LYP D3BJ def2-SVP TightSCF TightOpt defGrid4 Freq/MP2 def2-SVP TightSCF/g' temp
    #sed -i 's/r2scan-3c TightOpt Freq/MP2 def2-SVP TightSCF/g' temp
    ## remove 1st line from temp
    tail -n +2 temp > temp1
    echo "$header" > temp
    cat temp1 >> temp
    sed -i "s/SVP/$Basis_set/g" temp
    #echo "$(head -n1 temp)"

    sed -i 's/maxcore 2000/maxcore 10000/g' temp
    sed -i '6,$d' temp
    echo "*xyzfile 0 1 "$basename.xyz" " >> temp
    cp $basename.xyz ../$OUTPUT_DIR/.
    mv temp ../$OUTPUT_DIR/$file	
    rm -f temp1
    done
cd ..
fi

#----------------------------------------------------
# Submit jobs
#----------------------------------------------------

if ((submit_jobs));then
    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "Directory $OUTPUT_DIR not found."
        exit 1
    fi
    
    cd $OUTPUT_DIR
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
fi
