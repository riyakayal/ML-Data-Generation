#!/bin/bash

# Author: Riya Kayal
# Created: 05/03/2025



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
OUTPUT_DIR="mp2_$Basis_set"


## write results in a csv file

    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "Directory $OUTPUT_DIR not found."
        exit 1
    fi
    
    cd $OUTPUT_DIR
    rm -f log
    echo "Molecule,    E(HF) (au),   E(MP2) (au), ΔE(MP2 - HF) (au), dipole moment (au) " > mp2

    for file in *.inp; do
	basename="${file%.inp}"
	outfile="${basename}.mpi8.out"
        [ -e "$outfile" ] || continue  # skip if no output file found
       	
	if  ! grep -q "****ORCA TERMINATED NORMALLY****" "$outfile"; then
           echo " $file did not run properly" >> log
           continue
        fi


	echo "current file: $outfile"
        HFEnergy=$(awk '/Total Energy       :/{printf "%.12g\n", $4}' $outfile)
        MP2Energy=$(awk '/FINAL SINGLE POINT ENERGY/{printf "%.12g\n", $5}' $outfile)	
	DipoleMoment=$(awk '/Magnitude \(a.u.\)/{printf "%.8g\n", $4}' $outfile)	
        deltaEnergy=`awk "BEGIN{ print $MP2Energy - $HFEnergy}"`
        echo "$basename, $HFEnergy, $MP2Energy, $deltaEnergy, $DipoleMoment" >> mp2

      done
      cp mp2 ../mp2_$Basis_set.csv
      cp log ../.
    cd ../

echo "------------------------------------------------"   
echo "Results are written in "mp2_$Basis_set.csv"" 
echo "------------------------------------------------"   
