#!/bin/bash


# Author: Riya Kayal
# Created: 05/03/2025

## change basis set name as needed, e.g. "SVP", "TZVP", "QZVP" 
OUTPUT_DIR="orca_opt"



    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "Directory $OUTPUT_DIR not found."
        exit 1
    fi
    
    cd $OUTPUT_DIR
    rm -f log

    ## write results in a csv file
    echo "Molecule Id,  E(opt) (au), ZPE (au),Thermal Correction to Enthalpy(au), Gibbs Free Energy (au), Entropy(au)" > data

    for file in *.inp; do
	basename="${file%.inp}"
	outfile="${basename}.mpi8.out"
        [ -e "$outfile" ] || continue  # skip if no output file found
       	
	if  ! grep -q "****ORCA TERMINATED NORMALLY****" "$outfile"; then
           echo " $file did not run properly" >> log
           continue
        fi
	
	if  grep -q "The optimization did not converge but reached the maximum" "$outfile"; then
           echo " $file maximum optimization cycly reached" >> log
           continue
        fi
	
	if  grep -q "***imaginary mode***" "$outfile"; then
           echo " $file imaginary frequency found" >> log
           continue
        fi
	
	echo "current file: $outfile"
        
        ## E_DFT will be an array
	E_DFT=$(awk '/FINAL SINGLE POINT ENERGY/{printf "%.12g\n", $5}' $outfile)
        arrEnergy=(${E_DFT// / })

	ZPE=$(awk '/Zero point energy                .../{printf "%.12g\n", $5}' $outfile)	
	ThermalCorr=$(awk '/Thermal Enthalpy correction       .../{printf "%.12g\n", $5}' $outfile)	
	GibbsEnergy=$(awk '/Final Gibbs free energy/{printf "%.12g\n", $6}' $outfile)	

	Entropy=$(awk '/Final entropy term/{printf "%.12g\n", $5}' $outfile)	
        echo "$basename,  ${arrEnergy[-1]},  $ZPE, $ThermalCorr,  $GibbsEnergy, $Entropy" >> data

      done
      mv data ../data_thermochemistry.csv
      cp log ../.
    cd ../

echo "--------------------------------------------------"   
echo "Results are written in data_thermochemistry.csv"
echo "--------------------------------------------------"   
