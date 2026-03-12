#!/bin/bash

# Author: Riya Kayal
# Created: 10/03/2025

let run_opt=1,analyse_opt_result=0,run_dz=0
let run_tz=0,run_qz=0,analyse_mp2_result=0

#==============================================================================
#  I. SMILES --> XYZ --> DFT/semi-empirical optimisation (ORCA)+ FREQ 
#==============================================================================
if ((run_opt)); then
    if [ ! -d "output" ]; then
	    mkdir output
    fi	    
    cd output	
    ## 1. create a csv file with 2 columns: peptide id & composition
    cp ../scripts/gp_parallel.py .
    python  gp_parallel.py &
    wait

    cp ../scripts/create_orca_inputs.py .
    ## 3. create ORCA input files for optimisation + freq calculation
    python create_orca_inputs.py &
    wait

    cp ../scripts/submit_opt_jobs.sh .
    ## 4. submit jobs on the queue
    ./submit_opt_jobs.sh &
    wait

    rm -f *py *sh

    ## wait till all optimisation jobs are finished. Adjust the time
    ## of sleep as needed.
    # sleep 300
    
fi

#==============================================================================
#  II. Optimised Geometry: Imaginary FREQ check + get Thermochemistry data
#==============================================================================
if ((analyse_opt_result));then

    cd output
    cp ../scripts/check_imag_freq_parallel.py .    
    ## 5. check imaginary frequencies
    python check_imag_freq_parallel.py orca_opt
    wait
   
    cp ../scripts/thermochemistry.sh .
    cp ../scripts/add_AA_seq.sh .
    ## 6. Grep thermochemistry data from DFT optimised output files
    ./thermochemistry.sh &
    wait
    ./add_AA_seq.sh data_thermochemistry.csv peptides.csv &
    wait

    rm -f *py *sh
    cd ../
fi
#exit()
#==============================================================================
#  III. Run MP2 at different basis sets
#==============================================================================
    if ((run_dz));then
    cd output
    cp ../scripts/create_submit_mp2_input.sh .
    ## 7. create & submit MP2 ORCA input files for def2-SVP
    ./create_submit_mp2_input.sh -a SVP &
    wait
    cd ../
    fi

    if ((run_tz));then
    cd output
    cp ../scripts/create_submit_mp2_input.sh .
    ## 8. create & submit MP2 ORCA input files for def2-TZVP
    ./create_submit_mp2_input.sh -a TZVP &
    wait
    cd ../
    fi
    
    if ((run_qz));then
    cd output
    cp ../scripts/create_submit_mp2_input.sh .
    ## 9. create & submit MP2 ORCA input files for def2-QZVP
    ./create_submit_mp2_input.sh -a QZVP &
    wait
    cd ../
    fi

#==============================================================================
#  IV. Extract HF, MP2 energies, dipole moment 
#==============================================================================
if ((analyse_mp2_result));then
    cd output
    cp ../scripts/grep_mp2.sh .
    cp ../scripts/add_AA_seq.sh .
    ## 10. grep MP2 data from ORCA output files for def2-SVP
    ./grep_mp2.sh -a SVP &
    wait
    ./add_AA_seq.sh mp2_SVP.csv peptides.csv &
    wait

    ## 11. grep MP2 data from ORCA output files for def2-TZVP
    # ./grep_mp2.sh -a TZVP &
    # wait
    # ./add_AA_seq.sh mp2_TZVP.csv peptides.csv &
    wait

    ## 12. grep MP2 data from ORCA output files for def2-QZVP
    #./grep_mp2.sh -a QZVP &
    #wait
    #./add_AA_seq.sh mp2_QZVP.csv peptides.csv &
    # wait
	rm -f *sh
fi
