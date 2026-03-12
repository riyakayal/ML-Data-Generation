# Example 01: Fulll HPC-scale pipelene for Thermochmistry and MP2 Energies for polypeptides 

To run,
```
./run_mp2.sh &
```
##### Follow these steps sequentially.
#### 1. Run with setting run_opt=1 and all other flags at the top at 0 --> Create Peptide + Optimise with ORCA
#### 2. Run with setting analyse_opt_result=1 and all other flags at 0 --> Generate Thermochemistry data
#### 3. Run with setting run_dz=1 and all other flags at 0 --> Run mp2 input files
#### 4. Run with setting analyse_mp2_result=1 and all other flags at 0 --> Generate MP2 data
