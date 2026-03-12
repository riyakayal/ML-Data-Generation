#!/bin/bash

# Author: Riya Kayal
# Created: 21/11/2025

## Adds the AA sequence after molecule Id/CID/molecule name. 
## The AA sequence is present in the MolId_SeqFile.

# sort by Id after Mol_
 # cp data_thermochemistry.csv temp
  #sort -t_ -k2 -n temp > energy_sorted.csv

## EnergyFile will be the thermochemistry data file (data_thermochemistry.csv) or 
## MP2 data files
## i.e. mp2_SVP.csv, mp2_TZVP.svp or mp2_QZVP.csv.
EnergyFile="$1"
header=$(head -n 1 $EnergyFile) 

## insert the "AA sequence" col name after 1st comma
header="${header/,/,AA sequence,}"
MolId_SeqFile="$2"
#MolId_SeqFile="peptides.csv"

#Sort energy file (skip header)
(head -n 1 $EnergyFile && tail -n +2 $EnergyFile | sort -t, -k1,1) > energy_sorted.csv

# Sort seq file (skip header)
(head -n 1 $MolId_SeqFile && tail -n +2 $MolId_SeqFile | sort -t, -k1,1) > smiles_sorted.csv

# Remove headers temporarily
tail -n +2 energy_sorted.csv > energy_noheader.csv
tail -n +2 smiles_sorted.csv > smiles_noheader.csv

# Join on first column
join -t, -1 1 -2 1 energy_noheader.csv smiles_noheader.csv > joined_tmp.csv
cat -A joined_tmp.csv | head

# Reorder last column to 2nd
awk -F, 'BEGIN{OFS=","} {
    out = $1 "," $NF
    for (i=2; i<NF; i++) out = out "," $i
    print out
}'  joined_tmp.csv > joined_reordered.csv

# Add proper header
echo "$header" > final.csv
cat joined_reordered.csv >> final.csv

#Remove the "^M" after smiles
tr -d '\r' < final.csv > tmp && mv tmp final.csv

# Check if ^M still present or not
cat -A final.csv | head


# Sort by molecule Id, not alphabetically
sort -t_ -k2 -n final.csv > temp
mv temp seq_$EnergyFile
rm -f final.csv

rm -f *sorted*csv *joined*csv *noheader*csv
