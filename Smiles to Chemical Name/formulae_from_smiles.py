#!/usr/bin/env python3


# Author: Riya Kayal
# Created: 12/03/2025

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

# Function to generate a name from SMILES using RDKit
def get_molecule_name(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # Try to get a common name if available, else return the SMILES as fallback
        try:
            name = rdMolDescriptors.CalcMolFormula(mol)
            return name
        except:
            return "Unknown"
    else:
        return "Invalid SMILES"

# Input CSV file path
input_file = 'molecules.csv'  # Replace with your file path
input_file = 'input.csv'  # Replace with your file path

# Read the input CSV file
df = pd.read_csv(input_file)

# Apply the function to get the molecule names
df['Name'] = df['SMILES'].apply(get_molecule_name)

# Write the output to a new CSV file
output_file = 'output.csv'  # Name for the output file
df.to_csv(output_file, index=False)

print(f"Conversion complete. Output file saved as {output_file}")
