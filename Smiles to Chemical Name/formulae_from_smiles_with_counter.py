#!/usr/bin/env python3


# Author: Riya Kayal
# Created: 14/03/2025

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import time
import logging
import csv
import multiprocessing

# Setup logging
logging.basicConfig(
    filename='molecule_formula.log',
    filemode='w',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

# Function to generate a formula from SMILES using RDKit
def get_molecule_name(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Generate chemical formula
            formula = rdMolDescriptors.CalcMolFormula(mol)
            return formula
        else:
            return "Invalid SMILES"
    except Exception as e:
        logging.error(f"Error processing SMILES '{smiles}': {e}")
        return "Unknown"

# Main execution
if __name__ == "__main__":
    start_total = time.time()

    input_file = 'input.csv'  # Replace with your file path
    output_file = 'output.csv'
    summary_file = 'summary.csv'

    df = pd.read_csv(input_file)
    smiles_list = df['SMILES'].tolist()
    num_entries = len(smiles_list)
    cpu_cores = multiprocessing.cpu_count()

    logging.info(f"Starting processing of {num_entries} molecules")

    names = []
    status_list = []

    for i, smi in enumerate(smiles_list, start=1):
        name = get_molecule_name(smi)
        names.append(name)

        # Determine status
        if name in ["Unknown", "Invalid SMILES"]:
            status = name
        else:
            status = "Converted"
        status_list.append(status)

        logging.info(f"{i}/{num_entries} - SMILES: {smi} → Formula: {name}")
        print(f"Processed {i}/{num_entries} molecules", end='\r')

    print()

    df['Name'] = names

    # Write the output CSV
    df.to_csv(output_file, index=False, quoting=csv.QUOTE_MINIMAL)

    # Write the summary CSV
    summary_df = pd.DataFrame({
        'SMILES': smiles_list,
        'Formula': names,
        'Status': status_list
    })
    summary_df.to_csv(summary_file, index=False, quoting=csv.QUOTE_MINIMAL)

    # Summary statistics
    total_elapsed = time.time() - start_total
    correct = sum(1 for s in status_list if s == "Converted")
    unknown = sum(1 for s in status_list if s == "Unknown")
    invalid = sum(1 for s in status_list if s == "Invalid SMILES")

    label_width = 30
    label_width = 30  # adjust width for alignment
    print()
    print("=======================================================")
    print("        Molecule Conversion Summary (Serial)")
    print("=======================================================")
    print(f"{'Number of CPU cores detected':<{label_width}} : {cpu_cores}")
    print(f"{'Total entries processed':<{label_width}} : {num_entries}")
    print(f"{'Successfully converted':<{label_width}} : {correct}")
    print(f"{'Unknown conversions':<{label_width}} : {unknown}")
    print(f"{'Invalid SMILES':<{label_width}} : {invalid}")
    print(f"{'Total processing time (s)':<{label_width}} : {total_elapsed:.2f}")
    print("=======================================================")
    print(f"Output file saved as: {output_file}")
    print(f"Summary file saved as: {summary_file}")
    print(f"Log file saved as: molecule_formula.log")
