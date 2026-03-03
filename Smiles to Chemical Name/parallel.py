#!/usr/bin/env python3


# Author: Riya Kayal
# Created: 14/03/2025

import pandas as pd
from rdkit import Chem
import pubchempy as pcp
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import time
import logging
import multiprocessing

# Setup logging
logging.basicConfig(
    filename='molecule_fetch.log',
    filemode='w',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

# Function to fetch molecule name with retries
def get_molecule_name(smiles, max_retries=3, delay=0.2):
    mol_name = "Unknown"
    for attempt in range(max_retries):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                compound = pcp.get_compounds(smiles, 'smiles')
                if compound:
                    name = compound[0].synonyms[0]
                    mol_name = name.replace('"', '').strip().title()
            else:
                mol_name = "Invalid Smiles"
            if mol_name != "Unknown":
                break  # successful fetch
        except Exception as e:
            logging.warning(f"Attempt {attempt+1} failed for SMILES '{smiles}': {e}")
        time.sleep(delay)  # short delay before retry
    logging.info(f"Final result for SMILES '{smiles}': '{mol_name}'")
    return mol_name

# Main execution
if __name__ == "__main__":
    start_total = time.time()
    input_file = 'input.csv'
    output_file = 'output1.csv'

    df = pd.read_csv(input_file)
    smiles_list = df['SMILES'].tolist()
    num_entries = len(smiles_list)
    cpu_cores = multiprocessing.cpu_count()
    num_threads = min(5, cpu_cores)  # keep threads low to avoid API throttling

    logging.info(f"Processing {num_entries} molecules with {num_threads} threads")

    # Placeholder for results
    names = [None] * num_entries

    # Parallel processing with retries
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        future_to_index = {executor.submit(get_molecule_name, smi): i for i, smi in enumerate(smiles_list)}
        for future in as_completed(future_to_index):
            i = future_to_index[future]
            try:
                names[i] = future.result()
            except Exception as e:
                names[i] = "Unknown"
                logging.error(f"Unexpected error for index {i}, SMILES '{smiles_list[i]}': {e}")

    df['Name'] = names

    # Write output CSV
    df.to_csv(output_file, index=False, quoting=csv.QUOTE_MINIMAL)

    # Summary statistics
    total_elapsed = time.time() - start_total
    correct = sum(1 for n in names if n not in ["Unknown", "Invalid Smiles"])
    unknown = sum(1 for n in names if n in ["Unknown", "Invalid Smiles"])

    # Formatted final report with aligned colons
    label_width = 30  # adjust width for alignment
    print()
    print("=======================================================")
    print("        Molecule Conversion Summary (Parallel)")
    print("=======================================================")
    print(f"{'Number of CPU cores detected':<{label_width}} : {cpu_cores}")
    print(f"{'Total entries processed':<{label_width}} : {num_entries}")
    print(f"{'Successfully converted':<{label_width}} : {correct}")
    print(f"{'Unknown/Failed conversions':<{label_width}} : {unknown}")
    print(f"{'Total processing time (s)':<{label_width}} : {total_elapsed:.2f}")
    print("=======================================================")
    print(f"Output file saved as: {output_file}")
