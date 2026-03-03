#!/usr/bin/env python3


# Author: Riya Kayal
# Created: 12/03/2025


import pandas as pd
from rdkit import Chem
import requests
import time
import logging
import multiprocessing
import csv

# Setup logging
logging.basicConfig(
    filename='molecule_fetch_pubchem.log',
    filemode='w',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

def smiles_to_inchikey(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    inchi = Chem.MolToInchi(mol)
    inchikey = Chem.InchiToInchiKey(inchi)
    return inchikey

def fetch_pubchem_name(inchikey: str, max_retries=3, delay=0.2):
    mol_name = "Unknown"
    for attempt in range(max_retries):
        try:
            url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
            resp = requests.get(url_cid)
            if resp.status_code != 200:
                time.sleep(delay)
                continue
            data = resp.json()
            if 'IdentifierList' not in data or 'CID' not in data['IdentifierList']:
                continue
            cid = data['IdentifierList']['CID'][0]

            url_data = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
            resp2 = requests.get(url_data)
            if resp2.status_code != 200:
                continue
            compound_data = resp2.json()

            try:
                props = compound_data['PC_Compounds'][0]['props']
                for p in props:
                    urn = p.get('urn', {})
                    if urn.get('label') == 'IUPAC Name':
                        mol_name = p.get('value', {}).get('sval', 'Unknown').replace('"','').strip().title()
                        break
            except Exception:
                mol_name = "Unknown"

            if mol_name != "Unknown":
                break
        except Exception as e:
            logging.warning(f"Attempt {attempt+1} failed for InChIKey {inchikey}: {e}")
        time.sleep(delay)
    logging.info(f"Final result for InChIKey {inchikey}: {mol_name}")
    return mol_name

# Main script
if __name__ == "__main__":
    start_total = time.time()
    input_file = "input.csv"
    output_file = "output.csv"
    summary_file = "summary.csv"

    df = pd.read_csv(input_file)
    smiles_list = df['SMILES'].tolist()
    num_entries = len(smiles_list)
    cpu_cores = multiprocessing.cpu_count()

    logging.info(f"Starting processing of {num_entries} molecules")

    names = []
    status_list = []

    for i, smi in enumerate(smiles_list, start=1):
        inchikey = smiles_to_inchikey(smi)
        if inchikey:
            name = fetch_pubchem_name(inchikey)
        else:
            name = "Invalid SMILES"

        # Determine status
        if name == "Unknown":
            status = "Unknown"
        elif name == "Invalid SMILES":
            status = "Invalid SMILES"
        else:
            status = "Converted"

        names.append(name)
        status_list.append(status)

        print(f"Processed {i}/{num_entries} molecules", end='\r')

    print()
    df['Name'] = names

    # Save output CSV (names included)
    df.to_csv(output_file, index=False, quoting=csv.QUOTE_MINIMAL)

    # Save summary CSV
    summary_df = pd.DataFrame({
        'SMILES': smiles_list,
        'Name': names,
        'Status': status_list
    })
    summary_df.to_csv(summary_file, index=False, quoting=csv.QUOTE_MINIMAL)

    # Summary statistics
    total_elapsed = time.time() - start_total
    correct = sum(1 for s in status_list if s == "Converted")
    unknown = sum(1 for s in status_list if s == "Unknown")
    invalid = sum(1 for s in status_list if s == "Invalid SMILES")

    label_width = 30
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
