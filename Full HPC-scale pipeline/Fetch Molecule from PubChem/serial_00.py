
#!/usr/bin/env python3


# Author: Riya Kayal
# Created: 01/09/2025

import gzip
import pandas as pd
from rdkit import Chem
import time
import traceback

# ==========================================================
# CONFIG
# ==========================================================
N = 1000
INPUT_SDF = "Compound_048000001_048500000.sdf.gz"  # example FTP file
OUTPUT_FILE = "organic_molecules.csv"
LOG_FILE = "organic_molecules.log"

MIN_HEAVY = 2
MAX_HEAVY = 20
ALLOWED_ELEMENTS = {"C", "H", "N", "O", "F"}

# ==========================================================
# Molecule filter
# ==========================================================
def valid_molecule(mol):
    if mol is None:
        return False

    heavy = mol.GetNumHeavyAtoms()
    if not (MIN_HEAVY <= heavy <= MAX_HEAVY):
        return False

    elements = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if not elements.issubset(ALLOWED_ELEMENTS):
        return False

    return True

# ==========================================================
# MAIN
# ==========================================================
start_time = time.time()
collected = []
failures = []

print("Reading SDF and filtering molecules...")

with gzip.open(INPUT_SDF, "rb") as f:
    supplier = Chem.ForwardSDMolSupplier(f)

    for idx, mol in enumerate(supplier, 1):
        if len(collected) >= N:
            break

        try:
            if not valid_molecule(mol):
                failures.append(f"Skipped molecule at index {idx}: invalid element/atom count")
                continue

            smiles = Chem.MolToSmiles(mol, canonical=True)
            if not smiles:
                failures.append(f"Skipped molecule at index {idx}: missing SMILES")
                continue

            # Extract PubChem properties safely
            cid = mol.GetProp("PUBCHEM_COMPOUND_CID") if mol.HasProp("PUBCHEM_COMPOUND_CID") else ""
            iupac = mol.GetProp("PUBCHEM_IUPAC_NAME") if mol.HasProp("PUBCHEM_IUPAC_NAME") else ""
            formula = mol.GetProp("PUBCHEM_MOLECULAR_FORMULA") if mol.HasProp("PUBCHEM_MOLECULAR_FORMULA") else ""

            collected.append({
                "CID": cid,
                "IUPAC_Name": iupac,
                "Canonical_SMILES": smiles,
                "Molecular_Formula": formula
            })

            if len(collected) % 100 == 0:
                print(f"Collected {len(collected)} molecules...")

        except Exception as e:
            failures.append(f"Error at molecule index {idx}: {str(e)}\n{traceback.format_exc()}")

# ==========================================================
# SAVE CSV
# ==========================================================
df = pd.DataFrame(collected)
df.to_csv(OUTPUT_FILE, index=False)

# Save log
with open(LOG_FILE, "w") as logf:
    logf.write(f"Total failures/skipped entries: {len(failures)}\n\n")
    for entry in failures:
        logf.write(entry + "\n")

# ==========================================================
# SUMMARY
# ==========================================================
end_time = time.time()
elapsed = end_time - start_time

print("\n" + "="*70)
print("Summary of Molecule Extraction")
print("="*70)
print(f"    Input SDF file         : {INPUT_SDF}")
print(f"    Output CSV file        : {OUTPUT_FILE}")
print(f"    Log file               : {LOG_FILE}")
print(f"    Requested molecules    : {N}")
print(f"    Successfully collected : {len(collected)}")
print(f"    Failed/skipped entries : {len(failures)}")
print(f"    Elapsed time           : {elapsed:.2f} seconds")
print("="*70 + "\n")
