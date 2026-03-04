
#!/usr/bin/env python3


# Author: Riya Kayal
# Created: 01/09/2025


import gzip
import pandas as pd
from rdkit import Chem
import time
from itertools import islice

# ===========================
# CONFIG
# ===========================
N = 5000  # Number of molecules to collect
INPUT_SDF = "Compound_048000001_048500000.sdf.gz"
OUTPUT_FILE = "organic_molecules.csv"
LOG_FILE = "organic_molecules.log"

MIN_HEAVY = 2
MAX_HEAVY = 60
ALLOWED_ELEMENTS = {"C", "H", "N", "O", "F"}

BATCH_PROGRESS = 100  # Print progress every N molecules

# ===========================
# Helper: molecule valid check
# ===========================
def valid_molecule(mol):
    if mol is None:
        return False, "None molecule"

    heavy = mol.GetNumHeavyAtoms()
    if heavy < MIN_HEAVY or heavy > MAX_HEAVY:
        return False, f"heavy atom count {heavy}"

    elements = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if not elements.issubset(ALLOWED_ELEMENTS):
        return False, f"elements {elements}"

    smiles = Chem.MolToSmiles(mol, canonical=True)
    if not smiles:
        return False, "missing SMILES"

    return True, smiles

# ===========================
# MAIN
# ===========================
start_time = time.time()
collected = []
seen_smiles = set()
failures = []

print("Streaming SDF and extracting molecules...")

with gzip.open(INPUT_SDF, "rb") as f:
    supplier = Chem.ForwardSDMolSupplier(f)
    
    for idx, mol in enumerate(supplier, 1):
        if len(collected) >= N:
            break

        try:
            valid, result = valid_molecule(mol)
            if not valid:
                failures.append(f"Skipped molecule {idx}: {result}")
                continue

            smiles = result
            if smiles in seen_smiles:
                failures.append(f"Skipped molecule {idx}: duplicate SMILES")
                continue

            # Extract properties
            cid = mol.GetProp("PUBCHEM_COMPOUND_CID") if mol.HasProp("PUBCHEM_COMPOUND_CID") else ""
            iupac = mol.GetProp("PUBCHEM_IUPAC_NAME") if mol.HasProp("PUBCHEM_IUPAC_NAME") else ""
            formula = mol.GetProp("PUBCHEM_MOLECULAR_FORMULA") if mol.HasProp("PUBCHEM_MOLECULAR_FORMULA") else ""

            collected.append({
                "CID": cid,
                "IUPAC_Name": iupac,
                "Canonical_SMILES": smiles,
                "Molecular_Formula": formula
            })
            seen_smiles.add(smiles)

            if len(collected) % BATCH_PROGRESS == 0:
                elapsed = time.time() - start_time
                print(f"Collected {len(collected)}/{N} molecules, elapsed time: {elapsed:.1f}s")

        except Exception as e:
            failures.append(f"Error at molecule {idx}: {str(e)}")

# ===========================
# Save CSV
# ===========================
df = pd.DataFrame(collected)
df.to_csv(OUTPUT_FILE, index=False)

# Save log
with open(LOG_FILE, "w") as logf:
    logf.write(f"Total failures/skipped entries: {len(failures)}\n\n")
    for entry in failures:
        logf.write(entry + "\n")

# ===========================
# Summary
# ===========================
end_time = time.time()
elapsed = end_time - start_time

print("\n" + "="*70)
print("Summary of Molecule Extraction (Serial Streaming + Deduplication)")
print("="*70)
print(f"    Input SDF file           : {INPUT_SDF}")
print(f"    Output CSV file          : {OUTPUT_FILE}")
print(f"    Log file                 : {LOG_FILE}")
print(f"    Requested molecules      : {N}")
print(f"    Successfully collected   : {len(collected)}")
print(f"    Failed/skipped entries   : {len(failures)}")
print(f"    Elapsed time             : {elapsed:.2f} seconds")
print("="*70 + "\n")
