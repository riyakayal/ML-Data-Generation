#!/usr/bin/env python3


# Author: Riya Kayal
# Created: 14/03/2025

import rdkit
from rdkit import Chem
import requests
import argparse

def smiles_to_inchikey(smiles: str) -> str:
    """Convert SMILES to InChIKey using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    inchi = Chem.MolToInchi(mol)
    inchikey = Chem.InchiToInchiKey(inchi)
    return inchikey

def fetch_pubchem_names(inchikey: str):
    """Given an InChIKey, query PubChem to get CID, IUPAC name, and synonyms."""
    url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
    resp = requests.get(url_cid)
    if resp.status_code != 200:
        return None

    data = resp.json()
    if 'IdentifierList' not in data or 'CID' not in data['IdentifierList']:
        return None

    cid = data['IdentifierList']['CID'][0]

    url_data = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
    resp2 = requests.get(url_data)
    if resp2.status_code != 200:
        return None

    compound_data = resp2.json()
    info = {'CID': cid}

    # Extract IUPAC name
    try:
        props = compound_data['PC_Compounds'][0]['props']
        for p in props:
            urn = p.get('urn', {})
            if urn.get('label') == 'IUPAC Name':
                info['IUPAC'] = p.get('value', {}).get('sval')
                break
    except Exception:
        info['IUPAC'] = None

    # Extract synonyms
    try:
        syns = compound_data['PC_Compounds'][0].get('synonyms', [])
        info['Synonyms'] = syns[:5] if isinstance(syns, list) else None
    except Exception:
        info['Synonyms'] = None

    return info

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch PubChem CID and names for a SMILES string")
    parser.add_argument("smiles", help="SMILES string of the molecule")
    args = parser.parse_args()

    smiles_input = args.smiles
    inchikey = smiles_to_inchikey(smiles_input)

    if not inchikey:
        print("Invalid SMILES string!")
    else:
        print(f"InChIKey: {inchikey}")
        result = fetch_pubchem_names(inchikey)
        if result:
            print(f"PubChem CID: {result.get('CID')}")
            iupac = result.get('IUPAC')
            if iupac:
                print(f"IUPAC Name  : {iupac}")
            syns = result.get('Synonyms')
            if syns:
                print("Synonyms     :")
                for name in syns:
                    print("  ", name)
        else:
            print("No PubChem record found.")
