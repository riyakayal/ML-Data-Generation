**Requirements:**
  ```
  pip install rdkit pubchempy
  ```

**A) Get the chemical name of a molecule from smiles**\
We can do it:
* using RDKit and pubchempy\
**scripts:** \
i) no CPU parallelisation: serial.py, serial_with_counter.py\
ii) CPU parallelisation enabled: parallel.py

* looking it up if it exists in the PubChem records.\
  Here, **smiles --> InChIKey --> PubChem lookup** is done which is more reliable than the above procedure.\
**script:** \
serial_with_pubchem.py


**B) Get the chemical formula of a molecule from smiles**\
**scripts:** \
formulae_from_smiles.py\
formulae_from_smiles_with_counter.py

**NOTE:**
* *the input file (input.csv) for all these will be in the format:*\
  \
Name,SMILES\
Mol_1,CC#N\
Mol_2,CC(=O)O\
...\
Mol_1000,CCCCC(=O)NCCC
* All the scripts for A) and B) can be run with
  ```
  python script_name.py
  ```

**C) To look up the IUPAC name of a single smiles string, e.g. "CCC"**\
use this command:
```
python fetch_pubchem_name.py CCC
```
