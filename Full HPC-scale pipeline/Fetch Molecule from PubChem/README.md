# Automatic stream-lined HPC-ready production level scripts to generate large molecular data

#### NOTE: 
Parallel version is not possible due to limitation in RDKit structure. RDKit molecules are not fully picklable —
when we submit an rdkit.Chem.Mol to a ProcessPoolExecutor, Python pickles it to send to another process. 
Serial streaming is already fast and memory-efficient. We can even add batching + deduplication, 
and it scales well to 20k–100k molecules on HPC.

**Requirements**
```
pip install rdkit-pypi pandas
```
We will use the PubChem database to fetch molecules. 
Navigate to 
```
https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF
```
Download an SDF file, e.g.
```
wget "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_048000001_048500000.sdf.gz"
```
### To run
```
python scriptname.py
```
We have clean, stable, HPC-safe script that:
* Uses PubChem FTP SDF dump
* Filters:\
Only C, H, N, O, F (customisable)
* Heavy atoms between 2 and 60  (customisable)

The output CSV has
```
CID
IUPAC_Name
Canonical_SMILES
Molecular_Formula
```
Like this:
```
CID,IUPAC_Name,Canonical_SMILES,Molecular_Formula
2244,propan-2-one,CC(=O)C,C3H6O
702,ethanol,CCO,C2H6O
...
```
**Basic Features:**
* CID + IUPAC + Canonical SMILES + Molecular Formula
* Strict filtering: only C,H,N,O,F and heavy atoms 2–60
* Failure logging: missing SMILES, invalid atoms, exceptions
* Nicely indented summary at the end
* Progress print every 100 molecules
* Stops at N molecules
* Memory-efficient streaming using ForwardSDMolSupplier
  i.e. never loads full SDF into memory
  


|  Scripts    | Features                                                     |
| ----------- | ---------------------------------------------------------------------------------- |
|   i) serial_00.py    | Basis features (listed above)   |
| ii) serial_01.py     |  Basic features + on screen timer/counter to track progress 
|                      |    + deduplication of repeated smiles    |
|iii) serial_02.py |  all of ii) + logging of repeated formulae & unique formulae |
|iv) serial_03.py | all of iii) + keep only only neutral, single-fragment molecules|

                        
                                                  
                                      
