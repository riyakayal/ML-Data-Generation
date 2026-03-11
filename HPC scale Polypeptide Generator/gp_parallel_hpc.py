import os
import sys
import random
import multiprocessing as mp
import csv
from itertools import product

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms

########################################
# Amino acid dictionary
########################################

AA_CODES = {
    "ala": "A",
    "gly": "G",
    "val": "V",
    "leu": "L",
    "ser": "S",
    "thr": "T",
    "asp": "D",
    "glu": "E"
}

AA_KEYS = list(AA_CODES.keys())

########################################
# Build peptide
########################################

def build_peptide(seq):

    try:
        seq_letters = "".join([AA_CODES[a] for a in seq])

        mol = Chem.MolFromSequence(seq_letters)

        if mol is None:
            return None

        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)

        return mol

    except:
        return None


########################################
# Exact atom counting
########################################

def atom_count_ok(mol, max_atoms):

    return mol.GetNumAtoms() <= max_atoms


########################################
# Conformer generation
########################################

def embed_and_optimize(mol):

    try:

        params = AllChem.ETKDGv3()

        if AllChem.EmbedMolecule(mol, params) != 0:
            return None, None

        AllChem.UFFOptimizeMolecule(mol, maxIters=500)

        ff = AllChem.UFFGetMoleculeForceField(mol)

        energy = ff.CalcEnergy()

        return mol, energy

    except:
        return None, None


########################################
# Steric clash detection
########################################

def steric_clash(mol):

    conf = mol.GetConformer()
    n = mol.GetNumAtoms()

    for i in range(n):

        p1 = conf.GetAtomPosition(i)

        for j in range(i+1, n):

            if p1.Distance(conf.GetAtomPosition(j)) < 0.8:
                return False

    return True


########################################
# Bond validation
########################################

def bond_validation(mol):

    conf = mol.GetConformer()

    for bond in mol.GetBonds():

        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        d = conf.GetAtomPosition(i).Distance(
            conf.GetAtomPosition(j)
        )

        if d < 0.9 or d > 2.0:
            return False

    return True


########################################
# Ramachandran validation
########################################

def ramachandran_ok(mol):

    try:

        conf = mol.GetConformer()

        atoms = mol.GetAtoms()

        backbone = [
            a.GetIdx() for a in atoms
            if a.GetSymbol() in ["N", "C"]
        ]

        if len(backbone) < 4:
            return True

        phi = rdMolTransforms.GetDihedralDeg(
            conf,
            backbone[0],
            backbone[1],
            backbone[2],
            backbone[3],
        )

        psi = rdMolTransforms.GetDihedralDeg(
            conf,
            backbone[1],
            backbone[2],
            backbone[3],
            backbone[4] if len(backbone) > 4 else backbone[3],
        )

        # Allowed α/β regions
        if (-160 < phi < -40 and -80 < psi < 50) or \
           (-180 < phi < -60 and 90 < psi < 180):

            return True

    except:
        return False

    return False


########################################
# XYZ writer
########################################

def write_xyz(mol, energy, path):

    conf = mol.GetConformer()

    with open(path, "w") as f:

        f.write(f"{mol.GetNumAtoms()}\n")
        f.write(f"Energy {energy:.6f}\n")

        for atom in mol.GetAtoms():

            pos = conf.GetAtomPosition(atom.GetIdx())

            f.write(
                f"{atom.GetSymbol()} "
                f"{pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"
            )


########################################
# Worker generator
########################################

def worker(max_atoms):

    while True:

        length = random.randint(2,6)

        seq = [random.choice(AA_KEYS) for _ in range(length)]

        comp = "-".join(seq)

        mol = build_peptide(seq)

        if mol is None:
            continue

        if not atom_count_ok(mol, max_atoms):
            continue

        mol, energy = embed_and_optimize(mol)

        if mol is None:
            continue

        if not steric_clash(mol):
            continue

        if not bond_validation(mol):
            continue

        if not ramachandran_ok(mol):
            continue

        return comp, mol, energy

########################################
# Estimate maximum possible peptides
########################################

def estimate_max_possible(max_atoms):

    possible = set()

    # explore sequence space up to length 6
    for length in range(2,7):

        for seq in product(AA_KEYS, repeat=length):

            comp = "-".join(seq)

            mol = build_peptide(seq)

            if mol is None:
                continue

            if mol.GetNumAtoms() <= max_atoms:
                possible.add(comp)

    return len(possible)

########################################
# Main
########################################

def main():

    if len(sys.argv) != 3:
        print("Usage: python generate_peptides_hpc.py max_atoms n_peptides")
        sys.exit()
    
    max_atoms = int(sys.argv[1])
    target = int(sys.argv[2])

    max_possible = estimate_max_possible(max_atoms)
    print("\nMaximum possible peptides for this atom limit:", max_possible)

    if target > max_possible:
        print("Requested peptides exceeds maximum possible.")
        print("Capping requested number to:", max_possible)
        target = max_possible

    os.makedirs("xyz_files", exist_ok=True)

    pool = mp.Pool(mp.cpu_count())

    seen = set()

    csv_file = open("peptides.csv","w",newline="")
    writer = csv.writer(csv_file)
    writer.writerow(["peptide_id","composition"])

    count = 0

    while count < target:

        results = pool.map(worker,[max_atoms]*mp.cpu_count())

        for comp, mol, energy in results:

            if comp in seen:
                continue

            seen.add(comp)

            count += 1

            pid = f"pp_{count}"

            writer.writerow([pid,comp])

            xyz_path = f"xyz_files/{pid}.xyz"

            write_xyz(mol,energy,xyz_path)

            if count % 100 == 0:
                print("Generated:",count)

            if count >= target:
                break

    csv_file.close()

    print("\nGeneration finished")
    print("Peptides generated:",count)
    print("XYZ files folder: xyz_files/")
    print("CSV dataset: peptides.csv")
    print("CPU cores used:",mp.cpu_count())


########################################

if __name__ == "__main__":
    main()
