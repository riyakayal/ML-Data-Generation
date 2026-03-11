
# Author: Riya Kayal
# Created: 05/11/2025

import os
import sys
import random
import multiprocessing as mp
from itertools import product

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms

###############################################
# Amino acids
###############################################

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

# approximate atom counts (with H)
AA_ATOMS = {
    "ala": 13,
    "gly": 10,
    "val": 19,
    "leu": 22,
    "ser": 14,
    "thr": 17,
    "asp": 16,
    "glu": 19,
}

###############################################
# Estimate maximum peptides possible
###############################################

def estimate_max_peptides(max_atoms):

    possible = set()
    max_len = 6

    for L in range(2, max_len + 1):

        for seq in product(AA_KEYS, repeat=L):

            atoms = sum(AA_ATOMS[a] for a in seq) - (L - 1) * 2

            if atoms <= max_atoms:
                comp = "-".join(seq)
                possible.add(comp)

    return len(possible)

###############################################
# Peptide builder
###############################################

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

###############################################
# Conformer generation + energy
###############################################

def embed_and_optimize(mol):

    try:

        params = AllChem.ETKDGv3()

        res = AllChem.EmbedMolecule(mol, params)

        if res != 0:
            return None, None

        AllChem.UFFOptimizeMolecule(mol, maxIters=500)

        ff = AllChem.UFFGetMoleculeForceField(mol)

        energy = ff.CalcEnergy()

        return mol, energy

    except:

        return None, None

###############################################
# Steric clash detection
###############################################

def steric_clash_check(mol):

    conf = mol.GetConformer()
    n = mol.GetNumAtoms()

    for i in range(n):

        p1 = conf.GetAtomPosition(i)

        for j in range(i + 1, n):

            p2 = conf.GetAtomPosition(j)

            if p1.Distance(p2) < 0.8:
                return False

    return True

###############################################
# Bond length validation
###############################################

def bond_length_validation(mol):

    conf = mol.GetConformer()

    for bond in mol.GetBonds():

        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        p1 = conf.GetAtomPosition(i)
        p2 = conf.GetAtomPosition(j)

        d = p1.Distance(p2)

        if d < 0.9 or d > 2.0:
            return False

    return True

###############################################
# Ramachandran validation
###############################################

def ramachandran_validation(mol):

    try:

        conf = mol.GetConformer()

        atoms = mol.GetAtoms()

        backbone = [a.GetIdx() for a in atoms if a.GetSymbol() in ["N", "C"]]

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

        if -180 <= phi <= 180 and -180 <= psi <= 180:
            return True

    except:
        return False

    return False

###############################################
# Write XYZ
###############################################

def write_xyz(mol, energy, fname):

    conf = mol.GetConformer()

    with open(fname, "w") as f:

        f.write(f"{mol.GetNumAtoms()}\n")
        f.write(f"Energy {energy:.6f}\n")

        for atom in mol.GetAtoms():

            pos = conf.GetAtomPosition(atom.GetIdx())

            f.write(
                f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"
            )

###############################################
# Single peptide generation
###############################################

def generate_single(args):

    max_atoms, existing = args

    while True:

        length = random.randint(2, 6)

        seq = [random.choice(AA_KEYS) for _ in range(length)]

        comp = "-".join(seq)

        if comp in existing:
            continue

        mol = build_peptide(seq)

        if mol is None:
            continue

        if mol.GetNumAtoms() > max_atoms:
            continue

        mol, energy = embed_and_optimize(mol)

        if mol is None:
            continue

        if not steric_clash_check(mol):
            continue

        if not bond_length_validation(mol):
            continue

        if not ramachandran_validation(mol):
            continue

        return comp, mol, energy

###############################################
# Main
###############################################

def main():

    if len(sys.argv) != 3:
        print("Usage: python generate_peptides.py max_atoms n_peptides")
        sys.exit()

    max_atoms = int(sys.argv[1])
    n_peptides = int(sys.argv[2])

    max_possible = estimate_max_peptides(max_atoms)

    if n_peptides > max_possible:

        print(f"\nRequested peptides: {n_peptides}")
        print(f"Maximum possible under {max_atoms} atoms: {max_possible}")
        print("Capping requested peptides.\n")

        n_peptides = max_possible

    os.makedirs("xyz_files", exist_ok=True)

    manager = mp.Manager()

    compositions = manager.list()

    pool = mp.Pool(mp.cpu_count())

    results = []

    while len(results) < n_peptides:

        jobs = [
            (max_atoms, list(compositions))
            for _ in range(mp.cpu_count())
        ]

        outputs = pool.map(generate_single, jobs)

        for comp, mol, energy in outputs:

            if comp in compositions:
                continue

            compositions.append(comp)

            pid = f"pp_{len(compositions)}"

            xyz_path = f"xyz_files/{pid}.xyz"

            write_xyz(mol, energy, xyz_path)

            results.append((pid, comp))

            if len(results) >= n_peptides:
                break

    df = pd.DataFrame(results, columns=["peptide_id", "composition"])

    df.to_csv("peptides.csv", index=False)

    print("\n==============================")
    print("Peptide Generation Summary")
    print("==============================")
    print(f"Max atoms allowed : {max_atoms}")
    print(f"Peptides requested: {sys.argv[2]}")
    print(f"Peptides generated: {len(results)}")
    print(f"Maximum possible  : {max_possible}")
    print(f"CPU cores used    : {mp.cpu_count()}")
    print("XYZ folder        : xyz_files/")
    print("CSV file          : peptides.csv")
    print("==============================")

###############################################

if __name__ == "__main__":
    main()
