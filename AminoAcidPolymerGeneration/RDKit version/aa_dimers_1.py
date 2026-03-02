
# Author: Riya Kayal
# Created: 10/10/2025

# Generate all possible dimers from the 20 natural amino acids
# FASTA sequence using RDKit. Unlike the basic version, this 
# version is chemically more rigorous. Here, we implement
# N–H···O = 180°
# H···O = 2.0 Å
# Rotational scan around hydrogen bond axis
# Best steric separation 


import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

HBOND_DISTANCE = 2.0
CLASH_THRESHOLD = 1.6  # Å minimum interatomic distance
ROTATION_STEPS = 24    # 15 degree increments

# Natural amino acids
AA_SMILES = {
    "ALA": "NCC(=O)O",
    "ARG": "NC(C(=O)O)CCCNC(N)=N",
    "ASN": "NC(C(=O)O)CC(=O)N",
    "ASP": "NC(C(=O)O)CC(=O)O",
    "CYS": "NC(C(=O)O)CS",
    "GLN": "NC(C(=O)O)CCC(=O)N",
    "GLU": "NC(C(=O)O)CCC(=O)O",
    "GLY": "NCC(=O)O",
    "HIS": "NC(C(=O)O)CC1=CN=CN1",
    "ILE": "NC(C(=O)O)C(C)CC",
    "LEU": "NC(C(=O)O)CC(C)C",
    "LYS": "NC(C(=O)O)CCCCN",
    "MET": "NC(C(=O)O)CCSC",
    "PHE": "NC(C(=O)O)CC1=CC=CC=C1",
    "PRO": "N1CCC(C1)C(=O)O",
    "SER": "NC(C(=O)O)CO",
    "THR": "NC(C(=O)O)C(O)C",
    "TRP": "NC(C(=O)O)CC1=CNC2=CC=CC=C12",
    "TYR": "NC(C(=O)O)CC1=CC=C(O)C=C1",
    "VAL": "NC(C(=O)O)C(C)C"
}


def embed(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def coords(mol):
    conf = mol.GetConformer()
    return np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])


def set_coords(mol, arr):
    conf = mol.GetConformer()
    for i, p in enumerate(arr):
        conf.SetAtomPosition(i, p)


def find_N_H_O(mol):
    N_idx = H_idx = O_idx = None

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "N":
            for n in atom.GetNeighbors():
                if n.GetSymbol() == "H":
                    N_idx = atom.GetIdx()
                    H_idx = n.GetIdx()
                    break

        if atom.GetSymbol() == "O":
            if atom.GetDegree() == 1:
                O_idx = atom.GetIdx()

    return N_idx, H_idx, O_idx


def align_linear_hbond(m1, m2):
    c1 = coords(m1)
    c2 = coords(m2)

    N_idx, H_idx, _ = find_N_H_O(m1)
    _, _, O_idx = find_N_H_O(m2)

    N = c1[N_idx]
    H = c1[H_idx]
    O = c2[O_idx]

    # N-H vector
    NH = H - N
    NH = NH / np.linalg.norm(NH)

    # Target O position for 180° geometry
    target_O = H + NH * HBOND_DISTANCE

    translation = target_O - O
    c2 += translation

    set_coords(m2, c2)

    return m1, m2


def rotate_about_axis(points, axis_point, axis_dir, angle):
    axis_dir = axis_dir / np.linalg.norm(axis_dir)
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)

    R = np.array([
        [cos_a + axis_dir[0]**2*(1-cos_a),
         axis_dir[0]*axis_dir[1]*(1-cos_a) - axis_dir[2]*sin_a,
         axis_dir[0]*axis_dir[2]*(1-cos_a) + axis_dir[1]*sin_a],

        [axis_dir[1]*axis_dir[0]*(1-cos_a) + axis_dir[2]*sin_a,
         cos_a + axis_dir[1]**2*(1-cos_a),
         axis_dir[1]*axis_dir[2]*(1-cos_a) - axis_dir[0]*sin_a],

        [axis_dir[2]*axis_dir[0]*(1-cos_a) - axis_dir[1]*sin_a,
         axis_dir[2]*axis_dir[1]*(1-cos_a) + axis_dir[0]*sin_a,
         cos_a + axis_dir[2]**2*(1-cos_a)]
    ])

    shifted = points - axis_point
    rotated = shifted @ R.T
    return rotated + axis_point


def steric_score(c1, c2):
    min_dist = np.inf
    for p1 in c1:
        for p2 in c2:
            d = np.linalg.norm(p1 - p2)
            min_dist = min(min_dist, d)
    return min_dist


def optimize_rotation(m1, m2):
    c1 = coords(m1)
    c2 = coords(m2)

    _, H_idx, _ = find_N_H_O(m1)
    axis_point = c1[H_idx]
    axis_dir = c1[H_idx] - c1[find_N_H_O(m1)[0]]
    axis_dir = axis_dir / np.linalg.norm(axis_dir)

    best_score = -np.inf
    best_coords = c2.copy()

    for i in range(ROTATION_STEPS):
        angle = 2 * np.pi * i / ROTATION_STEPS
        rotated = rotate_about_axis(c2, axis_point, axis_dir, angle)

        score = steric_score(c1, rotated)

        if score > best_score:
            best_score = score
            best_coords = rotated

    set_coords(m2, best_coords)

    return m1, m2


def combine_and_save(m1, m2, filename):
    combo = Chem.CombineMols(m1, m2)
    conf = combo.GetConformer()

    with open(filename, "w") as f:
        f.write(f"{combo.GetNumAtoms()}\n\n")
        for atom in combo.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n")


def generate_all(output_dir="dimers_1"):
    os.makedirs(output_dir, exist_ok=True)

    for name1, smi1 in AA_SMILES.items():
        for name2, smi2 in AA_SMILES.items():

            m1 = embed(smi1)
            m2 = embed(smi2)

            m1, m2 = align_linear_hbond(m1, m2)
            m1, m2 = optimize_rotation(m1, m2)

            filename = os.path.join(output_dir, f"{name1}_{name2}.xyz")
            combine_and_save(m1, m2, filename)

            print(f"Generated {name1}-{name2}")


if __name__ == "__main__":
    generate_all()
