
# Author: Riya Kayal
# Created: 10/10/2025

# Generate all possible dimers from the 20 natural amino acids
# FASTA sequence using RDKit.

import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# 20 natural amino acids (side chain simplified SMILES)
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

HBOND_DISTANCE = 2.0  # Å


def embed_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def get_atom_coords(mol):
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return coords


def set_atom_coords(mol, coords):
    conf = mol.GetConformer()
    for i, pos in enumerate(coords):
        conf.SetAtomPosition(i, pos)


def find_backbone_atoms(mol):
    """
    Returns:
    N atom index
    One H attached to N
    Carbonyl oxygen index
    """
    N_idx = None
    H_idx = None
    O_idx = None

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "N":
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == "H":
                    N_idx = atom.GetIdx()
                    H_idx = neighbor.GetIdx()
                    break
        if atom.GetSymbol() == "O":
            if atom.GetDegree() == 1:  # likely carbonyl O
                O_idx = atom.GetIdx()

    return N_idx, H_idx, O_idx


def align_for_hbond(mol1, mol2):
    """
    Move mol2 such that:
    H (from mol1 NH2) is at HBOND_DISTANCE from carbonyl O of mol2
    and aligned linearly
    """

    coords1 = get_atom_coords(mol1)
    coords2 = get_atom_coords(mol2)

    N_idx, H_idx, _ = find_backbone_atoms(mol1)
    _, _, O_idx = find_backbone_atoms(mol2)

    H_pos = coords1[H_idx]
    N_pos = coords1[N_idx]
    O_pos = coords2[O_idx]

    # Direction of N-H bond
    NH_vec = H_pos - N_pos
    NH_vec = NH_vec / np.linalg.norm(NH_vec)

    # Desired O position = H + HBOND_DISTANCE * NH direction
    target_O = H_pos + NH_vec * HBOND_DISTANCE

    # Translation vector
    translation = target_O - O_pos

    coords2 += translation

    set_atom_coords(mol2, coords2)

    return mol1, mol2


def combine_molecules(mol1, mol2):
    combo = Chem.CombineMols(mol1, mol2)
    combo = Chem.AddHs(combo, addCoords=True)
    return combo


def save_xyz(mol, filename):
    conf = mol.GetConformer()
    with open(filename, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n")


def generate_all_dimers(output_dir="dimers_0"):
    os.makedirs(output_dir, exist_ok=True)

    for aa1_name, aa1_smiles in AA_SMILES.items():
        for aa2_name, aa2_smiles in AA_SMILES.items():

            mol1 = embed_molecule(aa1_smiles)
            mol2 = embed_molecule(aa2_smiles)

            mol1, mol2 = align_for_hbond(mol1, mol2)
            combo = combine_molecules(mol1, mol2)

            filename = os.path.join(output_dir, f"{aa1_name}_{aa2_name}.xyz")
            save_xyz(combo, filename)

            print(f"Generated {aa1_name}-{aa2_name}")


if __name__ == "__main__":
    generate_all_dimers()
