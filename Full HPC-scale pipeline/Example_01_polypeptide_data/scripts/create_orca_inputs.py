#!/usr/bin/env python3


# Author: Riya Kayal
# Created: 05/03/2025


import os
import glob

# === USER SETTINGS ===
XYZ_DIR = "xyz_files"
ORCA_INPUT_DIR = "orca_opt"
NPROCS = 8           # number of cores per ORCA job
DEFAULT_CHARGE = 0
DEFAULT_MULTIPLICITY = 1

os.makedirs(ORCA_INPUT_DIR, exist_ok=True)


def read_xyz_coordinates(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    natoms = int(lines[0].strip())
    coords = lines[2:2 + natoms]  # skip first two lines
    return coords


def generate_orca_input(xyz_path):
    filename = os.path.basename(xyz_path)
    name = os.path.splitext(filename)[0]

    coords = read_xyz_coordinates(xyz_path)

    inp_path = os.path.join(ORCA_INPUT_DIR, f"{name}.inp")

    with open(inp_path, "w") as f:
        ## write input here
        #f.write("! B3LYP D3BJ def2-SVP TightSCF TightOpt defGrid4 Freq\n\n")
        f.write("! r2scan-3c TightOpt  Freq\n\n")
        f.write("%maxcore 2000\n")  # memory per core in MB
        f.write(f"%pal nprocs {NPROCS} end\n\n")
        f.write(f"* xyz {DEFAULT_CHARGE} {DEFAULT_MULTIPLICITY}\n")

        for line in coords:
            f.write(line)

        f.write("*\n")

    print(f"Generated {inp_path}")


def main():
    xyz_files = glob.glob(os.path.join(XYZ_DIR, "*.xyz"))

    if not xyz_files:
        print("No XYZ files found.")
        return

    print(f"Found {len(xyz_files)} XYZ files.")
    print("Generating ORCA input files...\n")

    for xyz in xyz_files:
        generate_orca_input(xyz)

    print("\nDone.")


if __name__ == "__main__":
    main()
