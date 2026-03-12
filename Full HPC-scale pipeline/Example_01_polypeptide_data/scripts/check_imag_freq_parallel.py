#!/usr/bin/env python3


# Author: Riya Kayal
# Created: 04/03/2025


## This a parallelized ORCA imaginary-frequency checker designed for thousands (or tens of thousands) of output files.

## Recursively scans a target directory
## Uses multiprocessing (all CPU cores by default)
## Robustly parses ORCA frequency blocks
## Ignores small numerical noise (default −10 cm⁻¹)

## Produces:
# Full CSV report
# Separate file lists (OK / IMAGINARY / NO_FREQ)
# Summary timing + throughput statistics

## Usage: python check_imag_freq_parallel.py /path/to/orca_outputs
## OR, Optional: specify number of cores
## python check_imag_freq_parallel.py /path/to/orca_outputs 32


import os
import re
import sys
import csv
import time
from multiprocessing import Pool, cpu_count

#THRESHOLD = -0.99   # for tighter setting, ORCA does not write imag mode if freq > 0.99 
THRESHOLD = -10.00
REPORT = "imaginary_frequency_report.csv"


def extract_frequencies(filepath):
    freqs = []
    in_block = False

    try:
        with open(filepath, "r", errors="ignore") as f:
            for line in f:
                if "VIBRATIONAL FREQUENCIES" in line:
                    in_block = True
                    continue

                if in_block:
                    if "cm" in line:
                        numbers = re.findall(r"-?\d+\.\d+", line)
                        freqs.extend(float(n) for n in numbers)
    except Exception:
        return filepath, None, "READ_ERROR"

    if not freqs:
        return filepath, None, "NO_FREQ_SECTION"

    minfreq = min(freqs)

    if minfreq < THRESHOLD:
        return filepath, minfreq, "IMAGINARY"
    else:
        return filepath, minfreq, "OK"


def collect_out_files(directory):
    out_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".out"):
                out_files.append(os.path.join(root, file))
    return out_files


def main():

    if len(sys.argv) < 2:
        print("Usage: python check_imag_freq_parallel.py <output_directory> [n_cores]")
        sys.exit(1)

    target_dir = sys.argv[1]

    if not os.path.isdir(target_dir):
        print(f"Error: {target_dir} is not a valid directory.")
        sys.exit(1)

    if len(sys.argv) == 3:
        ncores = int(sys.argv[2])
    else:
        ncores = cpu_count()

    print(f"Scanning directory: {target_dir}")
    print(f"Using {ncores} CPU cores")

    start_time = time.time()

    out_files = collect_out_files(target_dir)
    total = len(out_files)

    print(f"Total .out files found: {total}")

    with Pool(ncores) as pool:
        results = pool.map(extract_frequencies, out_files)

    passed = failed = nofreq = readerr = 0

    with open(REPORT, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["File", "Min_Frequency_cm-1", "Status"])

        for file, minfreq, status in results:
            if status == "OK":
                passed += 1
                writer.writerow([file, f"{minfreq:.4f}", status])
            elif status == "IMAGINARY":
                failed += 1
                writer.writerow([file, f"{minfreq:.4f}", status])
            elif status == "NO_FREQ_SECTION":
                nofreq += 1
                writer.writerow([file, "NA", status])
            else:
                readerr += 1
                writer.writerow([file, "NA", status])

    # Write separate job lists
    with open("failed_imaginary.txt", "w") as f:
        for file, _, status in results:
            if status == "IMAGINARY":
                f.write(file + "\n")

    with open("no_freq_section.txt", "w") as f:
        for file, _, status in results:
            if status == "NO_FREQ_SECTION":
                f.write(file + "\n")

    end_time = time.time()
    elapsed = end_time - start_time

    print("\n----------------------------------------------------------------------------------------")
    print(f"Total files           : {total}")
    print(f"Passed                : {passed}")
    print(f"Imaginary             : {failed}")
    print(f"No freq section       : {nofreq}")
    print(f"Read errors           : {readerr}")
    print(f"Time elapsed (s)      : {elapsed:.2f}")
    if elapsed > 0:
        print(f"Files per second      : {total/elapsed:.2f}")
    print(f"Report file           : {os.path.abspath(REPORT)}")
    print("-----------------------------------------------------------------------------------------\n")


if __name__ == "__main__":
    main()
