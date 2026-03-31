#!/usr/bin/env python3

import argparse
import gzip
from tqdm import tqdm

def open_file(filename, mode="rt"):
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode)

def count_lines(filename):
    print("Counting total lines (for progress bar)...")
    with open_file(filename, "rt") as f:
        return sum(1 for _ in f) - 1  # minus header

def main():
    parser = argparse.ArgumentParser(description="Filter GWAS by rsID list (HapMap etc.) with progress bar")
    parser.add_argument("--gwas", required=True, help="GWAS input file (TSV or TSV.gz)")
    parser.add_argument("--rsid-list", required=True, help="File with rsIDs (one per line)")
    parser.add_argument("--rsid-col", required=True, help="Column name in GWAS containing rsIDs")
    parser.add_argument("--out", required=True, help="Output filtered GWAS file (can be .gz)")
    args = parser.parse_args()

    # Load rsID set
    print("Loading rsID list...")
    with open(args.rsid_list) as f:
        rsid_set = set(line.strip().split()[0] for line in f)

    print(f"Loaded {len(rsid_set):,} rsIDs")

    total_lines = count_lines(args.gwas)

    kept = 0
    removed = 0

    with open_file(args.gwas, "rt") as infile, open_file(args.out, "wt") as outfile:

        header = infile.readline().strip().split("\t")

        if args.rsid_col not in header:
            raise ValueError(f"Column '{args.rsid_col}' not found in GWAS header")

        rsid_index = header.index(args.rsid_col)

        outfile.write("\t".join(header) + "\n")

        for line in tqdm(infile, total=total_lines, desc="Filtering variants"):
            fields = line.strip().split("\t")

            if len(fields) <= rsid_index:
                removed += 1
                continue

            if fields[rsid_index] in rsid_set:
                outfile.write(line)
                kept += 1
            else:
                removed += 1

    print("\nFiltering complete")
    print(f"Total variants processed: {total_lines:,}")
    print(f"Kept variants: {kept:,}")
    print(f"Removed variants: {removed:,}")
    print(f"Retention rate: {100 * kept / total_lines:.2f}%")

if __name__ == "__main__":
    main()