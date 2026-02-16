#!/usr/bin/env python3
"""
Convert HLAHD typing output (final.result.txt) to Hapster haplotype CSV format.

HLAHD format (tab-separated):
    A    HLA-A*02:01:01:01    HLA-A*11:01:01:01

Hapster format (headerless CSV, one allele per line):
    A_02_01_01_01
    A_11_01_01_01
"""

import argparse
import os
import re
import sys


def parse_hlahd_allele(allele_str):
    """
    Parse an HLAHD allele string like 'HLA-A*02:01:01:01' into Hapster format 'A_02_01_01_01'.

    Returns None if the allele is 'Not typed' or '-'.
    """
    allele_str = allele_str.strip()
    if allele_str in ("Not typed", "-", ""):
        return None

    # Match pattern: HLA-X*FF:FF:FF:FF (2-4 field resolution)
    match = re.match(r"HLA-([A-Z]+\d?)\*(\d+):(\d+)(?::(\d+))?(?::(\d+))?", allele_str)
    if not match:
        print(f"WARNING: Could not parse allele: {allele_str}", file=sys.stderr)
        return None

    gene = match.group(1)
    fields = [match.group(i) for i in range(2, 6) if match.group(i) is not None]
    return f"{gene}_{'_'.join(fields)}"


def parse_hlahd_results(hlahd_file, genes):
    """
    Parse HLAHD final.result.txt and return a dict of gene -> [allele1, allele2].
    """
    results = {}
    with open(hlahd_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue

            gene = parts[0].strip()
            if gene not in genes:
                continue

            allele1 = parse_hlahd_allele(parts[1])
            allele2 = parse_hlahd_allele(parts[2])

            if allele1 is None and allele2 is None:
                print(f"WARNING: No valid alleles for gene {gene}", file=sys.stderr)
                continue

            # Handle homozygous: duplicate the single allele
            if allele1 is None:
                allele1 = allele2
            if allele2 is None:
                allele2 = allele1

            # Sort alphabetically (Hapster convention)
            alleles = sorted([allele1, allele2])
            results[gene] = alleles

    return results


def write_hapster_csv(results, output_dir):
    """
    Write per-gene CSVs and a consolidated CSV in Hapster format.
    """
    os.makedirs(output_dir, exist_ok=True)

    all_alleles = []

    for gene in sorted(results.keys()):
        alleles = results[gene]

        # Per-gene CSV
        gene_file = os.path.join(output_dir, f"{gene}_haplotypes.csv")
        with open(gene_file, "w") as f:
            for allele in alleles:
                f.write(f"{allele}\n")

        all_alleles.extend(alleles)

    # Consolidated CSV
    consolidated_file = os.path.join(output_dir, "haplotypes.csv")
    with open(consolidated_file, "w") as f:
        for allele in all_alleles:
            f.write(f"{allele}\n")

    return consolidated_file


def main():
    parser = argparse.ArgumentParser(
        description="Convert HLAHD typing output to Hapster haplotype CSV format"
    )
    parser.add_argument(
        "--hlahd_results",
        required=True,
        help="Path to HLAHD final.result.txt file",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for Hapster haplotype CSVs",
    )
    parser.add_argument(
        "--genes",
        default="A,B,C",
        help="Comma-separated list of HLA genes to convert (default: A,B,C)",
    )
    args = parser.parse_args()

    genes = [g.strip() for g in args.genes.split(",")]

    results = parse_hlahd_results(args.hlahd_results, genes)

    if not results:
        print("ERROR: No valid HLA alleles found in input", file=sys.stderr)
        sys.exit(1)

    consolidated = write_hapster_csv(results, args.output_dir)

    print(f"Converted {len(results)} genes to Hapster format")
    for gene, alleles in sorted(results.items()):
        print(f"  {gene}: {', '.join(alleles)}")
    print(f"Consolidated file: {consolidated}")


if __name__ == "__main__":
    main()
