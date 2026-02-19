#!/usr/bin/env python3
"""
Download IMGT/HLA genomic sequences and build a bgzipped FASTA reference
with pipeline-compatible headers.

Header transformation:
  >HLA:HLA00001 A*01:01:01:01 3503 bp  →  >A_01_01_01_01

Output:
  assets/imgt_hla_gen.fasta.gz   (bgzipped)
  assets/imgt_hla_gen.fasta.gz.fai
  assets/imgt_hla_gen.fasta.gz.gzi
"""

import argparse
import io
import os
import re
import subprocess
import sys
import tempfile
import zipfile
from urllib.request import urlopen

IMGT_URL = (
    "https://github.com/ANHIG/IMGTHLA/raw/Latest/fasta/hla_gen.fasta.zip"
)


def parse_imgt_header(header_line):
    """Extract allele name from IMGT header and convert to pipeline format.

    Example:
        '>HLA:HLA00001 A*01:01:01:01 3503 bp'
        → 'A_01_01_01_01'
    """
    # Second field is the allele name, e.g. 'A*01:01:01:01'
    parts = header_line.lstrip(">").split()
    if len(parts) < 2:
        return None
    allele = parts[1]
    # Strip 'HLA-' prefix if present (e.g. 'HLA-DRB1*01:01:01:01')
    allele = re.sub(r"^HLA-", "", allele)
    # Replace * and : with _
    allele = allele.replace("*", "_").replace(":", "_")
    return allele


def extract_from_zip(zip_path):
    """Extract FASTA content from a local or downloaded zip file."""
    with zipfile.ZipFile(zip_path) as zf:
        fasta_names = [n for n in zf.namelist() if n.endswith(".fasta")]
        if not fasta_names:
            sys.exit("ERROR: No .fasta file found in zip archive")
        content = zf.read(fasta_names[0]).decode("utf-8")
    return content


def download_and_extract(url):
    """Download zip from URL and return the FASTA content as a string."""
    print(f"Downloading {url} ...")
    response = urlopen(url)
    zip_bytes = response.read()
    print(f"Downloaded {len(zip_bytes)} bytes")
    return extract_from_zip(io.BytesIO(zip_bytes))


def reformat_fasta(raw_fasta):
    """Reformat IMGT FASTA to pipeline-compatible headers.

    Returns list of (header, sequence) tuples.
    """
    records = []
    current_header = None
    seq_lines = []

    for line in raw_fasta.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(seq_lines)))
            current_header = parse_imgt_header(line)
            seq_lines = []
        else:
            seq_lines.append(line)

    if current_header is not None:
        records.append((current_header, "".join(seq_lines)))

    return records


def main():
    parser = argparse.ArgumentParser(
        description="Build IMGT/HLA reference FASTA for hlasomatic pipeline"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default=os.path.join(os.path.dirname(__file__), "..", "assets"),
        help="Output directory (default: ../assets relative to this script)",
    )
    parser.add_argument(
        "--url",
        default=IMGT_URL,
        help="URL to IMGT hla_gen.fasta.zip",
    )
    parser.add_argument(
        "--local-zip",
        default=None,
        help="Path to a pre-downloaded hla_gen.fasta.zip (skips download)",
    )
    args = parser.parse_args()

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    out_gz = os.path.join(outdir, "imgt_hla_gen.fasta.gz")

    # Load FASTA from local zip or download
    if args.local_zip:
        print(f"Reading local zip: {args.local_zip}")
        raw_fasta = extract_from_zip(args.local_zip)
    else:
        raw_fasta = download_and_extract(args.url)
    records = reformat_fasta(raw_fasta)
    print(f"Parsed {len(records)} allele records")

    # Collect gene names for summary
    genes = sorted(set(r[0].split("_")[0] for r in records))
    print(f"Genes present ({len(genes)}): {', '.join(genes)}")

    # Write to a temp FASTA, then bgzip + index
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, dir=outdir
    ) as tmp:
        tmp_path = tmp.name
        for header, seq in records:
            tmp.write(f">{header}\n")
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                tmp.write(seq[i : i + 80] + "\n")

    print(f"Wrote {len(records)} records to {tmp_path}")

    # bgzip
    print("Running bgzip ...")
    subprocess.run(["bgzip", "-c", tmp_path], stdout=open(out_gz, "wb"), check=True)
    os.unlink(tmp_path)

    # samtools faidx
    print("Running samtools faidx ...")
    subprocess.run(["samtools", "faidx", out_gz], check=True)

    print(f"\nOutput files:")
    print(f"  {out_gz}")
    print(f"  {out_gz}.fai")
    print(f"  {out_gz}.gzi")


if __name__ == "__main__":
    main()
