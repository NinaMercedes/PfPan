#!/usr/bin/env python
import sys
import os
import argparse
import subprocess

# Taking Pacbio reference, lab strains and clinical samples from Pf3K:
# https://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/PF3K/ReferenceGenomes_Version1/GENOMES/
# These samples are long-read assemblies from pacbio data
# Run in the directory you need the files

def wget_files(file_path):
    # Parse the index page for .fasta.gz hrefs, then download each directly
    parse_cmd = (
        f"wget -q -O - {file_path} | "
        "grep -o 'href=\"[^\"]*\\.fasta\\.gz\"' | "
        "grep -v '\\.fai' | "
        "sed 's/href=\"//;s/\"//'"
    )
    result = subprocess.run(parse_cmd, shell=True, check=True, capture_output=True, text=True)
    filenames = result.stdout.strip().split("\n")

    if not filenames or filenames == [""]:
        raise RuntimeError("No .fasta.gz files found at the provided URL.")

    for fname in filenames:
        url = file_path.rstrip("/") + "/" + fname
        print(f"Downloading: {url}")
        subprocess.run(f"wget -nd {url}", shell=True, check=True)

def main(args):
    wget_files(args.file_path)

parser = argparse.ArgumentParser(
    description='Get fasta files from ftp',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument('--file_path', help='ftp file_path', required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)