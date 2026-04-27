#! /usr/bin/env python

import sys
import os
import argparse
import subprocess
import itertools
import pandas as pd
import fastq2matrix as fm


# Taking Pacbio reference, lab strains and clinical samples from Pf3K:
# https://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/PF3K/ReferenceGenomes_Version1/GENOMES/
# These samples are long-read assemblies from pac bio data
# Run in the directory you need the files

#Functions
def wget_files(file_path):
    wget_cmd = "wget -nd -np -r " + file_path
    fm.run_cmd(wget_cmd)
    

def main(args):
  wget_files(args.file_path)

parser = argparse.ArgumentParser(description='Get fasta files from ftp',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--file_path',help='ftp file_path',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)