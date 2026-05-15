#!/usr/bin/env python
import sys
import os
import argparse
import subprocess

# Taking list of fasta files and making the sequence file list required for cactus-pangenome

def write_seq(file_name):
    ls_cmd = "ls *.fasta.gz > fasta_list.txt"
    subprocess.run(ls_cmd, shell=True, check=True)

    newfile = open(file_name, 'w')
    fasta_list = open('fasta_list.txt', 'r')
    lines = fasta_list.readlines()
    for line in lines:
        strain = line.strip().split('.')[0]
        path = os.getcwd() + "/" + line.strip()
        newfile.write("{} {}\n".format(strain, path))
    newfile.close()

def main(args):
    write_seq(args.file_name)

parser = argparse.ArgumentParser(
    description='Get seq file for cactus pangenome',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument('--file_name', help='file name to output for cactus pangenome seq file', required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)