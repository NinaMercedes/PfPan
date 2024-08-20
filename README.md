# Pangenomes: The Plasmodium falciparum pangenome
This repository contains code to perform pangenomics research. Currently the code is implemented predominantly using python scripts. 

## Conda environment and installing packages
The first step is to set up a conda environment with all the python libraries and modules required to make and analyse a pangenome. 
```
conda create -n pangenomes
### First we need to install fastq2matrix and some dependencies
git clone https://github.com/LSHTMPathogenSeqLab/fastq2matrix.git
conda install python=3.8 bwa samtools bcftools parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree
cd fastq2matrix
python setup.py install
### Install cactus- you may wish to check installation using the test data provided by the cactus team
cd -
git clone https://github.com/ComparativeGenomicsToolkit/cactus.git --recursive
cd cactus
python3 -m pip install -U setuptools pip wheel
python3 -m pip install -U .
python3 -m pip install -U -r ./toil-requirement.txt
### Instal vg
conda install bioconda::vg
### Install panacus
conda install -c conda-forge -c bioconda panacus
### Clone this repository and move to the pangenomes directory
cd ~
git clone https://github.com/NinaMercedes/pangenomes.git
cd pangenomes
```

## Optional Stage: Retrieving assemblies
Assemblies were retrieved from Malariagen Pf3k database. These can obtained using the code below. 
```
### In the pangenomes directory or navigate to where you would like the assemblies to be stored.
python "./pipeline/wget_pf3k_fasta.py" --file_path https://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/PF3K/ReferenceGenomes_Version1/GENOMES/
```

## Step 1. Make a seq file
The seq file is a two column text file used by cactus pangenome that maps the sample names to their paths. This could be produced manually, one column with sample names and one column with the paths (spaced separately). Here a python script has been written to automate this (anything before the first '.' has been used as the sample name). 
```
### Again we chose to do this in the pangenomes directory.
python "./pipeline/make_seqtxt.py" --file_name pf3k_seq.txt
```

### Step 2. Make a pangenome graph (or graphs!)
After making the seq text file, we are now ready to make the pangenome graph. This github link has further details on how this is erformed as well as some of the downstream analysis (https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/sa_refgraph_hackathon_2023.md#short-read-mapping). Although a pangenome is described as somewhat reference-free, we are still required to input some of these sequences 





