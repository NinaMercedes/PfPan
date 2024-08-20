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

##


