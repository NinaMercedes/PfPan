# Pangenomes: The Plasmodium falciparum pangenome
This repository contains code to perform pangenomics research. Currently the code is implemented predominantly through python. The first step is to set up a conda environment with all the python libraries and modules required to make and analyse a pangenome. 

```
conda create -n pangenomes
### First we need to install fastq2matrix and some dependencies
git clone https://github.com/LSHTMPathogenSeqLab/fastq2matrix.git
conda install python=3.8 bwa samtools bcftools parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree
cd fastq2matrix
python setup.py install
### Install cactus
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
```
## Optional Stage: Retrieving assemblies
Assemblies were retrieved from Malariagen Pf3k database. These 



