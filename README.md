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

## Step 2. Make a pangenome graph (or graphs!)
After making the seq text file, we are now ready to make the pangenome graph. This github link has further details on how this is erformed as well as some of the downstream analysis (https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/sa_refgraph_hackathon_2023.md#short-read-mapping). Although a pangenome is described as somewhat reference-free, we are still required to input a reference sequence to project variants on. This can help us to interpret variants called or genotyped from the graph. Minigraph-Cactus uses the reference sequence as a backbone and the choice of reference can influence the topology of the graph. Minigraph-Cactus does permit the use of more than one reference. We can explore the use of both Pf3D7 (West Africa origin) and the PfDd2 (Indochina origin) reference strains to assess the impact.
```
python "./pipeline/make_pangenome.py" --prefix pf_pan_v1 --seqfile pf3k_seq.txt --ref_name Pf3D7
python "./pipeline/make_pangenome.py" --prefix pf_pan_v2 --seqfile pf3k_seq.txt --ref_name PfDd2
### New versions without PfSD01 (see below for explanation)
python "./pipeline/make_pangenome.py" --prefix pf_pan_v1b --seqfile pf3k_seq.txt --ref_name Pf3D7
python "./pipeline/make_pangenome.py" --prefix pf_pan_v2b --seqfile pf3k_seq.txt --ref_name PfDd2 

```
## Step 3. Get pangenome statistics
We can get pangenome statistics from the use of the cactus-mingraph (e.g. "pf_pan_v1.stats.tgz"), this includes sample, paths and graph statistics. Using the code below we will generate a 'basic_stats.txt' file which will give us the number of nodes, edges and total sequence length. We will also look at the contig sizes to check the length of each haplotype in the chromosome "./pf_pan_v1_Pf3D7_pan/chrom-subproblems/contig_sizes.tsv". We can see there is a contig size of '0' for sample PfSD01- this sample was actually assembled using a different assembler so we can exclude this from our pangenome and make a new set (v1b and v2b) (Otto et al., 2018). We do this by rerunning the above but change the same of the prefix and we also manually remove PfSD01 from the seq file. Panacus is also a useful package to explore pangenome statistics, including a coverage histogram, pangenome growth statistics and path-/group-resolved coverage table. Panacus relies on 'countable' features including nodes, edges and base pairs. Coverage is defined as the number of distinct paths including the countable. Meanwhile, 'quorum' is the proportion of paths in which one of the countable features needs to be present to be considered part of the core. Essentially using this tool we can estimate the number of 'common' and 'core' bases in the pangenome and look at how the addition of samples affects its growth using growth curves and statistics. This should generate a '*.basic_statistics.txt' and 'histgrowth' node tsv and pdf files.

```
cd pf_pan_v1b_Pf3D7_pan
python ~/pangenomes/pipeline/pangenome_stats.py --prefix pf_pan_v1b --ref_name Pf3D7 --graph_gfa pf_pan_v1b.gfa.gz --graph_gbz pf_pan_v1b.gbz
cd ~/pangenomes/pf_pan_v2b_PfDd2_pan
python ~/pangenomes/pipeline/pangenome_stats.py --prefix pf_pan_v2b --ref_name PfDd2 --graph_gfa pf_pan_v2b.gfa.gz --graph_gbz pf_pan_v2b.gbz

```

