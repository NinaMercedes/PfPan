# Pangenomes: The *Plasmodium falciparum* pangenome
This repository contains code to perform pangenomics research. Currently the code is implemented predominantly using python and shell scripts. Depending on usability, these scripts may be updated in the future and built into Nextflow or Snakemake pipelines (suggestions welcome!).

## Conda environment and installing packages
The first step is to set up a conda environment with all the python libraries and modules required to make and analyse a pangenome. 
```
conda create -n cactus python=3.8
conda activate pangenomes
### First we need to install fastq2matrix and some dependencies
#git clone https://github.com/LSHTMPathogenSeqLab/fastq2matrix.git
conda config --set channel_priority flexible    
conda install -c bioconda bwa samtools bcftools parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree
cd fastq2matrix
python setup.py install
### Install cactus- you may wish to check installation using the test data provided by the cactus team
python3 -m pip install virtualenv
### please see these notes here to install from bianry.
https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/BIN-INSTALL.md
### Instal vg
conda install bioconda::vg
### Install panacus
conda install bioconda::panacus
conda install matplotlib
### Clone this repository and move to the pangenomes directory
cd ~
git clone https://github.com/NinaMercedes/pangenomes.git
cd pangenomes
```

## Optional Stage: Retrieving assemblies
Assemblies were retrieved from Malariagen Pf3k database. These can obtained using the code below. 
```
### In the pangenomes directory or navigate to where you would like the assemblies to be stored.
mkdir genomes
cd genomes
python "./pipeline/wget_pf3k_fasta.py" --file_path https://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/PF3K/ReferenceGenomes_Version1/GENOMES/
```

## Step 1. Make a seq file
The seq file is a two column text file used by cactus pangenome that maps the sample names to their paths. This could be produced manually, one column with sample names and one column with the paths (spaced separately). Here a python script has been written to automate this (anything before the first '.' has been used as the sample name). 
```
### Again we chose to do this in the pangenomes directory.
python "./pipeline/make_seqtxt.py" --file_name pf3k_seq.txt
```

## Step 2. Make a pangenome graph (or graphs!)
After making the seq text file, we are now ready to make the pangenome graph. This github link has further details on how this is performed as well as some of the downstream analysis (https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/sa_refgraph_hackathon_2023.md#short-read-mapping). Although a pangenome is described as somewhat reference-free, we are still required to input a reference sequence to project variants on. This can help us to interpret variants called or genotyped from the graph. Minigraph-Cactus uses the reference sequence as a backbone and the choice of reference can influence the topology of the graph. Minigraph-Cactus does permit the use of more than one reference. We can explore the use of both Pf3D7 (West Africa origin) and the PfDd2 (Indochina origin) reference strains to assess the impact. Note in version 2 we exclude SD01, ML01 and TG01 due to contig size inconsistencies and mixed infections. 
```
source "/mnt/storage13/nbillows/pangenome/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate" 
cactus-pangenome ./PfPan_Pf3D7_js pf3k_seq_v2.txt --outDir ./PfPan_Pf3D7_pan  --outName PfPan --reference Pf3D7  --filter 2  --haplo --giraffe clip filter --viz --odgi --chrom-vg clip filter --chrom-og --gbz clip filter full   --gfa clip full --vcf --vcfReference Pf3D7 --logFile ./PfPan_Pf3D7_log.log --workDir ./  --consCores 8 --mgMemory  128Gi
source deactivate
```

# Cactus Pangenome Pipeline Outputs

## Command

```bash
cactus-pangenome ./PfPan_Pf3D7_js pf3k_seq_v2.txt \
  --outDir ./PfPan_Pf3D7_pan \
  --outName PfPan \
  --reference Pf3D7 \
  --filter 2 \
  --haplo \
  --giraffe clip filter \
  --viz --odgi \
  --chrom-vg clip filter \
  --chrom-og \
  --gbz clip filter full \
  --gfa clip full \
  --vcf --vcfReference Pf3D7 \
  --logFile ./PfPan_Pf3D7_log.log \
  --workDir ./ \
  --consCores 8 \
  --mgMemory 128Gi \
  chrom-alignments
```

---

## Output Files

| File | Format | Description |
|------|--------|-------------|
| `PfPan.gbz` | GBZ | Clipped pangenome graph (primary, used for most downstream tools) |
| `PfPan.full.gbz` | GBZ | Full (unclipped) pangenome graph in GBZ format |
| `PfPan.d2.gbz` | GBZ | Decorated/distance-indexed graph for Giraffe mapping |
| `PfPan.gfa.gz` | GFA | Clipped pangenome graph in GFA format |
| `PfPan.full.gfa.gz` | GFA | Full pangenome graph in GFA format |
| `PfPan.sv.gfa.gz` | GFA | Structural variant subgraph in GFA format |
| `PfPan.sv.gfa.fa.gz` | FASTA | Sequences associated with the SV subgraph |
| `PfPan.full.hal` | HAL | Full whole-genome alignment in HAL format (Cactus native) |
| `PfPan.full.og` | OG | Full graph in ODGI format |
| `PfPan.full.snarls` | Snarls | Snarl decomposition of the full graph |
| `PfPan.snarls` | Snarls | Snarl decomposition of the clipped graph |
| `PfPan.dist` | Dist | Distance index for the clipped graph |
| `PfPan.d2.dist` | Dist | Distance index for the Giraffe-ready graph |
| `PfPan.d2.min` | Min | Minimiser index for Giraffe read mapping |
| `PfPan.d2.snarls` | Snarls | Snarl index for the Giraffe-ready graph |
| `PfPan.min` | Min | Minimiser index for the clipped graph |
| `PfPan.hapl` | Hapl | Haplotype index for `vg haplotypes` / personalised graph construction |
| `PfPan.gaf.gz` | GAF | Read-to-graph alignments (graph alignment format) |
| `PfPan.paf` | PAF | Pairwise sequence alignments of input assemblies to the graph |
| `PfPan.paf.unfiltered.gz` | PAF | Unfiltered version of the pairwise alignments |
| `PfPan.paf.filter.log` | Log | Log of filtering applied to PAF alignments |
| `PfPan.vcf.gz` | VCF | Filtered variant calls against `Pf3D7` reference |
| `PfPan.vcf.gz.tbi` | TBI | Tabix index for `PfPan.vcf.gz` |
| `PfPan.raw.vcf.gz` | VCF | Raw (unfiltered) variant calls |
| `PfPan.raw.vcf.gz.tbi` | TBI | Tabix index for `PfPan.raw.vcf.gz` |
| `PfPan.viz` | Dir/PNG | Visual representation of the pangenome graph (from `--viz`) |
| `PfPan.stats.tgz` | TGZ | Archive of graph statistics |
| `PfPan.chroms` | Dir | Per-chromosome graph files |
| `chrom-alignments/` | Dir | Per-chromosome alignment files |
| `chrom-subproblems/` | Dir | Intermediate files from per-chromosome cactus alignment |
| `pf3k_seq_v2.txt` | TXT | Input sequence file listing assemblies |


## Step 3. Get pangenome statistics
#######EDIT
We can get pangenome statistics from the use of the cactus-mingraph (e.g. "pf_pan_v1.stats.tgz"), this includes sample, paths and graph statistics. Using the code below we will generate a '...basic_statistics.txt' file which will give us the number of nodes, edges and total sequence length. We will also look at the contig sizes to check the length of each haplotype in the chromosome "./pf_pan_v1_Pf3D7_pan/chrom-subproblems/contig_sizes.tsv". We can see there is a contig size of '0' for sample PfSD01- this sample was actually assembled using a different assembler so we can exclude this from our pangenome and make a new set (v1b and v2b) (Otto et al., 2018). We do this by rerunning the above but change the same of the prefix and we also manually remove PfSD01 from the seq file. Panacus is also a useful package to explore pangenome statistics, including a coverage histogram, pangenome growth statistics and path-/group-resolved coverage table. Panacus relies on 'countable' features including nodes, edges and base pairs. Coverage is defined as the number of distinct paths including the countable. Meanwhile, 'quorum' is the proportion of paths in which one of the countable features needs to be present to be considered part of the core. Essentially using this tool we can estimate the number of 'common' and 'core' bases in the pangenome and look at how the addition of samples affects its growth using growth curves and statistics. The following code should generate a '*.basic_statistics.txt' and 'histgrowth' node tsv and pdf files.

```
###NOOOO
"/mnt/storage13/nbillows/pangenome/final_no_mix/PfPan_Pf3D7_pan/"
python ~/pangenomes/pipeline/pangenome_stats.py --prefix PfPan --ref_name Pf3D7 --graph_gfa PfPan.gfa.gz --graph_gbz PfPan.gbz

```




## Step 4. Mapping short reads, make gaf and genotype
This is quite computationally intensive and requires the input of fastq files. Within the python script the short read fastq files are generated from bam/cram files that have been stored on the server (due to space requirements). This process is optional. However, if running straight from fastq, please use the naming format: 'sample_name.R1.fastq.gz and sample_name.R2.fastq.gz". There is also an option to remove fastq files (default: False) to save space. This is not reccomended if you would like to continue using your fastq files.

```
cd "/mnt/storage13/nbillows/pangenome/genomes/pf_pan_2024_Pf3D7_pan/"
mkdir output
mkdir fastqs

#### make BAM
## get paths from full graph
#vg paths -x ./pf_pan_2024.full.gbz -S Pf3D7 -L > Pf3D7.paths.txt
### surject to bam using d2 graph
vg surject -x ./pf_pan_2024.full.gbz -G ./output/ERR012227.gaf.gz --interleaved -F Pf3D7.paths.txt -b -N ERR012227 -R 'ID:1 LB:lib1 SM:./output/ERR012227 PL:illumina PU:unit1' > ./output/ERR012227.bam
#### change minimiser
vg minimizer -p -t 16 -k 11 -w 5 -d pf_pan_2024.dist -g pf_pan_2024.gbwt -o pf_pan_2024.min pf_pan_2024.gbz

#vg paths -x  ./pf_pan_2024.gbz z -S Pf3D7 -F > ./Pf3D7_graph.fa
#samtools faidx ./Pf3D7_graph.fa
#gatk CreateSequenceDictionary --REFERENCE Pf3D7_graph.fa


```
The BAM files produced here can be utilised for de novo variant calling (different to genotyping) using linear variant callers such as GATK, DeepVariant (Docker required, currently trained on human data), or your variant caller of choice! Could be more to come on this (watch this space!). 

### Optional Stage: Surject to BAM- illumina only


