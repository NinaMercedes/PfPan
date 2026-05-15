export OUTDIR="sv_results"
export MIN_SVLEN=2000

bash PfPan_sv_pipeline.sh \
  -v "/mnt/storage13/nbillows/pangenome/final_no_mix/analyse/PfPan.norm.svinfo.vcf.gz" \
  -b Core_genome_Pf3D7_v3_ext.bed \
  -r Pfalciparum.genome.fasta \
  -g Pfalciparum.genome.modified.new.gff3 \
  -o sv_results \
  -t 16 \
  -m 2000