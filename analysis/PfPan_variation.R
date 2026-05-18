############################################################
# PfPan_variation.R
# Variation analysis of the PfPan pangenome
#
# Sections:
#   1.  Load data
#   2.  Variant processing and classification
#   3.  GRanges builder
#   4.  Optional BED filter
#   5.  Summary exports
#   6.  Core variant plots
#   7.  GOI analysis
#   8.  HRP analysis
#   9.  Complex loci analysis
#   10. Load genomic region BED files
#   11. Map variants to regions
#   12. Regional plots
#   13. SV genomic distribution (chromosomal density)
#   14. SV fold enrichment
#   15. SV enrichment BED export
############################################################

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(scales)

############################################################
# CONFIGURATION
############################################################

VARIANTS_FILE <- "PfPan.all_variants_info_GT.tsv"
SAMPLES_FILE  <- "samples.txt"
GFF_FILE      <- "Pfalciparum.genome.modified.new.gff3"
GOI_BED       <- "genes_of_interest.bed"
HRP_BED       <- "hrp.bed"

pals <- c("#2BAE84","#3366CC","#8153A6","#E87DBF","#FF7033",
          "#F4A736","#D6A419","#3FB1C2","#8ACB4A","#A5426D","#6EC4E8")

############################################################
# 1. LOAD DATA
############################################################

variants <- read.delim(VARIANTS_FILE, header = FALSE)
samples  <- read.table(SAMPLES_FILE)

colnames(variants) <- c("CHR","POS","REF","ALT","QUAL",
                        "AC","AN","AF","NS",
                        "SVTYPE","SVLEN",
                        samples$V1)

sample_cols       <- samples$V1
variants$SVLEN    <- as.numeric(as.character(variants$SVLEN))

############################################################
# 2. VARIANT PROCESSING AND CLASSIFICATION
############################################################

process_variants <- function(df,
                             sample_cols,
                             missing_threshold = 0.9,
                             filter_NS         = FALSE,
                             remove_fixed      = FALSE,
                             max_len           = 10000) {

  df <- df %>%
    mutate(
      LEN    = nchar(ALT) - nchar(REF),
      ABSLEN = abs(LEN),
      n_missing    = rowSums(across(all_of(sample_cols),
                                    ~ . %in% c(".", NA))),
      missing_frac = n_missing / length(sample_cols)
    )

  df <- df %>%
    filter(missing_frac <= missing_threshold) %>%
    filter(AF != 0)

  if (filter_NS)    df <- df %>% filter(NS != 1)
  if (remove_fixed) df <- df %>% filter(AF != 1)

  df <- df %>% filter(ABSLEN < max_len)

  df <- df %>%
    mutate(
      TYPE = case_when(
        LEN == 0             ~ "SNP",
        LEN >= 50             ~ "Structural Insertion",
        LEN > 0 & LEN < 50  ~ "Small Insertion",
        LEN <= -50            ~ "Structural Deletion",
        LEN < 0 & LEN > -50 ~ "Small Deletion",
        TRUE ~ NA_character_
      ),
      TYPE = factor(TYPE, levels = c("SNP",
                                     "Small Insertion",
                                     "Small Deletion",
                                     "Structural Insertion",
                                     "Structural Deletion"))
    )

  return(df)
}

variants_clean <- process_variants(
  variants,
  sample_cols       = sample_cols,
  missing_threshold = 0.9
)

############################################################
# 3. GRANGES BUILDER
############################################################

make_granges <- function(df) {
  end_pos <- ifelse(df$LEN < 0,
                    df$POS + abs(df$LEN),
                    df$POS)
  GRanges(seqnames = df$CHR,
          ranges   = IRanges(start = df$POS, end = end_pos),
          LEN      = df$LEN)
}

############################################################
# 4. OPTIONAL BED FILTER
############################################################

filter_by_bed <- function(df, bed_file) {
  bed <- read.delim(bed_file, header = FALSE)
  colnames(bed) <- c("CHR","START","END","NAME")
  var_gr <- make_granges(df)
  bed_gr <- GRanges(seqnames = bed$CHR,
                    ranges   = IRanges(start = bed$START, end = bed$END))
  hits <- findOverlaps(var_gr, bed_gr)
  df[unique(queryHits(hits)), ]
}

# Example usage:
# variants_clean <- filter_by_bed(variants_clean, "regions.bed")

############################################################
# 5. SUMMARY EXPORTS
############################################################

write_variant_summaries <- function(df, prefix) {

  write_csv(df %>% count(TYPE),
            paste0(prefix, "_variant_type_counts.csv"))

  write_csv(df %>%
              group_by(TYPE) %>%
              summarise(mean_AF   = mean(AF),
                        median_AF = median(AF),
                        n         = n()),
            paste0(prefix, "_AF_summary.csv"))

  write_csv(df %>%
              filter(abs(LEN) > 50) %>%
              mutate(SVTYPE = ifelse(LEN > 0, "Insertion", "Deletion")) %>%
              count(SVTYPE),
            paste0(prefix, "_SV_counts.csv"))

  # Additional manuscript statistics
  singletons <- df %>% group_by(TYPE, AC) %>% summarise(n = n(), .groups = "drop")
  write_csv(singletons, paste0(prefix, "_singletons.csv"))

  for (vtype in c("Small Insertion","Small Deletion",
                   "Structural Insertion","Structural Deletion")) {
    slug <- gsub(" ", "_", vtype)
    stats <- df %>%
      filter(TYPE == vtype) %>%
      summarise(median_LEN = median(LEN, na.rm = TRUE),
                min_LEN    = min(LEN,    na.rm = TRUE),
                max_LEN    = max(LEN,    na.rm = TRUE))
    write_csv(stats, paste0(prefix, "_", slug, "_length_stats.csv"))
  }
}

write_variant_summaries(variants_clean, "all_variants")

############################################################
# 6. CORE VARIANT PLOTS
############################################################

# Variant types
ggsave("variant_types.tiff",
       ggplot(variants_clean, aes(TYPE, fill = TYPE)) +
         geom_bar() +
         geom_text(stat = "count",
                   aes(label = after_stat(count)),
                   vjust = -0.3, size = 6) +
         scale_fill_manual(values = setNames(pals[1:5],
                                             levels(variants_clean$TYPE))) +
         theme_minimal(base_size = 20) +
         theme(axis.text.x   = element_text(angle = 45, hjust = 1),
               legend.position = "none"),
       width = 10, height = 10, dpi = 300, device = "tiff")

# SV length distribution
sv_data <- variants_clean %>%
  filter(abs(LEN) > 50) %>%
  mutate(SVTYPE = ifelse(LEN > 0, "Insertion", "Deletion"),
         absLEN = abs(LEN))

ggsave("sv_length_distribution.tiff",
       ggplot(sv_data, aes(absLEN, fill = SVTYPE)) +
         geom_histogram(bins = 40, colour = "white", linewidth = 0.2) +
         scale_fill_manual(values = c("Insertion" = pals[1],
                                      "Deletion"  = pals[5])) +
         theme_minimal(base_size = 20) +
         theme(legend.position = "bottom"),
       width = 10, height = 8, dpi = 300, device = "tiff")

# AF distribution
ggsave("AF_distribution.tiff",
       ggplot(variants_clean, aes(AF, fill = TYPE)) +
         geom_histogram(binwidth = 0.05, colour = "white", linewidth = 0.2) +
         scale_fill_manual(values = setNames(pals[1:5],
                                             levels(variants_clean$TYPE))) +
         theme_minimal(base_size = 20) +
         theme(legend.position = "bottom"),
       width = 12, height = 6, dpi = 300, device = "tiff")

############################################################
# 7. GOI ANALYSIS
############################################################

analyse_GOI <- function(df, gene_bed, prefix,
                        promoter_upstream   = 2000,
                        promoter_downstream = 500) {

  genes <- read.delim(gene_bed, header = FALSE)
  colnames(genes) <- c("CHR","START","END","GENE")

  gene_gr <- GRanges(seqnames = genes$CHR,
                     ranges   = IRanges(start = genes$START,
                                        end   = genes$END),
                     gene     = genes$GENE)

  upstream_gr   <- flank(gene_gr, width = promoter_upstream,  start = TRUE)
  downstream_gr <- flank(gene_gr, width = promoter_downstream, start = FALSE)

  gene_gr_extended       <- punion(punion(gene_gr, upstream_gr), downstream_gr)
  gene_gr_extended$gene  <- gene_gr$gene

  var_gr <- make_granges(df)
  hits   <- findOverlaps(var_gr, gene_gr_extended)

  GOI_df <- df[queryHits(hits), ] %>%
    mutate(gene = gene_gr_extended$gene[subjectHits(hits)])

  write_csv(GOI_df %>% count(gene),
            paste0(prefix, "_GOI_counts.csv"))

  sum_goi <- GOI_df %>%
    select(CHR, POS, gene, TYPE) %>%
    unique() %>%
    group_by(gene, TYPE) %>%
    summarise(n = n(), .groups = "drop")
  write_csv(sum_goi, paste0(prefix, "_summary.csv"))

  return(GOI_df)
}

goi_df <- analyse_GOI(variants_clean,
                      GOI_BED, "GOI",
                      promoter_upstream   = 2000,
                      promoter_downstream = 500)

############################################################
# 8. HRP ANALYSIS
############################################################

hrp_df <- analyse_GOI(variants_clean,
                      HRP_BED, "HRP",
                      promoter_upstream   = 100,
                      promoter_downstream = 100)

############################################################
# 9. COMPLEX LOCI ANALYSIS
############################################################

analyse_complex_loci <- function(df,
                                 sample_cols,
                                 gff_file,
                                 window_size  = 200,
                                 buffer       = 200,
                                 min_alleles  = 2,
                                 min_span     = 1000,
                                 prefix       = "complex") {

  # Prepare variants
  df <- df %>%
    rowwise() %>%
    mutate(
      all_alleles = list(c(str_split(ALT, ",")[[1]], REF)),
      END         = ifelse(LEN < 0, POS + abs(LEN), POS + LEN)
    ) %>%
    ungroup() %>%
    mutate(
      n_site_alleles = map_int(all_alleles, length),
      n_missing      = rowSums(across(all_of(sample_cols),
                                      ~ . %in% c(".", NA))),
      missing_frac   = n_missing / length(sample_cols)
    ) %>%
    arrange(CHR, POS)

  # Cluster positions into windows
  df <- df %>%
    group_by(CHR) %>%
    mutate(cluster = cumsum(c(1, diff(POS) > window_size))) %>%
    ungroup()

  # Aggregate clusters
  complex_loci <- df %>%
    group_by(CHR, cluster) %>%
    summarise(
      n_alleles    = sum(n_site_alleles, na.rm = TRUE),
      start        = min(POS),
      end          = max(END),
      span         = end - start,
      mean_missing = mean(missing_frac, na.rm = TRUE),
      n_sites      = n(),
      .groups      = "drop"
    ) %>%
    filter(n_alleles >= min_alleles, span >= min_span) %>%
    mutate(start_buffered = start - buffer,
           end_buffered   = end   + buffer)

  # GRanges for complex loci
  complex_gr <- GRanges(
    seqnames     = complex_loci$CHR,
    ranges       = IRanges(start = complex_loci$start_buffered,
                           end   = complex_loci$end_buffered),
    n_alleles    = complex_loci$n_alleles,
    mean_missing = complex_loci$mean_missing
  )

  # Read GFF and extract genes
  gff <- read.delim(gff_file, header = FALSE, sep = "\t", comment.char = "#")
  colnames(gff) <- c("seqid","source","type","start","end",
                     "score","strand","phase","attributes")

  genes <- gff %>%
    filter(type == "gene") %>%
    mutate(
      gene_id          = sub(".*ID=([^;]+);.*",   "\\1", attributes),
      gene_name        = sub(".*Name=([^;]+);?.*","\\1", attributes),
      start_buffered   = pmax(start - buffer, 1),
      end_buffered     = end + buffer
    ) %>%
    select(seqid, start_buffered, end_buffered, strand, gene_id, gene_name)

  gene_gr <- GRanges(
    seqnames  = genes$seqid,
    ranges    = IRanges(start = genes$start_buffered,
                        end   = genes$end_buffered),
    gene_id   = genes$gene_id,
    gene_name = genes$gene_name
  )

  # Find overlaps
  hits <- findOverlaps(complex_gr, gene_gr)

  complex_with_genes <- data.frame(
    complex_loci[queryHits(hits), ],
    gene_id   = gene_gr$gene_id[subjectHits(hits)],
    gene_name = gene_gr$gene_name[subjectHits(hits)]
  )

  complex_with_genes$gene_name <- ifelse(
    grepl("^ID=", complex_with_genes$gene_name), NA,
    complex_with_genes$gene_name
  )

  complex_with_genes_collapsed <- complex_with_genes %>%
    group_by(CHR, cluster, start, end, span, n_alleles, mean_missing, n_sites) %>%
    summarise(
      gene_ids   = paste(unique(gene_id),                    collapse = ","),
      gene_names = paste(unique(gene_name[!is.na(gene_name)]), collapse = ","),
      .groups    = "drop"
    )

  write_csv(complex_with_genes_collapsed, paste0(prefix, "_summary.csv"))
  return(complex_with_genes_collapsed)
}

analyse_complex_loci(
  variants_clean %>% filter(abs(LEN) > 1000),
  sample_cols,
  GFF_FILE,
  window_size = 200,
  buffer      = 500,
  min_alleles = 2,
  min_span    = 1000,
  prefix      = "pf_complex"
)

############################################################
# 10. LOAD GENOMIC REGION BED FILES
############################################################

region_files <- list(
  Centromere                 = "centromere.bed",
  Internal_Hypervariable     = "internal_hypervar.bed",
  Subtelomeric_Hypervariable = "subtelomere_hypervar.bed",
  Core_Genome                = "Core_genome_Pf3D7_v3_ext.bed",
  Subtelomeric_Repeat        = "subtelomere_repeat.bed",
  Coding_Regions             = "Pf3D7_R.coding.regions.bed"
)

load_region_granges <- function(bed_file, region_name) {
  bed <- read.delim(bed_file, header = FALSE)
  colnames(bed)[1:3] <- c("CHR","START","END")
  GRanges(seqnames = bed$CHR,
          ranges   = IRanges(start = bed$START, end = bed$END),
          region   = region_name)
}

region_gr_list <- lapply(names(region_files), function(n) {
  load_region_granges(region_files[[n]], n)
})
names(region_gr_list) <- names(region_files)

############################################################
# 11. MAP VARIANTS TO REGIONS
############################################################

variant_gr <- make_granges(variants_clean)

assign_variants_to_regions <- function(var_gr, region_gr_list) {
  results <- list()
  for (region_name in names(region_gr_list)) {
    hits <- findOverlaps(var_gr, region_gr_list[[region_name]])
    if (length(hits) > 0) {
      df_region          <- variants_clean[queryHits(hits), ]
      df_region$REGION   <- region_name
      results[[region_name]] <- df_region
    }
  }
  bind_rows(results)
}

variants_by_region <- assign_variants_to_regions(variant_gr, region_gr_list) %>%
  mutate(signed_LEN = LEN)

# Summary table
type_counts_summary <- variants_by_region %>%
  group_by(REGION, TYPE) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(REGION) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

write_csv(type_counts_summary, "variant_type_counts_by_region_summary.csv")

cat("Variants per region:\n")
print(type_counts_summary %>% group_by(REGION) %>% summarise(n = sum(n)))

############################################################
# 12. REGIONAL PLOTS
############################################################

# Signed SV length by region
ggsave("multipanel_signed_SV_by_region.tiff",
       ggplot(variants_by_region %>% filter(TYPE != "SNP"),
              aes(x = signed_LEN, fill = TYPE)) +
         geom_histogram(bins = 50, color = "black") +
         facet_wrap(~ REGION, scales = "free_y") +
         scale_fill_manual(values = setNames(pals[1:5],
                                             levels(variants_clean$TYPE))) +
         theme_minimal(base_size = 18) +
         labs(title = "Signed SV length distribution by genome region",
              x     = "SV length (deletions negative, insertions positive)",
              y     = "Count"),
       width = 14, height = 10, dpi = 300, device = "tiff")

# Variant type counts by region (stacked bar)
ggsave("multipanel_variant_types_counts_by_region.tiff",
       ggplot(type_counts_summary,
              aes(x = gsub("_", " ", REGION), y = n, fill = TYPE)) +
         geom_bar(stat = "identity", position = "stack") +
         scale_fill_manual(values = setNames(pals[1:5],
                                             levels(variants_clean$TYPE))) +
         scale_y_continuous(labels = comma_format()) +
         theme_minimal(base_size = 16) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
         labs(title = "Variant type counts by genome region",
              y = "Count of variants", x = ""),
       width = 14, height = 10, dpi = 300, device = "tiff")

# Variant types faceted by region
ggsave("multipanel_variant_types_by_region.tiff",
       ggplot(variants_by_region, aes(TYPE, fill = TYPE)) +
         geom_bar() +
         facet_wrap(~ REGION, scales = "free_y") +
         scale_fill_manual(values = setNames(pals[1:5],
                                             levels(variants_clean$TYPE))) +
         theme_minimal(base_size = 18) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
         labs(title = "Variant types by genome region"),
       width = 14, height = 10, dpi = 300, device = "tiff")

# AF distribution by region
ggsave("multipanel_AF_by_region.tiff",
       ggplot(variants_by_region, aes(AF, fill = TYPE)) +
         geom_histogram(binwidth = 0.05, color = "black") +
         facet_wrap(~ REGION) +
         scale_fill_manual(values = setNames(pals[1:5],
                                             levels(variants_clean$TYPE))) +
         theme_minimal(base_size = 18) +
         labs(title = "Allele frequency by genome region"),
       width = 14, height = 10, dpi = 300, device = "tiff")

############################################################
# 13. SV GENOMIC DISTRIBUTION (CHROMOSOMAL DENSITY)
############################################################

sv_chr <- variants %>%
  filter(SVTYPE %in% c("INS","DEL")) %>%
  mutate(
    SVLEN     = as.numeric(SVLEN),
    size_bin  = case_when(
      abs(SVLEN) < 50                          ~ "<50 bp",
      abs(SVLEN) >= 50  & abs(SVLEN) < 200    ~ "50-200 bp",
      abs(SVLEN) >= 200 & abs(SVLEN) < 1000   ~ "200-1000 bp",
      abs(SVLEN) >= 1000                       ~ ">1000 bp"
    ),
    size_bin  = factor(size_bin,
                       levels = c("<50 bp","50-200 bp","200-1000 bp",">1000 bp")),
    CHR_label = gsub("Pf3D7_","Chr ", gsub("_v3","", CHR))
  )

size_cols <- c("<50 bp"      = pals[1],
               "50-200 bp"   = pals[5],
               "200-1000 bp" = pals[6],
               ">1000 bp"    = pals[10])

ggsave("sv_density_by_chr.pdf",
       ggplot(sv_chr, aes(x = POS / 1e6, fill = size_bin)) +
         geom_histogram(binwidth = 0.05, colour = NA, alpha = 0.9) +
         facet_grid(SVTYPE ~ CHR_label, scales = "free_x", space = "free_x") +
         scale_fill_manual(values = size_cols, name = "SV size") +
         scale_x_continuous(labels = scales::comma_format(suffix = " Mb")) +
         labs(x = "Chromosomal position (Mb)", y = "Count",
              title = "Structural variant density across P. falciparum chromosomes") +
         theme_bw(base_size = 11) +
         theme(strip.background  = element_rect(fill = "grey92", colour = "grey60"),
               strip.text        = element_text(size = 8, face = "bold"),
               axis.text.x       = element_text(angle = 45, hjust = 1, size = 7),
               legend.position   = "bottom",
               panel.spacing.x   = unit(0.15, "lines"),
               panel.spacing.y   = unit(0.4,  "lines")),
       width = 18, height = 6)

ggsave("sv_density_by_chr.png",
       last_plot(), width = 18, height = 6, dpi = 300)

# SV length distribution (log scale)
sv_len <- variants %>%
  filter(SVTYPE %in% c("INS","DEL")) %>%
  mutate(SVLEN = as.numeric(SVLEN)) %>%
  filter(!is.na(SVLEN))

ggsave("sv_length_dist.png",
       ggplot(sv_len, aes(x = abs(SVLEN), fill = SVTYPE)) +
         geom_histogram(bins = 60, colour = NA, alpha = 0.9) +
         facet_wrap(~ SVTYPE, ncol = 1, scales = "free_y") +
         scale_fill_manual(values = c("INS" = pals[1], "DEL" = pals[5]),
                           guide  = "none") +
         scale_x_log10(breaks = c(1, 10, 50, 200, 500, 1000, 5000),
                       labels = c("1","10","50","200","500","1 kb","5 kb")) +
         annotation_logticks(sides = "b") +
         labs(x = "SV length (bp, log scale)", y = "Count") +
         theme_bw(base_size = 12) +
         theme(strip.background = element_rect(fill = "grey92", colour = "grey60"),
               strip.text       = element_text(size = 11, face = "bold"),
               panel.spacing    = unit(0.5, "lines")),
       width = 8, height = 6, dpi = 300)

############################################################
# 14. SV FOLD ENRICHMENT
############################################################

sv_enrich <- variants %>%
  filter(SVTYPE %in% c("INS","DEL")) %>%
  mutate(
    SVLEN     = as.numeric(SVLEN),
    window    = floor(POS / 10000) * 10000,
    CHR_label = gsub("Pf3D7_","Chr ", gsub("_v3","", CHR))
  ) %>%
  filter(!is.na(SVLEN)) %>%
  group_by(CHR_label, window, SVTYPE) %>%
  summarise(n_sv = n(), .groups = "drop") %>%
  group_by(CHR_label, SVTYPE) %>%
  mutate(mean_rate       = mean(n_sv),
         fold_enrichment = n_sv / mean_rate) %>%
  ungroup()

# Line plot
ggsave("sv_fold_enrichment.pdf",
       ggplot(sv_enrich, aes(x = window / 1e6, y = fold_enrichment,
                             colour = SVTYPE)) +
         geom_hline(yintercept = 1, linetype = "dashed",
                    colour = "grey50", linewidth = 0.4) +
         geom_line(linewidth = 0.4, alpha = 0.6) +
         geom_point(size = 0.8, alpha = 0.8) +
         facet_grid(SVTYPE ~ CHR_label, scales = "free_x", space = "free_x") +
         scale_colour_manual(values = c("INS" = pals[1], "DEL" = pals[5]),
                             guide  = "none") +
         scale_x_continuous(labels = scales::comma_format(suffix = "Mb")) +
         scale_y_continuous(trans   = "log2",
                            breaks  = c(0.25, 0.5, 1, 2, 4, 8),
                            labels  = c("0.25x","0.5x","1x","2x","4x","8x")) +
         labs(x = "Chromosomal position (Mb)",
              y = "Fold enrichment over chr mean (log2)",
              title = "SV fold enrichment per 10 kb window") +
         theme_bw(base_size = 10) +
         theme(strip.background = element_rect(fill = "grey92", colour = "grey60"),
               strip.text       = element_text(size = 7, face = "bold"),
               axis.text.x      = element_text(angle = 45, hjust = 1, size = 6),
               panel.spacing.x  = unit(0.1, "lines"),
               panel.spacing.y  = unit(0.4, "lines")),
       width = 20, height = 5)

ggsave("sv_fold_enrichment.png", last_plot(), width = 20, height = 5, dpi = 300)

# Area plot (per chromosome)
ggsave("sv_fold_enrichment_coverage.pdf",
       ggplot(sv_enrich, aes(x = window / 1e6, y = fold_enrichment,
                             fill = SVTYPE)) +
         geom_hline(yintercept = 1, linetype = "dashed",
                    colour = "grey60", linewidth = 0.3) +
         geom_area(alpha = 0.7, colour = NA) +
         geom_line(aes(colour = SVTYPE), linewidth = 0.3, alpha = 0.9) +
         facet_grid(CHR_label ~ SVTYPE, scales = "free_x", space = "free_x") +
         scale_fill_manual(values   = c("INS" = pals[1], "DEL" = pals[5]),
                           guide    = "none") +
         scale_colour_manual(values = c("INS" = pals[1], "DEL" = pals[5]),
                             guide  = "none") +
         scale_x_continuous(labels  = scales::comma_format(suffix = " Mb"),
                            expand  = c(0, 0)) +
         scale_y_continuous(breaks  = c(0, 1, 2, 4),
                            labels  = c("0x","1x","2x","4x"),
                            expand  = c(0, 0)) +
         labs(x = "Chromosomal position (Mb)",
              y = "Fold enrichment",
              title = "SV fold enrichment per 10 kb window") +
         theme_bw(base_size = 10) +
         theme(strip.background  = element_rect(fill = "grey92", colour = "grey60"),
               strip.text.y      = element_text(size = 7, face = "bold", angle = 0),
               strip.text.x      = element_text(size = 9, face = "bold"),
               axis.text.x       = element_text(size = 7),
               axis.text.y       = element_text(size = 6),
               panel.spacing.x   = unit(0.8, "lines"),
               panel.spacing.y   = unit(0.1, "lines"),
               panel.grid.minor  = element_blank(),
               panel.grid.major.x = element_blank()),
       width = 10, height = 14)

ggsave("sv_fold_enrichment_coverage.png",
       last_plot(), width = 10, height = 14, dpi = 300)

############################################################
# 15. SV ENRICHMENT BED EXPORT
############################################################

# Chromosome boundaries (10kb buffer each end)
chr_limits <- data.frame(
  chr_num = 1:14,
  chr_max = c(640851, 947102, 1067971, 1200490, 1343557,
              1418242, 1445207, 1472805, 1541735, 1687656,
              2038340, 2271494, 2925236, 3291936)
) %>%
  mutate(chr_start = 10000,
         chr_end   = chr_max - 10000)

sv_bed <- variants %>%
  filter(SVTYPE %in% c("INS","DEL")) %>%
  mutate(
    SVLEN   = abs(as.numeric(SVLEN)),
    window  = floor(POS / 100000) * 100000,
    chr_num = as.integer(gsub("Pf3D7_0*([0-9]+)_v3","\\1", CHR))
  ) %>%
  filter(!is.na(SVLEN)) %>%
  group_by(chr_num, window, SVTYPE) %>%
  summarise(n_sv = n(), .groups = "drop") %>%
  group_by(chr_num, SVTYPE) %>%
  mutate(mean_rate       = mean(n_sv),
         fold_enrichment = ifelse(mean_rate > 0, n_sv / mean_rate, 0)) %>%
  ungroup() %>%
  group_by(chr_num, window) %>%
  summarise(fold_enrichment = mean(fold_enrichment), .groups = "drop") %>%
  mutate(start = window, end = window + 100000) %>%
  left_join(chr_limits, by = "chr_num") %>%
  mutate(start = pmax(start, chr_start),
         end   = pmin(end,   chr_end)) %>%
  filter(start < end) %>%
  mutate(sv_class = case_when(
    fold_enrichment == 0   ~ "low",
    fold_enrichment <  1.0 ~ "low",
    fold_enrichment <  2.0 ~ "mid",
    fold_enrichment <  3.0 ~ "high",
    fold_enrichment >= 3.0 ~ "hotspot"
  ))

write.table(
  sv_bed %>%
    mutate(start = format(start, scientific = FALSE, trim = TRUE),
           end   = format(end,   scientific = FALSE, trim = TRUE),
           label = sv_class) %>%
    select(chr_num, start, end, label),
  "sv_enrichment.bed",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

idmap <- data.frame(
  label  = c("low",   "mid",    "high",   "hotspot"),
  short  = c("LOW",   "MED",    "HIGH",   "HOTSPOT"),
  colour = c("#FFF",  pals[8],  pals[3],  pals[10])
)

write.table(idmap, "sv_enrichment.idmap",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

############################################################
# END
############################################################
