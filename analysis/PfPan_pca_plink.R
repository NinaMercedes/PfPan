library(dplyr)
library(ggplot2)
library(ape)
library(patchwork)

# ---- Palette fix ---------------------------------------------------------
# Original palette (unchanged except one entry):
#   "#2BAE84","#3366CC","#8153A6","#E87DBF","#FF7033",
#   "#F4A736","#D6A419","#3FB1C2","#8ACB4A","#A5426D","#6EC4E8"
# Reviewer flagged the golds/yellows (#F4A736 and #D6A419) as too similar
# to distinguish. Fix: #D6A419 (mustard gold) replaced with #E63946 (red),
# a hue not otherwise used in the palette. Every other color unchanged.
pals <- c("#2BAE84", "#3366CC", "#8153A6", "#E87DBF", "#FF7033",
          "#F4A736", "#E63946", "#3FB1C2", "#8ACB4A", "#A5426D", "#6EC4E8")
shapes <- c(16, 17, 15, 18, 8, 4, 3, 7, 9, 10, 12)

workdir <- "."
prefix1 <- "pan_SR_GTK.bi"          # A: pangenome, all variants
prefix2 <- "norm_SR_GTK.bi"         # B: linear reference
prefix3 <- "pan_novel_SR_GTK.bi.filt"    # C: pangenome, corroborated novel insertions only
metadata <- "../metadata_high_quality.csv"

calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

load_dist_pca <- function(prefix, met) {
    dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
    id   <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))
    desc <- id %>% left_join(met, by = c("V1" = "sample"))

    dist_m <- as.matrix(dist)
    colnames(dist_m) <- desc$V1
    rownames(dist_m) <- desc$V1

    cmd  <- cmdscale(dist_m, k = 5, eig = TRUE, x.ret = TRUE)
    vars <- calc_variance_explained(cmd)

    df <- as.data.frame(cmd$points, stringsAsFactors = FALSE)
    df$country   <- gsub("_", " ", desc$Country)
    df$region    <- gsub("_", " ", desc$Region)
    df$sample_id <- rownames(dist_m)
    colnames(df) <- gsub("V", "D", colnames(df))

    list(df = df, vars = vars, dist_m = dist_m)
}

# METADATA
met <- read.csv(metadata, stringsAsFactors = FALSE, header = TRUE)

res1 <- load_dist_pca(prefix1, met)
res2 <- load_dist_pca(prefix2, met)
res3 <- load_dist_pca(prefix3, met)

df1 <- res1$df; vars1 <- res1$vars
df2 <- res2$df; vars2 <- res2$vars
df3 <- res3$df; vars3 <- res3$vars

# ---- Axis orientation (reviewer: Figure 4) --------------------------------
# The sign of each principal coordinate is arbitrary: cmdscale() can return
# an axis mirrored relative to another analysis of very similar data, which
# makes panels look different when they are not. Each axis of panels B and C
# is therefore oriented to match panel A, using the samples shared between
# panels: if an axis correlates negatively with the same axis in A, it is
# multiplied by -1. This changes only the direction of the axis, never the
# distances between samples or the variance explained.
align_axes <- function(df, ref, axes = c("D1", "D2"), label = "") {
    shared <- inner_join(df  %>% select(sample_id, all_of(axes)),
                         ref %>% select(sample_id, all_of(axes)),
                         by = "sample_id", suffix = c("", "_ref"))
    cat(sprintf("--- Orientation of panel %s relative to panel A (%d shared samples) ---\n",
                label, nrow(shared)))
    # Full correlation matrix, so a PC1 <-> PC2 swap between analyses is visible
    print(round(cor(shared[, axes], shared[, paste0(axes, "_ref")]), 2))
    for (ax in axes) {
        r <- cor(shared[[ax]], shared[[paste0(ax, "_ref")]])
        if (r < 0) {
            df[[ax]] <- -df[[ax]]
            cat(sprintf("  %s flipped (r with panel A = %.2f)\n", ax, r))
        } else {
            cat(sprintf("  %s kept     (r with panel A = %.2f)\n", ax, r))
        }
    }
    cat("\n")
    df
}

df2 <- align_axes(df2, df1, label = "B")
df3 <- align_axes(df3, df1, label = "C")

# Manual override, if you prefer to choose orientations by eye instead:
# df2$D1 <- -df2$D1   # mirror panel B left-right
# df3$D2 <- -df3$D2   # mirror panel C top-bottom

color_by <- "region"

# ---- Outlier / unusual clustering check ----------------------------------
# For each sample, compare its PC1 rank in the novel-only PCA (df3) against
# its PC1 rank in the full pangenome PCA (df1). A sample whose rank shifts
# drastically suggests its "novel" sequence content behaves unusually
# relative to the rest of the population.
# (Runs after axis alignment, so rank shifts are not inflated by a mirrored PC1.)
rank_shift <- df1 %>%
    select(sample_id, D1_full = D1) %>%
    inner_join(df3 %>% select(sample_id, D1_novel = D1), by = "sample_id") %>%
    mutate(
        rank_full  = rank(D1_full),
        rank_novel = rank(D1_novel),
        rank_shift = abs(rank_full - rank_novel)
    ) %>%
    arrange(desc(rank_shift))

cat("=== Samples with the largest PC1 rank shift (full vs. novel-only PCA) ===\n")
print(head(rank_shift, 10))
write.csv(rank_shift, "pca_novel_vs_full_rank_shift.csv", row.names = FALSE)
cat("Full table written to: pca_novel_vs_full_rank_shift.csv\n\n")

# ---- SHARED THEME ---------------------------------------------------------
base_theme <- theme_classic(base_size = 22) +
    theme(
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold", size = 22),
        legend.text      = element_text(size = 20),
        axis.title       = element_text(face = "bold", size = 22),
        axis.text        = element_text(size = 19),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.title       = element_text(face = "bold", size = 28)
    )

make_panel <- function(df, vars, title_letter) {
    ggplot(df, aes(x = D1, y = D2, fill = !!sym(color_by))) +
        geom_point(shape = 21, size = 2.3, alpha = 0.75,
                   colour = "grey20", stroke = 0.25) +
        scale_fill_manual(values = pals, name = "Region") +
        # Larger text makes the 11-region legend wider than the figure on one
        # line, so wrap it to two rows and enlarge the legend keys to match.
        guides(fill = guide_legend(nrow = 2,
                                   override.aes = list(size = 5, alpha = 1))) +
        labs(
            title = title_letter,
            x = paste0("PC1 (", vars["PC1"], "%)"),
            y = paste0("PC2 (", vars["PC2"], "%)")
        ) +
        base_theme
}

pA <- make_panel(df1, vars1, "A")
pB <- make_panel(df2, vars2, "B")
pC <- make_panel(df3, vars3, "C")

# ---- COMBINE AND SAVE ------------------------------------------------------
combined <- pA + pB + pC +
    plot_layout(guides = "collect", nrow = 1) &
    theme(legend.position = "bottom")

ggsave(
    plot     = combined,
    filename = "PCA_panelled_ABC.tiff",
    device   = "tiff",
    width    = 22,
    height   = 8,   # was 7; extra room for the larger text and two-row legend
    dpi      = 300,
    bg       = "white"
)

cat("Saved: PCA_panelled_ABC.tiff\n")
cat(sprintf("Panel A (full pangenome)      variance PC1/PC2: %.1f%% / %.1f%%\n", vars1["PC1"], vars1["PC2"]))
cat(sprintf("Panel B (linear)              variance PC1/PC2: %.1f%% / %.1f%%\n", vars2["PC1"], vars2["PC2"]))
cat(sprintf("Panel C (corroborated novel)  variance PC1/PC2: %.1f%% / %.1f%%\n", vars3["PC1"], vars3["PC2"]))
