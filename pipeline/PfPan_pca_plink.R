library(dplyr)
library(ggplot2)
library(ape)
library(patchwork)

pals <- c("#2BAE84", "#3366CC", "#8153A6", "#E87DBF", "#FF7033",
          "#F4A736", "#D6A419", "#3FB1C2", "#8ACB4A", "#A5426D", "#6EC4E8")

workdir <- "."
prefix1 <- "pan_SR_GTK.bi"
prefix2 <- "norm_SR_GTK.bi"
metadata <- "../metadata_high_quality.csv"

calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

# METADATA
met <- read.csv(metadata, stringsAsFactors = FALSE, header = TRUE)

# DISTANCE MATRICES
dist1 <- read.table(file.path(workdir, paste0(prefix1, ".dist")), header = FALSE)
id1   <- read.table(file.path(workdir, paste0(prefix1, ".dist.id")))
dist2 <- read.table(file.path(workdir, paste0(prefix2, ".dist")), header = FALSE)
id2   <- read.table(file.path(workdir, paste0(prefix2, ".dist.id")))

desc1 <- id1 %>% left_join(met, by = c("V1" = "sample"))
desc2 <- id2 %>% left_join(met, by = c("V1" = "sample"))

dist_m1 <- as.matrix(dist1)
colnames(dist_m1) <- desc1$V1
rownames(dist_m1) <- desc1$V1

dist_m2 <- as.matrix(dist2)
colnames(dist_m2) <- desc2$V1
rownames(dist_m2) <- desc2$V1

# PCA
cmd1 <- cmdscale(dist_m1, k = 5, eig = TRUE, x.ret = TRUE)
vars1 <- calc_variance_explained(cmd1)

cmd2 <- cmdscale(dist_m2, k = 5, eig = TRUE, x.ret = TRUE)
vars2 <- calc_variance_explained(cmd2)

# BUILD DATAFRAMES - fixed: cmd1/cmd2 respectively, and rownames from dist_m
color_by <- "region"

df1 <- as.data.frame(cmd1$points, stringsAsFactors = FALSE)
df1$country   <- gsub("_", " ", desc1$Country)
df1$region    <- gsub("_", " ", desc1$Region)
df1$sample_id <- rownames(dist_m1)
colnames(df1) <- gsub("V", "D", colnames(df1))

df2 <- as.data.frame(cmd2$points, stringsAsFactors = FALSE)
df2$country   <- gsub("_", " ", desc2$Country)
df2$region    <- gsub("_", " ", desc2$Region)
df2$sample_id <- rownames(dist_m2)
colnames(df2) <- gsub("V", "D", colnames(df2))

# SHARED THEME
base_theme <- theme_classic(base_size = 16) +
    theme(
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold"),
        legend.text      = element_text(size = 14),
        axis.title       = element_text(face = "bold"),
        axis.text        = element_text(size = 14),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.title       = element_text(face = "bold", size = 18)
    )

# PANEL A - prefix1
pA <- ggplot(df1, aes(x = D1, y = D2, color = !!sym(color_by))) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = pals, name = "Region") +
    labs(
        title = "A",
        x = paste0("PC1 (", vars1["PC1"], "%)"),
        y = paste0("PC2 (", vars1["PC2"], "%)")
    ) +
    base_theme

# PANEL B - prefix2
pB <- ggplot(df2, aes(x = D1, y = D2, color = !!sym(color_by))) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = pals, name = "Region") +
    labs(
        title = "B",
        x = paste0("PC1 (", vars2["PC1"], "%)"),
        y = paste0("PC2 (", vars2["PC2"], "%)")
    ) +
    base_theme

# COMBINE AND SAVE
combined <- pA + pB +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

ggsave(
    plot     = combined,
    filename = "PCA_panelled_AB.tiff",
    device   = "tiff",
    width    = 16,
    height   = 7,
    dpi      = 300,
    bg       = "white"
)