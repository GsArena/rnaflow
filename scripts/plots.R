#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(tidyverse)
  library(optparse)
})

# ---------------------------
# Arguments
# ---------------------------
option_list <- list(
  make_option("--de", type = "character", help = "DESeq2 results CSV"),
  make_option("--norm", type = "character", help = "Normalized counts CSV"),
  make_option("--meta", type = "character", help = "Sample metadata CSV"),
  make_option("--outdir", type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Load data
# ---------------------------
de <- read.csv(opt$de)
norm <- read.csv(opt$norm, row.names = 1)
meta <- read.csv(opt$meta)

# ---------------------------
# Volcano plot
# ---------------------------
de$significance <- "NS"
de$significance[!is.na(de$padj) & de$padj < 0.05 & abs(de$log2FoldChange) >= 1] <- "DEG"

p_volcano <- ggplot(de, aes(log2FoldChange, -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_classic() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value")

ggsave(file.path(opt$outdir, "volcano_plot.png"),
       p_volcano, width = 7, height = 6, dpi = 300)

# ---------------------------
# PCA plot
# ---------------------------
log_norm <- log2(norm + 1)
pca <- prcomp(t(log_norm), scale. = TRUE)

pca_df <- as.data.frame(pca$x)
pca_df$sample <- rownames(pca_df)
pca_df <- left_join(pca_df, meta, by = c("sample" = "sample"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title = "PCA of Samples")

ggsave(file.path(opt$outdir, "pca_plot.png"),
       p_pca, width = 7, height = 6, dpi = 300)

# ---------------------------
# Heatmap (Top 20 DEGs)
# ---------------------------
top_genes <- de %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(20) %>%
  pull(gene_id)

heat <- log_norm[top_genes, ]

ann <- meta %>%
  column_to_rownames("sample") %>%
  select(condition)

png(file.path(opt$outdir, "heatmap_top20.png"),
    width = 1800, height = 2000, res = 300)

pheatmap(heat,
         annotation_col = ann,
         scale = "row",
         fontsize_row = 8,
         main = "Top 20 Differentially Expressed Genes")

dev.off()

