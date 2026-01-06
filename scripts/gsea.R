#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(tidyverse)
  library(optparse)
})

# ---------------------------
# Arguments
# ---------------------------
option_list <- list(
  make_option("--de", type = "character", help = "Full DESeq2 results CSV"),
  make_option("--outdir", type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Load DE results
# ---------------------------
de <- read.csv(opt$de)

de <- de %>%
  filter(!is.na(log2FoldChange), !is.na(pvalue))

# Ranking metric
de$rank <- sign(de$log2FoldChange) * -log10(de$pvalue)

# Ensembl → Entrez
gene_map <- bitr(
  de$gene_id,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

de <- inner_join(de, gene_map, by = c("gene_id" = "ENSEMBL"))

gene_list <- de$rank
names(gene_list) <- de$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# ---------------------------
# GSEA GO BP
# ---------------------------
gsea_go <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

write.csv(as.data.frame(gsea_go),
          file = file.path(opt$outdir, "GSEA_GO_BP_results.csv"),
          row.names = FALSE)

png(file.path(opt$outdir, "GSEA_GO_BP_dotplot.png"),
    width = 1800, height = 1600, res = 300)

dotplot(gsea_go, showCategory = 15) +
  ggtitle("GSEA – GO Biological Process")

dev.off()

# ---------------------------
# GSEA KEGG
# ---------------------------
gsea_kegg <- gseKEGG(
  geneList     = gene_list,
  organism     = "hsa",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

write.csv(as.data.frame(gsea_kegg),
          file = file.path(opt$outdir, "GSEA_KEGG_results.csv"),
          row.names = FALSE)

png(file.path(opt$outdir, "GSEA_KEGG_dotplot.png"),
    width = 1800, height = 1600, res = 300)

dotplot(gsea_kegg, showCategory = 15) +
  ggtitle("GSEA – KEGG Pathways")

dev.off()
