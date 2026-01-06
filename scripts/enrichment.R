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
  make_option("--deg", type = "character", help = "Significant DEG CSV"),
  make_option("--outdir", type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Load DEG list
# ---------------------------
deg <- read.csv(opt$deg)

# Convert Ensembl → Entrez
gene_map <- bitr(
  deg$gene_id,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

entrez <- unique(gene_map$ENTREZID)

# ---------------------------
# GO enrichment (BP)
# ---------------------------
ego <- enrichGO(
  gene          = entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(as.data.frame(ego),
          file = file.path(opt$outdir, "GO_BP_enrichment.csv"),
          row.names = FALSE)

png(file.path(opt$outdir, "GO_BP_dotplot.png"),
    width = 1800, height = 1600, res = 300)

dotplot(ego, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment")

dev.off()

# ---------------------------
# KEGG enrichment
# ---------------------------
ekegg <- enrichKEGG(
  gene         = entrez,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

write.csv(as.data.frame(ekegg),
          file = file.path(opt$outdir, "KEGG_enrichment.csv"),
          row.names = FALSE)

png(file.path(opt$outdir, "KEGG_dotplot.png"),
    width = 1800, height = 1600, res = 300)

dotplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment")

dev.off()
