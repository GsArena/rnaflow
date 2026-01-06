#!/usr/bin/env Rscript
library(DESeq2)
library(optparse)
library(tidyverse)

option_list <- list(
  make_option("--counts", type="character"),
  make_option("--samples", type="character"),
  make_option("--contrast", type="character", nargs=2),
  make_option("--out", type="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

counts <- read.delim(opt$counts, comment.char="#")
rownames(counts) <- counts$Geneid
counts <- counts[,7:ncol(counts)]

meta <- read.csv(opt$samples)
meta <- meta[match(colnames(counts), meta$sample),]
rownames(meta) <- meta$sample

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10,]
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition", opt$contrast[1], opt$contrast[2]))
res <- lfcShrink(dds, coef=2, res=res)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

write.csv(res_df, opt$out, row.names=FALSE)
write.csv(counts(dds, normalized=TRUE),
          "results/de/normalized_counts.csv")
