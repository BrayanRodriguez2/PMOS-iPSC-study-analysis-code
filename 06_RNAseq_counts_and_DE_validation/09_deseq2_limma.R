#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Release file: 06_RNAseq_counts_and_DE_validation/09_deseq2_limma.R   (original step script: 67_deseq2_limma.R)
# Produces: DESeq2 and limma-voom on raw counts and comparison with the OLS.
# Run: called by the SLURM script of the same number: Rscript 06_RNAseq_counts_and_DE_validation/09_deseq2_limma.R <args>
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
# Step 67a - methodological cross-check on RAW counts (our featureCounts counts, stranded): DESeq2 (negative binomial, Wald test) and limma-voom, PMOS vs control,
# with and without log2(XIST + 1) as covariate. Same gene set as the OLS replication of step 65 (the 16,495 symbols expressed in our normalised counts), counts summed by
# gene symbol (as in the normalised table). Size factors are estimated by DESeq2 (median of ratios) on ALL genes and then used for the tested subset; the limma-voom
# library sizes use edgeR TMM normalisation factors computed on all genes.
# Usage: Rscript 06_RNAseq_counts_and_DE_validation/09_deseq2_limma.R <counts_raw.tsv.gz (step 64)> <donor_meta.tsv> <de_ols.tsv.gz of ours_s2 (step 65)> <size_factors.tsv (step 64)> <out_dir> [exclude_line]
suppressPackageStartupMessages({ library(DESeq2); library(limma); library(edgeR) })
args <- commandArgs(trailingOnly = TRUE)
cnt <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE); meta <- read.delim(args[2], row.names = 1)
ols <- read.delim(args[3], stringsAsFactors = FALSE); sf_py <- read.delim(args[4], row.names = 1); out <- args[5]
excl <- if (length(args) >= 6) args[6] else NA
dir.create(out, showWarnings = FALSE, recursive = TRUE)
lines <- setdiff(rownames(meta), excl)
m <- as.matrix(cnt[, rownames(meta)]); sym <- cnt$gene_name
m <- rowsum(m, group = sym)                                        # symbol-level raw counts (sum over gene IDs sharing a symbol)
group <- factor(meta[rownames(meta), "group"], levels = c("CTRL", "PMOS")); names(group) <- rownames(meta)
dds_all <- DESeqDataSetFromMatrix(round(m), data.frame(group = group, row.names = colnames(m)), ~ group)
dds_all <- estimateSizeFactors(dds_all); sf <- sizeFactors(dds_all)
cat("size factors DESeq2 (R) vs ours (python): max |diff| =", max(abs(sf - sf_py[names(sf), 1])), "\n")
expressed <- intersect(ols$gene, rownames(m)); cat("genes tested:", length(expressed), "of", nrow(ols), "expressed in the OLS table\n")
nrm <- t(t(m) / sf); xist <- log2(nrm["XIST", ] + 1)
keep <- colnames(m) %in% lines
mm <- m[expressed, keep]; grp <- group[keep]; xi <- xist[keep]
res <- list()
for (adj in c(FALSE, TRUE)) {
  lab <- if (adj) "group_xist" else "group"
  cd <- data.frame(group = grp, xist = xi, row.names = colnames(mm))
  dds <- DESeqDataSetFromMatrix(round(mm), cd, if (adj) ~ xist + group else ~ group); sizeFactors(dds) <- sf[colnames(mm)]
  dds <- DESeq(dds, quiet = TRUE); r <- results(dds, contrast = c("group", "PMOS", "CTRL"))
  write.table(data.frame(gene = rownames(r), log2FC = r$log2FoldChange, stat = r$stat, p = r$pvalue, padj = r$padj, baseMean = r$baseMean), file.path(out, paste0("deseq2_", lab, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("DESeq2 ~ %s: p < 0.05 %d | padj < 0.1 %d (of %d tested)\n", lab, sum(r$pvalue < 0.05, na.rm = TRUE), sum(r$padj < 0.1, na.rm = TRUE), nrow(r)))
  dge <- DGEList(counts = m[, keep]); dge <- calcNormFactors(dge); dge <- dge[expressed, , keep.lib.sizes = TRUE]
  design <- if (adj) model.matrix(~ xi + grp) else model.matrix(~ grp)
  v <- voom(dge, design, plot = FALSE); fit <- eBayes(lmFit(v, design)); tt <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
  write.table(data.frame(gene = rownames(tt), log2FC = tt$logFC, stat = tt$t, p = tt$P.Value, padj = tt$adj.P.Val), file.path(out, paste0("limma_voom_", lab, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("limma-voom ~ %s: p < 0.05 %d | adj.P < 0.1 %d (of %d tested)\n", lab, sum(tt$P.Value < 0.05), sum(tt$adj.P.Val < 0.1), nrow(tt)))
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
