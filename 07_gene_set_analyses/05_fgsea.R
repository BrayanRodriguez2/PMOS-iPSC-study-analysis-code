#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Release file: 07_gene_set_analyses/05_fgsea.R   (original step script: 74_fgsea.R)
# Produces: fgsea preranked GSEA (custom sets and MSigDB Hallmark).
# Run: called by the SLURM script of the same number: Rscript 07_gene_set_analyses/05_fgsea.R <args>
# Paths are variables from config/config.sh; individual-level data are NOT part of this repository.
# ---------------------------------------------------------------------------
# Step 74 - preranked GSEA with fgsea (fgseaMultilevel) on the rankings of step 73, for the custom gene sets (GWAS consensus 33, literature 60, Min et al. if the lists exist)
# and the MSigDB Hallmark collection (50 sets; msigdbr, human symbols). minSize 5, maxSize 500 (counted on the genes present in the ranking), set.seed(2026), serial (no parallel
# RNG differences), scoreType "std". padj = Benjamini-Hochberg WITHIN each collection (custom sets: within the collection they belong to; Hallmark: over its 50 sets).
# A positive NES = the set is enriched among the genes HIGHER in PMOS than in control lines.
# Usage: Rscript 07_gene_set_analyses/05_fgsea.R <out_dir (with rankings/ and gene_sets/)>
suppressPackageStartupMessages({ library(fgsea); library(msigdbr); library(ggplot2); library(data.table) })
out <- commandArgs(trailingOnly = TRUE)[1]
dir.create(file.path(out, "tables"), showWarnings = FALSE); dir.create(file.path(out, "figures"), showWarnings = FALSE)
cs <- read.delim(file.path(out, "gene_sets", "custom_gene_sets.tsv"), stringsAsFactors = FALSE)
custom <- split(cs, cs$set); coll_of <- tapply(cs$collection, cs$set, function(x) x[1])
pw_custom <- lapply(custom, function(d) unique(d$gene))
hm <- msigdbr(species = "Homo sapiens", category = "H")
pw_h <- split(hm$gene_symbol, hm$gs_name); pw_h <- lapply(pw_h, unique)
write.table(unique(data.frame(set = hm$gs_name, gene = hm$gene_symbol)), file.path(out, "gene_sets", "hallmark.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Hallmark sets:", length(pw_h), "| msigdbr", as.character(packageVersion("msigdbr")), "| fgsea", as.character(packageVersion("fgsea")), "\n")
rk_files <- sort(list.files(file.path(out, "rankings"), pattern = "\\.tsv$", full.names = TRUE)); all_res <- list()
for (f in rk_files) {
  rname <- sub("\\.tsv$", "", basename(f)); rk <- read.delim(f, stringsAsFactors = FALSE); stats <- setNames(rk$stat, rk$gene); stats <- stats[!is.na(stats)]
  for (coll in c("custom", "Hallmark")) {
    pw <- if (coll == "custom") pw_custom else pw_h
    set.seed(2026)
    res <- fgseaMultilevel(pathways = pw, stats = stats, minSize = 5, maxSize = 500, scoreType = "std", eps = 0, BPPARAM = BiocParallel::SerialParam())
    res <- as.data.frame(res); res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), "")
    res$n_in_set <- vapply(res$pathway, function(p) length(pw[[p]]), 1L)
    res$collection <- if (coll == "custom") unname(coll_of[res$pathway]) else "Hallmark"
    res$padj_BH_collection <- ave(res$pval, res$collection, FUN = function(p) p.adjust(p, "BH"))
    res$ranking <- rname
    res <- res[, c("ranking", "collection", "pathway", "n_in_set", "size", "ES", "NES", "pval", "padj_BH_collection", "leadingEdge")]; names(res)[5] <- "n_present_in_ranking"
    all_res[[paste(rname, coll)]] <- res
    # sets below minSize in the ranking are not tested: list them
    if (coll == "custom") { miss <- setdiff(names(pw), res$pathway); if (length(miss)) cat(rname, ": not tested (< 5 genes present):", paste(miss, collapse = ", "), "\n") }
    # figures
    if (coll == "custom" && rname %in% c("deseq2_wald", "ols_signed", "deseq2_wald_xist", "deseq2_wald_no350")) {
      for (p in res$pathway) {
        r <- res[res$pathway == p, ]
        g <- plotEnrichment(pw[[p]], stats) + labs(title = sprintf("%s | ranking: %s", p, rname), subtitle = sprintf("%d of %d genes present; NES %.2f; p %.3g; padj %.3g", r$n_present_in_ranking, r$n_in_set, r$NES, r$pval, r$padj_BH_collection))
        ggsave(file.path(out, "figures", sprintf("enrichment_%s_%s.png", rname, p)), g, width = 6.5, height = 4, dpi = 160)
      }
    }
    if (coll == "Hallmark" && rname %in% c("deseq2_wald", "ols_signed", "deseq2_wald_xist", "deseq2_wald_no350")) {
      d <- res[order(res$NES), ]; d$pathway <- factor(sub("HALLMARK_", "", d$pathway), levels = sub("HALLMARK_", "", d$pathway)); d$sig <- ifelse(d$padj_BH_collection < 0.05, "padj < 0.05", ifelse(d$pval < 0.05, "p < 0.05", "n.s."))
      g <- ggplot(d, aes(NES, pathway, fill = sig)) + geom_col() + scale_fill_manual(values = c("padj < 0.05" = "#c0392b", "p < 0.05" = "#e59866", "n.s." = "#bbbbbb")) +
        labs(title = sprintf("MSigDB Hallmark, ranking: %s (positive = higher in PMOS)", rname), x = "NES", y = NULL, fill = NULL) + theme_bw(base_size = 8)
      ggsave(file.path(out, "figures", sprintf("hallmark_NES_%s.png", rname)), g, width = 7, height = 8, dpi = 160)
    }
  }
}
R <- do.call(rbind, all_res); rownames(R) <- NULL
write.table(R, file.path(out, "tables", "fgsea_results_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("results:", nrow(R), "rows\n"); print(R[R$collection != "Hallmark", c("ranking", "pathway", "n_in_set", "n_present_in_ranking", "NES", "pval", "padj_BH_collection")], row.names = FALSE)
for (rn in unique(R$ranking)) { h <- R[R$ranking == rn & R$collection == "Hallmark", ]; cat(sprintf("Hallmark %s: sets with p < 0.05: %d | padj < 0.05: %d | best padj %.3f (%s)\n", rn, sum(h$pval < 0.05), sum(h$padj_BH_collection < 0.05), min(h$padj_BH_collection), h$pathway[which.min(h$padj_BH_collection)])) }
writeLines(capture.output(sessionInfo()), file.path(out, "tables", "sessionInfo.txt"))
