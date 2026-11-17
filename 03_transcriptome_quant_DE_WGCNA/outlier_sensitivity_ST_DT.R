#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_root <- if (length(script_args)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_args[[1]])), ".."), mustWork = FALSE)
} else {
  normalizePath(".", mustWork = FALSE)
}
BASE <- Sys.getenv("EU_ANALYSIS_DIR", file.path(repo_root, "output", "03_ms_trans_analysis"))


args <- commandArgs(trailingOnly = TRUE)
BETA_FIXED <- if (length(args) >= 1 && nzchar(args[1])) as.numeric(args[1]) else NA_real_
OUT <- if (length(args) >= 2) args[2] else
       file.path(BASE, "02_wgcna/ST_DT_outlier_removed")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUTLIER <- "SRR14975289"

design <- read.delim(file.path(BASE, "00_input/sample_design.tsv"),
                     stringsAsFactors = FALSE)
d <- read.delim(file.path(BASE, "02_wgcna/input_vst/ST_DT_vst.tsv"),
                row.names = 1, check.names = FALSE)
pass <- readLines(file.path(BASE, "02_wgcna/input_vst/ST_DT_tpm_pass.txt"))
tg <- intersect(rownames(d), pass[nzchar(pass)])
pd <- design[design$panel == "ST_DT", ]
pd <- pd[pd$sample != OUTLIER, ]
datExpr <- as.matrix(t(d[tg, pd$sample]))
cat("After outlier removal:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

sft <- pickSoftThreshold(datExpr, powerVector = 1:30,
                         networkType = "signed", verbose = 0)
beta <- if (!is.na(BETA_FIXED)) BETA_FIXED else sft$powerEstimate
cat(if (!is.na(BETA_FIXED)) "fixed" else "reselected", "beta =", beta,
    " (original 15-sample beta=12)\n")
write.table(data.frame(chosen_beta = beta, rule = ifelse(is.na(BETA_FIXED),
            "official reselection", "fixed beta=12 control"), sft$fitIndices, check.names = FALSE),
            file.path(OUT, "softthreshold.tsv"), sep = "\t", row.names = FALSE)

adj <- adjacency(datExpr, power = beta, type = "signed")
TOM <- TOMsimilarity(adj, TOMType = "signed"); dissTOM <- 1 - TOM
tree <- hclust(as.dist(dissTOM), method = "average")
lab <- cutreeDynamic(dendro = tree, distM = dissTOM, deepSplit = 2,
                     pamRespectsDendro = FALSE, minClusterSize = 20)
m <- mergeCloseModules(datExpr, labels2colors(lab), cutHeight = 0.25)
module_col <- m$colors; names(module_col) <- colnames(datExpr)
MEs <- moduleEigengenes(datExpr, colors = unname(module_col))$eigengenes
MEs <- orderMEs(MEs); rownames(MEs) <- rownames(datExpr)


tr <- data.frame(row.names = pd$sample)
for (cc in sort(unique(pd$condition)))
  tr[[cc]] <- as.integer(pd$condition == cc)
n <- nrow(pd)
dh <- setNames(rep(NA_real_, n), pd$sample)
sh <- dh
dh[pd$condition %in% c("CK","DT_48","DT_144")] <-
  c(CK=0, DT_48=48, DT_144=144)[pd$condition[pd$condition %in% c("CK","DT_48","DT_144")]]
sh[pd$condition %in% c("CK","ST_48","ST_144")] <-
  c(CK=0, ST_48=48, ST_144=144)[pd$condition[pd$condition %in% c("CK","ST_48","ST_144")]]
tr$drought_hours <- unname(dh); tr$salt_hours <- unname(sh)
tr <- tr[rownames(MEs), , drop = FALSE]

mod_names <- sub("^ME", "", colnames(MEs))
corM <- matrix(NA_real_, length(mod_names), ncol(tr),
               dimnames = list(mod_names, colnames(tr)))
pM <- corM
for (mi in seq_along(mod_names)) for (ti in seq_len(ncol(tr))) {
  xy <- data.frame(me = MEs[, mi], tr = tr[, ti])
  xy <- xy[stats::complete.cases(xy), ]
  if (nrow(xy) < 4 || sd(xy$tr) == 0) next
  ct <- suppressWarnings(cor.test(xy$me, xy$tr))
  corM[mi, ti] <- ct$estimate; pM[mi, ti] <- ct$p.value
}
padjM <- matrix(NA_real_, nrow(pM), ncol(pM), dimnames = dimnames(pM))
pv <- pM[!is.na(pM)]
if (length(pv) > 1) padjM[!is.na(pM)] <- p.adjust(pv, method = "BH")

write.table(round(corM, 3), file.path(OUT, "ME_trait_correlation.tsv"), sep = "\t", na = "")
write.table(signif(pM, 3), file.path(OUT, "ME_trait_pvalue.tsv"), sep = "\t", na = "")
write.table(signif(padjM, 3), file.path(OUT, "ME_trait_padj.tsv"), sep = "\t", na = "")
kme <- signedKME(datExpr, MEs)
own <- sapply(names(module_col), function(g) {
  col <- paste0("kME", module_col[[g]])
  if (col %in% colnames(kme)) kme[g, col] else NA
})
write.table(data.frame(transcript_id = names(module_col),
                       module = unname(module_col),
                       kME_own_module = round(unname(own), 3)),
            file.path(OUT, "gene_module_assignment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
msz <- table(module_col)
write.table(data.frame(module = names(msz), n_genes = as.integer(msz)),
            file.path(OUT, "module_sizes.tsv"), sep = "\t", row.names = FALSE)
write.table(data.frame(panel = "ST_DT_outlier_removed", n_samples = nrow(datExpr),
                       n_genes = ncol(datExpr), beta = beta,
                       mean_k_at_beta = round(sft$fitIndices$mean.k.[sft$fitIndices$Power == beta], 1),
                       beta_rule = "default powerEstimate", n_modules = length(msz),
                       n_sig_module_trait_FDR05 = sum(padjM < 0.05 & abs(corM) >= 0.6, na.rm = TRUE),
                       n_sample_outliers = 0, exploratory = "N"),
            file.path(OUT, "panel_meta.tsv"), sep = "\t", row.names = FALSE)


old_gm <- read.delim(file.path(BASE, "02_wgcna/ST_DT/gene_module_assignment.tsv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
old_mod <- setNames(old_gm$module, old_gm$transcript_id)
common <- intersect(names(old_mod), names(module_col))

old_c <- old_mod[common]; new_c <- module_col[common]
g_old <- split(names(old_c), old_c)
g_new <- split(names(new_c), new_c)

match_rows <- data.frame()
for (mo in names(g_old)) {
  jac <- sapply(g_new, function(gx) length(intersect(g_old[[mo]], gx)) /
                                      length(union(g_old[[mo]], gx)))
  best <- which.max(jac)
  match_rows <- rbind(match_rows, data.frame(
    module_15s = mo, n_genes_15s = length(g_old[[mo]]),
    best_match_14s = names(jac)[best], n_genes_14s = length(g_new[[best]]),
    jaccard = round(unname(jac[best]), 3)))
}
match_rows <- match_rows[order(-match_rows$jaccard), ]
write.table(match_rows, file.path(OUT, "module_matching_vs_15s.tsv"),
            sep = "\t", row.names = FALSE)

old_cor <- as.matrix(read.delim(file.path(BASE, "02_wgcna/ST_DT/ME_trait_correlation.tsv"),
                                row.names = 1, check.names = FALSE, na.strings = "NA"))
old_padj <- as.matrix(read.delim(file.path(BASE, "02_wgcna/ST_DT/ME_trait_padj.tsv"),
                                 row.names = 1, check.names = FALSE, na.strings = "NA"))
sig_old <- which(old_padj < 0.05 & abs(old_cor) >= 0.6, arr.ind = TRUE)
mmap <- setNames(match_rows$best_match_14s, match_rows$module_15s)
stability <- data.frame()
for (i in seq_len(nrow(sig_old))) {
  mo <- rownames(old_cor)[sig_old[i, 1]]
  trt <- colnames(old_cor)[sig_old[i, 2]]
  mn <- if (mo %in% names(mmap)) mmap[[mo]] else NA
  new_r <- if (!is.na(mn) && mn %in% rownames(corM) && trt %in% colnames(corM))
             corM[mn, trt] else NA
  new_f <- if (!is.na(mn) && mn %in% rownames(padjM) && trt %in% colnames(padjM))
             padjM[mn, trt] else NA
  stability <- rbind(stability, data.frame(
    trait = trt, module_15s = mo, r_15s = round(old_cor[mo, trt], 3),
    matched_module_14s = mn, r_14s = round(new_r, 3),
    FDR_14s = signif(new_f, 2),
    still_FDR_sig = ifelse(!is.na(new_f) & new_f < 0.05 & abs(new_r) >= 0.6, "Y", "N")))
}
write.table(stability, file.path(OUT, "association_stability_vs_15s.tsv"),
            sep = "\t", row.names = FALSE)
cat("\n== Outlier sensitivity summary ==\n")
cat("15-sample module count:", length(g_old), "; 14-sample module count:", length(g_new), "\n")
cat("Original FDR-significant associations:", nrow(stability), "; retained FDR-significant with |r|>=0.6 in 14 samples:",
    sum(stability$still_FDR_sig == "Y"), "\n")
cat("Median matched-module Jaccard:", round(median(match_rows$jaccard), 3), "\n")
