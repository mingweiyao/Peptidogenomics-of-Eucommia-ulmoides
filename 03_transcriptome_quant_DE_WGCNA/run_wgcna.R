#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(WGCNA); library(pheatmap)
})
options(stringsAsFactors = FALSE)

script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_root <- if (length(script_args)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_args[[1]])), ".."), mustWork = FALSE)
} else {
  normalizePath(".", mustWork = FALSE)
}
BASE <- Sys.getenv("EU_ANALYSIS_DIR", file.path(repo_root, "output", "03_ms_trans_analysis"))
VSTDIR <- file.path(BASE, "02_wgcna/input_vst")
design <- read.delim(file.path(BASE, "00_input/sample_design.tsv"),
                     stringsAsFactors = FALSE)

pick_beta <- function(datExpr) {
  sft <- pickSoftThreshold(datExpr, powerVector = 1:30,
                           networkType = "signed", verbose = 0)
  fi <- sft$fitIndices
  n_genes <- ncol(datExpr)
  k_at <- setNames(fi$mean.k., fi$Power)
  rule <- sprintf("Official powerEstimate (first beta=%s with sign(slope)*R2 > 0.85)",
                  sft$powerEstimate)
  beta <- sft$powerEstimate
  if (is.null(beta) || beta < 2 || k_at[as.character(beta)] > 0.1 * n_genes) {
    ok <- which(fi$Power >= 2 & fi$SFT.R.sq >= 0.8 & fi$slope < 0 &
                fi$mean.k. <= 0.1 * n_genes)
    if (length(ok)) {
      beta <- fi$Power[ok[1]]
      rule <- sprintf("Sparsity-adjusted rule: lowest beta >= 2 with R2 >= 0.8, slope < 0, and mean.k <= %.0f",
                      0.1 * n_genes)
    } else if (is.null(beta)) {
      ok2 <- which(fi$slope < 0)
      beta <- if (length(ok2)) fi$Power[ok2[which.max(fi$SFT.R.sq[ok2])]] else 6
      rule <- "No beta meets the criteria; select the negative-slope beta with highest R2"
    }
  }
  list(sft = sft, beta = beta, rule = rule,
       mean_k = unname(k_at[as.character(beta)]))
}

build_traits <- function(panel, pd) {
  n <- nrow(pd)
  conds <- sort(unique(pd$condition))
  tr <- data.frame(row.names = pd$sample)
  for (cc in conds) tr[[cc]] <- as.integer(pd$condition == cc)
  if (panel == "leaf_position") {
    ord <- c(top = 1, "close to tip" = 2, "middle and lower" = 3,
             lower = 4, skin = 5)
    tr$position_ord <- unname(ord[pd$condition])
  } else if (panel == "leaf_dev") {
    ord <- c("leaf buds" = 1, "growing leaves" = 2,
             "young leaves" = 3, "mature leaves" = 4)
    tr$dev_order <- unname(ord[pd$condition])
  } else if (panel == "ST_DT") {

    dh <- setNames(rep(NA_real_, n), pd$sample)
    sh <- dh
    dh[pd$condition %in% c("CK", "DT_48", "DT_144")] <-
      c(CK = 0, DT_48 = 48, DT_144 = 144)[pd$condition[pd$condition %in% c("CK","DT_48","DT_144")]]
    sh[pd$condition %in% c("CK", "ST_48", "ST_144")] <-
      c(CK = 0, ST_48 = 48, ST_144 = 144)[pd$condition[pd$condition %in% c("CK","ST_48","ST_144")]]
    tr$drought_hours <- unname(dh)
    tr$salt_hours <- unname(sh)
  }
  as.data.frame(tr)
}

modules_at_beta <- function(datExpr, beta, minClSize) {
  adj <- adjacency(datExpr, power = beta, type = "signed")
  TOM <- TOMsimilarity(adj, TOMType = "signed")
  dissTOM <- 1 - TOM
  tree <- hclust(as.dist(dissTOM), method = "average")
  lab <- cutreeDynamic(dendro = tree, distM = dissTOM, deepSplit = 2,
                       pamRespectsDendro = FALSE, minClusterSize = minClSize)
  m <- mergeCloseModules(datExpr, labels2colors(lab), cutHeight = 0.25)
  cols <- m$colors
  names(cols) <- colnames(datExpr)
  cols
}

run_panel <- function(panel, minClSize, exploratory) {
  d <- read.delim(file.path(VSTDIR, paste0(panel, "_vst.tsv")),
                  row.names = 1, check.names = FALSE)
  pass <- readLines(file.path(VSTDIR, paste0(panel, "_tpm_pass.txt")))
  tg <- intersect(rownames(d), pass[nzchar(pass)])
  d <- d[tg, ]
  pd <- design[design$panel == panel, ]
  datExpr <- as.matrix(t(d[, pd$sample]))
  cat("\n=====", panel, ":", nrow(datExpr), "samples x", ncol(datExpr),
      "genes", if (exploratory) "(exploratory)" else "", "=====\n")

  outd <- file.path(BASE, "02_wgcna", panel)
  dir.create(outd, recursive = TRUE, showWarnings = FALSE)


  pb <- pick_beta(datExpr)
  beta <- pb$beta


  sample_cor <- cor(t(datExpr))
  sample_adj <- (0.5 * (1 + sample_cor))^beta
  k_samp <- colSums(sample_adj, na.rm = TRUE)
  Z.k <- as.numeric(scale(k_samp))
  sample_qc <- data.frame(sample = rownames(datExpr),
                          condition = pd$condition,
                          k_sample = round(k_samp, 2), Z_k = round(Z.k, 2),
                          outlier = ifelse(abs(Z.k) > 2, "Y", "N"))
  stree <- hclust(as.dist(1 - sample_cor), method = "average")
  png(file.path(outd, "sample_tree.png"), width = 1200, height = 500)
  plot(stree, main = paste0(panel, " sample clustering"),
       sub = paste("outliers(|Z.k|>2):",
                   ifelse(any(sample_qc$outlier == "Y"),
                          paste(sample_qc$sample[sample_qc$outlier == "Y"],
                                collapse = ","), "none")),
       xlab = "")
  dev.off()
  cat("  beta =", beta, "| mean_k =", signif(pb$mean_k, 1),
      "| rule:", pb$rule, "\n")
  write.table(sample_qc, file.path(outd, "sample_qc.tsv"),
              sep = "\t", row.names = FALSE)

  st <- data.frame(chosen_beta = beta, selection_rule = pb$rule,
                   pb$sft$fitIndices, check.names = FALSE)
  write.table(st, file.path(outd, "softthreshold.tsv"),
              sep = "\t", row.names = FALSE)

  module_col <- modules_at_beta(datExpr, beta, minClSize)
  names(module_col) <- colnames(datExpr)
  MEs <- moduleEigengenes(datExpr, colors = unname(module_col))$eigengenes
  MEs <- orderMEs(MEs)
  rownames(MEs) <- rownames(datExpr)


  tr <- build_traits(panel, pd)
  tr <- tr[rownames(MEs), , drop = FALSE]
  mod_names <- sub("^ME", "", colnames(MEs))
  corM <- matrix(NA_real_, length(mod_names), ncol(tr),
                 dimnames = list(mod_names, colnames(tr)))
  pM <- corM; nM <- corM
  for (mi in seq_along(mod_names)) for (ti in seq_len(ncol(tr))) {
    xy <- data.frame(me = MEs[, mi], tr = tr[, ti])
    xy <- xy[stats::complete.cases(xy), ]
    if (nrow(xy) < 4 || sd(xy$tr) == 0) next
    ct <- suppressWarnings(cor.test(xy$me, xy$tr))
    corM[mi, ti] <- ct$estimate; pM[mi, ti] <- ct$p.value; nM[mi, ti] <- nrow(xy)
  }
  padjM <- matrix(NA_real_, nrow(pM), ncol(pM), dimnames = dimnames(pM))
  pv <- pM[!is.na(pM)]
  if (length(pv) > 1) padjM[!is.na(pM)] <- p.adjust(pv, method = "BH")

  write.table(round(corM, 3), file.path(outd, "ME_trait_correlation.tsv"),
              sep = "\t", na = "")
  write.table(signif(pM, 3), file.path(outd, "ME_trait_pvalue.tsv"),
              sep = "\t", na = "")
  write.table(signif(padjM, 3), file.path(outd, "ME_trait_padj.tsv"),
              sep = "\t", na = "")
  write.table(nM, file.path(outd, "ME_trait_n.tsv"), sep = "\t", na = "")

  kme <- signedKME(datExpr, MEs)
  own <- sapply(names(module_col), function(g) {
    col <- paste0("kME", module_col[[g]])
    if (col %in% colnames(kme)) kme[g, col] else NA
  })
  write.table(data.frame(transcript_id = names(module_col),
                         module = unname(module_col),
                         kME_own_module = round(unname(own), 3)),
              file.path(outd, "gene_module_assignment.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(sample = rownames(MEs), MEs, check.names = FALSE),
              file.path(outd, "module_eigengenes.csv"),
              sep = ",", quote = FALSE, row.names = FALSE)
  msz <- table(module_col)
  write.table(data.frame(module = names(msz), n_genes = as.integer(msz)),
              file.path(outd, "module_sizes.tsv"), sep = "\t", row.names = FALSE)

  sig <- which(padjM < 0.05 & abs(corM) >= 0.6, arr.ind = TRUE)


  sens <- data.frame()
  if (nrow(sig) > 0) {
    key_mods <- unique(rownames(sig))
    genes_of <- function(cols) split(names(cols), cols)
    g0 <- genes_of(module_col)
    for (b2 in c(max(1, beta - 2), beta + 2)) {
      if (b2 == beta) next
      col2 <- try(modules_at_beta(datExpr, b2, minClSize), silent = TRUE)
      if (inherits(col2, "try-error")) next
      g2 <- genes_of(col2)
      for (km in key_mods) {
        gk <- g0[[km]]
        jac <- sapply(g2, function(gx) length(intersect(gk, gx)) / length(union(gk, gx)))
        best <- which.max(jac)
        sens <- rbind(sens, data.frame(
          module = km, beta_used = b2, n_genes = length(gk),
          best_match_module = names(jac)[best],
          best_match_size = length(g2[[best]]),
          jaccard = round(unname(jac[best]), 3)))
      }
    }
  }
  if (nrow(sens)) write.table(sens, file.path(outd, "beta_sensitivity.tsv"),
                              sep = "\t", row.names = FALSE)


  write.table(data.frame(
      panel = panel, n_samples = nrow(datExpr), n_genes = ncol(datExpr),
      beta = beta, mean_k_at_beta = round(pb$mean_k, 1),
      beta_rule = pb$rule, n_modules = length(msz),
      n_sig_module_trait_FDR05 = nrow(sig),
      n_sample_outliers = sum(sample_qc$outlier == "Y"),
      exploratory = ifelse(exploratory, "Y", "N")),
    file.path(outd, "panel_meta.tsv"), sep = "\t", row.names = FALSE)

  corP <- corM; pP <- padjM
  corP[is.na(corP)] <- 0; pP[is.na(pP)] <- 1
  txt <- matrix(paste0(signif(corP, 2), "\n(FDR ", signif(pP, 1), ")"), nrow(corP))
  pheatmap(corP, display_numbers = txt, cluster_rows = TRUE,
           cluster_cols = FALSE, main = paste0(panel, ": ME-trait (r, FDR)"),
           filename = file.path(outd, "me_trait_heatmap.png"),
           width = 3 + 1.1 * ncol(corP), height = 2.5 + 0.55 * nrow(corP),
           fontsize_number = 7)
  cat("  Modules:", paste0(names(msz), "(", as.integer(msz), ")",
      collapse = " "), "\n")
  cat("  FDR<0.05 and  |r|>=0.6 module-trait associations:", nrow(sig),
      "; outlier samples:", sum(sample_qc$outlier == "Y"), "\n")
}

run_panel("ST_DT", 20, FALSE)
run_panel("leaf_position", 20, FALSE)
run_panel("tissue", 20, FALSE)
run_panel("leaf_dev", 20, TRUE)
run_panel("MT_melatonin", 15, TRUE)
cat("\nWGCNA v4 completed\n")
