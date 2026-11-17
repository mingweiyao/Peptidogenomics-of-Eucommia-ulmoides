#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))

script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_root <- if (length(script_args)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_args[[1]])), ".."), mustWork = FALSE)
} else {
  normalizePath(".", mustWork = FALSE)
}
BASE <- Sys.getenv("EU_ANALYSIS_DIR", file.path(repo_root, "output", "03_ms_trans_analysis"))
COUNTS <- Sys.getenv("EU_COUNTS_FILE", file.path(repo_root, "data", "rnaseq", "07new_transcript_quant", "transcript_count_matrix.csv"))
design <- read.delim(file.path(BASE, "00_input/sample_design.tsv"), stringsAsFactors = FALSE)
targets <- unique(readLines(file.path(BASE, "00_input/ncp_unique_transcripts.tsv")))
targets <- targets[nzchar(targets)]
ncp_map <- read.delim(file.path(BASE, "00_input/ncp_to_transcript.tsv"), stringsAsFactors = FALSE)
n_ncp <- table(ncp_map$transcript_id)

cnt <- read.csv(COUNTS, row.names = 1, check.names = FALSE)
stopifnot(all(design$sample %in% colnames(cnt)))
cat("Count matrix:", nrow(cnt), "x", ncol(cnt), "; target transcripts", length(targets), "\n")

dir.create(file.path(BASE, "02_wgcna/input_vst"), recursive = TRUE, showWarnings = FALSE)
san <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

panels <- list(
  list(key = "heat",          ref = "heat_CK",    extra = NULL),
  list(key = "CA_content",    ref = "low content", extra = NULL),
  list(key = "leaf_position", ref = "top",        extra = NULL),
  list(key = "ST_DT",         ref = "CK",         extra = NULL),
  list(key = "MT_melatonin",  ref = "CK",
       extra = list(alt = "Melatonin", base = "drought")),
  list(key = "tissue",        ref = "leaf",       extra = NULL),
  list(key = "leaf_dev",      ref = "leaf buds",  extra = NULL)
)

summary_rows <- list()
for (p in panels) {
  d <- design[design$panel == p$key, ]
  cd <- cnt[, d$sample]


  keep <- rowSums(cd >= 10) >= 2
  cd <- cd[keep, ]
  tg <- intersect(targets, rownames(cd))
  cat("====", p$key, ":", ncol(cd), "samples; after filtering:", nrow(cd), "genes;",
      "passing targets", length(tg), "/", length(targets), "====\n")

  lv <- unique(as.character(d$condition))
  lv <- c(p$ref, setdiff(sort(lv), p$ref))
  coldata <- data.frame(condition = factor(d$condition, levels = lv),
                        row.names = d$sample)
  dds <- DESeqDataSetFromMatrix(cd, coldata, ~ condition)
  dds <- DESeq(dds, quiet = TRUE)

  vsd <- vst(dds, nsub = min(1000, nrow(dds)), blind = FALSE)
  vm <- assay(vsd)[tg, ]
  write.table(data.frame(transcript_id = rownames(vm), vm, check.names = FALSE),
              file.path(BASE, "02_wgcna/input_vst", paste0(p$key, "_vst.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  outd <- file.path(BASE, "01_deseq2", p$key)
  dir.create(outd, recursive = TRUE, showWarnings = FALSE)

  comps <- lapply(setdiff(lv, p$ref), function(x) c(x, p$ref))
  if (!is.null(p$extra)) comps[[length(comps) + 1]] <- c(p$extra$alt, p$extra$base)
  nrep <- length(d$sample) / length(lv)

  for (cc in comps) {
    res <- results(dds, contrast = c("condition", cc[1], cc[2]))
    df <- data.frame(transcript_id = rownames(res), as.data.frame(res),
                     stringsAsFactors = FALSE, check.names = FALSE)
    df$locus <- sub("\\.[^.]+$", "", df$transcript_id)
    df$n_ncp <- as.integer(n_ncp[df$transcript_id]); df$n_ncp[is.na(df$n_ncp)] <- 0
    df <- df[tg, ]
    comp <- paste0(san(cc[1]), "_vs_", san(cc[2]))
    write.table(df, file.path(outd, paste0(comp, "_results.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    n_sig <- sum(df$padj < 0.05 & abs(df$log2FoldChange) >= 1, na.rm = TRUE)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      panel = p$key, comparison = comp, design = paste(lv, collapse = "/"),
      n_rep = nrep, n_target = nrow(df), n_sig = n_sig,
      n_up = sum(df$padj < 0.05 & df$log2FoldChange >= 1, na.rm = TRUE),
      n_down = sum(df$padj < 0.05 & df$log2FoldChange <= -1, na.rm = TRUE))
    cat(sprintf("  %-36s DEG: %3d (up %d / down %d)\n", comp, n_sig,
                summary_rows[[length(summary_rows)]]$n_up,
                summary_rows[[length(summary_rows)]]$n_down))
  }
}
write.table(do.call(rbind, summary_rows),
            file.path(BASE, "01_deseq2/deseq2_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("DESeq2 v2 completed (20 contrasts)\n")
