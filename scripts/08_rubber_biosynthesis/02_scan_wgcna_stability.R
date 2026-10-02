# =====================================================================
# S8-A: Rubber WGCNA parameter-sensitivity grid (power x mergeCutHeight)
#       -- R2-Major7 / supplementary table 8, paragraph 1
#
# Purpose:
#   Measure how sensitive the rubber candidate pairs are to the WGCNA parameters,
#   against the published results.
# Pipeline (= published WGCNA_yield.r):
#   rubber_count_5_fpkm.txt (23,163 x 15) -> MAD > 0 (21,303)
#   -> signed network via blockwiseModules(minModuleSize = 150, deepSplit = 2,
#      fixed seed)
#   -> pair EuNCP_trans with rubber genes within each module
#      (rubber_gca.txt, the original 52-entry list)
#   -> TOM threshold = min(quantile_q0.80, median + 2 x MAD)
#   -> Pearson r/p (n = 15) -> significant = r >= 0.7 & p <= 0.05
#   -> retention rate vs the published 85 NCPs / 301 pairs
#
# Inputs:
#   From REPRO_REJ: 02figure/figure5/rubber/WGCNA/rubber_count_5_fpkm.txt and
#                   02figure/figure5/rubber/WGCNA/rubber_gca.txt (original 52 entries)
#   From REPRO_ROOT: recons/wgcna_stability_package/published_85_NCP.txt and
#                    recons/wgcna_stability_package/published_301_pairs.csv
# Outputs:
#   output/S8_stability/out_rubber_grid/wgcna_stability_summary.csv
#   + per-parameter significant-pair CSVs sig_pairs_p<power>_mc<mergeCut>.csv
# Usage:
#   Rscript 02_scan_wgcna_stability.R
# Notes:
#   Must use the original 52-entry rubber_gca.txt (this script already points to
#   the archived original). Grids run earlier with the corrupted 75-entry list
#   (61-84/85 etc.) are void (RESULTS.md *** final verdict).
#   The 23k x 23k TOM requires >= 32 GB RAM; about 0.5-3 h per parameter set
#   (depends on threads).
#   Verification target: the baseline (22, 0.3) should fully reproduce 85/85 +
#   294/294; the full grid retains 56-85/85 (supplementary table 8, paragraph 1).
# =====================================================================

suppressMessages(library(WGCNA))
NTHREADS <- 8   # adjust per machine; set to 80 on an 80-thread device
allowWGCNAThreads(nThreads = NTHREADS)

REJ <- Sys.getenv("REPRO_REJ", unset = "./data/rejected_20260204")
REF <- file.path(Sys.getenv("REPRO_ROOT", unset = "."), "recons/wgcna_stability_package")
outRoot <- Sys.getenv("REPRO_OUT", unset = "./output")
OUT <- file.path(outRoot, "S8_stability", "out_rubber_grid")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

fpkm_file <- file.path(REJ, "02figure/figure5/rubber/WGCNA/rubber_count_5_fpkm.txt")
gca_file  <- file.path(REJ, "02figure/figure5/rubber/WGCNA/rubber_gca.txt")   # original 52 entries
pub85_file  <- file.path(REF, "published_85_NCP.txt")
pub301_file <- file.path(REF, "published_301_pairs.csv")

dataExpr0 <- read.table(fpkm_file, sep = "\t", row.names = 1, header = TRUE, quote = "", comment.char = "")
dataExpr0 <- as.data.frame(t(dataExpr0))
cat("input matrix (samples x features):", dim(dataExpr0), "\n")
m.mad <- apply(dataExpr0, 2, mad)
dataExpr <- dataExpr0[, m.mad > 0]
cat("features after MAD>0 filter:", ncol(dataExpr), "(expect 21,303)\n")

gca <- readLines(gca_file); pub85 <- readLines(pub85_file); pub301 <- read.csv(pub301_file)
cat("rubber_gca:", length(gca), "(expect 52) | published 85:", length(pub85), "| published pairs:", nrow(pub301), "\n")

SETS <- list(c(22, 0.3), c(20, 0.3), c(24, 0.3), c(22, 0.2), c(22, 0.4))

summary <- data.frame()
for (st in SETS) {
  P <- st[1]; MC <- st[2]
  tag <- sprintf("p%d_mc%s", P, gsub("\\.", "", as.character(MC)))
  cat("\n==== power =", P, "| mergeCutHeight =", MC, "====\n"); flush.console()
  t0 <- Sys.time()
  tmpTom <- file.path(tempdir(), paste0("TOM_", tag))
  net <- blockwiseModules(
    dataExpr, power = P, maxBlockSize = ncol(dataExpr),
    TOMType = "signed", networkType = "signed",
    minModuleSize = 150, mergeCutHeight = MC, deepSplit = 2,
    reassignThreshold = 0, numericLabels = TRUE, pamRespectsDendro = TRUE,
    saveTOMs = TRUE, saveTOMFileBase = tmpTom,
    verbose = 0, randomSeed = 20260906
  )
  colors <- labels2colors(net$colors)
  nmod <- length(unique(colors[colors != "grey"]))
  load(net$TOMFiles[1], verbose = FALSE)
  TOM <- as.matrix(TOM)
  colnames(TOM) <- rownames(TOM) <- colnames(dataExpr)
  novelIds <- grep("^EuNCP_trans", colnames(TOM), value = TRUE)
  gcaIds <- intersect(gca, colnames(TOM))
  sig_pairs <- data.frame()
  for (m in unique(net$colors)) {
    if (m == 0) next
    modGenes <- colnames(TOM)[net$colors == m]
    mod_novel <- intersect(novelIds, modGenes)
    mod_gca <- intersect(gcaIds, modGenes)
    if (length(mod_novel) == 0 || length(mod_gca) == 0) next
    modTOM <- TOM[modGenes, modGenes, drop = FALSE]
    glob <- quantile(modTOM, 0.8, na.rm = TRUE)
    rob <- median(modTOM, na.rm = TRUE) + 2 * mad(modTOM, na.rm = TRUE)
    thr <- min(glob, rob)
    for (nv in mod_novel) for (gn in mod_gca) {
      w <- TOM[nv, gn]
      if (w >= thr) {
        r <- cor(dataExpr[[nv]], dataExpr[[gn]], use = "p")
        p <- corPvalueStudent(r, nrow(dataExpr))
        if (r >= 0.7 && p <= 0.05) {
          sig_pairs <- rbind(sig_pairs, data.frame(novel = nv, annotated = gn,
                                                   r = r, p = p, module = labels2colors(m)))
        }
      }
    }
  }
  un <- unique(sig_pairs$novel)
  keep <- length(intersect(un, pub85))
  keep_pairs <- length(intersect(paste(sig_pairs$novel, sig_pairs$annotated),
                                 paste(pub301$novel, pub301$annotated)))
  rown <- data.frame(power = P, mergeCut = MC, n_modules = nmod,
                     sig_pairs = nrow(sig_pairs), unique_NCP = length(un),
                     retained_of_85 = keep, retained_pairs_of_301 = keep_pairs,
                     minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  summary <- rbind(summary, rown)
  write.csv(sig_pairs, file.path(OUT, paste0("sig_pairs_", tag, ".csv")), row.names = FALSE)
  cat(sprintf("modules=%d | sig pairs=%d | unique NCP=%d | retained of 85=%d | of 301 pairs=%d | %.1f min\n",
              nmod, nrow(sig_pairs), length(un), keep, keep_pairs, rown$minutes))
  flush.console()
  rm(TOM, net, sig_pairs); invisible(gc())
}
print(summary, row.names = FALSE)
write.csv(summary, file.path(OUT, "wgcna_stability_summary.csv"), row.names = FALSE)
cat("\nsaved:", file.path(OUT, "wgcna_stability_summary.csv"), "\n")
