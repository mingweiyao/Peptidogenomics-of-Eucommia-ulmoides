# =====================================================================
# S8-D v2: Stress WGCNA parameter-sensitivity grid (module-rule version,
#          anchorable to the archived baseline)
#
# Purpose:
#   Sweep power x mergeCutHeight for the stress WGCNA network under the
#   "trait-related module membership" candidate rule that was actually in effect
#   in the manuscript (reverse-confirmed from the archived heatmap PDF):
#     DT modules = modules whose eigengene is significantly correlated with CK or
#                  DT_144 (|r| > 0.5 and p < 0.05, corPvalueStudent)
#     ST modules = modules whose eigengene is significantly correlated with CK,
#                  ST_48 or ST_144
#   The baseline (16, 0.4) must exactly reproduce the archived 231/181 module-NCP
#   lists and the 30/25 candidates (script self-check baseline_pass). The four
#   variants report the true retention rate of the published candidates (X/30,
#   Y/25), not a variability range.
#   TOM is no longer saved/loaded (the rule is membership-based; TOM not needed).
#
# Inputs:
#   Anchor directory (default ./output/S8_stability) containing 6 CSVs whose
#   first column is named EuNCP_trans:
#     DE_fdr_DT.csv(46) DE_fdr_ST.csv(33) arch_mod_DT.csv(231) arch_mod_ST.csv(181)
#     cand_fdr_DT30.csv(30) cand_fdr_ST25.csv(25)
#   From REPRO_REJ: 02figure/figure6/WGCNA/dt_st_count_fpkm.txt and
#                   02figure/figure6/WGCNA/dt_st_trait.txt
# Outputs:
#   <REPRO_OUT>/S8_stability/out_stress_grid_v2/wgcna_stress_rule_stability.csv
#   + per-parameter candidate lists cand_{dt,st}_p<power>_mc<mergeCut>.csv
# Usage:
#   Rscript 02_scan_wgcna_stability.R [<anchor_dir>]
#   Env: REPRO_REJ = rejected-archive root; REPRO_OUT = output root
# Notes:
#   Differs from v1 (02_scan_wgcna_stability.R) only in the candidate rule above.
#   Requires >= 32 GB RAM; about 0.5-3 h per parameter set.
# =====================================================================

suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
LD <- if (length(args) >= 1) args[1] else file.path("./output/S8_stability")
REJ <- Sys.getenv("REPRO_REJ", unset = "./data/rejected_20260204")
outRoot <- Sys.getenv("REPRO_OUT", unset = "./output")
OUT <- file.path(outRoot, "S8_stability", "out_stress_grid_v2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) read.csv(file.path(LD, f), colClasses = "character")[[1]]
de_dt <- rd("DE_fdr_DT.csv");   de_st <- rd("DE_fdr_ST.csv")
arch_dt <- rd("arch_mod_DT.csv"); arch_st <- rd("arch_mod_ST.csv")
pub_dt <- rd("cand_fdr_DT30.csv"); pub_st <- rd("cand_fdr_ST25.csv")
cat("Input: FDR DE DT", length(de_dt), "| ST", length(de_st),
    "|| archived module NCP", length(arch_dt), "/", length(arch_st),
    "|| published candidates", length(pub_dt), "/", length(pub_st), "\n")

fpkm <- file.path(REJ, "02figure/figure6/WGCNA/dt_st_count_fpkm.txt")
traitf <- file.path(REJ, "02figure/figure6/WGCNA/dt_st_trait.txt")
dataExpr <- read.table(fpkm, sep = "\t", row.names = 1, header = TRUE, quote = "", comment.char = "")
dataExpr <- as.data.frame(t(dataExpr))
m.mad <- apply(dataExpr, 2, mad)
dataExpr <- dataExpr[, m.mad > 0]
gsg <- goodSamplesGenes(dataExpr, verbose = 0)
if (!gsg$allOK) dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
nSamples <- nrow(dataExpr)
trait <- read.table(traitf, sep = "\t", header = TRUE, row.names = 1)
trait <- trait[match(rownames(dataExpr), rownames(trait)), , drop = FALSE]
stopifnot(all(rownames(trait) == rownames(dataExpr)))
ncpAll <- grep("^EuNCP_trans", colnames(dataExpr), value = TRUE)
cat("Network input: samples", nSamples, "| features", ncol(dataExpr), "| NCP", length(ncpAll), "\n")

allowWGCNAThreads(nThreads = 8)

# Rule: module eigengene vs any trait in cols with |r| > 0.5 and p < 0.05 (corPvalueStudent)
selectModules <- function(MEs, cols) {
  mods <- setdiff(names(MEs), "MEgrey")
  hit <- character(0)
  for (m in mods) {
    for (cc in cols) {
      r <- cor(MEs[[m]], trait[[cc]], use = "p")
      p <- corPvalueStudent(r, nSamples)
      if (abs(r) > 0.5 && p < 0.05) { hit <- c(hit, sub("^ME", "", m)); break }
    }
  }
  sort(hit)
}
mods2ncp <- function(sel, colors) {
  if (!length(sel)) return(character(0))
  genes <- names(colors)[colors %in% sel]
  intersect(genes, ncpAll)
}

SETS <- list(c(16, 0.4), c(14, 0.4), c(18, 0.4), c(16, 0.3), c(16, 0.5))
summary <- data.frame()
for (st in SETS) {
  P <- st[1]; MC <- st[2]
  tag <- sprintf("p%d_mc%s", P, gsub("\\.", "", as.character(MC)))
  cat("\n==== power =", P, "| mergeCutHeight =", MC, "====\n"); flush.console()
  t0 <- Sys.time()
  net <- blockwiseModules(
    dataExpr, power = P, maxBlockSize = ncol(dataExpr),
    TOMType = "signed", networkType = "signed",
    minModuleSize = 200, mergeCutHeight = MC, deepSplit = 0,
    reassignThreshold = 0, numericLabels = TRUE, pamRespectsDendro = TRUE,
    saveTOMs = FALSE, verbose = 0, randomSeed = 20260916
  )
  colors <- labels2colors(net$colors); names(colors) <- colnames(dataExpr)
  MEs <- moduleEigengenes(dataExpr, colors)$eigengenes
  nmod <- length(setdiff(unique(colors), "grey"))
  sel_dt <- selectModules(MEs, c("CK", "DT_144"))
  sel_st <- selectModules(MEs, c("CK", "ST_48", "ST_144"))
  ncp_dt <- mods2ncp(sel_dt, colors); ncp_st <- mods2ncp(sel_st, colors)
  cand_dt <- intersect(ncp_dt, de_dt); cand_st <- intersect(ncp_st, de_st)
  ret_dt <- intersect(cand_dt, pub_dt); ret_st <- intersect(cand_st, pub_st)
  baseline_pass <- NA
  if (P == 16 && MC == 0.4) {
    ok1 <- setequal(ncp_dt, arch_dt); ok2 <- setequal(ncp_st, arch_st)
    ok3 <- setequal(cand_dt, pub_dt); ok4 <- setequal(cand_st, pub_st)
    baseline_pass <- ok1 && ok2 && ok3 && ok4
    cat("Baseline anchor self-check: DT module NCP == 231:", ok1, "| ST == 181:", ok2,
        "| candidates == 30:", ok3, "| == 25:", ok4, "=> PASS:", baseline_pass, "\n")
    if (!baseline_pass) cat("⚠ Baseline failed the anchor self-check - variant rows must not be read as retention rates; investigate first!\n")
  }
  write.csv(data.frame(EuNCP_trans = sort(cand_dt)),
            file.path(OUT, paste0("cand_dt_", tag, ".csv")), row.names = FALSE)
  write.csv(data.frame(EuNCP_trans = sort(cand_st)),
            file.path(OUT, paste0("cand_st_", tag, ".csv")), row.names = FALSE)
  rown <- data.frame(power = P, mergeCut = MC, n_modules = nmod,
                     sel_dt_modules = length(sel_dt), dt_ncp = length(ncp_dt),
                     dt_candidates = length(cand_dt), dt_retained_of30 = length(ret_dt),
                     sel_st_modules = length(sel_st), st_ncp = length(ncp_st),
                     st_candidates = length(cand_st), st_retained_of25 = length(ret_st),
                     baseline_pass = baseline_pass,
                     minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  summary <- rbind(summary, rown)
  cat(sprintf("modules=%d | DT: %d modules/%d NCP/%d candidates/%d retained | ST: %d modules/%d NCP/%d candidates/%d retained\n",
              nmod, length(sel_dt), length(ncp_dt), length(cand_dt), length(ret_dt),
              length(sel_st), length(ncp_st), length(cand_st), length(ret_st)))
  flush.console()
  invisible(gc())
}
print(summary, row.names = FALSE)
write.csv(summary, file.path(OUT, "wgcna_stress_rule_stability.csv"), row.names = FALSE)
cat("\nsaved:", file.path(OUT, "wgcna_stress_rule_stability.csv"), "\nDONE\n")
