#Rscript # S8-B: TOM-quantile sensitivity, corrected (valid 52-entry RB list + reproduced network) - supplementary table 8, part 4
#
# Rule = min(quantile_Q, median+2xMAD) + Pearson r>=0.7 & p<=0.05, Q in {0.80,0.85,0.90,0.95,0.97}
# Prerequisite: run 01_run_wgcna_network.R first (or have finalTOM-block.1.RData + hub_genes module member tables ready)
#
# Verification (valid list + reproduced network; settled in addendum A, 2026-09-09):
#   q0.80 (baseline) 294 pairs / 85 NCPs, 100% retained; q0.85 221/71 (83.5%); q0.90 160/52 (61.2%);
#   q0.95 107/40 (47.1%); q0.97 98/40 (47.1%); every stricter-threshold set nests perfectly within the published 85.
#   (The old numbers "q0.90: 213; q0.95: 114; q0.97: 62" were run with the invalid list and are void.)
#
# Usage: REPRO_REJ=./data/rejected_20260204 REPRO_OUT=./output REPRO_ROOT=./data Rscript 03_scan_tom_threshold.R
# Output: output/S8_stability/sensitivity_TOMquantile_correctList.csv
suppressMessages(library(WGCNA))
outRoot <- Sys.getenv("REPRO_OUT", unset = "./output")
setwd(file.path(outRoot, "S8_stability"))
dir.create(".", showWarnings = FALSE, recursive = TRUE)

REJ <- Sys.getenv("REPRO_REJ", unset = "./data/rejected_20260204")
RERUN <- file.path(outRoot, "S6_rubber_candidates", "wgcna_rerun")

load(file.path(RERUN, "finalTOM-block.1.RData")); TOM <- as.matrix(TOM)
dataExpr <- read.table(file.path(REJ, "02figure/figure5/rubber/WGCNA/rubber_count_5_fpkm.txt"),
                       sep="\t", row.names=1, header=TRUE, quote="", comment.char="")
dataExpr <- as.data.frame(t(dataExpr))
m.mad <- apply(dataExpr, 2, mad); dataExpr <- dataExpr[, m.mad > 0]
gsg <- goodSamplesGenes(dataExpr, verbose=0)
if (!gsg$allOK) dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
colnames(TOM) <- rownames(TOM) <- colnames(dataExpr)
gca <- intersect(readLines(file.path(REJ, "02figure/figure5/rubber/WGCNA/rubber_gca.txt")), colnames(dataExpr))

# Module member tables come from the 01_run_wgcna_network.R rerun Step-8 by-product hub_genes/genes_<color>_trait_LF.txt (all genes of the 12 modules)
mod_colors <- c("turquoise","blue","brown","yellow","green","red","black","pink","magenta","purple","greenyellow","tan")
colors <- unlist(lapply(mod_colors, function(cc)
  setNames(rep(cc, length(readLines(file.path(RERUN, "hub_genes", paste0("genes_", cc, "_trait_LF.txt"))))),
           readLines(file.path(RERUN, "hub_genes", paste0("genes_", cc, "_trait_LF.txt"))))))
stopifnot(length(colors) == 21079)   # 21,303 - 224 grey

novel <- grep("^EuNCP_trans", colnames(TOM), value=TRUE)
pub85 <- readLines(file.path(Sys.getenv("REPRO_ROOT", unset = "."), "recons/wgcna_stability_package/published_85_NCP.txt"))
out <- data.frame(); summary <- data.frame()
for (Q in c(0.80, 0.85, 0.90, 0.95, 0.97)) {
  sig <- data.frame()
  for (m in unique(colors)) {
    gm <- names(colors)[colors == m]
    mn <- intersect(novel, gm); mg <- intersect(gca, gm)
    if (!length(mn) || !length(mg)) next
    mt <- TOM[gm, gm]
    thr <- min(quantile(mt, Q, na.rm=TRUE), median(mt, na.rm=TRUE) + 2*mad(mt, na.rm=TRUE))
    for (nv in mn) for (gn in mg) {
      w <- TOM[nv, gn]
      if (w >= thr) {
        r <- cor(dataExpr[[nv]], dataExpr[[gn]], use="p")
        p <- corPvalueStudent(r, nrow(dataExpr))
        if (r >= 0.7 && p <= 0.05) sig <- rbind(sig, data.frame(q=Q, novel=nv, annotated=gn, module=m))
      }
    }
  }
  un <- unique(sig$novel)
  cat(sprintf("q=%.2f: sig_pairs=%d unique_NCP=%d retained_of_published85=%d (%.1f%%)\n",
              Q, nrow(sig), length(un), length(intersect(un, pub85)), 100*length(intersect(un, pub85))/85))
  summary <- rbind(summary, data.frame(q=Q, sig_pairs=nrow(sig), unique_NCP=length(un),
                                       retained_of_85=length(intersect(un, pub85))))
  out <- rbind(out, sig)
}
write.csv(out, "sensitivity_TOMquantile_correctList.csv", row.names=FALSE)
write.csv(summary, "sensitivity_TOMquantile_summary.csv", row.names=FALSE)
cat("saved: sensitivity_TOMquantile_correctList.csv / _summary.csv\n")
