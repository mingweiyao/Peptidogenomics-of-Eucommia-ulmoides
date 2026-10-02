# =====================================================================
# Reproduce the drought/salt-stress WGCNA network (parameters identical,
# item by item, to the original manuscript script new_code/WGCNA_dt_st.r).
#
# Purpose:
#   Rebuild the signed WGCNA network on the drought-salt expression matrix
#   and export per-gene module assignments.
# Inputs:
#   dt_st_count_fpkm.txt (22,286 rows = 22,009 evm + 277 EuNCP_trans; the
#   quantification universe after the row-sum >= 10 filter)
# Outputs:
#   <outDir>/gene_module_colors.txt, finalTOM-block.1.RData, net.rds
# Usage:
#   Rscript 01_run_wgcna_network.R <expr.txt> <outDir> <softPower> <minModuleSize> <mergeCut> <deepSplit>
#   Manuscript parameters (stress network):
#     Rscript 01_run_wgcna_network.R ./data/rejected_20260204/02figure/figure6/WGCNA/dt_st_count_fpkm.txt out_stress 16 200 0.4 0
# Notes:
#   Expected: 15 samples, 21,608 genes (after the MAD > 0 filter). Module
#   assignment on a re-run may legitimately differ from the original run (the
#   module membership table behind the original 9 modules and the "high-TOM NCP"
#   231/181 lists was never archived, so 231/181 cannot be re-derived digit by
#   digit; to reproduce 35/27 exactly, use the archived-list path in
#   s5c_candidates.py).
#   The 23k-scale TOM requires >= 32 GB RAM; do not run on machines with less
#   memory than that.
# =====================================================================
suppressMessages({library(WGCNA)})
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
exprMat <- args[1]; outputDir <- args[2]
softPower <- as.numeric(args[3]); minModuleSize <- as.numeric(args[4])
mergeCut <- as.numeric(args[5]); deepSplit <- as.numeric(args[6])
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=TRUE, quote="", comment.char="")
dataExpr <- as.data.frame(t(dataExpr))
geneType <- ifelse(grepl("^evm", colnames(dataExpr)), "annotated", "novel")
m.mad <- apply(dataExpr, 2, mad)
dataExpr <- dataExpr[, m.mad > 0]
gsg <- goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK) dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
nSamples <- nrow(dataExpr)
cat("samples:", nSamples, "genes:", ncol(dataExpr), "\n")

enableWGCNAThreads()
net <- blockwiseModules(
  dataExpr, power = softPower, maxBlockSize = ncol(dataExpr),
  TOMType = "signed", networkType = "signed",
  minModuleSize = minModuleSize, mergeCutHeight = mergeCut, deepSplit = deepSplit,
  reassignThreshold = 0, numericLabels = TRUE, pamRespectsDendro = TRUE,
  saveTOMs = TRUE, saveTOMFileBase = file.path(outputDir, "finalTOM"),
  verbose = 3, maxPOutliers = 1
)
moduleColors <- labels2colors(net$colors)
write.table(data.frame(Gene = colnames(dataExpr), Module = moduleColors),
            file.path(outputDir, "gene_module_colors.txt"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Number of modules:", length(table(moduleColors)), "\n")
print(table(moduleColors))
saveRDS(net, file.path(outputDir, "net.rds"))
cat("DONE\n")
