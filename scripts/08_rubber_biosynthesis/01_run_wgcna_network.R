# =====================================================================
# Integrated WGCNA analysis script (rubber network rerun)
#
# Purpose:
#   Rebuild the published rubber WGCNA network (softPower 22, signed network,
#   minModuleSize 150, mergeCutHeight 0.3, deepSplit 2), produce the standard QC
#   plots, and derive hub genes plus novel-annotated gene pairs per module.
# Inputs:
#   <REJ>/02figure/figure5/rubber/WGCNA/rubber_count_5_fpkm.txt
#   <REJ>/02figure/figure5/rubber/WGCNA/rubber_count_5_fpkm_trait.txt
#   <REJ>/02figure/figure5/rubber/WGCNA/rubber_gca.txt
# Outputs:
#   <outRoot>/S6_rubber_candidates/wgcna_rerun/ (plots, module info, hub-gene
#   lists, module-specific novel-annotated pairs, Pearson correlations)
# Usage:
#   Rscript 01_run_wgcna_network.R
#   Env: REPRO_REJ = rejected-archive root; REPRO_OUT = output root
# Notes:
#   Steps 12/13 screen novel-annotated gene pairs per module; Step 14 computes
#   Pearson r/p for the pairs and plots the network heatmap.
# =====================================================================
library(WGCNA)
library(reshape2)
library(dplyr)
library(stringr)
options(stringsAsFactors = FALSE)

# ==== Step 1: Data preparation ====
REJ <- Sys.getenv("REPRO_REJ", unset = "./data/rejected_20260204")
exprMat <- file.path(REJ, "02figure/figure5/rubber/WGCNA/rubber_count_5_fpkm.txt")
outRoot <- Sys.getenv("REPRO_OUT", unset = "./output")
outputDir <- file.path(outRoot, "S6_rubber_candidates", "wgcna_rerun")
dir.create(outputDir, showWarnings = FALSE)
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=TRUE, quote="", comment.char="")
dataExpr <- as.data.frame(t(dataExpr))
allGenes <- colnames(dataExpr)
geneType <- ifelse(grepl("^evm", allGenes), "annotated", "novel")
names(geneType) <- allGenes
novelGenes <- names(geneType)[geneType == "novel"]
annotatedGenes <- names(geneType)[geneType == "annotated"]

# Filter out genes with invariant expression
m.mad <- apply(dataExpr, 2, mad)
dataExpr <- dataExpr[, m.mad > 0]
novelGenes <- intersect(names(geneType)[geneType == "novel"], colnames(dataExpr))
annotatedGenes <- intersect(names(geneType)[geneType == "annotated"], colnames(dataExpr))

# Check sample/gene quality
gsg <- goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK) {
  dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nSamples = nrow(dataExpr)

# ==== Step 2: Sample clustering plot ====
pdf(file.path(outputDir, "sample_clustering.pdf"),  width = 15, height = 8)
sampleTree <- hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

# ==== Step 3: Soft-threshold selection ====
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
pdf(file.path(outputDir, "soft_thresholding.pdf"))
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
     main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.8, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

# ==== Step 4: Build network and identify modules ====
# softPower <- sft$powerEstimate
softPower <- 22
enableWGCNAThreads()
net <- blockwiseModules(
  dataExpr,
  power = softPower,
  maxBlockSize = ncol(dataExpr),
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 150,
  mergeCutHeight = 0.3,
  deepSplit = 2,
  reassignThreshold = 0,
  numericLabels = TRUE,
  pamRespectsDendro = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(outputDir, "finalTOM"),
  verbose = 3
)
moduleColors <- labels2colors(net$colors)

# ==== Step 5: Module dendrogram ====
pdf(file.path(outputDir, "module_dendrogram.pdf"))
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
table(moduleColors)

# ==== Step 6: Module correlations ====
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf(file.path(outputDir, "module_correlation.pdf"))
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
dev.off()
trait <- file.path(REJ, "02figure/figure5/rubber/WGCNA/rubber_count_5_fpkm_trait.txt")
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}
MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
pdf(file.path(outputDir, "module_correlation_tissues.pdf"))
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
dev.off()

# ==== Step 7: Attach phenotype data ====
modTraitCor = cor(MEs_col, traitData, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf(file.path(outputDir, "module_trait_relationships.pdf"), width = 10, height = 16)
par(mar = c(6, 12, 3, 2))
labeledHeatmap(
  Matrix = modTraitCor,
  xLabels = colnames(traitData),
  yLabels = names(MEs_col),
  ySymbols = names(MEs_col),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.lab = 1.0,
  cex.text = 0.9,
  zlim = c(-1,1),
  main = "Module-trait relationships"
)
dev.off()

# ==== Step 8: hub gene====
hubDir <- file.path(outputDir, "hub_genes")
if (!dir.exists(hubDir)) dir.create(hubDir)
moduleList <- unique(moduleColors)
moduleList <- moduleList[moduleList != "grey"]
for (traitName in colnames(traitData)) {
  traitVec <- traitData[[traitName]]
  for (module in moduleList) {
    moduleGenes <- moduleColors == module
    geneNames <- colnames(dataExpr)
    hubGenes <- geneNames[moduleGenes]
    if (length(hubGenes) > 0) {
      outFile <- file.path(hubDir, paste0("genes_", module, "_trait_", traitName, ".txt"))
      write.table(hubGenes, file = outFile, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }
}

# ==== Step 9: hub gene====
MM.threshold <- 0.8
hubDir <- file.path(outputDir, "hub_genes_by_MM")
if (!dir.exists(hubDir)) dir.create(hubDir, recursive = TRUE)
moduleList <- unique(moduleColors)
moduleList <- moduleList[moduleList != "grey"]
for (module in moduleList) {
  moduleGenes <- moduleColors == module
  MEname <- paste0("ME", module)
  MM <- as.numeric(cor(dataExpr, MEs_col[, MEname], use = "p"))
  names(MM) <- colnames(dataExpr)
  isHub <- abs(MM) >= MM.threshold
  hubGenes <- names(MM)[isHub]
  if (length(hubGenes) > 0) {
    hubGenes <- hubGenes[order(-abs(MM[hubGenes]))]
    resultDF <- data.frame(
      Gene = hubGenes,
      Module = module,
      MM = MM[hubGenes],
      stringsAsFactors = FALSE
    )
    outFile <- file.path(hubDir, paste0("hubGenes_", module, ".txt"))
    write.table(resultDF, file = outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# ==== Step 10: hub gene====
MM.threshold <- 0.8
GS.threshold <- 0.3
hubDir <- file.path(outputDir, "hub_genes_by_MM_GS")
if (!dir.exists(hubDir)) dir.create(hubDir)
moduleList <- unique(moduleColors)
moduleList <- moduleList[moduleList != "grey"]
for (traitName in colnames(traitData)) {
  traitVec <- traitData[[traitName]]
  for (module in moduleList) {
    moduleGenes <- moduleColors == module
    geneNames <- colnames(dataExpr)
    GS <- as.numeric(cor(dataExpr, traitVec, use = "p"))
    names(GS) <- geneNames
    MEname <- paste0("ME", module)
    MM <- as.numeric(cor(dataExpr, MEs_col[, MEname], use = "p"))
    names(MM) <- geneNames
    isHub <- moduleGenes & abs(MM) > MM.threshold & abs(GS) > GS.threshold
    hubGenes <- geneNames[isHub]
    if (length(hubGenes) > 0) {
      outFile <- file.path(hubDir, paste0("hubGenes_", module, "_trait_", traitName, ".txt"))
      write.table(hubGenes, file = outFile, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }
}

# ==== Step 11: Export module information ====
gene_module <- data.frame(Gene = names(net$colors),
                          Module = moduleColors[names(net$colors)],
                          Type = geneType[names(net$colors)])
write.table(gene_module, file.path(outputDir, "gene_module_info.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ==== Step 12: Screen novel-annotated gene pairs per module ====
load(net$TOMFiles[1], verbose = TRUE)
TOM <- as.matrix(TOM)
colnames(TOM) <- rownames(TOM) <- colnames(dataExpr)
moduleOutputDir <- file.path(outputDir, "module_specific_results")
dir.create(moduleOutputDir, showWarnings = FALSE, recursive = TRUE)
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
for (mod in unique(moduleLabels)) {
  if (mod == 0) next
  modColor <- labels2colors(mod)
  modGenes <- colnames(TOM)[moduleLabels == mod]
  mod_novel <- intersect(novelGenes, modGenes)
  mod_annotated <- intersect(annotatedGenes, modGenes)
  if (length(mod_novel) == 0 || length(mod_annotated) == 0) next
  modTOM <- TOM[modGenes, modGenes]
  if (length(modGenes) > 100) {
    global_threshold <- quantile(modTOM, 0.8, na.rm = TRUE)
  } else {
    global_threshold <- quantile(modTOM, 0.8, na.rm = TRUE)
  }
  median_tom <- median(modTOM, na.rm = TRUE)
  mad_tom <- mad(modTOM, na.rm = TRUE)
  robust_threshold <- median_tom + 2 * mad_tom
  min_threshold <- 0
  pairs <- expand.grid(
    novel = mod_novel,
    annotated = mod_annotated,
    stringsAsFactors = FALSE
  )
  pairs$weight <- apply(pairs, 1, function(x) {
    TOM[x["novel"], x["annotated"]]
  })
  high_conf_pairs <- pairs %>%
    filter(
      (weight >= global_threshold | weight >= robust_threshold) &
      weight >= min_threshold
    ) %>%
    mutate(module = modColor)
  if (nrow(high_conf_pairs) > 0) {
    mod_final_pairs <- high_conf_pairs %>%
      group_by(novel) %>%
      filter(n() >= 1) %>%
      ungroup()
    write.table(
      mod_final_pairs,
      file.path(moduleOutputDir, paste0(modColor, "_high_conf_pairs.txt")),
      sep = "\t", 
      row.names = FALSE, 
      quote = FALSE
    )
  }
}
threshold_summary <- data.frame(
  Module = labels2colors(sort(unique(moduleLabels[moduleLabels != 0]))),
  GlobalThreshold = sapply(sort(unique(moduleLabels[moduleLabels != 0])), function(m) {
    modTOM <- TOM[moduleLabels == m, moduleLabels == m]
    if (length(modTOM) > 100) quantile(modTOM, 0.85, na.rm = TRUE)
    else quantile(modTOM, 0.8, na.rm = TRUE)
  }),
  RobustThreshold = sapply(sort(unique(moduleLabels[moduleLabels != 0])), function(m) {
    modTOM <- TOM[moduleLabels == m, moduleLabels == m]
    median(modTOM, na.rm = TRUE) + 2 * mad(modTOM, na.rm = TRUE)
  })
)
write.csv(threshold_summary, file.path(outputDir, "module_threshold_summary.csv"))

# ==== Step 13: Screen novel-annotated gene pairs per module (annotated set = rubber_gca.txt) ====
known_annotated_genes_file <- file.path(REJ, "02figure/figure5/rubber/WGCNA/rubber_gca.txt")
known_annotated_genes <- readLines(known_annotated_genes_file)
load(net$TOMFiles[1], verbose = TRUE)
TOM <- as.matrix(TOM)
colnames(TOM) <- rownames(TOM) <- colnames(dataExpr)
moduleOutputDir <- file.path(outputDir, "module_rubber_gca")
dir.create(moduleOutputDir, showWarnings = FALSE, recursive = TRUE)
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
all_high_conf_pairs <- data.frame()
for (mod in unique(moduleLabels)) {
  if (mod == 0) next
  modColor <- labels2colors(mod)
  modGenes <- colnames(TOM)[moduleLabels == mod]
  mod_annotated <- intersect(known_annotated_genes, modGenes)
  mod_novel <- intersect(novelGenes, modGenes)
  if (length(mod_novel) == 0 || length(mod_annotated) == 0) next
  modTOM <- TOM[modGenes, modGenes]
  global_threshold <- quantile(modTOM, 0.8, na.rm = TRUE)
  median_tom <- median(modTOM, na.rm = TRUE)
  mad_tom <- mad(modTOM, na.rm = TRUE)
  robust_threshold <- median_tom + 2 * mad_tom
  min_threshold <- 0
  pairs <- expand.grid(
    novel = mod_novel,
    annotated = mod_annotated,
    stringsAsFactors = FALSE
  )
  pairs$weight <- apply(pairs, 1, function(x) {
    TOM[x["novel"], x["annotated"]]
  })
  high_conf_pairs <- pairs %>%
    filter(
      (weight >= global_threshold | weight >= robust_threshold) &
        weight >= min_threshold
    ) %>%
    mutate(module = modColor)
  if (nrow(high_conf_pairs) > 0) {
    mod_final_pairs <- high_conf_pairs %>%
      group_by(novel) %>%
      filter(n() >= 1) %>%
      ungroup()
    write.table(
      mod_final_pairs,
      file.path(moduleOutputDir, paste0(modColor, "_high_conf_pairs.txt")),
      sep = "\t", 
      row.names = FALSE, 
      quote = FALSE
    )
  }
}

# ==== Step 14: Correlation computation ====
pairDir <- file.path(outputDir, "module_rubber_gca")   # output directory of Step 13
outCsv  <- file.path(outputDir, "rubber_gca_pairs_with_pearson.csv")
pairFiles <- list.files(pairDir, pattern = "_high_conf_pairs\\.txt$", full.names = TRUE)
if (length(pairFiles) == 0) {
  stop(paste("No pair files found in:", pairDir))
}
allPairs <- do.call(rbind, lapply(pairFiles, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
  df$source_file <- basename(f)
  df
}))
# Defensive check: make sure the required columns exist
required_cols <- c("novel", "annotated", "weight", "module")
missing_cols <- setdiff(required_cols, colnames(allPairs))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns in pair files:", paste(missing_cols, collapse = ", ")))
}
# Deduplicate (the same pair may be duplicated when merging)
allPairs <- allPairs %>%
  distinct(novel, annotated, module, weight, .keep_all = TRUE)
# Number of samples (used for p-values)
nSamples <- nrow(dataExpr)
# Compute Pearson correlations
calc_one <- function(novel_id, annotated_id) {
  # Return NA if either gene is absent from the expression matrix
  if (!(novel_id %in% colnames(dataExpr)) || !(annotated_id %in% colnames(dataExpr))) {
    return(c(r = NA_real_, p = NA_real_))
  }
  x <- as.numeric(dataExpr[, novel_id])
  y <- as.numeric(dataExpr[, annotated_id])
  # Correlation is undefined for constant vectors
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    return(c(r = NA_real_, p = NA_real_))
  }
  r <- cor(x, y, use = "pairwise.complete.obs", method = "pearson")
  p <- corPvalueStudent(r, nSamples)
  c(r = r, p = p)
}
rp_mat <- mapply(calc_one, allPairs$novel, allPairs$annotated)
allPairs$pearson_r <- as.numeric(rp_mat["r", ])
allPairs$pearson_p <- as.numeric(rp_mat["p", ])
# Add a direction label column (optional)
allPairs$direction <- ifelse(is.na(allPairs$pearson_r), NA,
                             ifelse(allPairs$pearson_r >= 0, "positive", "negative"))
# Sort by module, annotated gene, TOM descending
allPairs <- allPairs %>%
  arrange(module, annotated, desc(weight))
write.csv(allPairs, outCsv, row.names = FALSE)
message("✅ Pearson results written to: ", outCsv)


# ==== Step 14: Visualize the gene network ====
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA
pdf(file.path(outputDir, "TOM_plot.pdf"))
TOMplot(plotTOM, net$dendrograms, moduleColors,
        main = "Network heatmap plot, all genes")
dev.off()

message("✅ Analysis complete: all outputs are saved under the 'WGCNA_Results' folder")
