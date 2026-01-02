# ==========================================================
#                  WGCNA COMPLETE SCRIPT
# ==========================================================

library(WGCNA)
library(dplyr)
library(stringr)
options(stringsAsFactors = FALSE)
disableWGCNAThreads()

exprMat <- "/data/Eu/WGCNA/deseq_dt_st.csv"
traitFile <- "/data/Eu/WGCNA/DT_ST_trait.txt"
outputDir <- "/data/Eu/WGCNA/dt_st_results"
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
expr_raw <- read.csv(exprMat, row.names = 1, check.names = FALSE)
dataExpr <- as.matrix(expr_raw)
storage.mode(dataExpr) <- "numeric"

# 过滤表达恒定的基因
mad_val <- apply(dataExpr, 1, mad, constant = 1)
selected_genes <- rownames(dataExpr)[mad_val >= 10]
dataExpr2 <- dataExpr[selected_genes, , drop = FALSE]
message("过滤完成：保留基因数 = ", nrow(dataExpr2))

# 检查样本/基因质量
datExpr <- t(dataExpr2)
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
nSamples <- nrow(datExpr)
geneNames <- colnames(datExpr)

# ==== Step 2: 样本聚类图 ====
pdf(file.path(outputDir, "sample_clustering.pdf"),  width = 15, height = 8)
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

# ==== Step 3: 选择软阈值 ====
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
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

# ==== Step 4: 构建网络并识别模块 ====
fit <- sft$fitIndices
targetR2 <- 0.8
softPower <- fit$Power[which.max(fit$SFT.R.sq >= targetR2)]
softPower
enableWGCNAThreads()
net <- blockwiseModules(
  datExpr,
  power = softPower,
  maxBlockSize = ncol(datExpr),
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 100,
  mergeCutHeight = 0.4,
  deepSplit = 2,
  reassignThreshold = 0,
  numericLabels = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(outputDir, "finalTOM"),
  pamRespectsDendro = TRUE,
  verbose = 3
)
moduleColors <- net$colors

# ==== Step 5: 模块树状图 ====
pdf(file.path(outputDir, "module_dendrogram.pdf"))
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
table(moduleColors)

# ==== Step 5.1: 计算模块内部连通度 kWithin ==== 
kWithinRes <- intramodularConnectivity.fromExpr(
  datExpr,
  moduleColors,
  power = softPower
)
kWithinRes$Gene   <- rownames(kWithinRes)
kWithinRes$Module <- moduleColors[kWithinRes$Gene]
write.table(
  kWithinRes,
  file = file.path(outputDir, "kWithin_all_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==== Step 6: 模块相关性 ====
MEs0 <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)
pdf(file.path(outputDir, "module_correlation.pdf"))
plotEigengeneNetworks(
  MEs,
  "Eigengene adjacency heatmap",
  marDendro = c(3,3,2,4),
  marHeatmap = c(3,4,2,2),
  plotDendrograms = TRUE
)
dev.off()

# ==== Step 7: 模块–性状关联 ====
traitData <- read.table(traitFile, sep = "\t", header = TRUE, row.names = 1)
traitData <- traitData[match(rownames(datExpr), rownames(traitData)), , drop = FALSE]
traitDummy <- model.matrix(~ Group, data = traitData)[, -1]
modTraitCor <- cor(MEs, traitDummy, use = "p")
modTraitP   <- corPvalueStudent(modTraitCor, nSamples)
textMatrix <- paste(signif(modTraitCor, 2),
                    "\n(", signif(modTraitP, 1), ")",
                    sep = "")
dim(textMatrix) <- dim(modTraitCor)
pdf(file.path(outputDir, "module_trait_relationships.pdf"), width = 10, height = 8)
labeledHeatmap(
  Matrix = modTraitCor,
  xLabels = colnames(traitDummy),
  yLabels = colnames(MEs),
  ySymbols = colnames(MEs),
  textMatrix = textMatrix,
  colors = blueWhiteRed(50),
  zlim = c(-1, 1),
  main = "Module–trait relationships"
)
dev.off()

# ==== Step 7.1: 模块 Eigengene 表达趋势图 ====
pdf(file.path(outputDir, "module_Eigengene_profiles.pdf"), width = 10, height = 8)
par(mfrow = c(ceiling(length(colnames(MEs))/3), 3), mar = c(5, 5, 3, 2))
sampleNames <- rownames(datExpr)
for (MEname in colnames(MEs)) {
  plot(
    x    = 1:nSamples,
    y    = MEs[, MEname],
    type = "b",
    pch  = 16,
    xlab = "Sample",
    ylab = "Eigengene expression",
    main = MEname,
    xaxt = "n"
  )
  axis(1, at = 1:nSamples, labels = sampleNames, las = 2, cex.axis = 0.6)
}
dev.off()

# ==== Step 8: Module gene ====
moduleDir_all <- file.path(outputDir, "module_all_genes")
dir.create(moduleDir_all, showWarnings = FALSE)
moduleList <- unique(moduleColors)
moduleList <- moduleList[moduleList != "grey"]
for (module in moduleList) {
  genes_in_module <- geneNames[moduleColors == module]  
  outFile <- file.path(moduleDir_all, paste0("allGenes_", module, ".txt"))
  write.table(genes_in_module, outFile,
              sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ==== Step 9: Hub gene by MM + GS (module significance) ====
MM_all <- as.data.frame(cor(datExpr, MEs, use = "p"))
colnames(MM_all) <- colnames(MEs)
rownames(MM_all) <- geneNames
MM.threshold <- 0.8
GS.threshold <- 0.2
hubDirMG <- file.path(outputDir, "hub_genes_MM_GS")
dir.create(hubDirMG, showWarnings = FALSE)
for (traitName in colnames(traitData)) {  
  GS_all <- as.numeric(cor(datExpr, traitData[[traitName]], use="p"))
  names(GS_all) <- geneNames  
  for (module in moduleList) {    
    MEname <- paste0("ME", module)
    genes_in_module <- geneNames[moduleColors == module]    
    MM_vec <- MM_all[genes_in_module, MEname]
    GS_vec <- GS_all[genes_in_module]    
    isHub <- abs(MM_vec) >= MM.threshold & abs(GS_vec) >= GS.threshold
    hubGenes <- genes_in_module[isHub]    
    if (length(hubGenes) > 0) {
      resultDF <- data.frame(
        Gene = hubGenes,
        Module = module,
        MM = MM_vec[hubGenes],
        GS = GS_vec[hubGenes]
      )
      outFile <- file.path(
        hubDirMG,
        paste0("hubGenes_", module, "_trait_", traitName, ".txt")
      )
      write.table(resultDF, outFile, sep="\t",
                  quote=FALSE, row.names=FALSE)
    }
  }
}

# ==== Step 10: 输出模块信息 ====
gene_module <- data.frame(Gene = geneNames, Module = moduleColors)
write.table(gene_module, file.path(outputDir, "gene_module_info.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
