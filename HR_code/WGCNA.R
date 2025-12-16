# ==========================================================
#                  WGCNA COMPLETE SCRIPT
# ==========================================================

library(WGCNA)
library(dplyr)
library(stringr)
options(stringsAsFactors = FALSE)
disableWGCNAThreads()

exprMat <- "/data/Eu/WGCNA/dt.csv"
traitFile <- "/data/Eu/WGCNA/dt_trait.txt"
outputDir <- "/data/Eu/WGCNA/dt_results"
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
expr_raw <- read.csv(exprMat, row.names = 1, check.names = FALSE)
expr_log <- log2(expr_raw + 0.1)
expr_z <- t(scale(t(expr_log))) 
dataExpr <- as.data.frame(t(expr_z))

# 过滤表达恒定的基因
mad_val <- apply(dataExpr, 2, mad)
N <- min(10000, length(mad_val))
selected_genes <- names(sort(mad_val, decreasing = TRUE))[1:N]
dataExpr <- dataExpr[, selected_genes]
message("过滤完成：保留基因数 = ", ncol(dataExpr))

allGenes <- colnames(dataExpr)
geneType <- ifelse(grepl("^evm", allGenes), "annotated", "novel")
names(geneType) <- allGenes
novelGenes     <- names(geneType)[geneType == "novel"]
annotatedGenes <- names(geneType)[geneType == "annotated"]

# 检查样本/基因质量
gsg <- goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK) {
  dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nSamples <- nrow(dataExpr)
geneNames <- colnames(dataExpr)

# ==== Step 2: 样本聚类图 ====
pdf(file.path(outputDir, "sample_clustering.pdf"),  width = 15, height = 8)
sampleTree <- hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

# ==== Step 3: 选择软阈值 ====
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

# ==== Step 4: 构建网络并识别模块 ====
fit <- sft$fitIndices
targetR2 <- 0.80
softPower <- fit$Power[which.max(fit$SFT.R.sq >= targetR2)]
softPower
enableWGCNAThreads()
net <- blockwiseModules(
  dataExpr,
  power = softPower,
  maxBlockSize = ncol(dataExpr),
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 150,
  mergeCutHeight = 0.4,
  deepSplit = 1,
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
  dataExpr,
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
MEs0 <- moduleEigengenes(dataExpr, colors = moduleColors)$eigengenes
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
traitData <- traitData[match(rownames(dataExpr), rownames(traitData)), , drop = FALSE]
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
sampleNames <- rownames(dataExpr)
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
MM_all <- as.data.frame(cor(dataExpr, MEs, use = "p"))
colnames(MM_all) <- colnames(MEs)
rownames(MM_all) <- geneNames
MM.threshold <- 0.8
GS.threshold <- 0.3
hubDirMG <- file.path(outputDir, "hub_genes_MM_GS")
dir.create(hubDirMG, showWarnings = FALSE)
for (traitName in colnames(traitData)) {  
  GS_all <- as.numeric(cor(dataExpr, traitData[[traitName]], use="p"))
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

# ==== Step 11: 分模块筛选 novel–annotated 高共表达基因对 ====
load(net$TOMFiles[1], verbose = TRUE)
blockGenes <- colnames(dataExpr)[net$blockGenes[[1]]]
TOM <- as.matrix(TOM)
rownames(TOM) <- colnames(TOM) <- blockGenes
moduleOutputDir <- file.path(outputDir, "module_specific_results")
dir.create(moduleOutputDir, showWarnings = FALSE, recursive = TRUE)
moduleLabels <- moduleColors
uniqueModules <- sort(unique(moduleLabels))
uniqueModules <- uniqueModules[uniqueModules != "grey"]
for (modColor in uniqueModules) {
  modGenes <- names(moduleLabels)[moduleLabels == modColor]
  modGenes <- intersect(modGenes, colnames(TOM))
  mod_novel     <- intersect(novelGenes,     modGenes)
  mod_annotated <- intersect(annotatedGenes, modGenes)
  if (length(mod_novel) == 0 || length(mod_annotated) == 0) next
  modTOM <- TOM[modGenes, modGenes]
  if (length(modGenes) > 100) {
    global_threshold <- quantile(modTOM, 0.8, na.rm = TRUE)
  } else {
    global_threshold <- quantile(modTOM, 0.8, na.rm = TRUE)
  }
  median_tom <- median(modTOM, na.rm = TRUE)
  mad_tom    <- mad(modTOM, na.rm = TRUE)
  robust_threshold <- median_tom + 2 * mad_tom
  min_threshold <- 0
  pairs <- expand.grid(
    novel     = mod_novel,
    annotated = mod_annotated,
    stringsAsFactors = FALSE
  )
  pairs$weight <- mapply(function(nv, an) TOM[nv, an],
                         pairs$novel, pairs$annotated)
  high_conf_pairs <- pairs %>%
    dplyr::filter(
      (weight >= global_threshold | weight >= robust_threshold) &
      weight >= min_threshold
    ) %>%
    dplyr::mutate(module = modColor)
  if (nrow(high_conf_pairs) > 0) {
    mod_final_pairs <- high_conf_pairs %>%
      dplyr::group_by(novel) %>%
      dplyr::filter(n() >= ifelse(length(mod_annotated) > 50, 5, 5)) %>%
      dplyr::ungroup()    
    if (nrow(mod_final_pairs) > 0) {
      outFile <- file.path(moduleOutputDir,
                           paste0(modColor, "_high_conf_pairs.txt"))
      write.table(
        mod_final_pairs,
        file      = outFile,
        sep       = "\t",
        row.names = FALSE,
        quote     = FALSE
      )
    }
  }
}
threshold_summary <- data.frame(
  Module = uniqueModules,
  GlobalThreshold = sapply(uniqueModules, function(modColor) {
    modGenes <- names(moduleLabels)[moduleLabels == modColor]
    modGenes <- intersect(modGenes, colnames(TOM))
    modTOM   <- TOM[modGenes, modGenes]
    if (length(modTOM) > 100) quantile(modTOM, 0.85, na.rm = TRUE)
    else quantile(modTOM, 0.8, na.rm = TRUE)
  }),
  RobustThreshold = sapply(uniqueModules, function(modColor) {
    modGenes <- names(moduleLabels)[moduleLabels == modColor]
    modGenes <- intersect(modGenes, colnames(TOM))
    modTOM   <- TOM[modGenes, modGenes]
    median(modTOM, na.rm = TRUE) + 2 * mad(modTOM, na.rm = TRUE)
  })
)
write.csv(
  threshold_summary,
  file.path(outputDir, "module_threshold_summary.csv"),
  row.names = FALSE
)

# ==== Step 12: 导出所有涉及 novel 基因的边 (Cytoscape edge list) ====
novel_in_TOM <- intersect(novelGenes, colnames(TOM))
if (length(novel_in_TOM) > 0) {
  TOM_novel <- TOM[novel_in_TOM, , drop = FALSE]
  tomThreshold_novel <- quantile(TOM_novel, 0.9, na.rm = TRUE)
  sel <- which(TOM_novel >= tomThreshold_novel, arr.ind = TRUE)
  sources <- rownames(TOM_novel)[sel[, 1]]
  targets <- colnames(TOM_novel)[sel[, 2]]
  weights <- TOM_novel[sel]
  edges_novel <- data.frame(
    source = sources,
    target = targets,
    weight = weights,
    stringsAsFactors = FALSE
  )
  edges_novel_unique <- edges_novel
  idx_both_novel <- edges_novel$source %in% novel_in_TOM &
                    edges_novel$target %in% novel_in_TOM
  if (any(idx_both_novel)) {
    tmp <- edges_novel[idx_both_novel, ]
    st_sorted <- t(apply(tmp[, c("source", "target")], 1, sort))
    tmp$source <- st_sorted[,1]
    tmp$target <- st_sorted[,2]
    tmp <- tmp[!duplicated(tmp[, c("source", "target")]), ]
    edges_novel_unique <- rbind(
      edges_novel[!idx_both_novel, ],
      tmp
    )
  }
  novel_edge_file <- file.path(outputDir, "Cytoscape_edges_novel_all.txt")
  write.table(
    edges_novel_unique,
    file      = novel_edge_file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )  
  message("✅ 已输出所有涉及 novel 基因的边文件：", novel_edge_file)
} else {
  message("⚠ TOM 中没有匹配到 novel 基因（novel_in_TOM 为空），未生成 novel edge list。")
}

# ==== Step 13: 可视化基因网络====
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA
pdf(file.path(outputDir, "TOM_plot.pdf"))
TOMplot(plotTOM, net$dendrograms, moduleColors,
        main = "Network heatmap plot, all genes")
dev.off()

message("✅ WGCNA 完成！结果目录： ", outputDir)