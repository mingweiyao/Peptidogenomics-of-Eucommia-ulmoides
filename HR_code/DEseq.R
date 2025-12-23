# ========= counts -> DESeq2 DEG with comparison list =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly=TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl")
if (!requireNamespace("apeglm", quietly=TRUE)) BiocManager::install("apeglm")
if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("pheatmap", quietly=TRUE)) install.packages("pheatmap")
if (!requireNamespace("ggrepel", quietly=TRUE)) install.packages("ggrepel")
library(DESeq2)
library(readxl)
library(apeglm)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# ---------- 1. 输入文件 ----------
file_xlsx <- "F:/Eu_peptido/prepare Horiticulture plant journal/figure6/yield_deseq.xlsx"
out_dir <- "F:/Eu_peptido/prepare Horiticulture plant journal/figure6/gene_dt"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 2. 读 counts ----------
counts_df <- read_excel(file_xlsx, sheet = "expression")
counts_df <- as.data.frame(counts_df)

rownames(counts_df) <- counts_df[,1]
counts_df <- counts_df[,-1]

counts <- as.matrix(counts_df)
counts <- round(counts)
mode(counts) <- "integer"

# ---------- 3. 读分组表 ----------
coldata <- read_excel(file_xlsx, sheet = "group")
coldata <- as.data.frame(coldata)

# 要求列名至少包含：Sample, Group
stopifnot(all(c("Sample","Group") %in% colnames(coldata)))

rownames(coldata) <- coldata$Sample
coldata <- coldata[colnames(counts), , drop=FALSE]
coldata$Group <- factor(coldata$Group)
coldata$Group <- relevel(coldata$Group, ref = "CK")

# ---------- 4. 构建 DESeq2 对象 ----------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ Group
)

# ---------- 5. 过滤低表达 ----------
grp <- colData(dds)$Group
keep <- apply(counts(dds), 1, function(x) {
  any(tapply(x > 10, grp, sum) >= 2)
})
dds <- dds[keep, ]

# ---------- 6. 跑 DESeq ----------
dds <- DESeq(dds)

# 7.1 VST 转换（盲目 = FALSE 保留分组信息，做 PCA 更合理）
vsd <- vst(dds, blind = FALSE)

# 7.2 PCA 图 ---------------------------------------------------------------
pca_data <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  ggtitle("PCA of VST-transformed counts")

ggsave(
  filename = file.path(out_dir, "PCA_Group.pdf"),
  plot     = p_pca,
  width    = 7,
  height   = 6,
  dpi      = 300
)

# 7.3 样本间距离热图 -------------------------------------------------------
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

annotation_col <- as.data.frame(colData(vsd)[, "Group", drop = FALSE])

pheatmap(
  mat                      = sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col           = annotation_col,
  main                     = "Sample-to-sample distances",
  filename                 = file.path(out_dir, "Sample_distance_heatmap.pdf")
)

norm_counts <- counts(dds, normalized = TRUE)
norm_df <- as.data.frame(norm_counts)
norm_df$GeneID <- rownames(norm_df)
norm_df <- norm_df[, c("GeneID", setdiff(colnames(norm_df), "GeneID"))]

write.csv(
  norm_df,
  file      = file.path(out_dir, "normalized_counts.csv"),
  row.names = FALSE
)

# ---------- 9. 读比较列表 ----------
comp_df <- read_excel(file_xlsx, sheet = "compare")
comp_df$Contrast <- paste0(comp_df$Group1, "_vs_", comp_df$Group2)

stopifnot(all(c("Group1","Group2") %in% colnames(comp_df)))

# ---------- 8. 按比较表逐个输出 DEG ----------
all_groups <- levels(coldata$Group)

for (i in seq_len(nrow(comp_df))) {
  g1 <- as.character(comp_df$Group1[i])
  g2 <- as.character(comp_df$Group2[i])
  cname <- as.character(comp_df$Contrast[i])
  message("正在比较：", cname, "  [", g1, " vs ", g2, "]")
  coef_name <- paste0("Group_", g1, "_vs_CK")
  
  # LFC shrink（推荐）
  res_shrink <- lfcShrink(dds, coef = coef_name, type="apeglm")
  
  res_df <- as.data.frame(res_shrink)
  res_df$GeneID <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df <- res_df[order(res_df$padj), ]
  
  # 输出全量结果
  out_all <- file.path(out_dir, paste0(cname, "_all.csv"))
  write.csv(res_df, out_all, row.names = FALSE)
  
  # 输出显著 DEG（阈值可自己改）
  deg_df <- subset(res_df, padj < 0.05 & abs(log2FoldChange) >= 1)
  out_sig <- file.path(out_dir, paste0(cname, "_sig.csv"))
  write.csv(deg_df, out_sig, row.names = FALSE)
  
  # MA 图（使用 shrink 后的结果）-----------------------------------------
  ma_file <- file.path(out_dir, paste0(cname, "_MA.pdf"))
  pdf(ma_file, width = 7, height = 7)
  plotMA(res_shrink, main = paste0("MA plot (", cname, ")"))
  abline(h = c(-1, 1), col = "gray", lty = 2)
  dev.off()
  
  ## ---------------- 火山图（Volcano plot） ----------------
  # 为了避免 padj=0 导致 log10 出 Inf，这里加一个极小值
  res_df$padj[res_df$padj == 0] <- .Machine$double.xmin
  
  # 定义显著性和上下调分组（和你筛 DEG 的阈值保持一致：padj<0.05 & |log2FC|>=2）
  res_df$Regulation <- "NS"
  res_df$Regulation[res_df$padj < 0.05 & res_df$log2FoldChange >= 1]  <- "Up"
  res_df$Regulation[res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Down"
  
  # 基础火山图
  p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("Down" = "blue", "NS" = "grey70", "Up" = "red")) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    xlab("log2(Fold Change)") +
    ylab("-log10(adjusted p-value)") +
    ggtitle(paste0("Volcano plot (", cname, ")")) +
    theme_bw(base_size = 14)
  
  # 如果装了 ggrepel，可以标注显著 DEG 中 |log2FC| 最大的前 10 个
  if ("ggrepel" %in% loadedNamespaces()) {
    library(ggrepel)
    top_lbl <- res_df[res_df$Regulation != "NS", ]
    if (nrow(top_lbl) > 0) {
      top_lbl$absLFC <- abs(top_lbl$log2FoldChange)
      top_lbl <- top_lbl[order(-top_lbl$absLFC, top_lbl$padj), , drop = FALSE] 
      top_lbl <- head(top_lbl, 10)  # 取最显著的前 10 个
      p_vol <- p_vol +
        geom_text_repel(data = top_lbl,
                        aes(label = GeneID),
                        size = 3,
                        max.overlaps = 20)
    }
  }
  vol_file <- file.path(out_dir, paste0(cname, "_Volcano.pdf"))
  ggsave(vol_file, p_vol, width = 7, height = 6, dpi = 300)
  
  message("完成：", cname,
          "  all=", nrow(res_df),
          "  sig=", nrow(deg_df))
}

cat("全部比较完成！结果输出到：", out_dir, "\n")
