# WGCNA 1.73).

# hclust(average) → cutreeDynamic(deepSplit=2, pamRespectsDendro=FALSE,
# minClusterSize=20) → mergeCloseModules(cutHeight=0.25).

suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_root <- if (length(script_args)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_args[[1]])), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(".", mustWork = FALSE)
}
repo_root <- Sys.getenv("EU_REPO_ROOT", repo_root)
BASE <- file.path(repo_root, "output", "04_article", "tmp", "r5wgcna")
ART <- file.path(BASE, "out")
VSTDIR <- BASE
BETA <- 12
MINCL <- 20

d <- read.delim(file.path(BASE, "ST_DT_vst.tsv"), row.names = 1, check.names = FALSE)
pass <- readLines(file.path(BASE, "ST_DT_tpm_pass.txt"))
tg <- intersect(rownames(d), pass[nzchar(pass)])
d <- d[tg, ]
design <- read.delim(file.path(BASE, "sample_design.tsv"), check.names = FALSE)
pd <- design[design$panel == "ST_DT", ]
datExpr <- as.matrix(t(d[, pd$sample]))
cat("datExpr:", nrow(datExpr), "samples x", ncol(datExpr), "transcripts\n")

adj <- adjacency(datExpr, power = BETA, type = "signed")
TOM <- TOMsimilarity(adj, TOMType = "signed")
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
lab <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,
                     pamRespectsDendro = FALSE, minClusterSize = MINCL)
m <- mergeCloseModules(datExpr, labels2colors(lab), cutHeight = 0.25)
moduleColors <- m$colors
names(moduleColors) <- colnames(datExpr)


arch <- read.delim(file.path(BASE, "gene_module_assignment.tsv"),
                   check.names = FALSE)
stopifnot(setequal(names(moduleColors), arch$transcript_id))
same <- unname(moduleColors[arch$transcript_id]) == arch$module
cat("module assignment identical to archive:", sum(same), "/", nrow(arch),
    sprintf("(%.2f%%)", 100 * mean(same)), "\n")


saveRDS(list(geneTree = geneTree, moduleColors = moduleColors,
             dissTOM_dim = dim(dissTOM), identical_to_archive = all(same)),
        file.path(ART, "data", "fig6C1_geneTree_moduleColors.rds"))
write.table(data.frame(transcript_id = names(moduleColors),
                       leaf_order = geneTree$order,
                       module = unname(moduleColors)),
            file.path(ART, "data", "fig6C1_leaf_order_modules.tsv"),
            sep = "\t", row.names = FALSE)


pdf(file.path(ART, "figures", "fig6C1_dendro_R.pdf"),
    width = 2.96, height = 2.66, onefile = TRUE)
par(mar = c(0.4, 2.4, 0.6, 0.3), cex = 0.4)
plotDendroAndColors(geneTree, moduleColors, "Modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "", sub = "", cex.colorLabels = 0.55,
                    cex.dendroLabels = 0.3, marAll = c(0.4, 2.4, 0.6, 0.3))
dev.off()
cat("saved fig6C1_dendro_R.pdf\n")
cat("DONE R\n")
