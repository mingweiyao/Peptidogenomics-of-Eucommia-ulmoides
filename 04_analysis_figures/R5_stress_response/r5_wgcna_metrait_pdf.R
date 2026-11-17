suppressMessages(library(pheatmap))

script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_root <- if (length(script_args)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_args[[1]])), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(".", mustWork = FALSE)
}
repo_root <- Sys.getenv("EU_REPO_ROOT", repo_root)
IN <- file.path(repo_root, "output", "04_article", "tmp", "r5wgcna", "mt")
OUT <- IN

corM <- as.matrix(read.delim(file.path(IN, "ME_trait_correlation.tsv"),
                             row.names = 1, check.names = FALSE))
padjM <- as.matrix(read.delim(file.path(IN, "ME_trait_padj.tsv"),
                              row.names = 1, check.names = FALSE))
stopifnot(identical(dim(corM), dim(padjM)))

corP <- corM; pP <- padjM
corP[is.na(corP)] <- 0; pP[is.na(pP)] <- 1
txt <- matrix(paste0(signif(corP, 2), "\n(FDR ", signif(pP, 1), ")"), nrow(corP))

w <- 3 + 1.1 * ncol(corP); h <- 2.5 + 0.55 * nrow(corP)
pheatmap(corP, display_numbers = txt, cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "ST_DT: ME-trait (r, FDR)",
         filename = file.path(OUT, "me_trait_heatmap_ST_DT.pdf"),
         width = w, height = h, fontsize_number = 7)
pheatmap(corP, display_numbers = txt, cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "ST_DT: ME-trait (r, FDR)",
         filename = file.path(OUT, "me_trait_heatmap_ST_DT.png"),
         width = w, height = h, fontsize_number = 7, res = 300)
cat("pheatmap PDF/PNG written:", w, "x", h, "in\n")
