rm(list = ls())
# =============================== 构建杜仲GO富集专属OrgDb数据库 =======================================
# 1. 安装和加载包
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("AnnotationForge", "GO.db"))
suppressPackageStartupMessages({
  library(AnnotationForge)
  library(GO.db)
})
# 2. 加载数据
gene_go <- read.delim("/data/Eu/WGCNA/GO/drought/out.emapper.annotations.GO.txt", stringsAsFactors = FALSE)
colnames(gene_go) <- c("GID", "GO")
gene_kegg <- read.delim("/data/Eu/WGCNA/GO/drought/out.emapper.annotations.KEGG_Knum.txt", stringsAsFactors = FALSE)
colnames(gene_kegg) <- c("GID", "PATH")
# 3. 获取GO本体信息
unique_go_terms <- unique(gene_go$GO)
go_ontology <- select(GO.db, keys = unique_go_terms, 
                      columns = c("GOID", "ONTOLOGY", "TERM"))
gene_go_annot <- merge(gene_go, go_ontology, by.x = "GO", by.y = "GOID")
gene_go_annot$EVIDENCE <- "IEA"
# 4. 准备基因信息
gene_info <- data.frame(GID = unique(c(gene_go_annot$GID, gene_kegg$GID)),
                        GENENAME = unique(c(gene_go_annot$GID, gene_kegg$GID)),
                        stringsAsFactors = FALSE
)
# 5. 构建数据库
tax_id <- "43911"
valid_go <- gene_go_annot[gene_go_annot$ONTOLOGY %in% c("MF", "BP", "CC"), ]
go_data <- valid_go[, c("GID", "GO", "EVIDENCE")]
makeOrgPackage(
  gene_info = gene_info,
  go = go_data,
  kegg = gene_kegg,
  version = "0.1",
  maintainer = "Mingwei Yao <yaomingwei0318@163.com>",
  author = "Mingwei Yao",
  outputDir = ".",
  tax_id = 43911,
  genus = "Eucommia",
  species = "ulmoides",
  goTable = "go"
)
package_name <- paste0("org.", substr("Eucommia", 1, 1), gsub(" ", "", "ulmoides"), ".eg.db")
install.packages(package_name, repos = NULL)
library(package_name, character.only = TRUE)

install.packages("org.Eulmoides.eg.db", repos = NULL)
library(org.Eulmoides.eg.db)

# =============================== 杜仲富集分析 =======================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "ggplot2", "enrichplot", "dplyr", "Biostrings", "topGO", "Rgraphviz"))
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
  library(enrichplot)
  library(dplyr)
  library(Biostrings)
  library(topGO)
  library(Rgraphviz)
  library(purrr)
  library(tidyr)
})
# 1. 设置路径和参数
gene_pairs_file <- "/data/Eu/WGCNA/GO/test/yield_tom_all_info.txt"
output_dir <- "/data/Eu/WGCNA/GO/test/yield_all"
protein_file <- "/data/Eu/Eu_genome/GWHBISF00000000.Protein.faa"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# 2. 构建ID映射表
build_id_mapping <- function(protein_file) {
  if(!file.exists(protein_file)) {
    stop("Protein file not found: ", protein_file)
  }
  protein_seq <- readAAStringSet(protein_file)
  headers <- names(protein_seq)
  map_dfr(headers, ~{
    parts <- unlist(strsplit(., "\t"))
    if (length(parts) >= 2 && grepl("OriGeneID=", parts[length(parts)-1])) {
      data.frame(
        old_id = sub(".*OriGeneID=", "", parts[length(parts)-1]),
        new_id = sub("^>", "", parts[1]),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })
}
id_mapping <- build_id_mapping(protein_file)
# 3. ID转换函数
convert_annotated_ids <- function(ids) {
  ids[!grepl("^NCP", ids)] %>%
    {id_mapping$new_id[match(., id_mapping$old_id)]} %>%
    na.omit() %>%
    unique()
}
# 4. GO分析函数（修改后的版本）
perform_gene_analysis <- function(novel_gene, annotated_genes, output_path) {
  if(!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  ego <- tryCatch({
    enrichGO(
      gene = annotated_genes,
      OrgDb = org.Eulmoides.eg.db,
      keyType = "GID",
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.3,
      readable = FALSE
    )
  }, error = function(e) {
    message("Error in enrichGO for gene ", novel_gene, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(ego) && nrow(ego) > 0) {
    ego_result <- as.data.frame(ego)
    write.csv(ego_result, file.path(output_path, "full_go_results.csv"), row.names = FALSE)
    
    pdf(file.path(output_path, "go_plots.pdf"), width = 12, height = 10)
    try({
      print(barplot(ego, showCategory = 15, title = paste(novel_gene, "- GO Enrichment")))
      print(dotplot(ego, showCategory = 15))
      if (nrow(ego) >= 3) {
        print(emapplot(pairwise_termsim(ego), showCategory = 15))
        print(cnetplot(ego, showCategory = 8, foldChange = NULL))
      }
    }, silent = TRUE)
    dev.off()
    
    tryCatch({
      ego_result %>%
        group_by(ONTOLOGY) %>%
        arrange(p.adjust) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        dplyr::select(ONTOLOGY, Description, p.adjust, Count)  # 修改为Count而不是geneID
    }, error = function(e) {
      message("Error processing results for ", novel_gene, ": ", e$message)
      ego_result %>%
        group_by(ONTOLOGY) %>%
        arrange(p.adjust) %>%
        dplyr::filter(row_number() == 1) %>%
        ungroup() %>%
        dplyr::select(ONTOLOGY, Description, p.adjust, Count)  # 修改为Count而不是geneID
    })
  } else {
    message("No enrichment for ", novel_gene)
    NULL
  }
}
# 5. 主分析流程（修改后的版本）
analyze_novel_genes <- function() {
  if(!file.exists(gene_pairs_file)) {
    stop("Gene pairs file not found: ", gene_pairs_file)
  }
  
  first_line <- readLines(gene_pairs_file, n = 1)
  n_columns <- length(strsplit(first_line, "\t")[[1]])
  col_classes <- c("character", "character", rep("NULL", max(0, n_columns - 2)))
  gene_pairs <- read.table(gene_pairs_file, 
                           colClasses = col_classes,
                           header = FALSE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           quote = "",
                           comment.char = "",
                           fill = TRUE)
  gene_pairs <- gene_pairs[complete.cases(gene_pairs[,1:2]), ]
  colnames(gene_pairs) <- c("novel", "annotated")
  
  novel_stats <- gene_pairs %>%
    group_by(novel) %>%
    summarise(
      annotated_count = n_distinct(annotated),
      .groups = "drop"
    ) %>%
    filter(annotated_count >= 5)
  
  results <- map_dfr(novel_stats$novel, ~{
    novel_gene <- .x
    message("Processing novel gene: ", novel_gene)
    annotated_genes <- gene_pairs %>%
      filter(novel == novel_gene) %>%
      pull(annotated) %>%
      convert_annotated_ids()
    
    if (length(annotated_genes) < 5) {
      message("Skipped ", novel_gene, ": only ", length(annotated_genes), " genes after conversion")
      return(NULL)
    }
    
    gene_output_dir <- file.path(output_dir, "gene_results", novel_gene)
    top_funcs <- perform_gene_analysis(novel_gene, annotated_genes, gene_output_dir)
    
    if (is.null(top_funcs)) return(NULL)
    
    # 修改后的结果整理部分
    top_funcs %>%
      pivot_wider(
        names_from = ONTOLOGY,
        values_from = c(Description, p.adjust, Count),  # 使用Count而不是geneID
        names_glue = "{ONTOLOGY}_{.value}"
      ) %>%
      mutate(
        novel_gene = novel_gene,
        annotated_count = length(annotated_genes),
        .before = 1
      )
  })
  
  if (nrow(results) > 0) {
    # 确保列名一致
    expected_cols <- c("novel_gene", "annotated_count",
                       "BP_Description", "BP_p.adjust", "BP_Count",
                       "MF_Description", "MF_p.adjust", "MF_Count",
                       "CC_Description", "CC_p.adjust", "CC_Count")
    
    # 添加缺失的列（如果有）
    for(col in expected_cols) {
      if(!col %in% names(results)) {
        results[[col]] <- NA
      }
    }
    
    # 按预期顺序排列列
    results <- results[, expected_cols]
    
    write.csv(
      results, 
      file.path(output_dir, "novel_genes_function_predictions.csv"), 
      row.names = FALSE
    )
    
    sink(file.path(output_dir, "analysis_summary.txt"))
    cat("=== Novel Gene Function Prediction Summary ===\n")
    cat("Analysis Date:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")
    cat("Input Statistics:\n")
    cat("- Total novel genes:", length(unique(gene_pairs$novel)), "\n")
    cat("- Analyzed novel genes:", nrow(novel_stats), "\n\n")
    cat("Output Statistics:\n")
    cat("- Genes with functional predictions:", nrow(results), "\n")
    cat("- Average annotated genes per novel gene:", 
        round(mean(results$annotated_count), 1), "\n")
    sink()
    
    message("\nAnalysis successfully completed!")
    message("Main results saved to: ", file.path(output_dir, "novel_genes_function_predictions.csv"))
    message("Detailed results for each gene in: ", file.path(output_dir, "gene_results"))
  } else {
    warning("No valid results generated! Please check input data.")
  }
}

# 执行分析
analyze_novel_genes()