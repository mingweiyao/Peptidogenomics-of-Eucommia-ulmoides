library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(stringr)
library(RIdeogram)
library(data.table)

# ---------------- 图1：RIdeogram 染色体密度图 ----------------
karyo_file <- "G:/peptidegenomics/02figure/figure3/karyotype_use.tsv"
density_file <- "G:/peptidegenomics/02figure/figure3/material_gene_density_total.tsv"
label_file <- "G:/peptidegenomics/02figure/figure3/Pi_for_CE_and_CW.tsv"
karyotype <- read.table(karyo_file, sep = "\t", header = T, stringsAsFactors = F)
density <- read.table(density_file, sep = "\t", header = T, stringsAsFactors = F)
label <- read.table(label_file, sep = "\t", header = T, stringsAsFactors = F)
ideogram(karyotype = karyotype, overlaid = density, label = label, label_type = "line")
convertSVG("chromosome.svg", device = "png")
svg2pdf("chromosome.svg")

karyo_file <- "G:/peptidegenomics/02figure/figure5a/karyotype.tsv"
density_file <- "G:/peptidegenomics/02figure/figure5a/peptide_density_percentile.tsv"
karyotype <- read.table(karyo_file, sep = "\t", header = T, stringsAsFactors = F)
density <- read.table(density_file, sep = "\t", header = T, stringsAsFactors = F)
ideogram(karyotype = karyotype, overlaid = density)
convertSVG("chromosome.svg", device = "png")
svg2pdf("chromosome.svg")