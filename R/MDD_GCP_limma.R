# Calculate GCP log2 fold changes and differentially expressed proteins.
#
# The purpose of this script is to calculate log2 fold changes for all conditions
# compared to ctrl_0, using limma.
# This script is intended as the starting point for other scripts using RPPA
# L2FCs and or DE proteins

library(tidyverse)
library(cowplot)
library(limma)
library(biomaRt)
library(ComplexHeatmap)
# library(plotly)
# library(DESeq2)
# library(UpSetR)

GCPlevel3File <- "../GCP/Data/MDD_GCP_Level3.csv"
MDDannoFile <- "../Metadata/MDD_sample_annotations.csv"
#antibodyFile   <- "../GCP/Metadata/MDD_RPPA_antibodyAnnotations.csv"

logFC_threshold <- 0.5
pval_threshold  <- 0.01

outDirPlots <- "../plots/GCP_limma"
outDirData <- "../GCP/Data/DEProteins"

###############################################################################

GCP.mat <- read.table(GCPlevel3File,
                       header = TRUE,
                       sep = ",",
                       row.names = 1)
GCP.meta <- read.table(MDDannoFile,
                        header = TRUE,
                        sep = ",") %>% 
  filter(specimenID %in% colnames(GCP.mat)) %>% 
  dplyr::select(specimenID, specimenName, ligand, experimentalTimePoint, experimentalCondition, replicate) %>% 
  mutate(experimentalCondition = fct_inorder(factor(experimentalCondition)))


###############################################################################
# Creating design for comparisons.
# All conditions to be compared to ctrl_0.
design <- model.matrix(~experimentalCondition, GCP.meta)

lm <- lmFit(GCP.mat, design)
lm <- eBayes(lm)

res <- decideTests(lm, p.value = pval_threshold, lfc = logFC_threshold)

resTP <- matrix(res, nrow = nrow(res))
rownames(resTP) <- rownames(res)
# 
if (!dir.exists(outDirPlots)) {
  dir.create(outDirPlots)
}
# 
pdf(sprintf("%s/GCP_significantAnalytes.pdf", outDirPlots), height = 12, width = 16)
Heatmap(resTP, name = "significant", cluster_columns = FALSE,
        column_names_gp = gpar(fontsize = 9), show_row_names = TRUE,
        column_title = "GCP: significant analytes compared to ctrl_0")
dev.off()
