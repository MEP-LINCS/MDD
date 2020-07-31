# Calculate RPPA log2 fold changes and differentially expressed proteins.
#
# The purpose of this script is to calculate log2 fold changes for all conditions
# compared to ctrl_0, using limma.
# This script is intended as the starting point for other scripts using RPPA
# L2FCs and or DE proteins

library(tidyverse)
library(cowplot)
library(limma)
library(biomaRt)
library(plotly)
library(DESeq2)
library(UpSetR)

RPPAlevel3File <- "../RPPA/Data/MDD_RPPA_Level3.csv"
RPPAannoFile <- "../Metadata/MDD_sample_annotations.csv"
antibodyFile   <- "../RPPA/Metadata/MDD_RPPA_antibodyAnnotations.csv"

logFC_threshold <- 0.25
pval_threshold  <- 0.01

outDirPlots <- "../plots/RPPA_limma"
outDirData <- "../RPPA/Data/DEProteins"

###############################################################################

RPPA.mat <- read.table(RPPAlevel3File,
                       header = TRUE,
                       sep = ",",
                       row.names = 1)
RPPA.meta <- read.table(RPPAannoFile,
                        header = TRUE,
                        sep = ",") %>% 
  filter(specimenID %in% colnames(RPPA.mat)) %>% 
  dplyr::select(specimenID, specimenName, ligand, experimentalTimePoint, experimentalCondition, replicate) %>% 
  mutate(experimentalCondition = fct_inorder(as.factor(experimentalCondition)))

aTrppa <- 
  read.csv(antibodyFile, header = TRUE, stringsAsFactors = FALSE) %>% 
  mutate(antibody = str_remove(MDD, "-[:alnum:]-[:alnum:]$")) %>% 
  dplyr::rename(hgnc_symbol = Symbols) %>% 
  dplyr::select(antibody, hgnc_symbol, Sites)

###############################################################################
# Creating design for comparisons.
# All conditions to be compared to ctrl_0.
design <- model.matrix(~experimentalCondition, RPPA.meta)

RPPA.mat <- RPPA.mat[aTrppa %>% pull(antibody), ]

lm <- lmFit(RPPA.mat, design)
# lm <- eBayes(lm)
lm <- treat(lm, lfc = logFC_threshold)

res <- decideTests(lm, p.value = pval_threshold)

# resTP <- data.matrix(res)
# 
# if (!dir.exists(outDirPlots)) {
#   dir.create(outDirPlots)
# }
# 
# pdf(sprintf("%s/RPPA_significantProteins.pdf", outDirPlots), height = 12, width = 16)
# Heatmap(resTP, name = "significant", cluster_columns = FALSE, 
#         column_names_gp = gpar(fontsize = 9), show_row_names = FALSE, 
#         column_title = "RPPA: significant antibodies compared to ctrl_0")
# dev.off()

dir.create(outDirData)

write_csv(rownames_to_column(data.frame(res[, -1]), "antibody"), path = sprintf("%s/MDD_RPPA_significantAntibodies.csv", outDirData))
save(res, file = sprintf("%s/MDD_RPPA_limmaResults.Rdata", outDirData))
