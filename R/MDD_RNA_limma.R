# Calculate RNAseq log2 fold changes and differentially expressed proteins.
#
# The purpose of this script is to calculate log2 fold changes for all conditions
# compared to ctrl_0, using limma.
# This script is intended as the starting point for other scripts using RPPA
# L2FCs and or DE proteins

library(tidyverse)
library(cowplot)
library(limma)
library(ComplexHeatmap)


RNAseqlevel3File <- "../RNAseq/Data/MDD_RNAseq_Level3.csv"
MDDannoFile <- "../Metadata/MDD_sample_annotations.csv"

logFC_threshold <- 0.5
pval_threshold  <- 0.01

outDirPlots <- "../plots/RNAseq_limma"
outDirData <- "../RNAseq/Data/DEResults"

if (!dir.exists(outDirPlots)) {
  dir.create(outDirPlots)
}

if (!dir.exists(outDirData)) {
  dir.create(outDirData)
}

###############################################################################

RNAseq.mat <- read.table(RNAseqlevel3File,
                       header = TRUE,
                       sep = ",",
                       row.names = 1)
RNAseq.meta <- read.table(MDDannoFile,
                        header = TRUE,
                        sep = ",") %>% 
  filter(str_detect(experimentalTimePoint, "0|24|48"),
         specimenID %in% colnames(RNAseq.mat)) %>% 
  dplyr::select(specimenID, specimenName, ligand, experimentalTimePoint, experimentalCondition, replicate) %>% 
  mutate(experimentalCondition = fct_inorder(factor(experimentalCondition)))

RNAseq.mat <- RNAseq.mat %>%
  dplyr::select(all_of(unique(RNAseq.meta$specimenID)))

###############################################################################
# Creating design for comparisons.
# All conditions to be compared to ctrl_0.
#limma anlaysis expects log-expression values
design <- model.matrix(~experimentalCondition, RNAseq.meta)

lm <- lmFit(RNAseq.mat, design)
lm <- eBayes(lm)

tt <- topTable(lm, number = Inf)
write.csv(tt, sprintf("%s/RNAseq_DE_topTables.csv", outDirData))

res <- decideTests(lm, p.value = pval_threshold, lfc = logFC_threshold)

resTP <- matrix(res, nrow = nrow(res))[,-1]
rownames(resTP) <- rownames(res)
colnames(resTP) <- colnames(design)[-1] %>%
  str_remove("experimentalCondition")


pdf(sprintf("%s/RNAseq_significantAnalytes.pdf", outDirPlots), height = 12, width = 16)
hm <- Heatmap(resTP,
              name = "significant",
              cluster_columns = FALSE,
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 9),
              show_row_names = TRUE,
              column_title = "RNAseq: significant histone epigentic marks compared to ctrl_0")
draw(hm)
dev.off()
