# Calculate motif log2 fold changes and differentially expressed proteins.
#
# The purpose of this script is to calculate log2 fold changes for all conditions
# compared to ctrl_0, using limma.
# This script is intended as the starting point for other scripts using RPPA
# L2FCs and or DE proteins

library(tidyverse)
#library(cowplot)
library(limma)
library(ComplexHeatmap)


motiflevel3File <- "../ATACseq/Data/MDD_ATACseq_MotifFamilyScores.csv"
MDDannoFile <- "../Metadata/MDD_sample_annotations.csv"

logFC_threshold <- 0.5
pval_threshold  <- 0.01

outDirPlots <- "../plots/motif_DE_T0"
outDirData <- "../ATACseq/Data/DEResults"

if (!dir.exists(outDirPlots)) {
  dir.create(outDirPlots)
}

if (!dir.exists(outDirData)) {
  dir.create(outDirData)
}

###############################################################################

motif.mat <- read.table(motiflevel3File,
                       header = TRUE,
                       sep = ",",
                       row.names = 1)
motif.meta <- read.table(MDDannoFile,
                        header = TRUE,
                        sep = ",") %>% 
  filter(str_detect(experimentalTimePoint, "0|24|48"),
         specimenID %in% colnames(motif.mat)) %>% 
  dplyr::select(specimenID, specimenName, ligand, experimentalTimePoint, experimentalCondition, replicate) %>% 
  mutate(experimentalCondition = fct_inorder(factor(experimentalCondition)))

motif.mat <- motif.mat %>%
  dplyr::select(all_of(unique(motif.meta$specimenID)))

###############################################################################
# Creating design for comparisons.
# All conditions to be compared to ctrl_0.
design <- model.matrix(~experimentalCondition, motif.meta)

lm <- lmFit(motif.mat, design)
lm <- eBayes(lm)

tt <- topTable(lm, number = Inf)
write.csv(tt, sprintf("%s/motifs_DE_topTables.csv", outDirData))

res <- decideTests(lm, p.value = pval_threshold, lfc = logFC_threshold)

resTP <- matrix(res, nrow = nrow(res))[,-1]
rownames(resTP) <- rownames(res)
colnames(resTP) <- colnames(design)[-1] %>%
  str_remove("experimentalCondition")


pdf(sprintf("%s/motif_significantAnalytes.pdf", outDirPlots), height = 12, width = 16)
hm <- Heatmap(resTP,
              name = "significant",
              cluster_columns = FALSE,
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 9),
              show_row_names = TRUE,
              column_title = "motif: significant motifs compared to ctrl_0")
draw(hm)
dev.off()
