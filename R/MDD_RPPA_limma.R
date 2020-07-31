# Calculate RPPA log2 fold changes and differentially expressed proteins.
#
# The purpose of this script is to calculate log2 fold changes for all conditions
# compared to ctrl_0, using limma.
# This script is intended as the starting point for other scripts using RPPA
# L2FCs and or DE proteins

library(tidyverse)
#library(cowplot)
library(limma)
library(ComplexHeatmap)


RPPAlevel3File <- "../RPPA/Data/MDD_RPPA_Level3.csv"
MDDannoFile <- "../Metadata/MDD_sample_annotations.csv"

logFC_threshold <- 0.25
pval_threshold  <- 0.01

outDirPlots <- "../plots/RPPA_limma"
outDirData <- "../RPPA/Data/DEResults"

if (!dir.exists(outDirPlots)) {
  dir.create(outDirPlots)
}

if (!dir.exists(outDirData)) {
  dir.create(outDirData)
}

###############################################################################

RPPA.mat <- read.table(RPPAlevel3File,
                       header = TRUE,
                       sep = ",",
                       row.names = 1)
RPPA.meta <- read.table(MDDannoFile,
                        header = TRUE,
                        sep = ",") %>% 
  filter(str_detect(experimentalTimePoint, "0|24|48"),
         specimenID %in% colnames(RPPA.mat)) %>% 
  dplyr::select(specimenID, specimenName, ligand, experimentalTimePoint, experimentalCondition, replicate) %>% 
  mutate(experimentalCondition = fct_inorder(factor(experimentalCondition)))

RPPA.mat <- RPPA.mat %>%
  dplyr::select(all_of(unique(RPPA.meta$specimenID)))

###############################################################################
# Creating design for comparisons.
# All conditions to be compared to ctrl_0.
#limma anlaysis expects log-expression values
design <- model.matrix(~experimentalCondition, RPPA.meta)

lm <- lmFit(RPPA.mat, design)
lm <- treat(lm, lfc = logFC_threshold)

res <- decideTests(lm, p.value = pval_threshold)

lfc_values <- lm[["coefficients"]]
p_values <- lm[['p.value']]
tt <- topTreat(lm, number = Inf)

write.csv(tt, sprintf("%s/RPPA_DE_topTreat.csv", outDirData))

resTP <- matrix(res, nrow = nrow(res))[,-1]
rownames(resTP) <- rownames(res)
colnames(resTP) <- colnames(design)[-1] %>%
  str_remove("experimentalCondition")


pdf(sprintf("%s/RPPA_significantAnalytes.pdf", outDirPlots), height = 12, width = 16)
hm <- Heatmap(resTP,
              name = "significant",
              cluster_columns = FALSE,
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 9),
              show_row_names = TRUE,
              column_title = "RPPA: significant histone epigentic marks compared to ctrl_0")
draw(hm)
dev.off()
