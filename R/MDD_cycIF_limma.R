# Calculate cycIF log2 fold changes and differentially expressed proteins.
#
# The purpose of this script is to calculate log2 fold changes for all conditions
# compared to ctrl_0, using limma.
# This script is intended as the starting point for other scripts using RPPA
# L2FCs and or DE proteins

library(tidyverse)
library(cowplot)
library(limma)
library(ComplexHeatmap)


cycIFlevel3File <- "../cycIF/Data/MDD_cycIF_Level3.csv"
MDDannoFile <- "../Metadata/MDD_sample_annotations.csv"

logFC_threshold <- 0.5
pval_threshold  <- 0.01

outDirPlots <- "../plots/cycIF_limma"
outDirData <- "../cycIF/Data/DEResults"

###############################################################################
#Limit time points to 0, 24 and 48 hours
#Limit features to intensity means, nuc and cytoplasm areas and limited textures
cycIF.mat <- read.table(cycIFlevel3File,
                       header = TRUE,
                       sep = ",",
                       row.names = 1)
cycIF.mat <- cycIF.mat[!str_detect(rownames(cycIF.mat),"_med_|dna.*cytoplasm|laws|nucrin|centcyto|plasmem|none"),] 

cycIF.meta <- read.table(MDDannoFile,
                        header = TRUE,
                        sep = ",") %>% 
  filter(str_detect(experimentalTimePoint, "0|24|48"),
    specimenID %in% colnames(cycIF.mat)) %>% 
  dplyr::select(specimenID, specimenName, ligand, experimentalTimePoint, experimentalCondition, replicate) %>% 
  mutate(experimentalCondition = fct_inorder(factor(experimentalCondition)))

cycIF.mat <- cycIF.mat %>%
  dplyr::select(all_of(unique(cycIF.meta$specimenID)))

###############################################################################
# Creating design for comparisons.
# All conditions to be compared to ctrl_0.
design <- model.matrix(~experimentalCondition, cycIF.meta)

lm <- lmFit(cycIF.mat, design)
lm <- eBayes(lm)

tt <- topTable(lm, number = Inf)
write.csv(tt, sprintf("%s/cycIF_DE_topTables.csv", outDirData))

res <- decideTests(lm, p.value = pval_threshold, lfc = logFC_threshold)

resTP <- matrix(res, nrow = nrow(res))

rownames(resTP) <- rownames(res)
colnames(resTP) <- colnames(design) %>%
  str_remove("experimentalCondition") %>%
  str_replace("[(]Intercept[)]","CTRL")
# 
if (!dir.exists(outDirPlots)) {
  dir.create(outDirPlots)
}
# 
res <- pdf(sprintf("%s/cycIF_significantAnalytes.pdf", outDirPlots), height = 12, width = 16)
hm <- Heatmap(resTP,
        name = "significant",
        cluster_columns = FALSE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 6),
        show_row_names = TRUE,
        column_title = "cycIF: significant analytes compared to ctrl_0")
draw(hm)
res <- dev.off()

if (!dir.exists(outDirData)) {
  dir.create(outDirData)
}

detect_sig_feature <- function(x){
  res <- x %>%
    unlist %>%
    str_detect("1") %>%
    sum
  !res == 0
}

sig_features <- bind_cols(resTP, sig_feature = apply(resTP, 1, detect_sig_feature), feature = rownames(resTP)) %>%
  filter(!sig_feature == 0) %>%
  select(feature)


res <- write_csv(sig_features, sprintf("%s/cycIF_DE_features.csv", outDirData))


