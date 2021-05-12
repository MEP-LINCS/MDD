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

assays <- c("RPPA", "GCP")


logFC_threshold <- 0.25
pval_threshold  <- 0.01

outDirPlots <- paste0("../plots/",assay,"_limma")
outDirData <- paste0("../",assay,"/Data/DEResults")

if (!dir.exists(outDirPlots)) {
  dir.create(outDirPlots)
}

if (!dir.exists(outDirData)) {
  dir.create(outDirData)
}

get_assay_values <- function(assay){
  level3File <- paste0("../",assay,"/Data/MDD_",assay,"_Level3.csv")
  MDDannoFile <- "../Metadata/MDD_sample_annotations.csv"
  
  mat <- read.table(level3File,
                    header = TRUE,
                    sep = ",",
                    row.names = 1)
  
  meta <- read.table(MDDannoFile,
                     header = TRUE,
                     sep = ",") %>% 
    filter(str_detect(experimentalTimePoint, "0|24|48"),
           specimenID %in% colnames(mat)) %>% 
    dplyr::select(specimenID, specimenName, ligand, experimentalTimePoint, experimentalCondition, replicate) %>% 
    mutate(experimentalCondition = fct_inorder(factor(experimentalCondition)))
  
  mat <- mat %>%
    dplyr::select(all_of(unique(meta$specimenID)))
  return(list(mat = mat, meta = meta))
}


###############################################################################
# Creating design for comparisons.
# All conditions to be compared to ctrl_0.
#limma anlaysis expects log-expression values

foo <- lapply(assays, function(assay){
  
  assay_values <- get_assay_values(assay)
  
  design <- model.matrix(~experimentalCondition, assay_values[["meta"]])
  
  lm <- lmFit(assay_values[["mat"]], design)
  lm <- treat(lm, lfc = logFC_threshold)
  
  #Filter to significant condition rows
  sig_conditions <- decideTests(lm, p.value = pval_threshold)
  #Filter to non-zero condition rows
  sig_conditions <- sig_conditions@.Data[,-1]
  feature_index <- !rowSums(abs(sig_conditions)) == 0
  features <- rownames(sig_conditions)[feature_index]
  
 lfc_values <- lm[["coefficients"]][features,]
  p_values <- lm[['p.value']][features,]

  write.csv(lfc_values, sprintf("%s/%s_DE_lfc_values.csv", outDirData, assay))
  write.csv(p_values, sprintf("%s/%s_DE_p_values.csv", outDirData, assay))
  
  resTP <- matrix(res, nrow = nrow(res))[,-1]
  rownames(resTP) <- rownames(res)
  colnames(resTP) <- colnames(design)[-1] %>%
    str_remove("experimentalCondition")
  
  pdf(sprintf("%s/%s_significantAnalytes.pdf", outDirPlots, assay), height = 12, width = 16)
  hm <- Heatmap(resTP,
                name = "significant",
                cluster_columns = FALSE,
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 9),
                show_row_names = TRUE,
                column_title = paste0(assay, ": significant features compared to ctrl_0"))
  draw(hm)
  dev.off()
  
})
