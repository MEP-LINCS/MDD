# Calculate RPPA log2 fold changes and differentially expressed proteins.
#
# The purpose of this script is to calculate log2 fold changes for all conditions
# compared to ctrl_0, using limma.
# This script is intended as the starting point for other scripts using RPPA
# L2FCs and or DE proteins

library(tidyverse)
library(biomaRt)
library(limma)
library(ComplexHeatmap)

assays <- c("cycIF","RPPA", "GCP", "motifs", "RNAseq")
assay <- c("RNAseq")

logFC_threshold <- 0.25
pval_threshold  <- 0.01


get_assay_values <- function(assay){
  #browser()
  level3File <- paste0("../",assay,"/Data/MDD_",assay,"_Level3.csv")
  if(assay == "motifs") level3File <- "../ATACseq/Data/MDD_ATACseq_MotifFamilyScores.csv"
  MDDannoFile <- "../Metadata/MDD_sample_annotations.csv"
  
  mat_raw <- read.table(level3File,
                        header = TRUE,
                        sep = ",",
                        row.names = 1)
  row_sd <- apply(mat_raw, 1, sd, na.rm = TRUE)

  mat <- mat_raw[complete.cases(mat_raw) &!row_sd == 0,]
  
  meta <- read.table(MDDannoFile,
                     header = TRUE,
                     sep = ",") %>% 
    filter(str_detect(experimentalTimePoint, "0|24|48"),
           specimenID %in% colnames(mat)) %>% 
    dplyr::select(specimenID, specimenName, ligand, experimentalTimePoint, experimentalCondition, replicate) %>% 
    mutate(experimentalCondition = fct_inorder(factor(experimentalCondition)))
  
  if(assay == "cycIF"){
    mat <- mat[!str_detect(rownames(mat),"_med_|dna.*cytoplasm|laws|nucrin|centcyto|plasmem|none"),]
    mat <- log2(mat+.0001)
  }
  
  if(assay == "motifs"){
    mat_sd <- sd(unlist(mat), na.rm = TRUE)
    mat <- mat/mat_sd
  }
  
  if(assay == "RNAseq"){
    #Get annotations to convert from ensemble to HGNC
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                    dataset = "hsapiens_gene_ensembl",
                    host = "uswest.ensembl.org")
    
    annoTable <- getBM(attributes = c("ensembl_gene_id",
                                      "hgnc_symbol"),
                       mart = mart)
    mat <- mat %>%
      rownames_to_column("ensembl_gene_id") %>%
      # left_join(annoTable, by = "ensembl_gene_id") %>%
      # dplyr::select(-ensembl_gene_id) %>%
      group_by(ensembl_gene_id) %>%
      summarise(across(everything(),mean), .groups = "drop") %>%
      data.frame()
rownames(mat) <- mat$ensembl_gene_id

  }
  
  mat <- mat %>%
    dplyr::select(all_of(unique(meta$specimenID)))
  return(list(mat = mat, meta = meta))
}


###############################################################################
# Creating design for comparisons.
# All conditions to be compared to ctrl_0.
#limma anlaysis expects log-expression values

foo <- lapply(assays, function(assay){
 #browser() 
  assay_values <- get_assay_values(assay)
  
  design <- model.matrix(~experimentalCondition, assay_values[["meta"]])
  
  lm <- lmFit(assay_values[["mat"]], design)
  lm <- treat(lm, lfc = logFC_threshold)
  
  #Filter to significant condition rows
  sig_conditions <- decideTests(lm, p.value = pval_threshold)
  #Filter to non-zero condition rows
  sig_conditions <- sig_conditions@.Data[,-1]
  colnames(sig_conditions) <- str_remove(colnames(sig_conditions), "experimentalCondition") 
  feature_index <- !rowSums(abs(sig_conditions)) == 0
  features <- rownames(sig_conditions)[feature_index] 
  
  lfc_values <- lm[["coefficients"]][features,]
  p_values <- lm[['p.value']][features,]
  
  
  outDirPlots <- paste0("../plots/",assay,"_limma")
  outDirData <- paste0("../",assay,"/Data/DEResults")
  if(assay == "motifs") outDirData <- "../ATACseq/Data/DEResults"
  
  if (!dir.exists(outDirPlots)) {
    dir.create(outDirPlots)
  }
  
  if (!dir.exists(outDirData)) {
    dir.create(outDirData)
  }
  write.csv(lfc_values, sprintf("%s/%s_DE_lfc_values.csv", outDirData, assay))
  write.csv(p_values, sprintf("%s/%s_DE_p_values.csv", outDirData, assay))
  # 
  # resTP <- matrix(sig_conditions, nrow = nrow(sig_conditions))
  # rownames(resTP) <- rownames(res)
  # colnames(resTP) <- colnames(design)[-1] %>%
  #   str_remove("experimentalCondition")
  show_row_names <- FALSE
  if(nrow(sig_conditions) < 100) show_row_names <- TRUE
  pdf(sprintf("%s/%s_significantAnalytes.pdf", outDirPlots, assay), height = 12, width = 16)
  hm <- Heatmap(sig_conditions,
                name = "significant",
                cluster_columns = FALSE,
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 9),
                show_row_names = show_row_names,
                column_title = paste0(assay, ": significant features compared to ctrl_0"))
  draw(hm)
  dev.off()
  
})
