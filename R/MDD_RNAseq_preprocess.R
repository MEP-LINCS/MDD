#preprocess MDD RNAseq data
#convert ensembl names to HGNC
library(tidyverse)
library(biomaRt)


pk_preprocess_level3 <- function(df, type){
  # median summarize the replicates

  #median summarize the replicates of each feature's T0 values
  df_med <- df %>%
    group_by(feature, experimentalCondition) %>%
    summarise(value = median(value, na.rm = TRUE)) %>%
    ungroup() %>%
    rename(condition = experimentalCondition) %>%
    mutate(Type = type)
  return(df_med)
}

#####

#Get annotations to convert from ensemble to HGNC
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "uswest.ensembl.org")

annoTable <- getBM(attributes = c("ensembl_gene_id",
                                  "hgnc_symbol"),
                   mart = mart)

md <- read_csv("metadata/MDD_sample_annotations.csv")

RNA_values <- read_csv("RNAseq/Data/MDD_RNAseq_Level3.csv") %>%
  left_join(annoTable, by = "ensembl_gene_id") %>%
  dplyr::select(-ensembl_gene_id) %>%
  rename(feature = hgnc_symbol) %>%  
  gather(specimenID, value = value, -feature) %>%
  left_join(md, by = "specimenID") %>%
  dplyr::select(specimenID, specimenName, feature, value, experimentalCondition, ligand, experimentalTimePoint, replicate)

res <- write_csv(RNA_values, "RNAseq/Data/MDD_RNAseq_Level3_HGNC.csv")

RNA_med <- pk_preprocess_level3(df = RNA_values, type = "RNAseq") %>%
  write_csv("RNAseq/Data/MDD_RNAseq_med.csv")
