# This script generates a DESeq2 results table for the comparison of
# condition IFNG_48 vs ctrl_0.

library(biomaRt)
library(tidyverse)
library(DESeq2)

if(!grepl("R$", getwd())) {setwd("R")}

load("../RNAseq/Data/MDD_RNAseq_DESeqDataSet.Rdata")

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

annoTable <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   mart = mart)

IFNG48_RNA <-
  lfcShrink(dds, coef = 7, type = "apeglm")  %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(annoTable)

save(dds, IFNG48_RNA, 
     file = "../RNAseq/MCF10A_RNAseq_IFNG48_Results.Rdata")
