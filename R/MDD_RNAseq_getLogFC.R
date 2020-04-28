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

PBS48_RNA <-
  lfcShrink(dds, coef = 3, type = "apeglm")  %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(annoTable)

BMP248_RNA <-
  lfcShrink(dds, coef = 5, type = "apeglm")  %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(annoTable)

IFNG48_RNA <-
  lfcShrink(dds, coef = 7, type = "apeglm")  %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(annoTable)

TGFB48_RNA <-
  lfcShrink(dds, coef = 9, type = "apeglm")  %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(annoTable)

HGF48_RNA <-
  lfcShrink(dds, coef = 11, type = "apeglm")  %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(annoTable)

OSM48_RNA <-
  lfcShrink(dds, coef = 13, type = "apeglm")  %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(annoTable)

EGF48_RNA <-
  lfcShrink(dds, coef = 15, type = "apeglm")  %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(annoTable)

RNAseq_resultsTables <- list(
  "PBS_48" = PBS48_RNA,
  "BMP2_48" = BMP248_RNA,
  "IFNG_48" = IFNG48_RNA,
  "TGFB_48" = TGFB48_RNA,
  "HGF_48" = HGF48_RNA,
  "OSM_48" = OSM48_RNA,
  "EGF_48" = EGF48_RNA)

save(dds, RNAseq_resultsTables,
     file = "../RNAseq/MCF10A_RNAseq_DESeq2Results48.Rdata")
