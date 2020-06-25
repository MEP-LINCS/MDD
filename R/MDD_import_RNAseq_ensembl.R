# This script imports RNAseq data into the environment.
# Sample annotations are also imported.
# Complex Heatmap annotations are also created from the sample annotations.
library(tidyverse)
library(ComplexHeatmap)
library(biomaRt)

colScript          <- "MDD_importColors_pretty.R"
RNAseqMetadataFile <- "../RNAseq/Metadata/MDD_RNAseq_sampleMetadata.csv"
RNAseqL3DataFile   <- "../RNAseq/Data/MDD_RNAseq_Level3.csv"

###############################################################################

if(!grepl("R$", getwd())) {
  setwd("R")
}

source(colScript)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "uswest.ensembl.org")

at.RNA <- getBM(attributes = c("ensembl_gene_id",
                               "hgnc_symbol",
                               "gene_biotype"), 
                mart = mart)

##############################################################################
# Importing sample annotations
sa.RNA.L3 <-
  read.csv(RNAseqMetadataFile, 
           stringsAsFactors = FALSE) %>% 
  mutate(ligand = fct_inorder(ligand),
         experimentalTimePoint = fct_inorder(as.factor(experimentalTimePoint)),
         experimentalCondition = fct_inorder(as.factor(experimentalCondition))) %>% 
  filter(RNAseq_QCpass) %>% 
  dplyr::select(-contains("RNA")) %>% 
  dplyr::select(-contains("sequencing")) %>% 
  dplyr::select(-contains("File")) %>% 
  dplyr::rename(Time = experimentalTimePoint,
                Ligand = ligand) %>% 
  # mutate(experimentalCondition = fct_relevel(as.factor(experimentalCondition),
  #                                            c("ctrl_0", "EGF_48", "IFNG_48"))) %>% 
  arrange(experimentalCondition) %>% 
  mutate(experimentalCondition = as.character(experimentalCondition)) %>% 
  mutate(Ligand = fct_recode(Ligand, 
                             "CTRL" = "ctrl",
                             "BMP2+EGF" = "BMP2",
                             "IFNG+EGF" = "IFNG",
                             "TGFB+EGF" = "TGFB"),
         experimentalCondition = fct_recode(experimentalCondition,
                                            "CTRL_0" = "ctrl_0")) %>% 
  mutate(Ligand = fct_relevel(Ligand, "CTRL", "PBS",
                              "HGF", "OSM", "EGF",
                              "BMP2+EGF", "IFNG+EGF", "TGFB+EGF")) %>% 
  arrange(Ligand, Time) %>% 
  mutate(experimentalCondition = fct_inorder(experimentalCondition))
col$Ligand <- col$Ligand[unique(as.character(sa.RNA.L3$Ligand))]

# Importing Level 3 data
RNAseqL3 <- 
  read.csv(RNAseqL3DataFile,
           header = TRUE, row.names = 1) %>% 
  data.matrix

exprFilt <- apply(RNAseqL3, 1, function(x) {x <- sum(x >= 0.5) >= 3}) %>% 
  data.frame(ensembl_gene_id = names(.),
             expressedInRNA = .,
             stringsAsFactors = FALSE)

RNAseqL3  <- RNAseqL3[, sa.RNA.L3$specimenID]

RNAseqL3Z <- 
  read.csv(RNAseqL3DataFile,
           header = TRUE, row.names = 1) %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  data.frame() %>% 
  .[, sa.RNA.L3$specimenID]

##############################################################################

# Creating Heatmap Annotations
ha.RNA.L3 <- HeatmapAnnotation(df = dplyr::select(sa.RNA.L3, 
                                                  Ligand, 
                                                  Time), col = col)

