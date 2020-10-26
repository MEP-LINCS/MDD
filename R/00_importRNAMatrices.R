#!/usr/bin/env Rscript
#
# Import MCF10A MDD RNASeq gene expression data
# author: Daniel Derrick
#
# The purpose of this script is to import gene expression matrices from the
# MCF10A MDD RNAseq dataset. The data imported by this script 
#
###############################################################################
# Setup
library(biomaRt)
library(tidyverse)

if (!exists("whichGenes")) {whichGenes <- "proteinCoding"}
if (!exists("identifier")) {identifier <- "hgnc_symbol"}
if (!exists("returnAsMatrix")) {returnAsMatrix <- FALSE}
if (!exists("makeHA")) {makeHA <- FALSE}
if (!exists("biotypeInAnnoTable")) {biotypeInAnnoTable <- FALSE}

if (!grepl("R$", getwd())) {setwd("R")}

inDirData      <- "../RNAseq/Data"
inDirMetadata  <- "../RNAseq/Metadata/"
colScript          <- "MDD_importColors_pretty.R"

###############################################################################
# Getting gene/transcript annotations

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

if (biotypeInAnnoTable) {
  annoTable <- getBM(attributes = c("ensembl_gene_id", 
                                    "gene_biotype",
                                    "hgnc_symbol"),
                     mart = mart)
} else {
  annoTable <- getBM(attributes = c("ensembl_gene_id", 
                                    "hgnc_symbol"),
                     mart = mart)
}

###############################################################################
# Creating color palette

source(colScript)                          
   


###############################################################################
# Importing sample annotations

saFile <- list.files(inDirMetadata, 
                     pattern = "sampleMetadata.csv", full.names = TRUE)
sampleAnno_MDD <- read.csv(saFile, stringsAsFactors = FALSE) %>% 
  dplyr::rename(Ligand = ligand,
                Time   = experimentalTimePoint) %>% 
  mutate(Ligand = as.factor(Ligand),
         experimentalCondition = as.factor(experimentalCondition)) %>% 
    mutate(Ligand = fct_recode(Ligand, 
                               "CTRL" = "ctrl",
                               "BMP2+EGF" = "BMP2",
                               "IFNG+EGF" = "IFNG",
                               "TGFB+EGF" = "TGFB"),
           experimentalCondition = fct_recode(experimentalCondition,
                                              "CTRL_0" = "ctrl_0")) %>% 
  mutate(Ligand = fct_relevel(Ligand, 
                              names(col_MDD$Ligand))) %>%
  arrange(Ligand, Time) %>% 
  mutate(experimentalCondition = fct_inorder(experimentalCondition))
  










###############################################################################
# Importing RNAseq matrices

TPMFile <- list.files(sprintf("%s/geneExpression/%s",
                              inDirData, whichGenes),
                      pattern = "_scaledTPM",
                      full.names = TRUE)
TPM_MDD    <- read.csv(TPMFile, stringsAsFactors = FALSE)


log2TPMFile <- list.files(sprintf("%s/geneExpression/%s",
                                  inDirData, whichGenes),
                          pattern = "log2scaledTPM",
                          full.names = TRUE)
log2TPM_MDD <- read.csv(log2TPMFile, stringsAsFactors = FALSE)

###############################################################################
# Selecting identifiers

if (identifier != "all") {
  
  rSums <- apply(TPM_MDD[, -c(1:3)], 1, sum)
  
  TPM_MDD <- 
    TPM_MDD[order(rSums, decreasing = TRUE), ] %>% 
    filter(!duplicated(identifier)) %>% 
    filter(identifier != "") %>% 
    dplyr::select(identifier, contains("sid"))
  
  log2TPM_MDD <- 
    log2TPM_MDD[order(rSums, decreasing = TRUE), ] %>% 
    filter(!duplicated(identifier)) %>% 
    filter(identifier != "") %>% 
    dplyr::select(identifier, contains("sid"))
}

###############################################################################
# Setting rownames and changing to matrix

if (returnAsMatrix) {
  TPM_MDD <- 
    TPM_MDD %>% 
    rownames_to_column(identifier) %>% 
    as.matrix()
  
  log2TPM_MDD <- 
    log2TPM_MDD %>% 
    rownames_to_column(identifier) %>% 
    as.matrix()
}

###############################################################################
# Making heatmap annotation

if (makeHA) {
  library(ComplexHeatmap)
  
  ha_MDD_L3 <- HeatmapAnnotation(df = sampleAnno_MDD %>% 
                                  dplyr::select(Ligand),
                                col = col_MDD)
  
  ha_MDD_L4 <- HeatmapAnnotation(df = sampleAnno_MDD %>% 
                                  dplyr::distinct(Ligand),
                                col = col_MDD)
  
  ha_MDD_L4T0 <- HeatmapAnnotation(df = sampleAnno_MDD %>% 
                                    dplyr::distinct(Ligand) %>% 
                                    filter(Ligand != "CTRL"),
                                  col = col_MDD)
}

