#!/usr/bin/env Rscript
#
# MCF10A RNASeq write gene expression matrices
# author: Daniel Derrick
#
# The purpose of this script is to write gene expression matrices,
# in various units (counts, Rlog-transformed counts, TM, 
# log2(TM + 1)), from the MDD RNAseq dataset.
#
###############################################################################
# Setup
library(tidyverse)
library(DESeq2)

sampleAnnoFile   <- "../RNAseq/Metadata/MDD_RNAseq_sampleMetadata.csv"
inDirData        <- "../RNAseq/Data"

outDirCounts_all             <- "../RNAseq/Data/geneExpression/allGenes"
outDirCounts_proteinCoding   <- "../RNAseq/Data/geneExpression/proteinCoding"

###############################################################################
# Loading DESeq2 Data Set and tximport objects
if (!grepl("R$", getwd())) {setwd("R")}

load(sprintf("%s/MDD_RNAseq_DESeqFull.Rdata", inDirData))

###############################################################################
# Getting counts, Rlog(counts), TPM, and log2(TPM + 1)
countsAll <- 
  dds %>% 
  counts %>% 
  data.frame() %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  left_join(annoTable) %>% 
  dplyr::select(ensembl_gene_id, hgnc_symbol, gene_biotype, everything())

rlogAll <- 
  rlog(dds) %>% 
  assay %>% 
  data.frame() %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  left_join(annoTable) %>% 
  dplyr::select(ensembl_gene_id, hgnc_symbol, gene_biotype, everything())

TPMAll <-
  txiTPM$counts %>% 
  data.frame() %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  left_join(annoTable) %>% 
  dplyr::select(ensembl_gene_id, hgnc_symbol, gene_biotype, everything())

log2TPMAll <-
  log2(txiTPM$counts + 1) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  left_join(annoTable) %>% 
  dplyr::select(ensembl_gene_id, hgnc_symbol, gene_biotype, everything())

# Writing gene expression matrices to outDirCounts_all directory
if (!dir.exists(outDirCounts_all)) {
  dir.create(outDirCounts_all, recursive = TRUE)
  }

countsAll %>% 
  write_csv(sprintf("%s/MDD_RNAseq_counts_allGenes.csv", outDirCounts_all))

rlogAll %>% 
  write_csv(sprintf("%s/MDD_RNAseq_Rlogcounts_allGenes.csv", outDirCounts_all))

TPMAll %>% 
  write_csv(sprintf("%s/MDD_RNAseq_scaledTPM_allGenes.csv", outDirCounts_all))

log2TPMAll %>% 
  write_csv(sprintf("%s/MDD_RNAseq_log2scaledTPM_allGenes.csv", outDirCounts_all))

###############################################################################
# Subsetting matrices to only include protein-coding genes

countsProteinCoding <- 
  countsAll %>% 
  filter(gene_biotype == "protein_coding")

rlogProteinCoding <- 
  rlogAll %>% 
  filter(gene_biotype == "protein_coding")

TPMProteinCoding <-
  TPMAll %>% 
  filter(gene_biotype == "protein_coding")

log2TPMProteinCoding <-
  log2TPMAll %>% 
  filter(gene_biotype == "protein_coding")

# Writing gene expression matrices to outDirCounts_proteinCoding directory
if (!dir.exists(outDirCounts_proteinCoding)) {
  dir.create(outDirCounts_proteinCoding, recursive = TRUE)
}

countsProteinCoding %>% 
  write_csv(sprintf("%s/MDD_RNAseq_counts_proteinCoding.csv",
                    outDirCounts_proteinCoding))

rlogProteinCoding %>% 
  write_csv(sprintf("%s/MDD_RNAseq_Rlogcounts_proteinCoding.csv", 
                    outDirCounts_proteinCoding))

TPMProteinCoding %>% 
  write_csv(sprintf("%s/MDD_RNAseq_scaledTPM_proteinCoding.csv",
                    outDirCounts_proteinCoding))

log2TPMProteinCoding %>% 
  write_csv(sprintf("%s/MDD_RNAseq_log2scaledTPM_proteinCoding.csv",
                    outDirCounts_proteinCoding))

