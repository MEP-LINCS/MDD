#!/usr/bin/env Rscript
#
# Summarize MCF10A RNASeq transcript abundances
# author: Daniel Derrick
#
# The purpose of this script is to import transcript abundances, summarize them
# to gene-level, and write the gene level summaries.
#
###############################################################################
# Setup
library(biomaRt)
library(tximport)
library(DESeq2)
library(rhdf5)
library(tidyverse)

annoFile   <- "../Metadata/MDD_sample_annotations.csv"
fileInfo   <- "../RNAseq/Metadata/MDD_RNAseq_sequencingInfo.csv"
outDirData     <- "../RNAseq/Data"
outDirMetadata <- "../RNAseq/Metadata"
###############################################################################
# Getting gene/transcript annotations

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

tx2g <- getBM(attributes = c("ensembl_transcript_id", 
                             "ensembl_gene_id"), 
              mart = mart)

annoTable <- getBM(attributes = c("ensembl_gene_id", 
                                  "hgnc_symbol",
                                  "gene_biotype"), 
                   mart = mart)

###############################################################################
# Getting sample annotations

if(!grepl("R$", getwd())) {setwd("R")}

sampleAnno <- 
  read.csv(annoFile, stringsAsFactors = FALSE) %>% 
  right_join(read.csv(fileInfo, stringsAsFactors = FALSE)) %>% 
  filter(specimenName != "BMP2_48_C1_C") %>% 
  mutate(abundanceFile = sprintf("../%s", abundanceFile)) %>% 
  mutate(experimentalCondition = fct_inorder(experimentalCondition))

##############################################################################
# Importing data from Kallisto .h5 files

files <- setNames(sampleAnno$abundanceFile,
                  sampleAnno$specimenID)

txi <- tximport(files,
                type = "kallisto", 
                tx2gene = tx2g,
                ignoreTxVersion = TRUE)

# txiTranscript <- tximport(files,
#                           type = "kallisto", 
#                           txOut = TRUE,
#                           tx2gene = tx2g,
#                           ignoreTxVersion = TRUE)

txiTPM        <- tximport(files,
                          type = "kallisto",
                          tx2gene = tx2g,
                          countsFromAbundance = "scaledTPM",
                          ignoreTxVersion = TRUE)

###############################################################################
# Creating DESeq2 data set from gene-summarized, non-TPM txi object
dds <- DESeqDataSetFromTximport(txi, sampleAnno, ~ experimentalCondition)
dds <- DESeq(dds)

# And comparing each experimentalCondition to the CTRL condition
# Slow step
res <- lapply(2:15, function(x) {
  X <- lfcShrink(dds, coef = x, method = "apeglm")
  X
})

names(res) <-
  resultsNames(dds)[-1] %>% 
  str_remove("ligand_") %>%
  str_remove("_vs_ctrl_0")

# Writing tximport, DESeq2, and annotation files to Rdata file
if (!dir.exists(outDirData)) {dir.create(outDirData)}

save(txi, txiTPM, 
     # txiTranscript,
     dds, sampleAnno, res,
     tx2g, annoTable,
     file = sprintf("%s/MDD_RNAseq_DESeqFull.Rdata", outDirData))

# Writing DESeq2 and annotation files to smaller Rdata file
save(dds, sampleAnno, res,
     tx2g, annoTable,
     file = sprintf("%s/MDD_RNAseq_DESeqSlim.Rdata", outDirData))
