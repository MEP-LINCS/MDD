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

collapseReps <- function(x, sA) {
  mat_medians <- sapply(unique(sA$experimentalCondition), function(X) {
    print(X)
    X <- enquo(X)
    
    Xids <-
      sA %>% 
      filter(experimentalCondition == !!X) %>% 
      pull(specimenID)
    
    print(Xids)
    Xmat <- x[, Xids]    # matrix of one condition's replicates
    Xmat_median <- apply(Xmat, 1, median)
    return(Xmat_median)
  })
  colnames(mat_medians) <- unique(sA$experimentalCondition)
  return(mat_medians)
}

###############################################################################

if(!grepl("R$", getwd())) {setwd("R")}

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = 'apr2019.archive.ensembl.org')

sampleAnno <- read.csv(annoFile, stringsAsFactors = FALSE)

# Creating metadata table, making paths relative to working directory
sampleAnnoRNAseq <-
  sampleAnno %>% 
  right_join(read.csv(fileInfo, stringsAsFactors = FALSE)) %>% 
  filter(specimenName != "BMP2_48_C1_C") %>% 
  mutate(abundanceFile = sprintf("../%s", abundanceFile))

# Constructing table to match transcripts to ensembl gene ids
tx2g <- getBM(attributes = c("ensembl_transcript_id", 
                             "ensembl_gene_id"), 
              mart = mart)

# Importing alignment files using tximport and summarizing to gene level
files <- setNames(sampleAnnoRNAseq$abundanceFile, sampleAnnoRNAseq$specimenID)
txi <- tximport(files,
                type = "kallisto", 
                tx2gene = tx2g,
                ignoreTxVersion = TRUE)

# making DESeqDataSet
dds <- DESeqDataSetFromTximport(txi, sampleAnnoRNAseq, ~ experimentalCondition)

dds <- DESeq(dds)

# Getting counts, log2(fpkm + 1), and median-summarized log2(fpkm + 1)
counts <- 
  dds %>% 
  counts %>% 
  data.frame() %>% 
  rownames_to_column(var = "ensembl_gene_id")

# for GEO submission:
# fpkm <- 
#   fpkm(dds) %>% 
#   data.frame() %>% 
#   rownames_to_column(var = "ensembl_gene_id")

log2fpkm <- 
  log2(fpkm(dds) + 1) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ensembl_gene_id")

log2fpkmMedian <- 
  log2(fpkm(dds) + 1) %>% 
  data.frame() %>% 
  collapseReps(sampleAnnoRNAseq) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ensembl_gene_id")

if (!dir.exists(outDirData)) {
  dir.create(outDirData, recursive = TRUE)
}

# Saving expression data to csv file
write_csv(counts,
          path = paste(outDirData, "MDD_RNAseq_Level1.csv", sep = "/"))

# for GEO submission:
# write_csv(counts,
#           path = paste(outDirData, "MDD_MCF10A_RNAseq_countsEnsembl.csv", sep = "/"))
# write_csv(fpkm,
#           path = paste(outDirData, "MDD_MCF10A_RNAseq_fpkmEnsembl.csv", sep = "/"))

write_csv(log2fpkm,
          path = paste(outDirData, "MDD_RNAseq_Level3.csv", sep = "/"))

write_csv(log2fpkmMedian,
          path = paste(outDirData, "MDD_RNAseq_Level4.csv", sep = "/"))

###############################################################################
# Creating table to match ensembl gene IDs to HUGO symbols
annoTable <- getBM(attributes = c("ensembl_gene_id",
                                  "hgnc_symbol"), 
                   mart = mart,
                   filters = "ensembl_gene_id",
                   values = rownames(dds),
                   uniqueRows = TRUE)  %>% 
  filter(ensembl_gene_id %in% rownames(dds)) %>% 
  filter(!duplicated(ensembl_gene_id)) %>% 
  mutate(hgnc_symbol = case_when(
    hgnc_symbol != "" ~ hgnc_symbol,
    hgnc_symbol == "" ~ "NA"
  ))

if (!dir.exists(outDirMetadata)) {dir.create(outDirMetadata, recursive = TRUE)}

write_csv(annoTable, 
          path = sprintf("%s/MDD_RNAseq_geneAnnotations.csv", outDirMetadata))

###############################################################################
# Creating sample metadata to write

sampleAnnoRNAseqFull <-
  sampleAnno %>% 
  right_join(read.csv(fileInfo, stringsAsFactors = FALSE)) %>% 
  dplyr::select(-directory,
                -contains("ATAC"),
                -contains("L1000"),
                -contains("RPPA"),
                -contains("GCP"),
                -contains("cycIF"),
                -contains("IF"),
                -MDACC_name)

write_csv(sampleAnnoRNAseqFull, 
          path = sprintf("%s/MDD_RNAseq_sampleMetadata.csv", outDirMetadata))

###############################################################################

save(dds, annoTable, sampleAnno, counts, log2fpkm, log2fpkm_EGFnorm,
     file = paste(outDirData, "MDD_RNAseq_DESeqDataSet.Rdata", sep = "/"))
