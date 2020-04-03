# This script creates chromatin accessibility matrices for the MCF10A MDD
# ATACseq data from per-sample alignment (.bam) and narrowPeak (.bed) files.

library(tidyverse)
library(biomaRt)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

if(!grepl("R$", getwd())) {
  setwd("R")
}

annoFile     <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"

outDirData     <- c("../ATACseq/Data")
outDirMetadata <- c("../ATACseq/Metadata")

##############################################################################

# Functions
makeDBSampleSheet <- function(sampleAnno) {
  # renames the columns of a sample annotation dataframe
  # to create a new dataframe that can be used as input to Diffbind.
  ss <- sampleAnno %>% 
    dplyr::rename("SampleID"   = specimenID,
                  "Collection" = collection,
                  "Time"       = experimentalTimePoint) %>% 
    mutate(Replicate           = as.factor(replicate),
           Factor              = as.factor(ligand),
           PeakCaller          = "narrow",
           Peaks               = paste0("../", Peaks),
           bamReads            = paste0("../", bamReads)
           )
  ss
}

mergePeaksandCount <- function(ssheet,
                               score = DBA_SCORE_TMM_READS_EFFECTIVE_CPM) {
  # constructs a consensus peakset based on samples in ssheet, removes
  # non-standard chromosomes, and counts reads in peaks using merged peakset
  
  filterPeaks <- function(peakset, filter) {
    # filters peakset (as a GRanges object) to only contain peaks within
    # chromosomes listed in "filter"
    
    peakset <- peakset[seqnames(peakset) %in% filter]
    seqlevels(peakset) <- filter
    peakset
  }
  
  filter <- sprintf("chr%s", c(seq(1, 22, 1), "X"))
  
  dob_temp          <- dba(sampleSheet = ssheet)
  consensus_peakset <- dba.peakset(dob_temp, bRetrieve = TRUE)
  consensus_peakset <- filterPeaks(consensus_peakset, filter)
  
  dob               <- dba.count(dob_temp, 
                                 peaks = consensus_peakset,
                                 score = score)
  
  return(dob)
}

makeMatrixFromDob <- function(x) {
  # extracts the peak accessibility matrix from a DiffBind object,
  # formats as a dataframe, and renames rows with genomic locus of peak.
  
  x <- dba.peakset(x, bRetrieve = TRUE)
  x <- data.frame(x)
  
  peakNames <- paste(x$seqnames, 
                     x$start,
                     x$end, 
                     sep = "_")
  
  x.mat      <- x[, -c(1:5)]
  rownames(x.mat) <- peakNames
  x.mat
}

collapseReps <- function(x, sA) {
  mat_medians <- sapply(unique(sA$experimentalCondition), function(X) {
    Xids <-
      sA %>% 
      filter(experimentalCondition == X) %>% 
      pull(specimenID)
    Xmat <- x[, Xids]    # matrix of one condition's replicates
    if (length(Xids) > 1) {
      Xmat_median <- apply(Xmat, 1, median)
    } else {
      Xmat_median <- Xmat
    }
    return(Xmat_median)
  })
  colnames(mat_medians) <- unique(sA$experimentalCondition)
  mat_medians <- data.frame(mat_medians)
  return(mat_medians)
}

###############################################################################

if(!grepl("R$", getwd())) {
  setwd("R")
}

# Getting annotations for annotating peaks with nearest gene, etc.
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = 'apr2019.archive.ensembl.org')

# Importing sample metadata
sampleAnno <- read.csv(annoFile, stringsAsFactors = FALSE)

# Creating DiffBind samplesheet. Low-quality samples are removed.
sSheet <- makeDBSampleSheet(filter(sampleAnno, ATACseq_QCpass))

# Getting counts and TMM-normalized expression
dobCounts   <- mergePeaksandCount(sSheet, 
                                   score = DBA_SCORE_READS)

dobTMM      <- mergePeaksandCount(sSheet, 
                                   score = DBA_SCORE_TMM_READS_EFFECTIVE)
consensusPeakset   <- dba.peakset(dobTMM, bRetrieve = TRUE)

MDD_ATACseq_Level1 <- makeMatrixFromDob(dobCounts)
MDD_ATACseq_Level3 <- log2(makeMatrixFromDob(dobTMM) + 1)
MDD_ATACseq_Level4 <- collapseReps(MDD_ATACseq_Level3, 
                                   sA = filter(sampleAnno, ATACseq_QCpass))

###############################################################################

# Annotating ATAC-seq peaks
tx2g <- getBM(attributes = c("ensembl_gene_id",
                             "entrezgene",
                             "hgnc_symbol"), 
              mart = mart) %>% 
  dplyr::rename(entrez_gene  = entrezgene) %>% 
  filter(!duplicated(entrez_gene))

peakAnnotations <- 
  annotatePeak(peak = consensusPeakset, TxDb = txdb, level = "gene") %>% 
  data.frame() %>% 
  mutate(entrez_gene = as.integer(geneId)) %>%
  dplyr::select(-geneId, -contains("sid")) %>% 
  left_join(tx2g) %>% 
  mutate(peak = paste(seqnames, start, end, sep = "_")) %>% 
  dplyr::select(peak, everything())

###############################################################################
# Writing accessibility matrices to csv files

if (!dir.exists(outDirData)) {
  dir.create(dir, recursive = TRUE)
}

MDD_ATACseq_Level1 %>% 
  rownames_to_column("peak") %>% 
  write_csv(path = sprintf("%s/%s",
                           outDirData, "MDD_ATACseq_Level1.csv"))

MDD_ATACseq_Level3 %>% 
  rownames_to_column("peak") %>% 
  write_csv(path = sprintf("%s/%s",
                           outDirData, "MDD_ATACseq_Level3.csv"))

MDD_ATACseq_Level4 %>% 
  data.frame() %>% 
  rownames_to_column("peak") %>% 
  write_csv(path = sprintf("%s/%s",
                           outDirData, "MDD_ATACseq_Level4.csv"))

# Writing peak metadata to csv files
if (!dir.exists(outDirMetadata)) {
  dir.create(dir, recursive = TRUE)
}

write_csv(peakAnnotations,
          path = sprintf("%s/%s",
                         outDirMetadata, "MDD_ATACseq_peakMetadata.csv"))

# Writing DiffBind object for differential expression analysis
MDD_ATACseq_dob <- dobTMM
save(MDD_ATACseq_dob, 
     file = sprintf("%s/%s",
                    outDirData, "MDD_ATACseq_DiffBindObject.Rdata"))


