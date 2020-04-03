# The purpose of this script is to compute differential chromatin accessibility
# signatures from the MCF10A ATAC-seq data.

library(DiffBind)
library(tidyverse)
library(ChIPseeker)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

annoFile <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"
dobFile  <- "../ATACseq/Data/MDD_ATACseq_DiffBindObject.Rdata"

##############################################################################

# Functions
makeMasks <- function(x, colname, sA) {
  colname <- enquo(colname)
  sA_t <- sA %>% 
    mutate(temp = (!!colname == x)) %>% 
    pull(temp)
  return(sA_t)
}

makeContrasts <- function(dob, x, ref) {
  x <- grep(ref, x, invert = TRUE, value = TRUE)
  for (i in x) {
    dob <- dba.contrast(dob, dob$masks[[i]], dob$masks[[ref]],
                        i, ref)
  }
  return(dob)
}

getSignatures <- function(x, threshold = .05) {
  signatureList <- GRangesList()
  for (i in 1:length(x$contrasts)) {
    signatureList[[i]] <- dba.report(x, contrast = i, 
                                     th = threshold, method = DBA_EDGER)
  }
  names(signatureList) <- sapply(x$contrasts, function(x) x[[1]])
  return(signatureList)
}

annotateSig <- function(sig) {
  sig <- annotatePeak(sig, TxDb = txdb, level = "gene")
  sig <- as.GRanges(sig)
  sig <-
    sig %>% 
    as.data.frame() %>% 
    dplyr::rename(entrezgene = geneId) %>%
    left_join(tx2g)
  return(invisible(sig))
}

###############################################################################

if(!grepl("R$", getwd())) {
  setwd("R")
}

load(dobFile)

sampleAnno <- 
  read.csv(annoFile, stringsAsFactors = FALSE) %>% 
  filter(ATACseq_QCpass)

eCs <- 
  sampleAnno %>% 
  pull(experimentalCondition) %>% 
  unique()

# Generating signatures with ctrl_0 as reference
dob.T0 <- MDD_ATACseq_dob

dob.T0$masks <- sapply(eCs, makeMasks,
                       colname = experimentalCondition, sA = sampleAnno,
                       simplify = FALSE)

dob.T0 <- makeContrasts(dob.T0, eCs, "ctrl_0")
dob.T0 <- dba.analyze(dob.T0, method = DBA_EDGER)

# Generating signatures with EGF timecourse as reference
# dob.TC <- MDD_ATACseq_dob
# 
# dob.TC$masks <- sapply(eCs, makeMasks,
#                        colname = experimentalCondition, sA = sampleAnno,
#                        simplify = FALSE)
# dob.TC <- makeContrasts(dob.TC, 
#                         x = grep("24", eCs, value = TRUE),
#                         ref = "EGF_24")
# dob.TC <- makeContrasts(dob.TC, 
#                         x = grep("48", eCs, value = TRUE),
#                         ref = "EGF_48")
# dob.TC <- dba.analyze(dob.TC, method = DBA_EDGER)

###############################################################################
# Getting differentially accessible regions as GRanges objects
signaturesT0fdr05  <- getSignatures(dob.T0, threshold = .05)
signaturesT0noFilt <- getSignatures(dob.T0, threshold = 1)

# signaturesTCfdr05  <- getSignatures(dob.TC, threshold = .05)
# signaturesTCnoFilt <- getSignatures(dob.TC, threshold = 1)

# ###############################################################################
# Annotating signatures
txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = 'http://apr2019.archive.ensembl.org')

tx2g <- getBM(attributes = c("ensembl_gene_id",
                             "entrezgene",
                             "hgnc_symbol"),
              mart = mart) %>%
  filter(!duplicated(entrezgene))
tx2g$entrezgene <-  as.character(tx2g$entrezgene)

signaturesT0fdr05Anno  <- lapply(signaturesT0fdr05, annotateSig)
signaturesT0noFiltAnno <- lapply(signaturesT0noFilt, annotateSig)

signaturesTCfdr05Anno  <- lapply(signaturesTCfdr05,  annotateSig)
signaturesTCnoFiltAnno <- lapply(signaturesTCnoFilt, annotateSig)

###############################################################################
# Saving signatures
dir <- "../ATACseq/Data/signatures/ctrl_0/fdr_05"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
mapply(function(sig, name) {
  sig <- data.frame(sig)
  filt <- apply(sig, 2, anyNA)
  sig <- sig[, !filt]
  
  write.csv(sig,
            row.names = FALSE,
            file = sprintf("%s/MDD_ATACseq_%s_Signature.csv", dir, name))
}, sig = signaturesT0fdr05, name = names(signaturesT0fdr05))
dir <- "../ATACseq/Data/signatures/ctrl_0/fdr_05/anno"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
mapply(function(sig, name) {
  sig <- data.frame(sig)
  
  write.csv(sig,
            row.names = FALSE,
            file = sprintf("%s/MDD_ATACseq_%s_SignatureAnno.csv", dir, name))
}, sig = signaturesT0fdr05Anno, name = names(signaturesT0fdr05Anno))


dir <- "../ATACseq/Data/signatures/ctrl_0/noFilt"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
mapply(function(sig, name) {
  sig <- data.frame(sig)
  filt <- apply(sig, 2, anyNA)
  sig <- sig[, !filt]
  
  write.csv(sig,
            row.names = FALSE,
            file = sprintf("%s/MDD_ATACseq_%s_Signature.csv", dir, name))
}, sig = signaturesT0noFilt, name = names(signaturesT0noFilt))
dir <- "../ATACseq/Data/signatures/ctrl_0/noFilt/anno"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
mapply(function(sig, name) {
  sig <- data.frame(sig)
  
  write.csv(sig,
            row.names = FALSE,
            file = sprintf("%s/MDD_ATACseq_%s_SignatureAnno.csv", dir, name))
}, sig = signaturesT0noFiltAnno, name = names(signaturesT0noFiltAnno))




dir <- "../ATACseq/Data/signatures/EGF_TC/fdr_05"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
mapply(function(sig, name) {
  sig <- data.frame(sig)
  filt <- apply(sig, 2, anyNA)
  sig <- sig[, !filt]
  
  write.csv(sig,
            row.names = FALSE,
            file = sprintf("%s/MDD_ATACseq_%s_Signature.csv", dir, name))
}, sig = signaturesTCfdr05, name = names(signaturesTCfdr05))
dir <- "../ATACseq/Data/signatures/EGF_TC/fdr_05/anno"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
mapply(function(sig, name) {
  sig <- data.frame(sig)
  
  write.csv(sig,
            row.names = FALSE,
            file = sprintf("%s/MDD_ATACseq_%s_SignatureAnno.csv", dir, name))
}, sig = signaturesTCfdr05Anno, name = names(signaturesTCfdr05Anno))


dir <- "../ATACseq/Data/signatures/EGF_TC/noFilt"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
mapply(function(sig, name) {
  sig <- data.frame(sig)
  filt <- apply(sig, 2, anyNA)
  sig <- sig[, !filt]
  
  write.csv(sig,
            row.names = FALSE,
            file = sprintf("%s/MDD_ATACseq_%s_Signature.csv", dir, name))
}, sig = signaturesTCnoFilt, name = names(signaturesTCnoFilt))
dir <- "../ATACseq/Data/signatures/EGF_TC/noFilt/anno"
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}
mapply(function(sig, name) {
  sig <- data.frame(sig)
  
  write.csv(sig,
            row.names = FALSE,
            file = sprintf("%s/MDD_ATACseq_%s_SignatureAnno.csv", dir, name))
}, sig = signaturesTCnoFiltAnno, name = names(signaturesTCnoFiltAnno))

