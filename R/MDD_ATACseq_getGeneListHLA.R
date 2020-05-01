# The purpose of this script is to get a list of ensembl_gene_ids that fall 
# within the HLA region of chr5. These will be removed from analyses
# involving ATAC-seq.
library(genomation)
library(biomaRt)
library(tidyverse)
library(cmapR)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(openxlsx)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)

theme_set(theme_cowplot())

RNAscript <- "MDD_import_RNAseq_ensembl.R"

geneFiltOutDir <- "../RNAseq/misc"

if(!grepl("R$", getwd())) {
  setwd("R")
}

###############################################################################

source(RNAscript)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",  
                host = "uswest.ensembl.org")
aT <-
  getBM(attributes = c("ensembl_gene_id", 
                       "ensembl_transcript_id"),
        mart = mart)

sampleAnno <- 
  read.csv(ATACseqSAFile, stringsAsFactors = FALSE) %>% 
  mutate(Peaks = paste0("../", Peaks)) %>% 
  mutate(bamReads = paste0("../", bamReads)) %>% 
  mutate(experimentalCondition = fct_relevel(as.factor(experimentalCondition),
                                             c("ctrl_0", "EGF_24", "EGF_48", "IFNG_24", "IFNG_48"))
  )

txFiltTable <- read.table(transcriptsFile,
                          header = TRUE,
                          stringsAsFactors = FALSE)

txdb <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb$tx_name <- str_extract(txdb$tx_name, "ENST[0-9]+")
txdb$ensembl_gene_id <- aT[match(txdb$tx_name, aT$ensembl_transcript_id), 1]

x <- 
  findOverlaps(GRanges("chr6:28510120-33480577"), txdb) %>% 
  data.frame %>% 
    .$subjectHits

geneFilt_HLA <-
  data.frame(ensembl_gene_id = unique(txdb[x, ]$ensembl_gene_id),
             HLA = TRUE)

if (!dir.exists(geneFiltOutDir)) {dir.create(geneFiltOutDir)}

write_csv(geneFilt_HLA, path = sprintf("%s/HLA_geneList.csv", geneFiltOutDir))
