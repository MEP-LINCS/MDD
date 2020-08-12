library(tidyverse)
library(DiffBind)
library(motifmatchr)
library(chromVAR)
library(chromVARmotifs)
library(Matrix)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)

set.seed(2019)

rdata     <- "/Users/derrickd/workspace/mcf10a_common_project/mcf10a_atacseq/archive/Data/rdata/dob_new.Rdata"
annoFile  <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"
RNAfile <- "../RNAseq/Data/MDD_RNAseq_Level3.csv"

# out
motifFamilies_out <- "../ATACseq/Data/motif/MDD_ATACseq_MotifFamilyScores.csv"
motifFamilies_out_L4 <- "../ATACseq/Data/motif/MDD_ATACseq_MotifFamilyScores_summarized.csv"
motifScores_out   <- "../ATACseq/Data/motif/MDD_ATACseq_MotifScores.csv"
motifAnno_out     <- "../ATACseq/Data/motif/MDD_ATACseq_MotifAnno_orig.csv"

###############################################################################
if( str_extract(getwd(), "[:alnum:]+$") != "R" ) {setwd("R")}

source("MDD_ATACseq_motifFunctions.R")
###############################################################################
# Loading ATACseq data
load(rdata)

# Loading motifs
motifs <- jaspar()

# Sample annotation file
sampleAnno <- 
  read.csv(annoFile, stringsAsFactors = FALSE) %>% 
  mutate(ligand = fct_inorder(ligand)) %>% 
  mutate(experimentalCondition = fct_inorder(experimentalCondition)) %>% 
  filter(ATACseq_QCpass)

# annotation table to map uniprot ids to ensembl gene ids/hugo symbols
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl", host = "uswest.ensembl.org")

annoTable <-
  getBM(attributes = c("ensembl_gene_id", 
                       "hgnc_symbol",
                       "uniprot_gn_id"),
        mart = mart)

##############################################################################
# Calculating motif scores
motifScoresObject <- calculateMotifScores(dob.counts, 
                                    motifs, sampleAnno, expect = FALSE,
                                    returnZ = TRUE, returnObject = TRUE)

motifScores <- motifScoresObject$scores

motifScores <- 
  motifScores %>%  
  data.frame %>%
  rownames_to_column("motif")

##############################################################################
# Creating motif annotations
# Motif families
mFam <-
  lapply(motifs@listData, function(x) {
    fam <- x@tags$family
    if(anyNA(fam)) {
      fam <- "None"
      }
    return(c(fam))
    }) %>% 
  stack
colnames(mFam) <- c("family", "motif")

# UniProt accession #
mAcc <-
  lapply(motifs@listData, function(x) {
    uniprot <- x@tags$acc
    if(anyNA(uniprot)) {
      uniprot <- "None"
    }
    return(c(uniprot))
  }) %>% 
  stack
colnames(mAcc) <- c("uniprot_gn_id", "motif")

motifAnno <- 
  mAcc %>% 
  left_join(mFam) %>% 
  dplyr::select(motif, family, uniprot_gn_id) %>% 
  left_join(annoTable)

write_csv(motifScores,
          path = motifScores_out)

write_csv(motifAnno, 
          path = motifAnno_out)

##############################################################################
# Summarizing motifs to families
motif_families <-
  motifAnno %>% 
  group_by(motif) %>% 
  left_join(motifAnno) %>% 
  dplyr::select(-uniprot_gn_id) %>%
  filter(!(grepl("var", motif))) %>% 
  mutate(family = str_trim(family, "right")) %>%  # cleaning names
  mutate(family = str_replace(family, "ü", "u")) %>% 
  left_join(motifScores) %>% 
  pivot_longer(contains("sid"), values_to = "motifEnrichment", names_to = "specimenID") %>% 
  group_by(motif) %>% 
  mutate(motifVar = var(motifEnrichment)) %>% 
  filter(!is.na(family)) %>%
  group_by(family, specimenID) %>% 
  summarize(motifEnrichment = mean(motifEnrichment)) %>% 
  data.frame() %>% 
  pivot_wider(names_from = specimenID,
              values_from = motifEnrichment)

motif_families <- motif_families[, 
                                 c("family", sampleAnno %>% filter(ATACseq_QCpass) %>% pull(specimenID))]
write_csv(motif_families, path = motifFamilies_out)

# Summarizing to experimental conditions

motif_families %>% 
  pivot_longer(-family, names_to = "specimenID", values_to = "motifEnrichment") %>% 
  left_join(sampleAnno) %>% 
  group_by(family, experimentalCondition) %>% 
  summarize(meanEnrichment = mean(motifEnrichment)) %>% 
  pivot_wider(names_from = experimentalCondition, values_from = meanEnrichment) %>% 
  write_csv(path = motifFamilies_out_L4)
  
