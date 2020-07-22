library(tidyverse)
library(biomaRt)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Gviz)
library(motifmatchr)

plottingFunFile <- "../ATACseq/scripts/plottingFunctions.R"
annoFile     <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"

###############################################################################
source(plottingFunFile)

# for annotating peaks with nearest gene, etc.
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "useast.ensembl.org")

sampleAnno <- read.csv(annoFile, stringsAsFactors = FALSE) %>% 
  filter(ATACseq_QCpass) %>% 
  filter(ligand %in% c("ctrl", "EGF", "IFNG")) %>% 
  filter(experimentalTimePoint %in% c(0, 48))

###############################################################################
# Creating DiffBind samplesheet. Low-quality samples are removed.
sSheet    <- makeDBSampleSheet(sampleAnno) %>% 
  mutate(Factor = fct_recode(Factor, CTRL = "ctrl_0", EGF = "EGF_48", "IFNG+EGF" = "IFNG_48")) %>% 
  mutate(Factor = fct_relevel(Factor, "CTRL", "EGF", "IFNG+EGF")) %>% 
  arrange(Factor)
dob_temp  <- dba(sampleSheet = sSheet)

###############################################################################
# Genes of interest
outDir <- "../rsync_folder/tracks/IFNG_induced/"

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

myGene <- "ACE2"
eg <- c("BMX")
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()


myGene <- "AIM2"
eg <- "IFI16"
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()

myGene <- "NLRP3"
eg <- "OR2B11"
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()

myGene <- "BATF2"
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno)
dev.off()

myGene <- "CXCL9"
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno)
dev.off()

myGene <- "CXCL10"
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno)
dev.off()

myGene <- "CXCL11"
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno)
dev.off()


###############################################################################

myGene <- "IDO1"
eg <- "IDO2"
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()

myGene <- "HIF1A"
eg <- c("HIF1A-AS1", "HIF1A-AS2", "HIF1A-AS3")
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()


myGene <- "NDRG1"
eg <- c("CCN4")
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()

myGene <- "ACE2"
eg <- c("BMX")
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp,
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()


myGene <- "IRF1"
eg <- "IRF1-AS1"
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 5)
geneToPlot(symbol = myGene, 
           buffer = c(10000),
           basicDob = dob_temp, 
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()

myGene <- "STAT1"
eg <- c("GLS", "STAT4")
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 8)
geneToPlot(symbol = myGene, 
           buffer = c(50000),
           basicDob = dob_temp,
           mart = mart,
           sA = sampleAnno,
           extraGenes = eg)

dev.off()  

myGene <- "STAT3"
eg     <- c("STAT5B", "STAT5A", "CAVIN1")
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 8)
geneToPlot(symbol = myGene, 
           buffer = c(50000),
           basicDob = dob_temp, 
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()  

myGene <- "CD274"
eg     <- c("PLGRKT", "PDCD1LG2")
pdf(sprintf("%s/MDD_ATACseq_%sTrack.pdf", outDir, myGene), 
    height = 5, width = 8)
geneToPlot(symbol = myGene, 
           buffer = c(50000),
           basicDob = dob_temp, 
           mart = mart, 
           sA = sampleAnno,
           extraGenes = eg)
dev.off()  

