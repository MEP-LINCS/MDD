# The purpose of this script is to create a paired ATAC-seq tornado plots
# and RNAseq heatmaps for genes in each integrated module. 

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

moduleDir <- "../misc/Ligand_tables"

hFScript  <- "MDD_RNAseq_heatmapFunctions.R"
RNAscript <- "MDD_import_RNAseq_ensembl.R"

transcriptsFile  <- "../misc/MDD_RNAseq_mostAbundantTranscripts.txt"
ATACseqSAFile    <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"
ReMapChIPSeqFile <- "../misc/CHEA3_ReMap_ChIP-seq.gmt"

promoterOutFile <- "~/Desktop/rsync_folder/allModule_promoterList.Rdata"
outDir <- "../plots/ATACseq_tornadoPlots/allModules"

###############################################################################

# Functions for ATAC-seq data
getProm <- function(txdb, annoTable, buffer = 1500, genesOnly = FALSE) {
  txdb <- keepStandardChromosomes(txdb)
  txdb <- dropSeqlevels(txdb, c("chrY", "chrM"))
  
  if (genesOnly) {
    txdb <- genes(txdb)
    prom <- promoters(txdb, upstream = buffer, downstream = buffer)
    prom$ensembl_gene_id <- aT[match(prom$gene_id, aT$entrezgene), 1]
  } else {
    prom <- promoters(txdb, upstream = buffer, downstream = buffer)
    prom$tx_name         <- str_remove(prom$tx_name, "[.][0-9]+")
    prom$ensembl_gene_id <- aT[match(prom$tx_name, aT$ensembl_transcript_id), 1]
  }
  
  filt <- is.na(prom$ensembl_gene_id)
  prom <- prom[!filt]
  
  return(prom)
}

geneGrangestoTSS <- function(x, buffer = 1000) {
  strands <- as.vector(strand(x))
  
  end(x[strands == "+"])   <- start(x[strands == "+"]) + buffer
  start(x[strands == "+"]) <- start(x[strands == "+"]) - buffer
  
  start(x[strands == "-"]) <- end(x[strands == "-"]) - buffer
  end(x[strands == "-"])   <- end(x[strands == "-"]) + buffer
  
  x
}

combineMat <- function(geneMat, combine) {
  dumSum <- function(x, y) {
    z <- x + y
    return(z)
  }
  
  geneMatAvg <-
    lapply(combine, function(x) {
      print(x)
      x_sum <-
        Reduce(dumSum, geneMat[x])
      x_avg <- (x_sum/length(x))
      return(x_avg)
    })
  
  names(geneMatAvg) <- names(combine)
  geneMatAvg <- as(geneMatAvg, "ScoreMatrixList")
  return(geneMatAvg)
}

###############################################################################

if(!grepl("R$", getwd())) {
  setwd("R")
}

source(hFScript)
source(RNAscript)

###############################################################################
# Step 0: 
# Import RNAseq L1 data. 
# Import ATACseq annotations. 
# Import sample annotations.
# Filter promoter dataset to relevant promoters only 
# (those in transcriptsFile)

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
  mutate(bamReads = paste0("../", bamReads))

txFiltTable <- read.table(transcriptsFile,
                          header = TRUE,
                          stringsAsFactors = FALSE)

prom <- getProm(TxDb.Hsapiens.UCSC.hg38.knownGene, annoTable)
prom <- prom[prom$tx_name %in% txFiltTable$ensembl_transcript_id]

#################################################
# Step 1: 
# Load in gene list.
# Annotate with promoter status.
# Annotate with IFNG_24 log2FC.
# Arrange in order of decreasing log2FC, so no subsequent ordering required.

myFiles        <- list.files(moduleDir, full.names = TRUE)
names(myFiles) <- str_replace(str_extract(myFiles, "cluster_[0-9]+"), 
                              "cluster", "module")
fT <- lapply(myFiles, read.csv,
             header = FALSE, stringsAsFactors = FALSE)

fT <-
  lapply(fT, function(x) {
    colnames(x) <- c("feature", sprintf("V%s", seq(2, 10, 1)), "type", "symbol")
    x
  }
  )

moduleRNAseqFeatures <- lapply(fT, function(x) {
  x <- x %>%
    filter(type == "RNA") %>%
    distinct(symbol) %>%
    pull(symbol)
})
moduleRNAseqFeatures <- moduleRNAseqFeatures[c(1, 5:12, 2:4)]

moduleList <- mapply(function(x, x_module) {
  x <- data.frame(hgnc_symbol = x,
                  module = x_module,
                  stringsAsFactors = FALSE) %>% 
    left_join(at.RNA) %>%
    mutate(transcriptInRNA = ensembl_gene_id %in% rownames(RNAseqL3),
           promoterInATAC  = ensembl_gene_id %in% prom$ensembl_gene_id) %>% 
    filter(transcriptInRNA,
           promoterInATAC)
  x
  }, x = moduleRNAseqFeatures,
     x_module = names(moduleRNAseqFeatures),
     SIMPLIFY = FALSE)

###############################################################################
# Step 2A: 
# Filter promoters to ensembl ids from IFNG gene list.
# Create HLA region grange and intersect with filtered promoters.
# Annotate cluster11_list with HLA status based on overlaps.
# Filter promoters to desired list (no HLA)
promList <- rep(list(prom), 12)
names(promList) <- names(moduleRNAseqFeatures)

promList <- mapply(function(promoters, moduleGenes) {
  promoters <- promoters[promoters$ensembl_gene_id %in% moduleGenes]
  promoters
  }, 
  prom = promList, 
  moduleGenes = lapply(moduleList, function(x) {x <- x %>% pull(ensembl_gene_id)}))

HLA_indices <- lapply(promList, function(promoters) {
  x <- 
    findOverlaps(GRanges("chr6:28510120-33480577"), promoters) %>% 
    data.frame %>% 
    .$subjectHits
    
  x
})

moduleListHLA <-
  mapply(function(moduleDF, promoters, HLA_ind) {
  moduleDF$promoterInHLA <- FALSE
  
  moduleDF$promoterInHLA[match(promoters$ensembl_gene_id[HLA_ind], moduleDF$ensembl_gene_id)] <- TRUE
  
  moduleDF <- filter(moduleDF, !promoterInHLA)
  
  return(moduleDF)
  }, 
  moduleDF = moduleList, promoters = promList, HLA_ind = HLA_indices, 
  SIMPLIFY = FALSE)

###############################################################################
# Step 2B:
# Import ReMap ChIP-seq data, subset to IRFs and STATs
# Combine all ChIP-seq targets into single vector
# Annotate cluster11_list with that vector
# Subset promoter to remove promoters that are not in ATACseq or are in HLA
# write to Rdata

remap_chipseq <- parse.gmt(ReMapChIPSeqFile)

irfs <- grep("IRF", names(remap_chipseq))
stats <- grep("STAT", names(remap_chipseq))[-c(3,5)]

remap_filt <- remap_chipseq[c(irfs, stats)]

remap_unique <- lapply(remap_filt, function(x) {
  x <- x$entry
  x
}) %>% 
  unlist %>% 
  unique

moduleListHLA <-
  lapply(moduleListHLA, function(x) {
    x <- x %>% 
      mutate(ChIP_target = hgnc_symbol %in% remap_unique)
    x
  })

promNoHLA <- mapply(function(promoters, moduleDF) {
  promoters <- promoters[match(moduleDF %>% 
                                 filter(promoterInATAC) %>% 
                                 filter(!promoterInHLA) %>% 
                                 pull(ensembl_gene_id),
                               promoters$ensembl_gene_id)]
  promoters
  }, 
  promoters = promList,
  moduleDF = moduleListHLA)

save(promNoHLA, sampleAnno, 
     file = promoterOutFile)

###############################################################################
# Step 3 - EXACLOUD
# Load filtered promoter list
# count from bam files
# name columns
# write to file

library(tidyverse)
library(genomation)

load("../rsync_folder/allModule_promoterList.Rdata")

geneMats <- lapply(promNoHLA, function(x) {
  x <- ScoreMatrixList(target = sampleAnno %>%
                         filter(ATACseq_QCpass) %>%
                         pull(bamReads),
                       windows = x, bin.num = 300,
                       rpm = TRUE, type = "bam", 
                       strand.aware = TRUE)
  x
  }
  )

geneMatsNamed <- lapply(geneMats, function(x) {
  names(x) <-
    sampleAnno %>%
    filter(ATACseq_QCpass) %>%
    pull(specimenName)
  x
})

save(geneMatsNamed, file = "../rsync_folder/allModule_promoterMat.Rdata")

###############################################################################
# Step 4
# Load quantified promoters
# combine conditions
# plot combined
load("~/Desktop/rsync_folder/allModule_promoterMat.Rdata")

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

temp <- 
  sampleAnno %>% 
  filter(ATACseq_QCpass) %>% 
  mutate(index = seq(1,35,1)) %>% 
  mutate(experimentalCondition = fct_inorder(as.factor(experimentalCondition))) %>% 
  dplyr::select(experimentalCondition, specimenName, index) %>% 
  group_split(experimentalCondition)

temp2 <- lapply(temp, function(x) {
  x <- x$index
  x
  }
  )

names(temp2) <- unique(sampleAnno$experimentalCondition)

geneMatsCombined <- lapply(geneMatsNamed,
                           combineMat, 
                           combine = temp2)

mapply(function(x, x_names) {
  png(sprintf("%s/%s/MDD_%s_genes_ATACseq_tornadoPlot.png", outDir, x_names, x_names),
      height = 12, width = 12, units = "in", res = 450)
  multiHeatMatrix(x,
                  winsorize = c(0, 97),
                  order = TRUE)
  dev.off()
  }, x = geneMatsCombined, x_names = names(geneMatsCombined))

###############################################################################
# Step 5
# find mean accessibility of each peak in gMat_combined
addAccOrderToModuleDF <- function(gmat, moduleDF) {
  gmat <- lapply(gmat, as, "matrix")
  gmat <- Reduce(cbind, gmat)
  
  sum_acc <- apply(gmat, 1, sum)
  order_acc <- order(sum_acc, decreasing = TRUE)
  
  orderedGenes   <- data.frame(ensembl_gene_id = moduleDF$ensembl_gene_id[order_acc],
                                   rank = seq(1, length(order_acc), 1),
                               stringsAsFactors = FALSE)
  
  moduleDF <- moduleDF %>% 
    left_join(orderedGenes) %>% 
    arrange(rank)
  
  moduleDF
}

moduleListHLAOrdered <- mapply(addAccOrderToModuleDF, 
               gmat = geneMatsCombined, moduleDF = moduleListHLA,
               SIMPLIFY = FALSE)

###############################################################################
# Step 6
# Subset RNAseq L3Z data to ensembl gene ids in Order of Promoters
# replace ensembl gene ids with hugo symbols
# heatmap without clustering

raList <- lapply(moduleListHLAOrdered, function(x) {
  ra <- rowAnnotation(df = x %>% 
                        filter(promoterInATAC) %>% 
                        filter(!promoterInHLA) %>% 
                        dplyr::select(ChIP_target) %>% 
                        dplyr::rename('IRF/STAT' = ChIP_target),
                      col = list('IRF/STAT' = c("TRUE" = "black",
                                                "FALSE" = "white")))
  ra
})

RNAseqL3Z_filtList <- lapply(moduleListHLAOrdered, function(x) {
  x <- RNAseqL3Z[x %>%
                   filter(promoterInATAC) %>%
                   filter(!promoterInHLA) %>%
                   pull(ensembl_gene_id), 
                 sa.RNA.L3 %>% 
                   pull(specimenID)]
  
  # rownames(x) <- at.RNA[match(rownames(x), at.RNA$ensembl_gene_id), ]$hgnc_symbol
  x
})

lapply(names(RNAseqL3Z_filtList), function(x) {dir.create(sprintf("%s/%s", outDir, x))})

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[1]), names(RNAseqL3Z_filtList[1])),
            height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[1]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[1]])
dev.off()

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[2]), names(RNAseqL3Z_filtList[2])),
            height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[2]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[2]])
dev.off()


pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[3]), names(RNAseqL3Z_filtList[3])),
    height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[3]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[3]])
dev.off()

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[4]), names(RNAseqL3Z_filtList[4])),
    height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[4]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[4]])
dev.off()

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[5]), names(RNAseqL3Z_filtList[5])),
    height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[5]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[5]])
dev.off()

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[6]), names(RNAseqL3Z_filtList[6])),
    height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[6]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[6]])
dev.off()

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[7]), names(RNAseqL3Z_filtList[7])),
    height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[7]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[7]])
dev.off()

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[8]), names(RNAseqL3Z_filtList[8])),
    height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[8]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[8]])
dev.off()

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[9]), names(RNAseqL3Z_filtList[9])),
    height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[9]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[9]])
dev.off()

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[10]), names(RNAseqL3Z_filtList[10])),
    height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[10]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[10]])
dev.off()

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[11]), names(RNAseqL3Z_filtList[11])),
    height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[11]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[11]])
dev.off()

pdf(sprintf("%s/%s/MDD_%s_genes_RNAseqHeatmap_tornadoOrder.pdf", outDir,
            names(RNAseqL3Z_filtList[12]), names(RNAseqL3Z_filtList[12])),
    height = 10, width = 12)
HeatmapL3C(RNAseqL3Z_filtList[[12]],
           clust = FALSE, cluster_rows = FALSE, left_annotation = raList[[12]])
dev.off()
