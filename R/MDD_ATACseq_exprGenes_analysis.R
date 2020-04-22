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

RNAscript <- "MDD_import_RNAseq_ensembl.R"

transcriptsFile <- "../misc/MDD_RNAseq_mostAbundantTranscripts.txt"
ATACseqSAFile   <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"

promoterOutFile <- "~/Desktop/rsync_folder/exprUnexpr_promoterList.Rdata"

# Functions for RNA-seq data plotting
HeatmapL3C <- function(mat, ca = ha.RNA.L3, featureName = "Z-score", showRow = FALSE, clust = FALSE, ...) {
  Heatmap(mat, 
          top_annotation = ca,
          cluster_columns = clust,
          show_row_names = showRow,
          cluster_row_slices = FALSE,
          row_names_gp = gpar(fontsize = 7),
          show_column_names = FALSE,
          name = featureName,
          col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
          ...)
}

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
                       "ensembl_transcript_id",
                       "gene_biotype"),
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

prom <- getProm(TxDb.Hsapiens.UCSC.hg38.knownGene, annoTable)
prom <- prom[prom$tx_name %in% txFiltTable$ensembl_transcript_id]

#################################################
# Step 1: 
# Load in gene list.
# Annotate with promoter status.
# Annotate with IFNG_24 log2FC.
# Arrange in order of decreasing log2FC, so no subsequent ordering required.

toPlotNoExpr <-
  exprFilt %>% 
  filter(!expressedInRNA) %>% 
  dplyr::select(ensembl_gene_id) %>% 
  left_join(at.RNA) %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>%     
  mutate(transcriptInRNA = ensembl_gene_id %in% rownames(RNAseqL3),
         promoterInATAC  = ensembl_gene_id %in% prom$ensembl_gene_id)

toPlotExpr <-
  exprFilt %>% 
  filter(expressedInRNA) %>% 
  dplyr::select(ensembl_gene_id) %>% 
  left_join(at.RNA) %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>%     
  mutate(transcriptInRNA = ensembl_gene_id %in% rownames(RNAseqL3),
         promoterInATAC  = ensembl_gene_id %in% prom$ensembl_gene_id)

tpList <- list("NoExpr" = toPlotNoExpr,
               "Expr"   = toPlotExpr)

###############################################################################
# Step 2A: 
# Filter promoters to ensembl ids from IFNG gene list.
# Create HLA region grange and intersect with filtered promoters.
# Annotate cluster11_list with HLA status based on overlaps.
# Filter promoters to desired list (no HLA)
promList <- rep(list(prom), 2)
names(promList) <- names(tpList)

promList <- mapply(function(promoters, moduleGenes) {
  promoters <- promoters[promoters$ensembl_gene_id %in% moduleGenes]
  promoters
}, 
prom = promList, 
moduleGenes = lapply(tpList, function(x) {x <- x %>% pull(ensembl_gene_id)}))

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
  moduleDF = tpList, promoters = promList, HLA_ind = HLA_indices, 
  SIMPLIFY = FALSE)

###############################################################################
# Step 2B:
# Import ReMap ChIP-seq data, subset to IRFs and STATs
# Combine all ChIP-seq targets into single vector
# Annotate cluster11_list with that vector
# Subset promoter to remove promoters that are not in ATACseq or are in HLA
# write to Rdata

remap_chipseq <- parse.gmt("/Users/derrickd/workspace/mcf10a_common_project/mcf10a_integrative_analysis/Data/ENCODE/CHEA3_GMT/CHEA3_ReMap_ChIP-seq.gmt")

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

load("../rsync_folder/exprUnexpr_promoterList.Rdata")

geneMats <- lapply(promNoHLA, function(x) {
  x <- ScoreMatrixList(target = sampleAnno %>%
                         filter(ATACseq_QCpass) %>%
                         filter(ligand %in% c("ctrl", "EGF", "IFNG")) %>%
                         filter(experimentalTimePoint %in% c(0, 48)) %>%
                         arrange(ligand) %>%
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
    filter(ligand %in% c("ctrl", "EGF", "IFNG")) %>%
    filter(experimentalTimePoint %in% c(0, 48)) %>%
    arrange(ligand) %>%
    pull(specimenName)
  x
})

save(geneMatsNamed, file = "../rsync_folder/exprUnexpr_promoterMat.Rdata")

###############################################################################
# Step 4
# Load quantified promoters
# combine conditions
# plot combined
load("~/Desktop/rsync_folder/exprUnexpr_promoterMat.Rdata")

dir <- "../plots/ATACseq_tornadoPlots/exprUnexpr/"

if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

geneMatsCombined <- lapply(geneMatsNamed,
                           combineMat, 
                           combine = list("CTRL" = 1:4,
                                          "EGF" = 5:8,
                                          "IFNG+EGF" = 9:10)
)

mapply(function(x, x_names) {
  png(sprintf("%s/MDD_%s_genes_ATACseq_tornadoPlot.png", dir, x_names),
      height = 12, width = 12, units = "in", res = 450)
  multiHeatMatrix(x,
                  winsorize = c(0, 97),
                  order = TRUE)
  dev.off()
}, x = geneMatsCombined, x_names = names(geneMatsCombined))

###############################################################################
# Step 5
# find mean accessibility of each peak in gMat_combined
addSumToModuleDF <- function(gmat, moduleDF) {
  gmat <- lapply(gmat, as, "matrix")
  gmat <- Reduce(cbind, gmat)
  
  sum_acc <- apply(gmat, 1, sum)
  head(sum_acc) %>% print
  order_acc <- order(sum_acc, decreasing = TRUE)
  head(sum_acc) %>% print
  
  orderedGenes   <- data.frame(ensembl_gene_id = moduleDF$ensembl_gene_id,
                               accSum = sum_acc,
                               accRank = order_acc,
                               stringsAsFactors = FALSE)
  orderedGenes
}

moduleListHLAOrdered <- mapply(addSumToModuleDF, 
                               gmat = geneMatsCombined, moduleDF = promNoHLA,
                               SIMPLIFY = FALSE)


geneMatAsMats <- lapply(geneMatsNamed$NoExpr, as, "matrix")

geneMatrixBig <- Reduce(cbind, geneMatAsMats)
rownames(geneMatrixBig) <- promNoHLA$NoExpr$ensembl_gene_id

geneAccSum  <- apply(geneMatrixBig, 1, sum)
at.RNA <-
  at.RNA %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE)

geneAcc <- data.frame(ensembl_gene_id = names(geneAccSum),
                      accessibility = geneAccSum) %>% 
  mutate(accessible = accessibility > 500) %>% 
  inner_join(at.RNA)

ggplot(geneAcc, aes(accessibility)) +
  geom_density() +
  scale_x_log10()


heatMatrix(geneMatsCombined$NoExpr$CTRL[geneAcc$accessible],
           order = TRUE)
heatMatrix(geneMatsCombined$NoExpr$CTRL[!geneAcc$accessible],
                order = TRUE)
