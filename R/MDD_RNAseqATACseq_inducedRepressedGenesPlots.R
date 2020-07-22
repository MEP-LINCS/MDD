# The purpose of this script is to create tornado plots and heat maps for 
# induced and repressed genes
library(tidyverse)
library(cowplot)

library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)

library(ComplexHeatmap)
library(circlize)
library(genomation)

# Scripts
RNAscript <- "MDD_import_RNAseq_ensembl.R"
heatmapFunctionScript <- "MDD_RNAseq_heatmapFunctions.R"

# Tables
ATACseqSAFile   <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"
# ATACseqAccTable <- "ATACseq/Metadata/MDD_ATACseq_promoterAccessibilityTable.csv"
transcriptsFile <- "../misc/MDD_RNAseq_mostAbundantTranscripts.txt"
combinedTable   <- "../RNAseq/misc/MDD_RNAseq_combinedTable.csv"

# Out directories
outDirPlots   <- "../plots/ATACseq_promoterAccessibility_inducedRepressed_final/"

# Functions for ATAC-seq and RNA-seq data
getProm <- function(txdb, annoTable, buffer = 1000, genesOnly = FALSE) {
  # This function returns a GRanges object of transcript promoter coordinates,
  # annotated with transcript id and ensembl gene id
  
  txdb <- keepStandardChromosomes(txdb)
  txdb <- dropSeqlevels(txdb, c("chrY", "chrM"))
  
  if (genesOnly) {
    txdb <- genes(txdb)
    prom <- promoters(txdb, upstream = buffer, downstream = buffer)
    prom$ensembl_gene_id <- annoTable[match(prom$gene_id, 
                                            annoTable$entrezgene), 1]
  } else {
    prom <- promoters(txdb, upstream = buffer, downstream = buffer)
    prom$tx_name         <- str_remove(prom$tx_name, "[.][0-9]+")
    prom$ensembl_gene_id <- annoTable[match(prom$tx_name, 
                                            annoTable$ensembl_transcript_id), ]$ensembl_gene_id
  }
  
  filt <- is.na(prom$ensembl_gene_id)
  prom <- prom[!filt]
  
  return(prom)
}

combineMat <- function(geneMat, combine) {
  # This function is used to combine accessibility matrices for multiple 
  # replicates of the same treatment, by taking the replicate mean.
  
  dumSum <- function(x, y) {
    z <- x + y
    return(z)
  }
  
  geneMatAvg <-
    lapply(combine, function(x) {
      # print(x)
      x_sum <-
        Reduce(dumSum, geneMat[x])
      x_avg <- (x_sum/length(x))
      return(x_avg)
    })
  
  names(geneMatAvg) <- names(combine)
  geneMatAvg <- as(geneMatAvg, "ScoreMatrixList")
  return(geneMatAvg)
}

makeTornadoPlots <- function(gmat, myBreaks = c(0, .1, .75, 2)) {
  # This function winsorizes scoreMatrix-type and creates ComplexHeatmap 
  # accessibility heatmaps ordered by decreasing accessibility.
  winsor <- function (x, fraction=.05) {
    if(length(fraction) != 1 || fraction < 0 ||
       fraction > 0.5) {
      stop("bad value for 'fraction'")
    }
    
    lim <- quantile(x, probs=c(fraction, 1-fraction))
    # x[ x < lim[1] ] <- lim[1]
    x[ x > lim[2] ] <- lim[2]
    return(x)
  }
  
  gmat <- lapply(gmat, as, "matrix")
  
  gmat <- lapply(gmat, winsor, fraction = .02)
  
  Heatmaps <- mapply(function(x, x_nm, b = myBreaks) {
    x <- Heatmap(x, show_row_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
                 name = "accessibility",
                 column_title = x_nm,
                 col = colorRamp2(breaks = b, c("darkblue", "blue", "yellow", "red")))
    x
  },
  x = gmat, x_nm = names(gmat))
  
  return(Heatmaps)
}

source(heatmapFunctionScript)

###############################################################################
# IMPORTING RNASEQ DATA
source(RNAscript)

###############################################################################
# IMPORTING ATACSEQ METADATA
sa.ATAC <- 
  read.csv(ATACseqSAFile, stringsAsFactors = FALSE) %>% 
  mutate(Peaks = paste0("../", Peaks)) %>%
  mutate(bamReads = paste0("../", bamReads)) %>%
  mutate(bigWig = str_replace(str_remove(bamReads, "ATACseq/Data/Level0/bam/"), ".bam", ".bw"),
         bwFile = str_replace(bamReads, "bam", "bigWig"))  %>% 
  mutate(bwFile = str_replace(bwFile, "bam", "bw")) %>% 
  mutate(experimentalCondition = fct_relevel(as.factor(experimentalCondition),
                                             c("ctrl_0", "EGF_24", "EGF_48", "IFNG_24", "IFNG_48"))
  )

txFiltTable <- read.table(transcriptsFile,
                          header = TRUE,
                          stringsAsFactors = FALSE)

prom <- getProm(TxDb.Hsapiens.UCSC.hg38.knownGene, txFiltTable)

###############################################################################
# Creating tornado plot for genes induced in IFNG_48
cTable <- 
  read.csv(combinedTable, stringsAsFactors = FALSE)

tpList <- list("Induced" = 
                 cTable %>% 
                 filter(experimentalCondition == "IFNG_48") %>% 
                 filter(Induced) %>% 
                 filter(geneType == "protein_coding") %>% 
                 filter(!HLAregion) %>%
                 arrange(desc(promCovTornado)) %>% 
                 dplyr::distinct(ensembl_gene_id),
               "Repressed"   = 
                 cTable %>% 
                 filter(experimentalCondition == "IFNG_48") %>% 
                 filter(Repressed) %>% 
                 filter(geneType == "protein_coding") %>% 
                 filter(!HLAregion) %>%
                 arrange(desc(promCovTornado)) %>% 
                 dplyr::distinct(ensembl_gene_id)
)

promList <- mapply(function(promoters, tpGenes) {
  promoters <- promoters[match(tpGenes, promoters$ensembl_gene_id)]
  promoters
  }, 
  promoters = setNames(rep(list(prom), length(tpList)),
                     names(tpList)), 
  tpGenes = lapply(tpList, function(x) {x <- x %>% pull(ensembl_gene_id)}))

sa.geneMats <- 
  sa.ATAC %>% 
  filter(ATACseq_QCpass) %>%
  filter(ligand %in% c("ctrl", "EGF", "IFNG")) %>%
  filter(experimentalTimePoint %in% c(0, 48)) %>%
  arrange(ligand)

geneMats <- lapply(promList, function(x) {
  x <- ScoreMatrixList(target = pull(sa.geneMats, bwFile),
                       windows = x, 
                       strand.aware = TRUE)
  names(x) <- sa.geneMats$specimenName
  return(x)
}
)

geneMatsCombined <- lapply(geneMats,
                           combineMat, 
                           combine = list("CTRL" = 1:4,
                                          "EGF" = 5:8,
                                          "IFNG+EGF" = 9:10)
)


tPlotsInduced <- makeTornadoPlots(geneMatsCombined$Induced)

if (!dir.exists(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PDF", outDirPlots))) {
  dir.create(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PDF", outDirPlots), recursive = TRUE)
}

pdf(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PDF/MDD_ATACseq_InducedIFNG_tornadoPlotCTRL.pdf", outDirPlots),
    height = 10, width = 3)
tPlotsInduced$CTRL
dev.off()

pdf(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PDF/MDD_ATACseq_InducedIFNG_tornadoPlotEGF48.pdf", outDirPlots),
    height = 10, width = 3)
tPlotsInduced$EGF
dev.off()


pdf(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PDF/MDD_ATACseq_InducedIFNG_tornadoPlotIFNG48.pdf", outDirPlots),
    height = 10, width = 3)
tPlotsInduced$`IFNG+EGF`
dev.off()

if (!dir.exists(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PNG_600ppi", outDirPlots))) {
  dir.create(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PNG_600ppi", outDirPlots), recursive = TRUE)
}

png(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PNG_600ppi/MDD_ATACseq_InducedIFNG_tornadoPlot_600ppi.png", outDirPlots),
    height = 10, width = 7, units = "in", res = 600)
tPlotsInduced$CTRL + tPlotsInduced$EGF + tPlotsInduced$`IFNG+EGF`
dev.off()

png(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PNG_600ppi/MDD_ATACseq_InducedIFNG_tornadoPlotCTRL_600ppi.png", outDirPlots),
    height = 10, width = 3, units = "in", res = 600)
tPlotsInduced$CTRL
dev.off()
png(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PNG_600ppi/MDD_ATACseq_InducedIFNG_tornadoPlotEGF48_600ppi.png", outDirPlots),
    height = 10, width = 3, units = "in", res = 600)
tPlotsInduced$EGF
dev.off()
png(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlot_PNG_600ppi/MDD_ATACseq_InducedIFNG_tornadoPlotIFNG48_600ppi.png", outDirPlots),
    height = 10, width = 3, units = "in", res = 600)
tPlotsInduced$`IFNG+EGF`
dev.off()

# png(sprintf("%s/MDD_ATACseq_InducedIFNG_tornadoPlotCTRL.png", outDirPlots),
#     height = 10, width = 3, units = "in", res = 300)
# HeatmapL3C(MDD_RNAseq_FPKM[tpList$Induced$ensembl_gene_id, ],
# cluster_rows = FALSE,
# featureName = "mean-centered\nlog2(fpkm + 1)")
# HeatmapL4C(MDD_RNAseq_FPKM_Medians[tpList$Induced$ensembl_gene_id, ],
# cluster_rows = FALSE,
# featureName = "mean-centered\nlog2(fpkm + 1)")

###############################################################################
# RNAseq Heatmaps of induced genes
inducedGenesAll <-
  cTable %>%
  filter(Induced) %>%
  filter(geneType == "protein_coding") %>%
  filter(!HLAregion) %>%
  dplyr::distinct(ensembl_gene_id) %>%
  pull(ensembl_gene_id)

inducedGenesIFNG48 <-
  cTable %>%
  filter(experimentalCondition == "IFNG_48") %>%
  filter(Induced) %>%
  filter(geneType == "protein_coding") %>%
  filter(!HLAregion) %>%
  arrange(desc(promCovTornado)) %>%
  dplyr::distinct(ensembl_gene_id) %>%
  pull(ensembl_gene_id)

inducedGenesNotIFNG <-
  setdiff(inducedGenesAll, inducedGenesIFNG48)

pdf(sprintf("%s/MDD_RNAseq_InducedHeatmap_allConditions.pdf",
            outDirPlots), height = 7, width = 7)
HeatmapL3(RNAseqL3[c(inducedGenesAll), ], bks = c(0, .5, 4))
HeatmapL3C(RNAseqL3Z[c(inducedGenesAll), ], featureName = "z-score")
dev.off()

pdf(sprintf("%s/MDD_RNAseq_InducedHeatmap_allConditions_IFNGseparate.pdf",
            outDirPlots), height = 7, width = 7)
HeatmapL3(RNAseqL3[c(inducedGenesIFNG48, inducedGenesNotIFNG), ],
          split = fct_inorder(
            as.factor(
              rep(c("IFNG_48", "Other"),
                  c(length(inducedGenesIFNG48),
                    length(inducedGenesNotIFNG))))),
        # cluster_row_slices = FALSE,
        bks = c(0, .5, 4))
HeatmapL3C(RNAseqL3Z[c(inducedGenesIFNG48, inducedGenesNotIFNG), ],
           featureName = "z-score",
           split = fct_inorder(
             as.factor(
               rep(c("IFNG_48", "Other"),
                   c(length(inducedGenesIFNG48),
                     length(inducedGenesNotIFNG)))))
           )
dev.off()


pdf(sprintf("%s/MDD_RNAseq_InducedHeatmap_allConditions_IFNGseparate_ordered.pdf",
            outDirPlots), height = 7, width = 7)
orderL3 <- hclust(dist(RNAseqL3[c(inducedGenesNotIFNG), ]))$order
HeatmapL3(RNAseqL3[c(inducedGenesIFNG48, inducedGenesNotIFNG[orderL3]), ],
          split = fct_inorder(
            as.factor(
              rep(c("IFNG_48", "Other"),
                  c(length(inducedGenesIFNG48),
                    length(inducedGenesNotIFNG))))),
          cluster_rows = FALSE,
          bks = c(0, .5, 4))

orderL3Z <- hclust(dist(RNAseqL3Z[c(inducedGenesNotIFNG), ]))$order
HeatmapL3C(RNAseqL3Z[c(inducedGenesIFNG48, inducedGenesNotIFNG[orderL3Z]), ],
           featureName = "z-score",
        split = fct_inorder(
          as.factor(
            rep(c("IFNG_48", "Other"),
                c(length(inducedGenesIFNG48),
                  length(inducedGenesNotIFNG)))))
        )
dev.off()


##############################################################################
# DEGEne Stackedbarplots

pdf(sprintf("%s/MDD_RNAseq_DEGeneStackedBarplot.pdf", outDirPlots),
    height = 6, width = 6)
cTable %>% 
  filter(DECategory != "None") %>% 
  mutate(experimentalCondition = fct_inorder(experimentalCondition)) %>% 
  ggplot(aes(experimentalCondition)) +
  geom_bar(aes(fill = DECategory)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # scale_color_manual(values = col$Ligand) +
  scale_fill_manual("Category", 
                    values = c("Induced" = "firebrick3",
                               "Repressed" = "steelblue3",
                               "DE" = "gray50"))
  # ggtitle("Differentially Expressed Protein-Coding Genes") 
dev.off()

# I think the math could be easy on this, but I don't have the raw numbers to try it myself. 
# Here is my thought - In comparing panels B and E - 
# is it possible to add up all the 'ligand' DE genes, induced genes, and repressed genes 
# and then compare these values to the respective values in B.
# This analysis should show that many of the DE gene changes are conserved across ligands, 
# whereas the induced are mostly ligand specific. In other words what % of DE genes show up
# in alteast two ligand conditions, and what % of induced genes show up in atleast two ligand conditions.
# Thoughts on this?

cTable %>% 
  filter(geneType == "protein_coding") %>% 
  filter(DECategory %in% c("DE", "Repressed")) %>% 
  dplyr::select(ensembl_gene_id, Ligand) %>% 
  distinct() %>% 
  group_by(ensembl_gene_id) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  ungroup %>% 
  mutate(ligandSpecific = n == 1) %>% 
  group_by(ligandSpecific) %>% 
  summarize(n = n())
1472 + 1115

1472/2587

## DE Genes - how many are ligand specific?
cTable %>% 
  filter(geneType == "protein_coding") %>% 
  filter(DECategory == "DE") %>% 
  dplyr::select(ensembl_gene_id, Ligand) %>% 
  distinct() %>% 
  group_by(ensembl_gene_id) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  ungroup %>% 
  mutate(ligandSpecific = n == 1) %>% 
  group_by(ligandSpecific) %>% 
  summarize(n = n())

1462/(1462 + 810)
1462 + 810
## Induced Genes - how many are ligand specific?
cTable %>% 
  filter(geneType == "protein_coding") %>% 
  filter(DECategory == "Induced") %>% 
  dplyr::select(ensembl_gene_id, Ligand) %>% 
  distinct() %>% 
  group_by(ensembl_gene_id) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  ungroup %>% 
  mutate(ligandSpecific = n == 1) %>% 
  group_by(ligandSpecific) %>% 
  summarize(n = n())

476/(476 + 149)

forDonut <-
  cTable %>% 
  filter(geneType == "protein_coding") %>%
  group_by(ensembl_gene_id) %>% 
  summarize(Category = case_when(any(Induced, na.rm = TRUE) ~ "Induced",
                                 any(Repressed, na.rm = TRUE) ~ "Repressed",
                                 any(DE, na.rm = TRUE) ~ "Differentially Expressed",
                                 !any(DE, na.rm = TRUE) & any(expressed, na.rm = TRUE) ~ "Expressed, no change",
                                 !any(expressed) ~ "Not Expressed")) %>% 
  mutate(Category = fct_relevel(as.factor(Category), c("Not Expressed",
                                                       "Expressed, no change",
                                                       "Differentially Expressed",
                                                       "Induced", "Repressed"))) %>% 
  ungroup %>% 
  group_by(Category) %>% 
  summarize(n = n()) %>% 
  ungroup %>% 
  mutate(fraction = n/sum(n))
forDonut$ymax = cumsum(forDonut$fraction)
forDonut$ymin = c(0, head(forDonut$ymax, n = -1))

lenGenes <- 
  cTable %>% 
  filter(geneType == "protein_coding") %>% 
  pull(ensembl_gene_id) %>% 
  unique %>% 
  length

a <- ggdonutchart(forDonut, "n", "n", 
                  fill = "Category") +
  ggtitle(sprintf("%s Protein-Coding Genes", lenGenes)) +
  scale_fill_manual(values = c("Not Expressed" = "gray30",
                               "Expressed, no change" = "gray95",
                               "Differentially Expressed" = "darkorchid2",
                               "Induced" = "firebrick3",
                               "Repressed" = "steelblue3"))
a <- ggpar(a, legend = "right")

if(!dir.exists(outDirPlots)) {dir.create(outDirPlots)}
pdf(sprintf("%s/MDD_RNAseq_geneCategoryDonutPlot.pdf", outDirPlots), height = 5, width = 6)
a
dev.off()
