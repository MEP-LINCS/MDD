## ----setup, message=FALSE, warning=FALSE---------------------------------
library(DESeq2)
library(tidyverse)
library(plotly)
library(biomaRt)
library(genomation)
library(cowplot)

transcriptsFile <- "../misc/MDD_RNAseq_mostAbundantTranscripts.txt"
RNAseqResultsFile <- "../RNAseq/MDD_RNAseq_shrunkenResults.Rdata"
ATACseqResultsFile <- "../ATACseq/MDD_ATACseqPromoter_shrunkenResults.Rdata"
combinedTableFile <- "../RNAseq/misc/MDD_RNAseq_combinedTable.csv"
# colScript <- "MDD_importColors_pretty.R"

promCovDir <- "../ATACseq/Data/promoterCoverage"
saFile     <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"
orderFile  <- "../ATACseq/Data/promoterCoverage/order.txt"

outDir <- "../plots/ATACseqAcc_vs_RNAseqExpr"

importRNAFile <- "MDD_import_RNAseq_ensembl.R"

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

aT <- getBM(mart = mart, 
            attributes = c("ensembl_gene_id",
                           "hgnc_symbol"))

theme_set(theme_cowplot())

###############################################################################
# source(colScript)
source(importRNAFile)

sa.ATAC.L3 <- read.csv(saFile,
                       header = TRUE,
                       stringsAsFactors = FALSE) %>% 
  mutate(specimenName = fct_inorder(as.factor(specimenName)))

###############################################################################

promCov <- readBed(sprintf("%s/promoterCov_01.bed", promCovDir))
promCov$score <- NULL
coverage_01 <- mcols(promCov)
coverage_02 <- read.table(sprintf("%s/coverages_02.tsv", promCovDir))
coverage_03 <- read.table(sprintf("%s/coverages_03.tsv", promCovDir))
coverage_04 <- read.table(sprintf("%s/coverages_04.tsv", promCovDir))

myOrder <- str_replace(readLines(sprintf("%s/order.txt", promCovDir)), "../", "ATACseq/")
mcolsAll <- bind_cols(data.frame(coverage_01), coverage_02, coverage_03, coverage_04)

# adding specimenIDs
colnames(mcolsAll) <- c("ensembl_transcript_id", 
                        sa.ATAC.L3 %>%
                          arrange(match(bamReads, myOrder)) %>% 
                          pull(specimenID)
                        )

# adding specimenIDs
MDD_ATACseq_promAcc <- mcolsAll[, c("ensembl_transcript_id",
                                    sa.ATAC.L3 %>% 
                                      filter(ATACseq_QCpass) %>%
                                      pull(specimenID))] %>% 
  column_to_rownames("ensembl_transcript_id") %>% 
  data.matrix


## ----filterToBestTranscripts---------------------------------------------
bestTranscripts <- read.table(transcriptsFile, header = TRUE) %>% 
  filter(!is.na(ensembl_gene_id))

MDD_ATACseq_promAcc <- MDD_ATACseq_promAcc[as.character(bestTranscripts$ensembl_transcript_id), ]
rownames(MDD_ATACseq_promAcc) <- bestTranscripts$ensembl_gene_id

## ----runDESeq2-----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(MDD_ATACseq_promAcc,
                              sa.ATAC.L3 %>% filter(ATACseq_QCpass), ~experimentalCondition)
dds <- DESeq(dds)

ATACseqPromotersL3 <- vst(dds) %>% 
  assay() %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id")

ATACseqL3Long <- 
  ATACseqPromotersL3 %>% 
  pivot_longer(-ensembl_gene_id, names_to = "specimenID", values_to = "accessibility") %>% 
  left_join(dplyr::select(sa.RNA.L3, specimenID, experimentalCondition, Ligand, Time))

###############################################################################

RNAseqL3Long <-
  RNAseqL3 %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  pivot_longer(-ensembl_gene_id, names_to = "specimenID", values_to = "expression") %>% 
  left_join(dplyr::select(sa.RNA.L3, specimenID, experimentalCondition, Ligand, Time))

###############################################################################
L3JoinedLong <-
  inner_join(ATACseqL3Long, RNAseqL3Long) 

L3JFilt <-
  L3JoinedLong %>% 
  filter(Ligand == "CTRL")

###############################################################################
cTable <- read.csv(combinedTableFile, stringsAsFactors = FALSE)

L3JFiltAnno <- 
  L3JFilt %>% 
  left_join(dplyr::select(cTable, experimentalCondition, ensembl_gene_id, expressed, accessible, geneType, HLAregion))

L3JFiltAnnoSummarized <-
  L3JFiltAnno %>% 
  group_by(ensembl_gene_id, expressed, accessible, geneType, HLAregion) %>% 
  summarize(meanExpr = median(expression),
            meanAcc = median(accessibility)) %>% 
  ungroup
###############################################################################

if(!dir.exists(outDir)) {dir.create(outDir)}
L3JFiltAnnoSummarized$accessible %>% anyNA

pdf(sprintf("%s/MDD_ATACseqTSS_vs_RNAseqFPKM_accLabels.pdf", outDir))
L3JFiltAnnoSummarized %>% 
  filter(geneType == "protein_coding") %>% 
  filter(!HLAregion) %>% 
  # slice(1:10) %>% 
  ggplot(aes(x = meanAcc, y = meanExpr, color = accessible)) +
  geom_point(alpha = .25, size = .75, shape = 16) +
  xlab("promoter accessibility\n(vst counts)") +
  ylab("gene expression\n(log2(fpkm + 1))") +
  ggtitle("TSS accessibililty vs gene expression", subtitle = "CTRL_0, protein-coding genes") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_color_discrete("Accessible")
dev.off()

pdf(sprintf("%s/MDD_ATACseqTSS_vs_RNAseqFPKM_noLabels.pdf", outDir))
L3JFiltAnnoSummarized %>% 
  filter(geneType == "protein_coding") %>% 
  filter(!HLAregion) %>% 
  ggplot(aes(x = meanAcc, y = meanExpr)) +
  geom_point(alpha = .25, size = .75, shape = 16) +
  xlab("promoter accessibility\n(vst counts)") +
  ylab("gene expression\n(log2(fpkm + 1))") +
  ggtitle("TSS accessibililty vs gene expression", subtitle = "CTRL_0, protein-coding genes") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
dev.off()
