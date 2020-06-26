library(DESeq2)
library(tidyverse)

windowWidth <- 5000

transcriptsFile <- "../misc/MDD_RNAseq_mostAbundantTranscripts.txt"
RNAseqResultsFile <- "../RNAseq/MDD_RNAseq_shrunkenResults.Rdata"
ATACseqResultsFile <- sprintf("../ATACseq/MDD_ATACseqGeneSum_%sbp_shrunkenResults.Rdata",
                              windowWidth)
ATACseqL1File <- "../ATACseq/data/MDD_ATACseq_Level1.csv"
ATACseqL1AnnoFile <- "../ATACseq/Metadata/MDD_ATACseq_peakMetadata.csv"
saFile     <- "../ATACseq/Metadata/MDD_ATACseq_sampleMetadata.csv"
colScript <- "MDD_importColors_pretty.R"
cTable <- "../RNAseq/misc/MDD_RNAseq_combinedTable.csv"
source(colScript)

outDirPlots <- "../plots/ATAC_RNA_scatterplots"

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

aT <- getBM(mart = mart, 
            attributes = c("ensembl_gene_id",
                           "entrezgene_id",
                           "hgnc_symbol"))
sampleAnno <- read.csv(saFile,
                       header = TRUE,
                       stringsAsFactors = FALSE) %>% 
  mutate(specimenName = fct_inorder(as.factor(specimenName)))


MDD_ATACseq_L1      <- read.csv(ATACseqL1File, 
                                header = TRUE, stringsAsFactors = FALSE)
MDD_ATACseq_L1_Anno <- read.csv(ATACseqL1AnnoFile, 
                                header = TRUE, stringsAsFactors = FALSE)
cTable <- read.csv(cTable, stringsAsFactors = FALSE)

MDD_ATACseq_L1_geneSum <-
  MDD_ATACseq_L1 %>% 
  left_join(dplyr::select(MDD_ATACseq_L1_Anno, peak, distanceToTSS)) %>% 
  filter(abs(distanceToTSS) <= windowWidth) %>% 
  dplyr::select(-distanceToTSS) %>% 
  pivot_longer(-peak, names_to = "specimenID", values_to = "accessibility") %>% 
  mutate(specimenID = fct_inorder(as.factor(specimenID))) %>% 
  left_join(dplyr::select(MDD_ATACseq_L1_Anno, peak, ensembl_gene_id)) %>% 
  dplyr::select(-peak) %>% 
  dplyr::filter(!is.na(ensembl_gene_id)) %>% 
  group_by(specimenID, ensembl_gene_id) %>% 
  summarize(sumAcc = sum(accessibility)) %>% 
  pivot_wider(values_from = "sumAcc", names_from = "specimenID") %>% 
  column_to_rownames("ensembl_gene_id") %>% 
  data.matrix

# dds <- DESeqDataSetFromMatrix(MDD_ATACseq_L1_geneSum,
#                               sampleAnno %>% filter(ATACseq_QCpass), ~experimentalCondition)
# dds <- DESeq(dds)
# 
# conditionsToTest <-
#   sampleAnno %>%
#   filter(ATACseq_QCpass) %>%
#   mutate(experimentalCondition = fct_inorder(as.factor(experimentalCondition))) %>%
#   filter(experimentalTimePoint != 0) %>%
#   group_by(experimentalCondition) %>%
#   summarize(n = n()) %>%
#   filter(n >= 2) %>%
#   pull(experimentalCondition) %>%
#   as.character()
# 
# resultsList <- lapply(conditionsToTest, function(x) {
#   print(x)
#   X <- lfcShrink(dds, contrast = c("experimentalCondition", x, "ctrl_0"))
#   X
#   })
# names(resultsList) <- conditionsToTest
# save(resultsList, file = ATACseqResultsFile)

load(ATACseqResultsFile)

resFixedATAC <- 
  mapply(function(x, x_nm) {
    x <- x %>% 
      data.frame %>% 
      rownames_to_column("ensembl_gene_id") %>% 
      mutate(experimentalCondition = x_nm) %>% 
      dplyr::select(ensembl_gene_id, experimentalCondition, padj, log2FoldChange) %>% 
      dplyr::rename(ATAC_l2fc = log2FoldChange,
                    ATAC_padj = padj)
    x
  }, x = resultsList, x_nm = names(resultsList), 
  SIMPLIFY = FALSE)

resFixedATACLong <-
  Reduce(bind_rows, resFixedATAC)

load(RNAseqResultsFile)

res_RNA_df <-
  mapply(function(x, x_nm) {
    X <- 
      x %>% 
      mutate(experimentalCondition = x_nm) %>% 
      dplyr::select(ensembl_gene_id, experimentalCondition,
                    padj, log2FoldChange) %>% 
      dplyr::rename(RNA_l2fc = log2FoldChange,
                    RNA_padj = padj)
    X
  }, 
  x = res_RNA_df, x_nm = names(res_RNA_df),
  SIMPLIFY = FALSE)

resRNALong <- Reduce(bind_rows, res_RNA_df)

resJoined <-
  inner_join(resFixedATACLong, resRNALong) %>% 
  left_join(sampleAnno %>% 
              distinct(experimentalCondition, .keep_all = TRUE) %>% 
              dplyr::select(experimentalCondition, ligand, experimentalTimePoint)
  ) %>%
  drop_na()


tableTp <- 
  resJoined %>% 
  filter(experimentalCondition == "IFNG_48") %>% 
  dplyr::select(ensembl_gene_id, experimentalCondition, contains("l2fc")) %>% 
  left_join(cTable)

pdf(sprintf("%s/MDD_ATACseqRNAseq_IFNG48_L2FC_Scatterplot.pdf", outDirPlots),
    height = 5, width = 6)
tableTp %>% 
  ggplot(aes(ATAC_l2fc, RNA_l2fc, col = fct_rev(as.factor(Induced)))) +
  geom_point(alpha = .5, size = 1, shape = 16) +
  geom_smooth(size = .75, method = "lm", color = "black") +
  scale_color_manual("Induced", 
                     values = c("TRUE" = "red", "FALSE" = "gray50")) +
  facet_wrap(~experimentalCondition) +
  xlab("Accessibility\nlog2(FC)") +
  ylab("Gene Expression\nlog2(FC)")
tableTp %>% 
  ggplot(aes(ATAC_l2fc, RNA_l2fc, col = fct_rev(as.factor(Induced)))) +
  geom_point(alpha = .5, size = 1, shape = 16) +
  geom_smooth(size = .75, method = "lm", color = "black") +
  scale_color_manual("Induced", 
                     values = c("TRUE" = "red", "FALSE" = "gray50")) +
  xlab("Accessibility\nlog2(FC)") +
  ylab("Gene Expression\nlog2(FC)")
dev.off()
