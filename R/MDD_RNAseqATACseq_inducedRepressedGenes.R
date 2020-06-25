# The purpose of this script is to create a large annotated file of 
# ensembl_gene_ids annotated with:
# 1. expression status: log2(fpkm + 1) ≥ 0.5
# 2. differential expression: induced, repressed, DE, DEcat
# 3. accessiblity: accessible, accessibleIFNG
# 4. promAcc: promoter coverage
# 5. HLA: HLAregion (TRUE/FALSE)

library(tidyverse)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Parameters
RNAseq_lfc_threshold  <- 1.5
RNAseq_padj_threshold <- 0.01

# Tables
HLA_filtFile <- "../RNAseq/misc/HLA_geneList.csv"
ATAC_accFile <- "../ATACseq/Metadata/MDD_ATACseq_promoterAccessibilityTable.csv"
ATAC_promFile <- "../ATACseq/Metadata/MDD_ATACseq_promoterAccessibilityScores.csv"

# Scripts
RNAscript <- "../R/MDD_import_RNAseq_ensembl.R"

# Data
RNAseqResultsFile <- "../RNAseq/MDD_RNAseq_shrunkenResults.Rdata"

# In directories
bwDir <- "../ATACseq/Data/Level0/bigWig"
promoterDir <- "../ATACseq/Data/misc"

# Out directories
outDirPlots   <- "../plots/RNAseqATACseq"
outDirExcel   <- "../RNAseq/misc/IFNG_inducedRepressedGenes"
exprFiltOutFile <- "../RNAseq/misc/MDD_RNAseq_exprFilt.csv"
# indRepFiltOutFile <- "../RNAseq/misc/MDD_RNAseq_inducedRepressed.csv"
everythingOutFile <- "../RNAseq/misc/MDD_RNAseq_combinedTable.csv"

###############################################################################
# STEP 1: CREATING EXPRESSION FILTER
source(RNAscript)

exprFilt <-
  RNAseqL3 %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  pivot_longer(-ensembl_gene_id, names_to = "specimenID", values_to = "log2fpkm") %>% 
  left_join(sa.RNA.L3, by = "specimenID") %>% 
  group_by(experimentalCondition, ensembl_gene_id) %>% 
  summarize(expressed = sum(log2fpkm >= 0.5) >= 2) %>% 
  ungroup %>% 
  left_join(distinct(sa.RNA.L3, experimentalCondition,
                     .keep_all = TRUE)) %>% 
  dplyr::select(-specimenID) %>%
  left_join(at.RNA) %>%                 # Add gene biotype information
  mutate(gene_biotype = as.factor(gene_biotype)) %>% 
  mutate(geneType = fct_collapse(gene_biotype,
                                 pseudogene = grep("pseudogene", levels(gene_biotype), 
                                                   value = TRUE),
                                 other = grep("_gene", levels(gene_biotype), 
                                              value = TRUE)))

###############################################################################
# STEP 2: CREATING DIFFERENTIAL EXPRESSION TABLE
load(RNAseqResultsFile)

deFilt <-
  Reduce(bind_rows, res_RNA_df) %>% 
  mutate(padj = replace_na(padj, 1)) %>% 
  dplyr::select(-hgnc_symbol)
  
deExprFilt <-
  exprFilt %>% 
  left_join(deFilt) %>% 
  left_join(exprFilt %>% 
              filter(experimentalCondition == "CTRL_0") %>% 
              dplyr::rename(expressedCtrl = expressed) %>% 
              dplyr::select(ensembl_gene_id, expressedCtrl),
            by = "ensembl_gene_id") %>% 
  # mutate(experimentalCondition = fct_inorder(as.factor(experimentalCondition))) %>% 
  mutate(DE = (expressed &
                 abs(log2FoldChange) >= RNAseq_lfc_threshold &
                 padj <= RNAseq_padj_threshold),
         Induced = (expressed & 
                      !expressedCtrl &
                      log2FoldChange >= RNAseq_lfc_threshold &
                      padj <= RNAseq_padj_threshold),
         Repressed = (!expressed &
                        expressedCtrl &
                        log2FoldChange <= -RNAseq_lfc_threshold &
                        padj <= RNAseq_padj_threshold)
  ) %>% 
  mutate(DE = DE | Repressed) %>% 
  mutate(DECategory = case_when(Induced ~ "Induced",
                                Repressed ~ "Repressed",
                                !Induced & !Repressed & DE ~ "DE",
                                !Induced & !Repressed & !DE ~ "None")) %>% 
  mutate(DE = replace_na(DE, FALSE),
         Induced = replace_na(Induced, FALSE),
         Repressed = replace_na(Repressed, FALSE),
         DECategory = replace_na(DECategory, "None")
         )

# STEP 3: ADD ACCESSIBILITY
accFilt <- read.csv(ATAC_accFile, stringsAsFactors = FALSE)
deExprFiltAcc <-
  deExprFilt %>% 
  left_join(accFilt) %>% 
  mutate(accessible = replace_na(accessible, FALSE))

# STEP 4: ADD PROM COVERAGE
promCov <- read.csv(ATAC_promFile, stringsAsFactors = FALSE) %>% 
  pivot_wider(names_from = experimentalCondition, values_from = promNorm) %>%
  dplyr::rename(promCovCtrl = ctrl_0,
                promCovEGF48 = EGF_48,
                promCovIFNG48 = IFNG_48) %>% 
  group_by(ensembl_gene_id) %>%
  mutate(promCovTornado = mean(c(promCovCtrl, promCovEGF48, promCovIFNG48)))

deExprFiltAccProm <-
  deExprFiltAcc %>% 
  left_join(promCov)

# STEP 5: ADD HLA
deExprFiltAccPromHLA <-
  deExprFiltAccProm %>% 
  left_join(read.csv(HLA_filtFile, stringsAsFactors = FALSE)) %>% 
  mutate(HLA = replace_na(HLA, FALSE)) %>% 
  dplyr::rename(HLAregion = HLA)

write_csv(deExprFiltAccPromHLA, path = everythingOutFile)

# ## ------------------------------------------------------------------------
# inducedGenes_IFNG_48 <- 
#   RNAseq_resExpr %>% 
#   filter(geneType == "protein_coding") %>% 
#   filter(experimentalCondition == "IFNG_48") %>% 
#   filter(Induced) %>% 
#   pull(ensembl_gene_id)
# 
# repressedGenes_IFNG_48 <- 
#   RNAseq_resExpr %>% 
#   filter(geneType == "protein_coding") %>% 
#   filter(experimentalCondition == "IFNG_48") %>% 
#   filter(Repressed) %>% 
#   pull(ensembl_gene_id)
# 
# 
# ## ------------------------------------------------------------------------
# inducedGenesTableIFNG_48 <-
#   RNAseq_resExpr %>%
#   filter(geneType == "protein_coding") %>% 
#   filter(experimentalCondition == "IFNG_48") %>%
#   filter(Induced) %>%
#   dplyr::select(experimentalCondition, hgnc_symbol, HLA, everything())
# 
# repressedGenesTableIFNG_48 <-
#   RNAseq_resExpr %>%
#   filter(geneType == "protein_coding") %>% 
#   filter(experimentalCondition == "IFNG_48") %>%
#   filter(Repressed) %>%
#   dplyr::select(experimentalCondition, hgnc_symbol, HLA, everything())
# 
# # if(!dir.exists(outDirExcel)) {dir.create(outDirExcel)}
# 
# write_csv(inducedGenesTableIFNG_48, path = sprintf("%s/MDD_RNAseq_InducedGenes_IFNG48.csv", outDirExcel))
# write_csv(repressedGenesTableIFNG_48, path = sprintf("%s/MDD_RNAseq_RepressedGenes_IFNG48.csv", outDirExcel))

# write.xlsx(inducedGenesTable, 
#            file = sprintf("%s/MDD_RNAseq_InducedRepressedGenes_IFNG48.xlsx", outDirExcel),
#            sheetName = "Induced", row.names = FALSE)
# write.xlsx(repressedGenesTable, 
#            file = sprintf("%s/MDD_RNAseq_InducedRepressedGenes_IFNG48.xlsx", outDirExcel),
#            sheetName = "Repressed", row.names = FALSE,
#            append = TRUE)
