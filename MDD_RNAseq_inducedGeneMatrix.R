# The purpose of this script is to create tornado plots and heat maps for 
# induced and repressed genes
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(genomation)

# Scripts
RNAscript <- "MDD_import_RNAseq_ensembl.R"
heatmapFunctionScript <- "MDD_RNAseq_heatmapFunctions.R"

# Tables
combinedTable   <- "../RNAseq/misc/MDD_RNAseq_combinedTable.csv"

# Out paths
inducedGeneMatrixOut <- "../RNAseq/misc/MDD_RNAseq_inducedGeneMatrix.csv"

###############################################################################
# loading table of induced genes
cTable <- 
  read.csv(combinedTable, stringsAsFactors = FALSE)

# formatting
inducedGeneMatrix <-
  cTable %>% 
  group_by(ensembl_gene_id) %>% 
  filter(any(Induced)) %>% 
  dplyr::select(ensembl_gene_id, hgnc_symbol, experimentalCondition, Induced) %>% 
  mutate(Induced = case_when(Induced ~ 1,
                             !Induced ~ 0)) %>%
  distinct %>% 
  ungroup %>% 
  pivot_wider(names_from = experimentalCondition, values_from = Induced)

# writing file
write_csv(inducedGeneMatrix, inducedGeneMatrixOut)
