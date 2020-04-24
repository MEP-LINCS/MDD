# The purpose of this script is to quickly plot and write tables
# of motif enrichments for specific TFs for a request from Laura

library(tidyverse)
library(cowplot)

theme_set(theme_cowplot())

source("R/MDD_importColors_default.R")

TFs <- c("FOXM1", "MYT1", "YBX1", "MYC", "POU5F1", "FOXO3", "STAT1", "SMAD3")

motifEnrichments <- read.csv("../mcf10a_common_project/mcf10a_atacseq/Data/motif/MDD_ATACseq_MotifScores.csv")
meAnno <- read.csv("../mcf10a_common_project/mcf10a_atacseq/Data/motif/MDD_ATACseq_MotifAnno.csv")
sampleAnno <- read.csv("../mcf10a_common_project/mcf10a_atacseq/Metadata/MDD_ATACseq_sampleMetadata.csv")
outDir <- "../plots/ATACseq_motifPlots_lauraEmail"

if (!dir.exists(outDir)) {dir.create(outDir)}

pdf(sprintf("%s/MDD_MotifEnrichment_FOXO3_MYC_OCT4.pdf", outDir), 
    height = 8, width = 12)
motifEnrichments %>% 
  left_join(meAnno) %>% 
  filter(hgnc_symbol %in% TFs) %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>% 
  dplyr::select(contains("sid"), hgnc_symbol) %>% 
  pivot_longer(-hgnc_symbol, names_to = "specimenID", values_to = "motifEnrichment") %>% 
  left_join(dplyr::select(sampleAnno, experimentalCondition,
                          ligand, experimentalTimePoint, specimenID, specimenName)) %>% 
  mutate(experimentalCondition = fct_inorder(as.factor(experimentalCondition)),
         ligand = fct_inorder(as.factor(ligand))) %>% 
  ggplot(aes(experimentalCondition, motifEnrichment)) +
  geom_boxplot(aes(color = ligand)) +
  geom_point(size = 2, shape = 21, aes(fill = ligand)) +
  facet_wrap(~hgnc_symbol, ncol = 2, scale = "free_y") +
  scale_color_manual(values = col$ligand) +
  scale_fill_manual(values = col$ligand) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

motifEnrichments %>% 
  left_join(meAnno) %>% 
  filter(hgnc_symbol %in% TFs) %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>% 
  dplyr::select(contains("sid"), hgnc_symbol) %>% 
  pivot_longer(-hgnc_symbol, names_to = "specimenID", values_to = "motifEnrichment") %>% 
  left_join(dplyr::select(sampleAnno,
                          specimenName, specimenID)) %>% 
  dplyr::select(-specimenID) %>% 
  mutate(motifEnrichment = round(motifEnrichment, digits = 5)) %>% 
  pivot_wider(names_from = specimenName, values_from = motifEnrichment) %>% 
  write_csv(path = sprintf("%s/MDD_MotifEnrichment_FOXO3_MYC_OCT4.csv", outDir))
