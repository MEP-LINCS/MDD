library(tidyverse)

hallmarkNESScoreFile <- "../RNAseq/misc/MDD_RNAseq_HallmarkNES.csv"
RNAscript <- "MDD_import_RNAseq_ensembl.R"
colorFile <- "MDD_importColors_pretty.R"
source(colorFile)
source(RNAscript)

outDirPlots <- "../plots/MDD_RNAseq_Hallmark_barplots"

###############################################################################
a <- read.csv(hallmarkNESScoreFile,
              stringsAsFactors = FALSE) %>% 
  dplyr::rename(hallmark = variable) %>% 
  filter(grepl("MYC", hallmark)) %>% 
  pivot_longer(-hallmark, names_to = "experimentalCondition") %>%
  left_join(sa.RNA.L3 %>% 
              distinct(experimentalCondition, .keep_all = TRUE) %>% 
              dplyr::select(experimentalCondition, Ligand)) %>% 
  mutate(experimentalCondition = as.factor(experimentalCondition)) %>% 
  mutate(experimentalCondition = fct_relevel(experimentalCondition,
                                             setdiff(levels(sa.RNA.L3$experimentalCondition), "CTRL_0"))
         ) %>% 
  arrange(experimentalCondition)

if(!dir.exists(outDirPlots)) {
  dir.create(outDirPlots)
}

pdf(sprintf("%s/MDD_RNAseq_Hallmark_MYC_TARGETS_V1_barplot.pdf", outDirPlots), height = 6, width = 6)
a %>% 
  filter(hallmark == "MYC_TARGETS_V1") %>% 
  ggplot(aes(experimentalCondition, value, fill = Ligand)) +
  geom_col() +
  scale_fill_manual(values = col$Ligand) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Normalized Enrichment Score")
dev.off()


  
pdf(sprintf("%s/MDD_RNAseq_Hallmark_MYC_TARGETS_V2_barplot.pdf", outDirPlots), height = 6, width = 6)
a %>% 
  filter(hallmark == "MYC_TARGETS_V2") %>% 
  ggplot(aes(experimentalCondition, value, fill = Ligand)) +
  geom_col() +
  scale_fill_manual(values = col$Ligand) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Normalized Enrichment Score")
dev.off()
