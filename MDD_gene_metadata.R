#MDD Integrated Analysis
#Mark Dane

library(tidyverse)
load("Data/selected_assay_pk_data_rr.rda")
#create a metadata file that annotates features with an associated gene

RPPA_genes <- read_csv("RPPA/Metadata/MDD_RPPA_antibodyAnnotations.csv") %>%
  select(MDD, Symbols) %>%
  rename(feature = MDD,
         symbol = Symbols) %>%
  mutate(feature = str_remove_all(feature,"-[RMG]-[VCQ]"),
         feature = paste0(feature,"_RPPA"))

cycIF_genes <- read_csv("cycIF/Data/cycIF_HUGO.csv") %>%
  rename(symbol = HUGOSymbol) %>%
  mutate(feature = str_replace(feature, "_.*","_cycIF"))

RNA_genes <- read_csv("RNAseq/Data/MDD_geneList_lfc15.csv") %>%
  rename(symbol = hgnc_symbol) %>%
  mutate(feature = paste0(symbol,"_RNA"))

ChEA3_genes <- TFs_selected_rr %>%
  select(feature) %>%
  mutate(symbol = str_remove(feature, "_TF")) %>%
  distinct()

feature_gene <- bind_rows(ChEA3_genes, cycIF_genes, RNA_genes, RPPA_genes) %>%
  distinct()

write_csv(feature_gene, "MDD_Integrated_feature_gene.csv")
                              