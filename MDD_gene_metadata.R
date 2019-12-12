#MDD Integrated Analysis
#Mark Dane

library(tidyverse)
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
         
# cycIF_genes <- read_csv("cycIF/Data/MDD_cycIF_Level3.csv") %>%
#   select(feature) %>%
#   mutate(feature = str_remove(feature, "_.*"),
#          symbol = str_replace(feature, "cateninbeta", "CTNNB1"),
#          symbol = str_replace(symbol, "cjun","JUN"),
#          symbol = str_replace(symbol, "cyclind1", "CCND1"),
#          symbol = str_replace(symbol, "cytokeratin18", "KRT18"),
#          symbol = str_replace(symbol, "cytokeratin7human", "KRT7"),
#          symbol = str_replace(symbol, "egfr", "EGFR"),
#          symbol = str_replace(symbol, "hes1", "HES1"),
#          symbol = str_replace(symbol, "ki67", "MKI67"),
#          symbol = str_replace(symbol, "lc3ab", "MAP1LC3A"),
#          symbol = str_replace(symbol, "met", "MET"),
#          symbol = str_replace(symbol, "ndg1pt346", "NDRG1"),
#          symbol = str_replace(symbol, "nfkbp65", "NFKB1"),
#          symbol = str_replace(symbol, "p21waf1cip1", "CDKN1A"),
#          symbol = str_replace(symbol, "pdl1", "CD274"),
#          symbol = str_replace(symbol, "s6", "RPS6"),
#          symbol = str_replace(symbol, "s6ps235s236", "RPS6"),
#          symbol = str_replace(symbol, "s6ps240244", "RPS6"),
#          symbol = str_replace(symbol, "stat1alphaisoform", "STAT1"),
#          symbol = str_replace(symbol, "stat1ps727", "STAT1"),
#          symbol = str_replace(symbol, "stat3", "STAT3"),
#          symbol = str_replace(symbol, "vimentin", "VIM"),
#          feature = paste0(feature, "_cycIF")) %>%
#   distinct()

RNA_genes <- read_csv("RNAseq/Data/MDD_geneList_lfc15.csv") %>%
  rename(symbol = hgnc_symbol) %>%
  mutate(feature = paste0(symbol,"_RNA"))

feature_gene <- bind_rows(RPPA_genes, cycIF_genes, RNA_genes)

foo <- selected_assay_pk_data_rr %>%
  left_join(feature_gene, by = "feature")
