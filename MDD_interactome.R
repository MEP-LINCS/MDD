library(tidyverse)

genes_RNAseq <- read_csv(file = "RNAseq/Data/MDD_geneList_lfc15.csv") %>%
  mutate(RNAseq = TRUE)

metadata <- read_tsv("HGNC_to_Ensembl.txt") %>%
  rename(Symbol = "Approved symbol",
         Ensembl_ID = "Ensembl gene ID",
         Gene_group_name = "Gene group name") %>%
  select(Symbol, Ensembl_ID, Gene_group_name)
metadata$Gene_group_name[is.na(metadata$Gene_group_name)] = "none"

Huri <- read_tsv("RNAseq/Data/HuRI.tsv")  %>%
  left_join(metadata, by = c("Ensembl_gene_id_a" = "Ensembl_ID")) %>%
  rename(Source_symbol = Symbol,
         Source_gene_group_name = Gene_group_name) %>%
  select(Source_symbol, Source_gene_group_name, Ensembl_gene_id_b) %>%
  left_join(metadata, by = c("Ensembl_gene_id_b" = "Ensembl_ID")) %>%
  rename(Target_symbol = Symbol,
         Target_gene_group_name = Gene_group_name) %>%
  select(-Ensembl_gene_id_b) %>%
  filter(!is.na(Source_symbol),
         !is.na(Target_symbol))

huri_RNA <- Huri %>%
  left_join(genes_RNAseq, by = c(Source_symbol = "hgnc_symbol")) %>%
  rename(Source_RNAseq = RNAseq) %>%
  left_join(genes_RNAseq, by = c(Target_symbol = "hgnc_symbol")) %>%
  rename(Target_RNAseq = RNAseq) %>%
  mutate(interaction = "pp")

huri_RNA$Source_RNAseq[is.na(huri_RNA$Source_RNAseq)] = FALSE
huri_RNA$Target_RNAseq[is.na(huri_RNA$Target_RNAseq)] = FALSE
huri_RNA$Source_symbol[is.na(huri_RNA$Source_symbol)] = "none"
huri_RNA$Target_symbol[is.na(huri_RNA$Target_symbol)] = "none"

RNAseq_not_in_HuRI <- setdiff(genes_RNAseq$hgnc_symbol, Huri$`Approved symbol`)
targets_not_nodes <- setdiff(Huri$Ensembl_gene_id_b, Huri$Ensembl_gene_id_a)
write_csv(huri_RNA, "MDD_RNAseq_huri.csv")
#   
# huri_RPPA <- read_csv("RPPA/Metadata/MDD_RPPA_antibodyAnnotations.csv") %>%
#   mutate(RPPA = TRUE) %>%
#   right_join(Huri, by = c("Symbols" = "Approved symbol")) %>%
#   select(matches("gene|^Symbols$|RPPA")) %>%
#   mutate(interaction = "pp") %>%
#   rename(Source = "Ensembl_gene_id_a",
#          Target = "Ensembl_gene_id_b",
#          Group = 'Gene group name')
# huri_RPPA$RPPA[is.na(huri_RPPA$RPPA)] = FALSE
# 
# write_csv(huri_RPPA, "MDD_RPPA_huri.csv")

huri_ligand_clusters <- read_csv("MDD_condition.csv") %>%
  rename(Symbol = feature) %>%
  select(Symbol, Cluster) %>%
  mutate(Symbol = str_remove(Symbol, "_.*")) %>%
  right_join(huri_RNA, by = c("Symbol" = "Source_symbol"))
huri_ligand_clusters$Cluster[is.na(huri_ligand_clusters$Cluster)] = 0


write_csv(huri_ligand_clusters, "MDD_huri_ligand_clusters.csv")

huri_cl4 <- read_csv("cluster_4_genes.csv",col_names = FALSE) %>%
  rename(Symbols = X1) %>%
  mutate(cl4 = TRUE) %>%
  right_join(Huri, by = c("Symbols" = "Approved symbol")) %>%
  select(matches("gene|^Symbols$|cl4")) %>%
  mutate(interaction = "pp") %>%
  rename(Source = "Ensembl_gene_id_a",
         Target = "Ensembl_gene_id_b",
         Group = 'Gene group name') %>%
  filter(!Source == Target)
huri_cl4$cl4[is.na(huri_cl4$cl4)] = FALSE
huri_cl4$Group[is.na(huri_cl4$Group)] = "none"
huri_cl4$Symbols[is.na(huri_cl4$Symbols)] = "none"

write_csv(huri_cl4, "MDD_cl4_huri.csv")

huri_cl7 <- read_csv("cluster_7_genes.csv",col_names = FALSE) %>%
  rename(Symbols = X1) %>%
  mutate(cl7 = TRUE) %>%
  right_join(Huri, by = c("Symbols" = "Approved symbol")) %>%
  select(matches("gene|^Symbols$|cl7")) %>%
  mutate(interaction = "pp") %>%
  rename(Source = "Ensembl_gene_id_a",
         Target = "Ensembl_gene_id_b",
         Group = 'Gene group name') %>%
  filter(!Source == Target)
huri_cl7$cl7[is.na(huri_cl7$cl7)] = FALSE
huri_cl7$Group[is.na(huri_cl7$Group)] = "none"
huri_cl7$Symbols[is.na(huri_cl7$Symbols)] = "none"

write_csv(huri_cl7, "MDD_cl7_huri.csv")