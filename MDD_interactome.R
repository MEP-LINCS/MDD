library(tidyverse)

genes_RNAseq <- read_csv(file = "RNAseq/Data/MDD_geneList_lfc15.csv") %>%
  mutate(RNAseq = TRUE)

metadata <- read_tsv("HGNC_to_Ensembl.txt")

Huri <- read_tsv("RNAseq/Data/HuRI.tsv")  %>%
  left_join(metadata, by = c("Ensembl_gene_id_a" = "Ensembl gene ID"))

huri_RNA <- Huri %>%
  left_join(genes_RNAseq, by = c("Approved symbol" = "hgnc_symbol")) %>%
  select(matches("gene|symbol$|RNAseq")) %>%
  mutate(interaction = "pp") %>%
  rename(Source = "Ensembl_gene_id_a",
         Target = "Ensembl_gene_id_b",
         Symbol = 'Approved symbol',
         Group = 'Gene group name') %>%
  filter(!Source == Target)

huri_RNA$RNAseq[is.na(huri_RNA$RNAseq)] = FALSE
huri_RNA$Group[is.na(huri_RNA$Group)] = "none"
huri_RNA$Symbol[is.na(huri_RNA$Symbol)] = "none"

write_csv(huri_RNA, "MDD_RNAseq_huri.csv")
  
huri_RPPA <- read_csv("RPPA/Metadata/MDD_RPPA_antibodyAnnotations.csv") %>%
  mutate(RPPA = TRUE) %>%
  right_join(Huri, by = c("Symbols" = "Approved symbol")) %>%
  select(matches("gene|^Symbols$|RPPA")) %>%
  mutate(interaction = "pp") %>%
  rename(Source = "Ensembl_gene_id_a",
         Target = "Ensembl_gene_id_b",
         Group = 'Gene group name')
huri_RPPA$RPPA[is.na(huri_RPPA$RPPA)] = FALSE

write_csv(huri_RPPA, "MDD_RPPA_huri.csv")

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