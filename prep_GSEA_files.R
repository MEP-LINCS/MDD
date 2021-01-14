library(tidyverse)
#load level 3 RNAseq data
l3 <- read_csv("RNAseq/Data/MDD_RNAseq_Level3_HGNC.csv") %>%
  select(specimenName, value, feature) %>%
  filter(str_detect(specimenName, "_[02]4?_"))

l3_wide <- l3 %>%
  pivot_wider(names_from = specimenName, values_from = value, values_fn = median) %>%
  mutate(DESCRIPTION = "") %>%
  rename(NAME = feature) %>%
  select(NAME, DESCRIPTION, everything())

#write out a GSEA compatible text file
res <- write_tsv(l3_wide, "tables/RNAseq_GSEA.txt")

#Prepare data for a CLS file
sample_names <- colnames(l3_wide)[str_detect(colnames(l3_wide),"_")] %>%
  str_remove("_C[12]_.*") %>%
  matrix(nrow = 1)

#write.table(sample_names,
            # file = "tables/RNAseq_GSEA.cls",
            # quote = FALSE,
            # sep = "\t",
            # row.names = FALSE,
            # col.names = FALSE)

#after running each ligand vs control, aggragate the results
fdr_thresh <- 0.2
gsea_reports <- dir("RNAseq/Data/GSEA/", full.names = TRUE)
gsea_report_tables <- map(gsea_reports,  read_tsv) %>%
  bind_rows(.id = "ligand") %>%
  janitor::clean_names() %>%
  mutate(ligand = gsea_reports[as.integer(ligand)],
         ligand = str_remove_all(ligand, ".*_for_|_[[:digit:]]*.tsv")) %>%
  select(ligand, name, size, es, nes, nom_p_val, fdr_q_val, fwer_p_val, rank_at_max) %>%
  filter(fdr_q_val <= fdr_thresh)

write_csv(gsea_report_tables, "tables/GSEA_report.csv")
