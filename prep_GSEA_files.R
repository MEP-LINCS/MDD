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

write.table(sample_names,
            file = "tables/RNAseq_GSEA.cls",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

#after running each lignad vs control, aggragte the results
