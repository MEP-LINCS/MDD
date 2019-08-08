library(scales)
library(RColorBrewer)
library(classInt)
library(tidyverse)

mddMetadata <- read_csv("../Metadata/MDD_sample_annotations.csv")

GCP_Data <- read_tsv("../GCP/Data/GCP_MCF10a_log2_noProbeMetadata.txt", skip = 2) %>%
  t() %>%
  as_tibble()
colnames(GCP_Data) <- GCP_Data[1,]
GCP_Data$pert_iname[GCP_Data$pert_time==0] <- "ctrl"


l3_long <- GCP_Data %>%
  slice(-1) %>%
  filter(!is.na(id)) %>%
  mutate(ligand = pert_iname,
         ligand = str_replace(ligand,"IFNg","IFNG"),
         ligand = str_replace(ligand,"TGFb","TGFB"),
         replicate = as.integer(pert_batch_internal_replicate),
         replicate = LETTERS[replicate],
         collection = "C2",
         experimentalTimePoint = pert_time,
         specimenName=paste(ligand, experimentalTimePoint, collection, replicate, sep="_")) %>%
 select(specimenName, matches("^H3K|^H4")) %>%
  mutate_at(vars(matches("^H3K|H4")), as.numeric) %>%
  inner_join(mddMetadata, by = "specimenName") %>%
  select(specimenID, matches("^H3K|H4")) %>%
  gather(key = histone, value = value, -specimenID)

l3 <- l3_long %>%
  select(specimenID, histone, value) %>%
  spread(key = specimenID, value = value) %>%
  mutate(histone = factor(histone, levels = unique(l3_long$histone))) %>%
  arrange(histone) %>%
  write_csv("../GCP/Data/MDD_GCP_Level3.csv")

#median summarize replicates by condition
l4 <- l3 %>% 
  gather(specimenID, value, -histone) %>%
  left_join(mddMetadata, by = "specimenID") %>%
  select(experimentalCondition, value, histone) %>%
  group_by(experimentalCondition, histone) %>%
  summarise(value = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  spread(key = experimentalCondition, value = value) %>%
  write_csv("../GCP/Data/MDD_GCP_Level4.csv")
