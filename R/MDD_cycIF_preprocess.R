library(tidyverse)

mddMetadata <- read_csv("../Metadata/MDD_sample_annotations.csv")

cycIF_Data <- read_csv("../cycIF/Data/biological_replicate_means.csv") %>%
  mutate(condition = str_remove(X1, "_[ABCDE]"),
         specimenName = paste0(condition, "_C3_",str_remove(X1,".*_"))) %>%
  select(-X1, -condition) %>%
  gather(key = feature, value = value, -specimenName) %>%
  inner_join(mddMetadata, by = "specimenName") %>%
  select(specimenID, feature, value) %>%
  spread(specimenID, value) %>%
  filter(!str_detect(feature, "background")) %>%
  write_csv("../cycIF/Data/MDD_cycIF_Level3.csv")

#get EGF values for each feature within each time point and replicate
egf <- cycIF_Data %>%
  gather(specimenID, value, -feature) %>%
  inner_join(mddMetadata, by = "specimenID") %>%
  select(specimenName, feature, value, ligand, experimentalTimePoint, replicate) %>%
  filter(str_detect(ligand, "EGF")) %>%
  select(feature, value, experimentalTimePoint, replicate) %>%
  rename(EGFValue = value)

l4 <- cycIF_Data %>%
  gather(specimenID, value, -feature) %>%
  inner_join(mddMetadata, by = "specimenID") %>%
  select(specimenName, feature, value, ligand, experimentalTimePoint, replicate) %>%
  left_join(egf,by=c("feature","replicate","experimentalTimePoint")) %>%
  filter(!str_detect(ligand, "ctrl")) %>%
  mutate(value = value/EGFValue,
         value = signif(value, 3)) %>%
  select(specimenName, feature, value) %>%
  inner_join(mddMetadata, by = "specimenName") %>%
  select(specimenID, feature, value) %>%
  spread(key = specimenID, value = value) %>%
  write_csv("../cycIF/Data/MDD_cycIF_Level4.csv")

