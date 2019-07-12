library(tidyverse)
library(sva)

mddMetadata <- read_csv("../Metadata/MDD_sample_annotations.csv")

cycIF_Data <- read_csv("../cycIF/Data/proper_FFC_biological_replicate_mean.csv") %>%
  mutate(condition = str_remove(X1, "_[ABCDE]"),
         specimenName = paste0(condition, "_C3_",str_remove(X1,".*_"))) %>%
  select(-X1, -condition) %>%
  gather(key = feature, value = value, -specimenName) %>%
  inner_join(mddMetadata, by = "specimenName") %>%
  select(specimenID, feature, value) %>%
  spread(specimenID, value) %>%
  filter(!str_detect(feature, "background")) 

batch <- cycIF_Data %>%
  gather(specimenID, value = "value", -feature) %>%
  left_join(mddMetadata, by="specimenID") %>%
  select(specimenID, replicate) %>%
  distinct %>%
  mutate(replicate = factor(replicate))
  
modcombat <- model.matrix(~replicate, data=batch)

cycIF_Data_combat <- ComBat(dat = select(cycIF_Data, -feature) %>% as.matrix,
                            batch = batch$replicate) %>%
  as.data.frame() %>%
  mutate(feature = cycIF_Data$feature) %>%
  write_csv("../cycIF/Data/MDD_cycIF_ComBat_Level3.csv")
