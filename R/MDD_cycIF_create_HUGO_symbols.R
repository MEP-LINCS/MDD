library(tidyverse)

cycIF_HUGO <- read_csv("cycIF_HUGO.csv")

l4_HUGO <- read_csv("../cycIF/Data/MDD_cycIF_Level4.csv") %>%
  select(feature) %>%
  mutate(cycIFName = str_remove(feature, "_.*")) %>%
  left_join(cycIF_HUGO) %>%
  select(-cycIFName)

write_csv(l4_HUGO, "cycIF_HUGO.csv")

