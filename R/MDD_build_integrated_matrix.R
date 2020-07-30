#Build the integrated analysis matrix

library(tidyverse)

assays <- c("GCP", "cycIF", "motifs")

read_tt <- function(x){
  
  if(x == "motifs"){
    data_type = "ATACseq"
  } else {
    data_type = x
  }
  foo <- paste0("../", data_type, "/Data/DEResults/",x,"_DE_topTables.csv") %>% 
    read_csv
}
T0_DE_integrated <- map(assays, read_tt) %>%
  bind_rows(.id = "Type") %>%
  mutate(Type = assays[as.integer(Type)]) %>%
  rename(feature = X1) %>%
  rename_with( ~gsub("experimentalCondition", "", .x, fixed = TRUE))

write_csv(T0_DE_integrated, "../integrated_analysis/T0_DE_integrated.csv")
