#Build the integrated analysis matrix

library(tidyverse)

assays <- c("GCP", "cycIF", "motifs", "RPPA", "RNAseq")

read_assay_values <- function(x, output_type = "lfc_z_scores"){
  #browser()
  if(x == "motifs"){
    data_type = "ATACseq"
  } else {
    data_type = x
  }
  foo <- paste0("../", data_type, "/Data/DEResults/",x,"_DE_",output_type,".csv") %>% 
    read_csv
}

output_type =  "lfc_z_scores"
T0_DE_integrated <- map(assays, read_assay_values, output_type) %>%
  bind_rows(.id = "Type") %>%
  mutate(Type = assays[as.integer(Type)]) %>%
  mutate(across(where(is.numeric), signif, digits = 4)) %>%
  rename(feature = X1) %>%
  rename_with( ~gsub("experimentalCondition", "", .x, fixed = TRUE))

write_csv(T0_DE_integrated, paste0("../integrated_analysis/integrated_T0_DE_",output_type,".csv"))

output_type =  "adj_p_values"
T0_DE_integrated <- map(assays, read_assay_values, output_type) %>%
  bind_rows(.id = "Type") %>%
  mutate(Type = assays[as.integer(Type)]) %>%
  mutate(across(where(is.numeric), signif, digits = 4))

write_csv(T0_DE_integrated, paste0("../integrated_analysis/integrated_T0_DE_",output_type,".csv"))

