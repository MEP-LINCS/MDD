#Build the integrated analysis matrix

library(tidyverse)

assays <- c("cycIF","GCP",  "motifs", "RNAseq", "RPPA")
read_assay_values <- function(x, output_type = "lfc_values"){
  #browser()
  if(x == "motifs"){
    data_type = "ATACseq"
  } else {
    data_type = x
  }
  foo <- paste0("../", data_type, "/Data/IntegratedResults/",x,"_int_",output_type,".csv") %>% 
    read_csv
}

output_type =  "lfc_values"
T0_DE_integrated <- map(assays, read_assay_values, output_type) %>%
  bind_rows(.id = "Type") %>%
  dplyr::select(-X1) %>%
  mutate(Type = assays[as.integer(Type)]) %>%
  mutate(across(where(is.numeric), signif, digits = 4))
  
write_csv(T0_DE_integrated, paste0("../integrated_analysis/integrated_matrix_",output_type,".csv"))

output_type =  "adj_p_values"
assays_with_p_values <- c("cycIF","GCP", "RNAseq", "RPPA")
T0_DE_integrated_p_values <- map(assays_with_p_values, read_assay_values, output_type) %>%
  bind_rows(.id = "Type") %>%
  mutate(Type = assays_with_p_values[as.integer(Type)]) %>%
  mutate(across(where(is.numeric), signif, digits = 4))

write_csv(T0_DE_integrated_p_values, paste0("../integrated_analysis/integrated_",output_type,".csv"))

