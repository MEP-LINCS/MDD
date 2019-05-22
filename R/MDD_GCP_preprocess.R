library(tidyverse)

#create EGF normalized level 4 data
#Reformat the file from the Broad to match the MDD format
#Select the histone names and sample identifiers then
#use the 20th row as the column names
GCP_level4_file <- read_csv("Data/GCP/MCF10a_GCP_data_with_OHSU_metadata_vocabulary.csv") %>%
  select(matches("histone|GMa")) %>%
  rename(histone = pr_gcp_histone_mark)
colnames(GCP_level4_file)[2:ncol(GCP_level4_file)] <- as.character(GCP_level4_file[20,2:ncol(GCP_level4_file)])

#Delete the first 23 rows since they contain unneeded metadata
GCP_l4 <-   slice(GCP_level4_file, 24:nrow(GCP_level4_file))

