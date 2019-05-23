library(scales)
library(RColorBrewer)
library(classInt)
library(tidyverse)

mddMetadata <- read_csv("../Metadata/MDD_sample_annotations.csv")

GCP_Data <- read_tsv("../GCP/Data/GCP_MCF10a_log2_noProbeMetadata.txt", skip = 2) %>%
  t() %>%
  as_tibble()
colnames(GCP_Data) <- GCP_Data[1,]

l2_long <- GCP_Data %>%
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
  select(specimenName, matches("^H3K")) %>%
  mutate_at(vars(matches("^H3K")), as.numeric) %>%
  inner_join(mddMetadata, by = "specimenName") %>%
  select(specimenID, matches("^H3K")) %>%
  gather(key = histone, value = value, -specimenID)

l2 <- l2_long %>%
  select(specimenID, histone, value) %>%
  spread(key = specimenID, value = value) %>%
  select(histone, str_sort(colnames(.), numeric=TRUE)) %>%
  write_csv("../GCP/Data/MDD_GCP_Level2.csv")

#create EGF normalized level 4 data
#Reformat the file from the Broad to match the MDD format
#Select the histone names and sample identifiers then
#use the 20th row as the column names
GCP_level4_file <- read_tsv("../GCP/Data/MCF10a_EGFtimepointnorm_RepRemoved_withPBS_OHSU_metadata.txt", skip = 2) 
GCP_level4_file[(GCP_level4_file$id == "sampleid"),] <- paste(GCP_level4_file[(GCP_level4_file$id == "condition"),],
                                                              GCP_level4_file[(GCP_level4_file$id == "collection"),],
                                                              GCP_level4_file[(GCP_level4_file$id == "replicate"),],
                                                              sep="_")
GCP_level4_file <- GCP_level4_file %>%
  select(matches("histone|GMa")) %>%
  rename(histone = pr_gcp_histone_mark)

colnames(GCP_level4_file)[2:ncol(GCP_level4_file)] <- GCP_level4_file[20,2:ncol(GCP_level4_file)] 

#Delete the first 23 rows since they contain unneeded metadata
#match with MDD metadata file to use specimenID values
l4 <- slice(GCP_level4_file, 24:nrow(GCP_level4_file)) %>%
  gather(specimenName, value = value, -histone) %>%
  mutate(specimenName = str_replace(specimenName,"IFNg","IFNG"),
         specimenName = str_replace(specimenName,"TGFb","TGFB")) %>%
  inner_join(mddMetadata) %>%
  select(specimenID, histone, value) %>%
  mutate(value = as.numeric(value)) %>%
  spread(key = specimenID, value = value) %>%
  filter(!str_detect(histone,"[.(]")) %>%
  write_csv("../GCP/Data/MDD_GCP_Level4.csv")

