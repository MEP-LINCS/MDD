library(scales)
library(RColorBrewer)
#library(MDDRPPA)
library(classInt)
library(tidyverse)


#Values used for out of bounds values
lowProb <- .02
upProb <- .98

#Filter threshold for Biomarker variances
varianceCut <- 0.04

#fold change threshold
fcThresh <- 1.5
z_score_thresh <- 1.5

#Create a function that takes in a vector and returns a same-length vector of z scores 
z_score <- function(x) (x-mean(x, na.rm = TRUE))/sqrt(var(x, na.rm = TRUE))


ProbeMetadata <- read_delim("../GCP/Metadata/ProbeMetadata.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

mddMetadata <- read_csv("../Metadata/MDD_sample_annotations.csv")

GCP_Data <- read_tsv("../GCP/Data/GCP_MCF10a_log2_noProbeMetadata.txt", skip = 2) %>%
  t() %>%
  as.tibble()
colnames(GCP_Data) <- GCP_Data[1,]

l2 <- GCP_Data %>%
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
  
  l2_synapse<- l2 %>%
  select(specimenID, histone, value) %>%
  spread(key = specimenID, value = value) %>%
  select(histone, str_sort(colnames(.), numeric=TRUE))

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
GCP_l4 <- slice(GCP_level4_file, 24:nrow(GCP_level4_file)) %>%
  gather(specimenName, value = value, -histone) %>%
  mutate(specimenName = str_replace(specimenName,"IFNg","IFNG"),
         specimenName = str_replace(specimenName,"TGFb","TGFB")) %>%
  inner_join(mddMetadata) %>%
  select(specimenID, histone, value) %>%
  mutate(value = as.numeric(value)) %>%
  spread(key = specimenID, value = value) %>%
  filter(!str_detect(histone,"[.(]"))

cycIF_Data <- read_csv("Data/biological_replicate_means.csv") %>%
  mutate(condition = str_remove(X1, "_[ABCDE]"),
         specimenName = paste0(condition, "_C3_",str_remove(X1,".*_"))) %>%
  select(-X1, -condition) %>%
  gather(key = feature, value = value, -specimenName) %>%
  inner_join(mddMetadata, by = "specimenName") %>%
  select(specimenID, feature, value) %>%
  spread(specimenID, value) %>%
  filter(!str_detect(feature, "background")) %>%
  write_csv("../cycIF/Data/MDD_cycIF_Level2.csv")

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
  write_csv("MDD_cycIF_Level4.csv")
