library(readxl)
suppressPackageStartupMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))

#
#' Read in and clean an RPPA data header
#'
#' Reads the first eight lines of the RPPA data file.
#'
#' @param fn The file path and name for an excel-formatted file
#' @return the header as a tibble
#'
#' @details The RPPA data is distributed in a multi-sheet excel file. Each sheet contains an eight row header
#'  of metadata for the antibodies. This function reads the header and
#'  formats it with the antibodies in rows and the metadata in columns.
#'
#' @export
#'
readRPPAHeader <- function(fn){
  header <- read_excel(fn,n_max = 8) %>%
    select(-starts_with("...")) %>%
    t() %>%
    as_tibble()
  colnames(header) <- header[1,]
  header <- header %>%
    slice(-1)
  header$ReplicateSet <- str_extract(fn, "_[[:alpha:]][.]") %>%
    str_extract("[[:alpha:]]")
  return(header)
}

#' Read in and clean RPPA data
#'
#' Read in RPPA data.
#'
#' @param fn The file path and name for an excel-formatted file
#' @return the file data as a tibble
#'
#' @details The RPPA data is distributed in a multi-sheet excel file. The antibody names
#' are in row 12 and the data for each antibody are in the rows below them. This function reads
#' the data and removes columns without names.
#'
#' @export
#'
readRPPAData <- function(fn, sheet="NormLinear"){
  df <- read_excel(fn,
                   sheet = sheet, col_names = FALSE,
                   skip = 11)
  colnames(df) <- read_excel(fn,
                             sheet = "NormLinear", col_names = FALSE,
                             skip = 9,n_max = 1)
  df <- df %>%
    select(-matches("^NA$"))
  return(df)
}

#' Calculate the CV of a numeric vector
#'
#' Calculate the coefficient of variation
#' @param x a numeric vector
#' @return the standard deviation divided by the mean
#' @export
#'
CV <- function(x){
  sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)
}


# The RPPA data is delivered in a multi-sheet excel file. This vignette starts by loading the header and data from the NormLinear sheet. 

# #temporary function for Ligand_2 metadata
#   label_ligand_2 <- function(x) {
#     ligand_2 <- rep("None", times = length(x))
#     ligand_2[x %in% c("BMP2", "IFNG", "TGFB")] <- "EGF"
#     return(ligand_2)
#   }
#   
l2Header <- lapply(dir("../RPPA/Data", pattern = "01_", full.names = TRUE), readRPPAHeader) %>%
  bind_rows() %>%
  spread(ReplicateSet,`Probabilities (QC Score)`) %>%
  mutate(SetAQCScore=signif(as.numeric(A),2),
         SetBQCScore=signif(as.numeric(B),2),
         SetCQCScore=signif(as.numeric(C),2)) %>%
  rename(Antibody=`Antibody Name in Heatmap`) %>%
  select(-A, -B, -C)

l2_data <- lapply(dir("../RPPA/Data", pattern = "01_", full.names = TRUE), readRPPAData) %>%
  bind_rows() %>%
  select(-matches("Order|Sample [(Source)(Name)(Type)]|Category_[23]|CF1&2 |^CF[12]$")) %>%
  gather(key = "antibody",value = "value",-Category_1, -`Sample description`) %>%
  mutate(antibody=str_remove(antibody,"-.-.$"),
         replicate=str_replace_all(Category_1,".*_",""),
         tempFeature=str_replace(`Sample description`,"null1","ctrl"),
         tempFeature=str_replace(tempFeature,"IFNg","IFNG"),
         tempFeature=str_replace(tempFeature,"TGFb","TGFB"),
         tempFeature=str_replace(tempFeature,"pbs","PBS"),
         ligand=str_extract(tempFeature,"ctrl|PBS|EGF|HGF|OSM|TGFB|IFNG|BMP2"),
         experimentalTimePoint=str_replace(tempFeature,".*_",""),
         experimentalTimePoint=str_replace(experimentalTimePoint,"^0",""),
         experimentalTimePoint=as.numeric(experimentalTimePoint),
         specimenName=paste(ligand, experimentalTimePoint, "C1", replicate, sep="_")) %>%
  select(specimenName, antibody, value) %>%
  mutate(value = log2(value)) 

mddMetadata <- read_csv("../Metadata/MDD_sample_annotations.csv")

l2 <- l2_data %>%
  left_join(mddMetadata, by = "specimenName")

l2_synapse<- l2 %>%
  select(specimenID, antibody, value) %>%
  spread(key = specimenID, value = value) %>%
  select(antibody, str_sort(colnames(.), numeric=TRUE)) %>%
  write_csv("../RPPA/Data/MDD_RPPA_level1.csv")



# Next we will add metadata about the antibodies.


RPPAMetadata <- read_csv("../RPPA/Metadata/RPPA_Metadata.csv",col_types = cols(
  MDD = col_character(),
  Symbols = col_character(),
  Sites = col_character(),
  Effect = col_character(),
  Pathway = col_character(),
  Protein = col_character(),
  Direction = col_integer(),
  RTK = col_logical(),
  IncycIF = col_logical(),
  Fig1 = col_logical(),
  Histone = col_logical(),
  CellCycle = col_logical()
)) %>%
  mutate(MDD = str_remove(MDD,"-.-.$"))



# The RPPA data is time course normalized by dividing each antibody by the corresponding time course values of EGF.  


egf <- l2 %>%
  filter(ligand=="EGF") %>%
  select(antibody,replicate,value,experimentalTimePoint) %>%
  rename(EGFValue=value)

l4 <- l2 %>%
  left_join(egf,by=c("antibody","replicate","experimentalTimePoint"))

l4 <- l4 %>%
  mutate(value = value - EGFValue,
         value = signif(value, 3))

l4_synapse<- l4 %>%
  select(specimenID, antibody, value) %>%
  spread(key = specimenID, value = value) %>%
  select(antibody, str_sort(colnames(.), numeric=TRUE),-sid1, -sid2, -sid3) %>%
  write_csv("../RPPA/Data/MDD_RPPA_level4.csv")




# The data is from three biological replicate experiments. Next, the replicates are summarized to the condition (ligand+time point) level.
 

###Data Descriptions
# **SampleDescription** - The experimental condition in the format ligand_timePoint\_replicate  
# **ProteinConcentration** - The protein concentration in the sample sent to MDACC  
# **Replicate** - The experiment includes triplicate samples gathered over three different weeks. `Replicate` values are only present in the Level 2 data. The Level 3 data is summarized across the replicates.  
# **CV** - Coefficient of variation of the replicates. Level 3 only.  
# **ValueVariance** - The variance of the time course values within a ligand treatment.  
# **Antibody** - The MDACC name of the antibody.  
# **Value** - The MDACC normalized linear value.  
# **Ligand** - The cells were incubated with one of six ligands or PBS. The time 0 samples have a ligand value of "ctrl".    
# **TimePoint** - Incubation time in hours. Values are 0, 1, 4, 8, 24 and 48.  
# **Ligand_TimePoint** - Concatenation of Ligand and TimePoint values with "_" separator.    
# **Symbols** - HUGO symbol for the protein associated with the antibody.   
# **Sites** - Phosphorylation site. Multiple sites are separated by a pipe symbol. If there are multiple symbols (genes) the phosphorylated sites for each will be separated by a space.  
# **Effect** - Used for CausalPath analysis    
# **Pathway** - Manually curated pathway assignment.    
# **Protein** - Redundant annotation to be depracated.   
# **Direction** - Used in pathway analysis.  
# **RTK** - Logical value for whether the antibody detects a Receptor Tyrosine Kinase.  
# **IncycIF** -  Logical value for whether the antibody overlaps with the cycIF antibodies.  
# **Fig1** -  Logical value for whether the antibody is included in figure 1. 
# **Antibody Origin** - MDACC annotation.  
# **Gene Name** - MDACC annotation.   
# **Validation Status** - MDACC annotation.  
# **Slide Id** - MDACC annotation.  
# **Slide no.** - MDACC annotation.   
# **Sample Type** - MDACC annotation.    
# **SetAQCScore** - MDACC QC value.  
# **SetBQCScore** - MDACC QC value.    
# **SetCQCScore** - MDACC QC value.   
# **T0Norm** - Time 0 normalized value calculated by dividing `Value` by the `Value` of `ctrl_0` of the same antibody in its replicate set.    
# **TCNorm**  - Time course normalized value calculated by dividing `Value` by either the PBS or EGF `Value` of the same antibody and same time point in its replicate set. BMP2, IFNG, and TGFB samples were incubated in the presence of EGF and are normalized against the EGF values. EGF, HGF and OSM samples were incubated in PBS and are normalized against their corresponding PBS values.  
# **PBSNorm**  - Time course normalized value calculated by dividing `Value` by the PBS `Value` of the same antibody and same time point in its replicate set.  
# **EGFNorm**  - Time course normalized value calculated by dividing `Value` by the egf `Value` of the same antibody and same time point in its replicate set.  