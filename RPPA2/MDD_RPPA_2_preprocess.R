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

#Loading the header and data from the raaw data NormLinear sheet. 
l3Header <- lapply(dir("./RPPA2/Data", pattern = "xlsx", full.names = TRUE), readRPPAHeader) %>%
  bind_rows() %>%
  spread(ReplicateSet,`Probabilities (QC Score)`) %>%
  mutate(SetAQCScore=signif(as.numeric(A),2),
         SetBQCScore=signif(as.numeric(B),2),
         SetCQCScore=signif(as.numeric(C),2)) %>%
  rename(Antibody=`Antibody Name in Heatmap`) %>%
  select(-A, -B, -C)

l3_data <- lapply(dir("../RPPA/Data", pattern = "01_", full.names = TRUE), readRPPAData) %>%
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

l3 <- l3_data %>%
  left_join(mddMetadata, by = "specimenName")

l3_synapse<- l3 %>%
  select(specimenID, antibody, value) %>%
  spread(key = specimenID, value = value) %>%
  select(antibody, str_sort(colnames(.), numeric=TRUE)) %>%
  write_csv("../RPPA/Data/MDD_RPPA_Level3.csv")

#median summarize replicates by condition
l4 <- l3 %>% 
  select(experimentalCondition, value, antibody) %>%
  group_by(experimentalCondition, antibody) %>%
  summarise(value = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  spread(key = experimentalCondition, value = value) %>%
  write_csv("../RPPA/Data/MDD_RPPA_Level4.csv")
