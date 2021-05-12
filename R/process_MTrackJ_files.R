library(tidyverse)
library(here)
library(rmarkdown)

#assume there is a subfolder named Data that contains data files
#with plateID_well alphanumeric_field_intensities.csv files

process_files <- function(path, plateID = NULL){
  out_filename <- str_replace(path,"Points","intensities")
  
  if(!dir.exists(paste0(here(),"/MTrackJ/",plateID, "/Reports"))) dir.create(paste0(here(),"/MTrackJ/",plateID, "/Reports"))
  
  report_name <- str_replace(path, "Data","Reports") %>%
    str_replace("_Points.csv",".html")
  
  render("MTrackJ/Manual_tracking_EDA.Rmd",
         output_file = report_name,
         output_format = "html_document",
         output_options = list(code_folding = "hide"))
}
 
#plateID <- "LI204701"
plateID <- "LI802303"

res <- dir(paste0("MTrackJ/",plateID,"/Data"), pattern = "Points.csv", full.names = TRUE) %>%
  here() %>%
  map(process_files, plateID = plateID)

render_summary_report <- function(plateID){
  report_name <- paste0(here(),"/MTrackJ/",plateID, "_tracking_summary.html")
  render("MTrackJ/Manual_tracking_summary.Rmd",
         output_file = report_name,
         output_format = "html_document",
         output_options = list(code_folding = "hide"))
}
render_summary_report(plateID)
