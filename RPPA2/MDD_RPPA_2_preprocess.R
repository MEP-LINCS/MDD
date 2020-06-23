library(readxl)
suppressPackageStartupMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))


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
l3_header <- dir("./RPPA2/Data", pattern = "xlsx", full.names = TRUE) %>%
  read_excel(n_max = 8) %>%
  select(-starts_with("...")) %>%
  t() %>%
  as_tibble()
colnames(l3_header) <- l3_header[1,]
l3_header <- janitor::clean_names(l3_header) %>%
  slice(-1) %>%
  mutate(probabilities_qc_score = as.numeric(probabilities_qc_score))

l3_data <- dir("./RPPA2/Data", pattern = "xlsx", full.names = TRUE) %>%
  read_excel(sheet =  "NormLinear", col_names = FALSE,
                   skip = 11)

colnames(l3_data) <-  dir("./RPPA2/Data", pattern = "xlsx", full.names = TRUE) %>%
  read_excel(sheet = "NormLinear", col_names = FALSE,skip = 9,n_max = 1) %>%
  slice(1)
colnames(l3_data) <- str_replace(colnames(l3_data), " ", "_")

write_csv(l3_header, "./RPPA2/Data/MDD_RPPA_2_header.csv")
write_csv(l3_data, "./RPPA2/Data/MDD_RPPA_2_data.csv")
