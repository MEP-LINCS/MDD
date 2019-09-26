library(tidyverse)
library(rrscale)

message("Starting rrscaling")
all_rr_values <- map(df_as_list, rrscaleDiagnostics)
