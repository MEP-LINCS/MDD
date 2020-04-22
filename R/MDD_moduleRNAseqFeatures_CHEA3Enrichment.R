# This script uses the CHEA3 API to test for TF enrichment in the RNAseq features from
# the MCF10A integrated clusters.
# Adapted from script from Mark Dane.

library(tidyverse)
library(xlsx)

fDir      <- "../misc/ligand_tables"
outDirAll  <- "../RNAseq/Data/Enrichment/IntegratedCluster_CHEA3/All_Libraries"
outDirReMap <- "../RNAseq/Data/Enrichment/IntegratedCluster_CHEA3/ReMap"

# Functions
getTFs <- function(genes, qname){
  library(httr)
  library(jsonlite)
  
  url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
  encode = "json"
  payload = list(query_name = qname, gene_set = genes)
  
  #POST to ChEA3 server
  response = POST(url = url, body = payload, encode = encode)
  json = content(response, "text")
  
  #results as list of R dataframes
  results <- try(jsonlite::fromJSON(json))
  return(results)
}

writeListToExcel <- function(x, fname) {
  nm <- names(x)
  
  write.xlsx(x[[1]], fname, nm[1], row.names = FALSE)
  
  mapply(write.xlsx, 
         x = x[-1], 
         sheetName = nm[-1],
         file = fname,
         row.names = FALSE,
         append = TRUE)
}

###############################################################################
# Read cluster feature files 

if(!grepl("R$", getwd())) { setwd("R") }

myFiles        <- list.files(fDir, full.names = TRUE)
names(myFiles) <- str_extract(myFiles, "cluster_[0-9]+") %>% 
  str_replace("cluster", "module")

featureTables <- lapply(myFiles, read.csv, 
                        header = FALSE, stringsAsFactors = FALSE)

featureTables <- 
  lapply(featureTables, function(x) {
    colnames(x) <- c("feature", sprintf("V%s", seq(2, 10, 1)), "type", "symbol")
    x
  })

clusterRNAseqGenes <- lapply(featureTables, function(x) {
  x <- x %>% 
    filter(type == "RNA") %>%
    distinct(symbol) %>% 
    pull(symbol)
})

# submit each gene list to CHEA3 API
clusterCHEA3List <- mapply(getTFs, 
                           clusterRNAseqGenes, 
                           names(clusterRNAseqGenes), SIMPLIFY = FALSE)



# Subsetting to ReMap results
CHEA3_remap <- lapply(clusterCHEA3List, function(x) {
  x <- x[["ReMap--ChIP-seq"]]
  x
  })

ord <- order(as.numeric(str_extract(names(CHEA3_remap), "[0-9]+")))
CHEA3_remap <- CHEA3_remap[ord]
CHEA3_remap_csv <- Reduce(bind_rows, CHEA3_remap)

# Save results for all libraries
if(!dir.exists(outDirAll)) {
  dir.create(outDirAll)
}

mapply(writeListToExcel, 
       x = clusterCHEA3List, 
       fname = sprintf("%s/MDD_%s_CHEA3.xlsx", 
                       outDirReMap, names(clusterCHEA3List)))

# Save results for remap
if(!dir.exists(outDirReMap)) {
  dir.create(outDirReMap)
}

# writing excel file
write.xlsx(CHEA3_remap$cluster_1, 
           file = sprintf("%s/MDD_Clusters_CHEA3_Remap.xlsx", outDirReMap),
           sheetName = "Cluster_1", row.names = FALSE)
lapply(CHEA3_remap[-1], function(x) {
  write.xlsx(x, 
             file = sprintf("%s/MDD_Clusters_CHEA3_Remap.xlsx", outDirReMap),
             append = TRUE,
             sheetName = unique(str_replace(x$`Query Name`,"^m", "M")))
})

# writing csv file
write_csv(CHEA3_remap_csv, "../Data/IntegratedCluster_CHEA3/ReMap/MDD_Clusters_CHEA3_Remap.csv")

