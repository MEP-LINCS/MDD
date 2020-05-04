# The purpose of this script is to run differential expression analyses
# for each condition vs ctrl_0. For each condition, the shrunken log2FoldChange
# estimates will then be used to rank genes. The top n (25-50) for each condition 
# will be selected. These lists will be combined (and any duplicates removed)
# to produce a list of genes for use in the "mega matrix" heatmap.

library(tidyverse)
library(DESeq2)

MDD_DESeqDataSetFile <- "../RNAseq/Data/MDD_RNAseq_DESeqDataSet.Rdata"
outDir <- "../RNAseq/Data"
if( str_extract(getwd(), "[:alnum:]+$") != "R" ) {setwd("R")}

##############################################################################

load(MDD_DESeqDataSetFile)

res <- list()

for (i in 2:15) {
  res[[i - 1]] <- lfcShrink(dds, coef = i, type = "apeglm")
}

names(res) <- str_extract(resultsNames(dds)[2:15], "[:alnum:]+_[248]{2}")

formatResults <- function(x, x_nm) {
  x <- 
    data.frame(x) %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    left_join(annoTable) %>% 
    filter(hgnc_symbol != "") %>% 
    filter(!duplicated(hgnc_symbol)) %>% 
    dplyr::select(hgnc_symbol, ensembl_gene_id, everything()) %>%
    mutate(experimentalCondition = x_nm) %>% 
    arrange(padj)
  x
}

res_RNA <- res
res_RNA_df  <- mapply(formatResults, 
                      x = res_RNA, x_nm = names(res_RNA), SIMPLIFY = FALSE)

save(res_RNA, res_RNA_df, 
     file = sprintf("%s/MDD_RNAseq_shrunkenResults.Rdata", outDir))
save(res_RNA_df, 
     file = sprintf("%s/MDD_RNAseq_shrunkenResults_df.Rdata", outDir))
