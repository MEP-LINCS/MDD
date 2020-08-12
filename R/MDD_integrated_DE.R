# Calculate RPPA log2 fold changes and differentially expressed proteins.
#
# The purpose of this script is to calculate log2 fold changes for all conditions
# compared to ctrl_0, using limma.
# This script is intended as the starting point for other scripts using RPPA
# L2FCs and or DE proteins

library(tidyverse)
library(limma)
library(ComplexHeatmap)

assays <- c("cycIF", "GCP", "motifs", "RNAseq", "RPPA")
#assay <- c("motifs")

get_assay_values <- function(assay){
  #browser()
  if(assay == "motifs"){
    level3File <- "../ATACseq/Data/MDD_ATACseq_MotifFamilyScores.csv"
  } else if(assay == "RNAseq"){
    level3File <- "../RNAseq/Data/MDD_RNAseq_DESeq2_Long.csv"
  } else if(assay == "RPPA"){
    level3File <- "../RPPA/Data/MDD_RPPA_limma_logFCTable.csv"
  }else {
    level3File <- paste0("../",assay,"/Data/MDD_",assay,"_Level3.csv")
  }
  
  MDDannoFile <- "../Metadata/MDD_sample_annotations.csv"
  
  if(assay == "RNAseq") {
    #passing through processed data
    mat <- read_csv(level3File,col_types = cols(
      experimentalCondition = col_character(),
      ensembl_gene_id = col_character(),
      hgnc_symbol = col_character(),
      baseMean = col_double(),
      log2FoldChange = col_double(),
      lfcSE = col_double(),
      pvalue = col_double(),
      padj = col_double()
    ))
    meta <- NULL
  } else if(assay == "RPPA"){
    mat <- read_csv(level3File,
                    col_types = cols(
                      antibody = col_character(),
                      experimentalCondition = col_character(),
                      coefficient = col_double(),
                      pvalue = col_double(),
                      adjusted_pvalue = col_double()
                    )) %>%
      filter(str_detect(experimentalCondition, "_24|_48"))
    
    meta <- NULL
  }
  
  else {
    mat_raw <- read.table(level3File,
                          header = TRUE,
                          sep = ",",
                          row.names = 1)
    row_sd <- apply(mat_raw, 1, sd, na.rm = TRUE)
    
    mat <- mat_raw[complete.cases(mat_raw) &!row_sd == 0,]
    
    meta <- read.table(MDDannoFile,
                       header = TRUE,
                       sep = ",") %>% 
      filter(str_detect(experimentalTimePoint, "0|24|48"),
             specimenID %in% colnames(mat)) %>% 
      dplyr::select(specimenID, specimenName, ligand, experimentalTimePoint, experimentalCondition, replicate) %>% 
      mutate(experimentalCondition = fct_inorder(factor(experimentalCondition)))
  }
  
  if(assay == "cycIF"){
    mat <- mat[!str_detect(rownames(mat),"_med_|dna.*cytoplasm|laws|nucrin|centcyto|plasmem|none"),]
    mat <- log2(mat+.0001)
  }
  
  if(!assay %in% c("RNAseq", "RPPA"))  mat <- mat %>%
    dplyr::select(all_of(unique(meta$specimenID)))
  
  return(list(mat = mat, meta = meta))
}

###############################################################################
# Creating design for comparisons.
# All conditions to be compared to ctrl_0.
#limma anlaysis expects log-expression values

combined_analysis <- lapply(assays, function(assay){
  #browser() 
  outDirPlots <- paste0("../plots/",assay,"_integrated")
  outDirData <- paste0("../",assay,"/Data/IntegratedResults")
  if(assay == "motifs") outDirData <- "../ATACseq/Data/IntegratedResults"
  
  if (!dir.exists(outDirPlots)) {
    dir.create(outDirPlots)
  }
  
  if (!dir.exists(outDirData)) {
    dir.create(outDirData)
  }
  
  assay_values <- get_assay_values(assay)
  
  if(assay == "motifs"){
    lfc_values <- assay_values$mat %>%
      rownames_to_column(var = "feature") %>%
      pivot_longer(cols = where(is.numeric), names_to = "specimenID") %>%
      left_join(assay_values$meta, by = "specimenID") %>%
      filter(str_detect(experimentalTimePoint, "24|48")) %>%
      group_by(feature, experimentalCondition) %>%
      summarise(value = mean(value), .groups = "drop") %>%
      pivot_wider(names_from = experimentalCondition, values_from = value)
    
    pdf(sprintf("%s/%s_int_values_hist.pdf", outDirPlots, assay), height = 12, width = 16)
    df <- lfc_values %>%
      as_tibble() %>%
      pivot_longer(cols = where(is.numeric))
    p <- ggplot(df, aes(x = value)) +
      geom_histogram(bins = 300) +
      #coord_cartesian(xlim = c(-6,6)) +
      labs(title = sprintf("%s deviation z scores histogram", assay)) +
      theme_bw()
    print(p)
    dev.off()
    
  } else if (assay == "RNAseq"){
    lfc_values <- assay_values$mat %>%
      filter(!is.na(hgnc_symbol)) %>%
      dplyr::select(experimentalCondition, hgnc_symbol, log2FoldChange) %>%
      pivot_wider(names_from = experimentalCondition, values_from = log2FoldChange) %>%
      drop_na() %>%
      rename(feature = hgnc_symbol)
    
    adjusted_p_values_df <- assay_values$mat %>%
      filter(hgnc_symbol %in% lfc_values$feature) %>%
      dplyr::select(experimentalCondition, hgnc_symbol, padj) %>%
      pivot_wider(names_from = experimentalCondition, values_from = padj) %>%
      #drop_na() %>%
      rename(feature = hgnc_symbol)
    
    lfc_mean <- lfc_values[,-1] %>%
      unlist() %>%
      mean
    
    lfc_sd <- lfc_values[,-1] %>%
      unlist() %>%
      sd
    
    lfc_values_z <- (lfc_values[,-1]-lfc_mean)/lfc_sd
  } else if (assay == "RPPA"){
    lfc_values <- assay_values$mat %>%
      dplyr::select(experimentalCondition, antibody, coefficient) %>%
      pivot_wider(names_from = experimentalCondition, values_from = coefficient) %>%
      drop_na() %>%
      rename(feature = antibody)
    
    adjusted_p_values_df <- assay_values$mat %>%
      dplyr::select(experimentalCondition, antibody, adjusted_pvalue) %>%
      pivot_wider(names_from = experimentalCondition, values_from = adjusted_pvalue) %>%
      #drop_na() %>%
      rename(feature = antibody)
    
    lfc_mean <- lfc_values[,-1] %>%
      unlist() %>%
      mean
    
    lfc_sd <- lfc_values[,-1] %>%
      unlist() %>%
      sd
    
    lfc_values_z <- (lfc_values[,-1]-lfc_mean)/lfc_sd
  } else {
    design <- model.matrix(~experimentalCondition, assay_values$meta)
    
    lm <- lmFit(assay_values$mat, design)
    lm <- eBayes(lm)
    adj_p_value_list <- lapply(colnames(lm$coefficients), function(condition){
      adj_p_values <- topTable(lm, n = Inf, coef = condition)["adj.P.Val"] %>%
        rownames_to_column("feature")
      return(adj_p_values)
    })
    adjusted_p_values_df <- bind_rows(adj_p_value_list, .id = "condition") %>%
      mutate(condition = colnames(lm$coefficients)[as.integer(condition)],
             condition = str_remove(condition, "experimentalCondition")) %>%
      dplyr::filter(!condition == '(Intercept)') %>%
      pivot_wider(names_from = condition, values_from = adj.P.Val)
    
    lfc_values <- lm[["coefficients"]] %>%
      data.frame() %>%
      rownames_to_column(var = "feature") %>%
      select(matches("experimentalCondition|feature")) %>%
      rename_with(~gsub("experimentalCondition", "", .x, fixed = TRUE))
    
    lfc_mean <- lfc_values[,-1] %>%
      unlist() %>%
      mean
    
    lfc_sd <- lfc_values[,-1] %>%
      unlist() %>%
      sd
    
    lfc_values_z <- (lfc_values[,-1]-lfc_mean)/lfc_sd
    
  }
  
  write.csv(lfc_values, sprintf("%s/%s_int_lfc_values.csv", outDirData, assay))
  #write.csv(lfc_values_z, sprintf("%s/%s_DE_lfc_z_scores.csv", outDirData, assay))
  if(!assay %in% c("motifs", "RNAseq", "RPPA")) {
    write_csv(adjusted_p_values_df, sprintf("%s/%s_int_adj_p_values.csv", outDirData, assay))
    
    show_row_names <- FALSE
    if(nrow(assay_values[['mat']]) < 100) show_row_names <- TRUE
    pdf(sprintf("%s/%s_input_values_hist.pdf", outDirPlots, assay), height = 12, width = 16)
    df <- assay_values[['mat']] %>%
      pivot_longer(cols = matches("sid"))
    p <- ggplot(df, aes(x = value)) +
      geom_histogram(bins = 300) +
      labs(title = sprintf("%s input values histogram", assay)) +
      theme_bw()
    print(p)
    dev.off()
    
    pdf(sprintf("%s/%s_int_values_z_hist.pdf", outDirPlots, assay), height = 12, width = 16)
    df <- lfc_values_z %>%
      as_tibble() %>%
      pivot_longer(cols = matches("_24|_48"))
    p <- ggplot(df, aes(x = value)) +
      geom_histogram(bins = 300) +
      coord_cartesian(xlim = c(-5,5)) +
      labs(title = sprintf("%s integrated Z scores histogram", assay)) +
      theme_bw()
    print(p)
    dev.off()
    
    
    pdf(sprintf("%s/%s_int_values_hist.pdf", outDirPlots, assay), height = 12, width = 16)
    df <- lfc_values %>%
      as_tibble() %>%
      pivot_longer(cols = matches("_24|_48"))
    p <- ggplot(df, aes(x = value)) +
      geom_histogram(bins = 300) +
      coord_cartesian(xlim = c(-6,6)) +
      labs(title = sprintf("%s integrated values histogram", assay)) +
      theme_bw()
    print(p)
    dev.off()
    
  } else {
    pdf(sprintf("%s/%s_int_values_hist.pdf", outDirPlots, assay), height = 12, width = 16)
    df <- lfc_values %>%
      as_tibble() %>%
      pivot_longer(cols = where(is.numeric))
    p <- ggplot(df, aes(x = value)) +
      geom_histogram(bins = 300) +
      #coord_cartesian(xlim = c(-6,6)) +
      labs(title = sprintf("%s input values histogram", assay)) +
      theme_bw()
    print(p)
    dev.off()
  }
  
})

