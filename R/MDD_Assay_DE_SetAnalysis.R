#import libraries 
library(readr)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(devtools)
library(ggpubr)
library(RColorBrewer)
library(rlang)
library(circlize)
`%notin%` <- Negate(`%in%`)

#set global parameters
ligand_cols <- c(
  "HGF" = "#80b1d3",
  "OSM" = "#fdb462",
  "EGF" = "#fb8072",
  "BMP2" = "#b3de69",
  "IFNG" = "#bebada",
  "TGFB" = "#ffd92f")
ligand_cols_shared <- c(
  "HGF" = "#80b1d3",
  "OSM" = "#fdb462",
  "EGF" = "#fb8072",
  "BMP2" = "#b3de69",
  "IFNG" = "#bebada",
  "TGFB" = "#ffd92f",
  "Shared" = "#808080")

#define custom functions to perform set analysis

#prep dataframes for set analysis
Set_Analysis_Prep <- function(DE_pval_path, DE_lfc_path, assay, pvalthresh = .05, lfcthresh = 1.5, SetType = 'Anyset', PBS=T){
  
  #import files
  t0_de_p<-read_csv(DE_pval_path)
  t0_de_lfc<-read_csv(DE_lfc_path)
  
  #format data files
  t0_de_p_long<-t0_de_p %>%
    gather(key = experimentalCondition, value=padj, -Type, -feature)
  
  t0_de_lfc_long<-t0_de_lfc %>%
    gather(key = experimentalCondition, value=log2FoldChange, -Type, -feature)
  
  #bind data together
  t0_de <- left_join(t0_de_lfc_long, t0_de_p_long)
  
  #filter for assay specific features
  t0_de_assay<-t0_de %>%
    filter(Type==assay)
  if(PBS==F){
    t0_de_assay<-t0_de %>%
      filter(Type==assay) %>%
      filter(experimentalCondition != 'PBS_24' & experimentalCondition != 'PBS_48')
  }
  #filter for p val and lfc threshold
  if(SetType=='Anyset'){ if(assay=='motifs'){
    t0_filtered<-t0_de_assay %>%
      filter(log2FoldChange > lfcthresh | log2FoldChange < (-lfcthresh))
  } else
    t0_filtered<-t0_de_assay %>%
      filter(padj < pvalthresh) %>%
      filter(log2FoldChange > lfcthresh | log2FoldChange < (-lfcthresh))
  } else if (SetType=='Upset'){
    t0_filtered<-t0_de_assay %>%
      filter(padj < pvalthresh) %>%
      filter(log2FoldChange > lfcthresh)
  } else if (SetType=='Downset'){
    t0_filtered<-t0_de_assay %>%
      filter(padj < pvalthresh) %>%
      filter(log2FoldChange < (-lfcthresh))
  } else {
    stop('Invalid SetType')
  }
  
  #Transform dataframe into long dataframe, indicating regulation by the time/condition
  t0_de_long <- t0_filtered %>%
    select(experimentalCondition, feature, log2FoldChange) %>%
    spread(key=feature, value = log2FoldChange) %>%
    mutate(ligand=gsub('_.*', '' , experimentalCondition)) %>%
    mutate(experimentalTimepoint = gsub('.*_', '' , experimentalCondition)) 
  
  
  return(t0_de_long)
  
}

#function to perform ligand dependent set analysis
Set_Analysis_Ligand_Unique<-function(prepped_DE_data, plot=F, unique=T, assay, SetType='Anyset'){
  
  #find features only regulated in one ligand condition
  ligand_summarized<-prepped_DE_data %>%
    group_by(ligand) %>%
    summarise_if(is.numeric, ~ median(., na.rm=TRUE))
  
  if(unique==T){
    ligand_unique_features<-data.frame(ligand_count=colSums(!is.na(ligand_summarized[-c(1)]))) %>%
      rownames_to_column('Feature') %>%
      filter(ligand_count<2)
  } else if(unique==F){
    ligand_unique_features<-data.frame(ligand_count=colSums(!is.na(ligand_summarized[-c(1)]))) %>%
      rownames_to_column('Feature') %>%
      filter(ligand_count>1)
  }
  
  #make dataframe of uniquely regulated features organized by ligand
  ligand_unique_frame <- ligand_summarized %>%
    select(ligand, ligand_unique_features$Feature) %>%
    gather(key=feature, value = 'Log2FoldChange', -ligand, na.rm = TRUE) %>%
    arrange(ligand, desc(Log2FoldChange))
  
  if(plot==T){
    
    t0_data<-ligand_summarized %>%
      gather(key=feature, value = 'Log2FoldChange', -ligand, na.rm = TRUE) %>%
      arrange(ligand, desc(Log2FoldChange))
    
    ligand_list<-unique(t0_data$ligand)
    
    Upset_list<-list()
    for(i in 1:length(ligand_list)){
      ligand_cur=ligand_list[i]
      
      ligand_dat <- t0_data %>%
        filter(ligand==ligand_cur)
      
      tmp<-list(ligand_dat$feature)
      Upset_list[ligand_cur]<-tmp
    }
    
    upset_plot<-make_comb_mat(Upset_list)
    ligand_unique_frame<- UpSet(upset_plot, column_title = paste0(assay,' ',SetType, ' Plot by Ligand'), row_names_gp = grid::gpar(fontsize = 20),comb_order = rev(order(comb_size(upset_plot))))
    
    
  }
  
  return(ligand_unique_frame)
  
}

#heatmap of unique features z-scored
unique_heatmap<-function(unique_features, condition_level, condition, assay, SetType = 'Anyset'){
  
  unique_set<-unique_features %>%
    filter(!!sym(condition_level) == condition)
  
  #import level 4 data
  lvl4<-read_csv(paste0('../',assay,'Data/MDD_',assay,'_Level4.csv'))
  
  if(assay=='RNAseq'){
    gene_md<-read_csv('../',assay,'Metadata/MDD_RNAseq_geneAnnotations.csv')
    
    
    RNAseq_md<-full_join(gene_md,lvl4) %>%
      select(-ensembl_gene_id) %>%
      rename(feature=hgnc_symbol)
    
    lvl4<-RNAseq_md
  } else if(assay=='RPPA'){
    lvl4<-lvl4 %>%
      rename(feature=antibody)
  } else if(assay=='GCP'){
    lvl4<-lvl4 %>%
      rename(feature=histone)
  } else if(assay=='motifs'){
    lvl4<-lvl4 %>%
      rename(feature=family)}
  
  lvl4_unique<-lvl4 %>%
    filter(feature %in% unique_set$feature) 
  
  unique_mat <-  lvl4_unique %>%
    select(-feature) %>%
    as.matrix() %>%
    t %>%
    scale %>%
    t %>%
    as.data.frame() %>%
    drop_na() %>%
    as.matrix()
  
  
  Heatmap(unique_mat, name = "Z score",cluster_columns = T, cluster_rows = T,
          row_gap = unit(2, "mm"), col=colorRamp2(c(-4, 0, 4), c("blue", "white", "red")), column_title = paste0(condition," Uniquely ",gsub('set','',SetType)  ," Regulated Features from ",condition_level,' level ', assay), width = unit(20, "cm"), row_title = NULL , show_row_names = FALSE, heatmap_legend_param = list(heatmap_legend_side = 'top', at = c(-4,-2, 0,2, 4)) 
  )
}

#heatmap of shared features z-scored
shared_heatmap<-function(shared_features, condition_level, assay, SetType='Anyset'){
  
  #import level 4 data
  lvl4<-read_csv(paste0('../',assay,'Data/MDD_',assay,'_Level4.csv'))
  
  if(assay=='RNAseq'){
    gene_md<-read_csv('../',assay,'Metadata/MDD_RNAseq_geneAnnotations.csv')
    
    RNAseq_md<-full_join(gene_md,lvl4) %>%
      select(-ensembl_gene_id) %>%
      rename(feature=hgnc_symbol)
    
    lvl4<-RNAseq_md
  } else if(assay=='RPPA'){
    lvl4<-lvl4 %>%
      rename(feature=antibody)
  } else if(assay=='motifs'){
    lvl4<-lvl4 %>%
      rename(feature=family)} else if(assay=='GCP'){
        lvl4<-lvl4 %>%
          rename(feature=histone)
      }
  
  lvl4_shared<-lvl4 %>%
    filter(feature %in% shared_features$feature) 
  
  shared_mat <-  lvl4_shared %>%
    select(-feature) %>%
    as.matrix() %>%
    t %>%
    scale %>%
    t
  
  
  rownames(shared_mat) <- lvl4_shared$feature
  
  Heatmap(shared_mat, name = "Z score",cluster_columns = T, cluster_rows = T, col=colorRamp2(c(-4, 0, 4), c("blue", "white", "red")), 
          row_gap = unit(2, "mm"), column_title = paste0("Shared ",gsub('set','',SetType) , " Regulated Features from ",condition_level,' level ', assay), width = unit(20, "cm"), row_title = NULL , show_row_names = FALSE, heatmap_legend_param = list(heatmap_legend_side = 'top', at = c(-4,-2, 0,2, 4)) 
  )
  
  
}

#heatmap of unique features by LFC
unique_heatmap_lfc<-function(unique_features, condition_level, condition, assay, SetType = 'Anyset', DE_lfc_path){
  
  
  unique_set<-unique_features %>%
    filter(!!sym(condition_level) == condition)
  
  #import lfc data
  lfc<-read_csv(DE_lfc_path)
  
  if(assay=='Integrated'){
    lfc_assay<-lfc
  }else{
    lfc_assay<-lfc %>%
      filter(Type==assay)
  }
  lfc_unique<-lfc_assay %>%
    filter(feature %in% unique_set$feature) 
  
  unique_mat <-  lfc_unique %>%
    select(-feature, -Type) %>%
    as.matrix()
  
  
  rownames(unique_mat) <- lfc_unique$feature
  
  Heatmap(unique_mat, name = "LFC",cluster_columns = F, cluster_rows = T, col=colorRamp2(c(-4, 0, 4), c("blue", "white", "red")), 
          row_gap = unit(2, "mm"), column_title = paste0(condition," Uniquely ",gsub('set','',SetType)  ," Regulated Features from ",condition_level,' level ', assay), width = unit(20, "cm"), row_title = NULL , show_row_names = FALSE, heatmap_legend_param = list(heatmap_legend_side = 'top', at = c(-4,-2, 0,2, 4)) 
  )
  
  
}

#heatmap of shared features by LFC
shared_heatmap_lfc<-function(shared_features, condition_level, assay, SetType='Anyset', DE_lfc_path){
  
  #import lfc data
  lfc<-read_csv(DE_lfc_path)
  
  if(assay=='Integrated'){
    lfc_assay<-lfc
  }else{
    lfc_assay<-lfc %>%
      filter(Type==assay)
  }
  lfc_shared<-lfc_assay %>%
    filter(feature %in% shared_features$feature) 
  
  shared_mat <-  lfc_shared %>%
    select(-feature, -Type) %>%
    as.matrix() 
  
  
  
  rownames(shared_mat) <- lfc_shared$feature
  
  Heatmap(shared_mat, name = "LFC",cluster_columns = F, cluster_rows = T, col=colorRamp2(c(-4, 0, 4), c("blue", "white", "red")), 
          row_gap = unit(2, "mm"), column_title = paste0("Shared ",gsub('set','',SetType) , " Regulated Features from ",condition_level,' level ', assay), width = unit(20, "cm"), row_title = NULL , show_row_names = FALSE, heatmap_legend_param = list(heatmap_legend_side = 'top', at = c(-4,-2, 0,2, 4)) 
  )
  
  
}

#bar graphs of feature counts
unique_bar<-function(unique_features, assay){
  
  ligand_counts<-unique_features %>%
    group_by(ligand) %>%
    summarize(count=n())
  
  ggplot(ligand_counts, aes(x=ligand, y=count, fill=ligand)) +
    geom_bar(stat='identity') + theme_bw() + labs(title = paste0(assay, ' unique features')) + theme(text=element_text(size=18),panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                     panel.grid.minor = element_blank(), axis.line = element_blank()) + scale_fill_manual(values=ligand_cols)
}

#function to perform set analysis on each assay dataset
perform_set_function <- function(DE_pval_path, DE_lfc_path, assay=c('cycIF', 'RPPA', 'RNAseq', 'motifs', 'GCP'), SetType = 'Anyset', PBS=F){
  
  #Perform set analysis on each assay dataset
  total_unique_lista<-data.frame()
  total_shared_lista<-data.frame()
  
  for(i in 1:length(assay)){
    cur_assay <- assay[i]
    
    if(cur_assay == 'motifs'){
      lfcthresh = 1.5
      pvalthresh = 0
    } else {
      lfcthresh = 1.5
      pvalthresh = .05
    }
    
    t0_anyset <- Set_Analysis_Prep(DE_pval_path = DE_pval_path, DE_lfc_path = DE_lfc_path, assay=cur_assay, pvalthresh = pvalthresh, lfcthresh = lfcthresh, SetType = SetType, PBS = PBS)
    assay_unique_features <- Set_Analysis_Ligand_Unique(t0_anyset, plot=F, assay = cur_assay, SetType = SetType)
    assay_shared_features <- Set_Analysis_Ligand_Unique(t0_anyset, plot=F, assay = cur_assay, unique=F, SetType = SetType)
    
    assay_unique_format<-assay_unique_features %>%
      mutate(Type=paste0(cur_assay))
    assay_shared_format<-assay_shared_features %>%
      mutate(Type=paste0(cur_assay))
    
    total_unique_lista<-rbind(total_unique_list, assay_unique_format)
    total_shared_lista<-rbind(total_shared_list, assay_shared_format)
  }
  
  #check for directionality across timepoints for each unique feature
  lfc_rr<-read_csv(DE_lfc_path)
  lfc_long_rr <-lfc_rr %>%
    gather(key = experimentalCondition, value=log2FoldChange, -Type, -feature) %>%
    mutate(ligand=gsub('_.*', '' , experimentalCondition)) 
  
  dir_check<-left_join(total_unique_list, lfc_long_rr) %>%
    group_by(ligand, feature) %>%
    mutate(pos=sum(log2FoldChange>0), neg=sum(log2FoldChange<0)) %>%
    mutate(Direction = ifelse(pos == 1, 'Mixed', ifelse(pos==2, 'Positive', "Negative"))) %>%
    mutate(Abs_Direction = ifelse(Log2FoldChange>0, 'Positive', 'Negative')) %>%
    ungroup()
  
  total_unique_list_dir<-dir_check %>%
    select(-experimentalCondition, -log2FoldChange, -pos, -neg) %>%
    unique()
  
  return(list(total_unique_list_dir, total_shared_list))
}

#run set analysis on assay datasets
set_by_ligand <-perform_set_function(DE_pval_path=('../Analysis/Integrated/integrated_adj_p_values.csv'), DE_lfc_path<-('../Analysis/Integrated/integrated_matrix_lfc_rr_noPBS.csv'), assay=c('cycIF', 'RPPA', 'RNAseq', 'motifs', 'GCP'), SetType = 'Anyset', PBS = F)

write_csv(set_by_ligand[[1]], '../Analysis/Integrated/integrated_unique_features_rr_noPBS2_absDir.csv')
write_csv(set_by_ligand[[2]], '../Analysis/Integrated/integrated_shared_features_rr_noPBS2.csv')

#integrated figures for total feature counts
unique_stacked<-set_by_ligand[[1]] %>%
  group_by(ligand) %>%
  summarize(unique=n())

shared_stacked<-set_by_ligand[[2]] %>%
  group_by(ligand) %>%
  summarize(shared=n())

stacked_combined<-full_join(unique_stacked, shared_stacked) %>%
  gather(key = 'class', value='count', -ligand) 

ggplot(stacked_combined, aes(alpha=class, fill= ligand, y=count, x=factor(ligand, levels = c('EGF', 'BMP2', 'HGF', 'OSM', 'TGFB', 'IFNG')))) + 
  geom_bar(position="stack", stat="identity") + theme_bw() + labs(title = 'Integrated unique and shared upregulated features') + labs(fill = "Ligand") + xlab('ligand') + theme(text=element_text(size=18),panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                panel.grid.minor = element_blank(), axis.line = element_blank()) + scale_fill_manual(values=ligand_cols) + scale_alpha_manual(values = c(shared = .5, unique = 1)) +
  labs(alpha = "Feature Class") +
  geom_text(aes(label=count), size=5, position = position_stack(vjust = 1))

