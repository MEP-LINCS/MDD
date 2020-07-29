suppressMessages(library(ComplexHeatmap))

#process all data from all assays
preprocess_level3 <- function(df, type){
  df_pp <- df %>%
    mutate(Type = type) %>%
    group_by(ligand, feature, Type, experimentalTimePoint) %>%
    summarise(value = median(value, na.rm = TRUE)) %>%
    ungroup()
  return(df_pp)
}

select_features <- function(df, var_quantile = 0){
  df_small <- df %>%
    group_by(feature) %>%
    mutate(feature_variance = var(value, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(feature_variance >= quantile(feature_variance,
                                        probs = var_quantile,
                                        na.rm = TRUE)) %>%
    dplyr::select(-feature_variance)
  return(df_small)
}

#' rrscale a numeric vector
rrscaleValues <- function(x, zeros = .001, ncores = 4, zscore_cutoff = Inf){
  x_rr <- as.matrix(x) %>%
    rrscale(zeros = zeros,  z = zscore_cutoff, ncores = ncores)
  return(x_rr)
}

getrrDetails <- function(x){
  T_name <- map(x, function(xx) xx[["T_name"]]) %>%
    unlist
  par_hat <- map(x, function(xx) xx[["par_hat"]]) %>%
    unlist
  df <- tibble(feature = names(x),
               par_hat = par_hat,
               T_name = T_name)
  return(df)
}

#' df a dataframe with feature, value, condition and Type columns. The values must
#' be all numeric and transformed to common scales

#prepare a datamatrix from a dataframe of data and metadta
prep_hm_matrix <- function(df, columns = "condition"){
  if(columns == "condition"){
    condition_order <-  c("ctrl_0",paste(rep(c("PBS", "HGF", "OSM", "EGF","BMP2", "IFNG", "TGFB"), each = 2), rep(c(24, 48),  times = 14),   sep = "_") )
    
    condition_order <- condition_order[condition_order %in% df$experimentalCondition]
    df_sp <- df %>%
      dplyr::select(feature, value, experimentalCondition) %>%
      spread(key = experimentalCondition, value = value, fill = median(df$value, na.rm = TRUE)) %>%
      dplyr::select(feature, condition_order)
    
    df_as_matrix <- df_sp %>%
      dplyr::select(-feature) %>%
      as.matrix()
    rownames(df_as_matrix) <- df_sp$feature
    return(df_as_matrix)
    
  } else if(columns == "ligand") {
    ligand_order <-  c("PBS", "HGF", "OSM", "EGF","BMP2", "IFNG", "TGFB")
    
    ligand_order <- ligand_order[ligand_order %in% df$ligand]
    df_sp <- df %>%
      filter(ligand %in% ligand_order) %>%
      #    mutate(feature = paste0(feature,"_",experimentalTimePoint)) %>%
      dplyr::select(feature, value, ligand) %>%
      spread(key = ligand, value = value, fill = median(df$value, na.rm = TRUE)) %>%
      dplyr::select(feature, ligand_order)
    
    df_as_matrix <- df_sp %>%
      dplyr::select(-feature) %>%
      as.matrix()
    rownames(df_as_matrix) <- df_sp$feature
    return(df_as_matrix)
    
  } else if(columns == "time") {
    #feature_ligand x time code
    time_order <-c(0,1,4,8,24,48)
    time_order <- time_order[time_order %in% df$experimentalTimePoint] %>%
      paste0("T",.)
    
    df_sp <- df %>%
      mutate(Time = paste0("T",experimentalTimePoint)) %>%
      filter(Time %in% time_order) %>%
      mutate(feature = paste0(feature,"_",ligand)) %>%
      dplyr::select(feature, value, Time) %>%
      spread(key = Time, value = value) %>%
      dplyr::select(feature, time_order)
    #removing median fill on missing data
    df_as_matrix <- df_sp %>%
      dplyr::select(-feature) %>%
      as.matrix()
    rownames(df_as_matrix) <- df_sp$feature
    return(df_as_matrix)
    
  } else stop("received columns value of ",columns, " which must be condition, ligand or time")
}

prep_hm_matrix_org <- function(df, columns = "condition"){
  if(columns == "condition"){
    condition_order <-  c("ctrl_0",paste(rep(c("PBS", "HGF", "OSM", "EGF","BMP2", "IFNG", "TGFB"), each = 2), rep(c(24, 48),  times = 14),   sep = "_") )
    
    condition_order <- condition_order[condition_order %in% df$experimentalCondition]
    df_sp <- df %>%
      dplyr::select(feature, value, experimentalCondition) %>%
      spread(key = experimentalCondition, value = value, fill = median(df$value, na.rm = TRUE)) %>%
      dplyr::select(feature, condition_order)
    
    df_as_matrix <- df_sp %>%
      dplyr::select(-feature) %>%
      as.matrix()
    rownames(df_as_matrix) <- df_sp$feature
    return(df_as_matrix)
    
  } else if(columns == "ligand") {
    ligand_order <-  c("PBS", "HGF", "OSM", "EGF","BMP2", "IFNG", "TGFB")
    
    ligand_order <- ligand_order[ligand_order %in% df$ligand]
    df_sp <- df %>%
      filter(ligand %in% ligand_order) %>%
      mutate(feature = paste0(feature,"_",experimentalTimePoint)) %>%
      dplyr::select(feature, value, ligand) %>%
      spread(key = ligand, value = value, fill = median(df$value, na.rm = TRUE)) %>%
      dplyr::select(feature, ligand_order)
    
    df_as_matrix <- df_sp %>%
      dplyr::select(-feature) %>%
      as.matrix()
    rownames(df_as_matrix) <- df_sp$feature
    return(df_as_matrix)
    
  } else if(columns == "time") {
    #feature_ligand x time code
    time_order <-c(0,1,4,8,24,48)
    time_order <- time_order[time_order %in% df$experimentalTimePoint] %>%
      paste0("T",.)
    
    df_sp <- df %>%
      mutate(Time = paste0("T",experimentalTimePoint)) %>%
      filter(Time %in% time_order) %>%
      mutate(feature = paste0(feature,"_",ligand)) %>%
      dplyr::select(feature, value, Time) %>%
      spread(key = Time, value = value) %>%
      dplyr::select(feature, time_order)
    #removing median fill on missing data
    df_as_matrix <- df_sp %>%
      dplyr::select(-feature) %>%
      as.matrix()
    rownames(df_as_matrix) <- df_sp$feature
    return(df_as_matrix)
    
  } else stop("received columns value of ",columns, " which must be condition, ligand or time")
}

prep_hm_annotations <- function(df, columns = "condition"){
  if(columns == "condition"){
    df_as_matrix <- prep_hm_matrix(df, columns = "condition")
    #create top annotations
    ann_nv_pairs <- df %>%
      dplyr::select(feature, Type) %>%
      distinct()
    
    ann_df <- tibble(feature = rownames(df_as_matrix)) %>%
      left_join(ann_nv_pairs, by = "feature") %>%
      dplyr::select(Type)
    return(ann_df)
  } else if(columns == "ligand") {
    df_as_matrix <- prep_hm_matrix(df, columns = "ligand")
    #create top annotations
    ann_nv_pairs <- df %>%
      dplyr::select(feature, Type, experimentalTimePoint) %>%
      mutate(feature = paste0(feature, "_", experimentalTimePoint)) %>%
      distinct()
    
    ann_df <- tibble(feature = rownames(df_as_matrix)) %>%
      left_join(ann_nv_pairs, by = "feature") %>%
      dplyr::select(Type, experimentalTimePoint) %>%
      rename(Time = experimentalTimePoint)
    return(ann_df)
  } else if(columns == "time") {
    df_as_matrix <- prep_hm_matrix(df, columns = "time")
    #create top annotations
    ann_nv_pairs <- df %>%
      dplyr::select(feature, Type, ligand) %>%
      mutate(feature = paste0(feature, "_", ligand)) %>%
      distinct()
    
    ann_df <- tibble(feature = rownames(df_as_matrix)) %>%
      left_join(ann_nv_pairs, by = "feature") %>%
      dplyr::select(Type, ligand)
    return(ann_df)
  } else stop("received columns value of ",columns, " which must be condition, ligand or time")
  
}

prep_hm_annotations_org <- function(df, columns = "condition"){
  if(columns == "condition"){
    df_as_matrix <- prep_hm_matrix(df, columns = "condition")
    #create top annotations
    ann_nv_pairs <- df %>%
      dplyr::select(feature, Type) %>%
      distinct()
    
    ann_df <- tibble(feature = rownames(df_as_matrix)) %>%
      left_join(ann_nv_pairs, by = "feature") %>%
      dplyr::select(Type)
    return(ann_df)
  } else if(columns == "ligand") {
    df_as_matrix <- prep_hm_matrix_org(df, columns = "ligand")
    #create top annotations
    ann_nv_pairs <- df %>%
      dplyr::select(feature, Type, experimentalTimePoint) %>%
      mutate(feature = paste0(feature, "_", experimentalTimePoint)) %>%
      distinct()
    
    ann_df <- tibble(feature = rownames(df_as_matrix)) %>%
      left_join(ann_nv_pairs, by = "feature") %>%
      dplyr::select(Type, experimentalTimePoint) %>%
      rename(Time = experimentalTimePoint)
    return(ann_df)
  } else if(columns == "time") {
    df_as_matrix <- prep_hm_matrix(df, columns = "time")
    #create top annotations
    ann_nv_pairs <- df %>%
      dplyr::select(feature, Type, ligand) %>%
      mutate(feature = paste0(feature, "_", ligand)) %>%
      distinct()
    
    ann_df <- tibble(feature = rownames(df_as_matrix)) %>%
      left_join(ann_nv_pairs, by = "feature") %>%
      dplyr::select(Type, ligand)
    return(ann_df)
  } else stop("received columns value of ",columns, " which must be condition, ligand or time")
  
}


get_iheatmap <- function(df_as_matrix, ...) {
  #Create the heatmap
  hm <- main_heatmap(data = df_as_matrix,
                     name = "score") 
}

format_hm <- function(hm, k = 6, cluster_method = "hclust", type_colors = NULL, groups = NULL, ...){
  hm_mod<- hm %>%
    # add_row_labels(font = list(size = 4),
    #                side = "right") %>%
    add_col_labels() %>%
    add_row_annotation(ann_df,
                       side = "left",
                       size = 0.05,
                       colors = list("Type" = type_colors))  %>%
    modify_layout(list(margin = list(r = 120)))
  if(!k==0){
    hm_mod <- hm_mod %>%
      add_row_clustering(name = "Cluster",
                         k = k,
                         method = cluster_method,
                         groups = groups,
                         colors = cluster_cols)
  }
  return(hm_mod)
}

#apply rrscale on all values of an assay
rrscale_assay <- function(df_long, zscore_cutoff = Inf){
  df_sp <- df_long %>%
    spread(feature, value)
  if("experimentalCondition" %in% colnames(df_sp)) {
    rr_objects <- df_sp %>%
      dplyr::select(-Type, -experimentalCondition) %>%
      as.matrix() %>%
      rrscale::rrscale(zeros = 0.01,z = zscore_cutoff)
    
    rr_mat <- rr_objects[["RR"]] %>%
      as_tibble()
    
    df_rr <- df_sp %>%
      dplyr::select(experimentalCondition, Type) %>%
      bind_cols(rr_mat) %>%
      gather("feature", "value", -Type, -experimentalCondition)
  } else {
    rr_objects <- df_sp %>%
      dplyr::select(-Type, -ligand, -experimentalTimePoint) %>%
      as.matrix() %>%
      rrscale::rrscale(zeros = 0.01,z = zscore_cutoff)
    
    rr_mat <- rr_objects[["RR"]] %>%
      as_tibble()
    
    df_rr <- df_sp %>%
      dplyr::select(ligand, experimentalTimePoint, Type) %>%
      bind_cols(rr_mat) %>%
      gather("feature", "value", -Type, -ligand, -experimentalTimePoint)
  }
  
  
  return(df_rr)
}

get_TFs <- function(sheet, dir_path, pattern) {
  #read in and combine MDD TF scores
  TFs_input <- dir(path = dir_path,
                   pattern = pattern,
                   full.names = TRUE) %>%
    map(readxl::read_excel, col_types = c("text"), sheet = sheet, .name_repair = "universal") %>%
    bind_rows() %>%
    mutate(Rank = as.numeric(Rank),
           experimentalTime = str_extract(Query.Name, "[24][48]"),
           experimentalTime = as.numeric(experimentalTime),
           ligand = str_remove(Query.Name, "ctrl vs "),
           ligand = str_remove(ligand, ",.*"),
           condition = paste0(ligand, "_", experimentalTime),
           value =(1633-Rank)/1632,
           feature = TF,
           Library_only = str_remove(Library, ",.*"),
           Library_only = str_replace(Library_only, " ", "--"))
}

hafun <- function(x, k, ...){
  hm <- get_iheatmap(x, assay_name = "Prior knowledge, rrscale-transformed", ...) %>%
    format_hm(k = k, cluster_method = cluster_method, ...)
  Cluster <- hm@plots@listData$score@data %>%
    as.data.frame() %>%
    mutate(cluster = hm@plots@listData$Cluster@data) %>%
    drop_na() %>%
    dplyr::select(cluster) %>%
    as.list()
}

plot_gap = function(x) {
  gstab = data.frame(x$Tab, k = seq_len(nrow(x$Tab)))
  ggplot(gstab, aes(k, gap)) + 
    geom_line() +
    geom_errorbar(aes(ymax = gap + SE.sim,
                      ymin = gap - SE.sim), width=0.1) +
    geom_point(size = 3, col=  "red") +
    labs(title = paste("Gap analysis to determine cluster number,",cluster_method))
}


plot_correlation <- function(df, md, assay_name, EGF_normed = TRUE,
                             ligand_cols = c( "PBS" = "#8dd3c7",
                                              "HGF" = "#80b1d3",
                                              "OSM" = "#fdb462",
                                              "EGF" = "#fb8072",
                                              "BMP2" = "#b3de69",
                                              "IFNG" = "#bebada",
                                              "TGFB" = "#ffd92f"),
                             ligand_2_cols = c("EGF" = "#ff0000",
                                               "None" = "#00000000")){
  #calculate correlations across the conditions and show in a heatmap
  
  #convert specimenID column names to specimanName, reorder delete EGF samples
  df <- prep_data(df, md) 
  
  title_suffix <- ""
  if(EGF_normed){
    title_suffix <- ", EGF Normalized"
    df <- df %>%
      dplyr::select(-matches("EGF"))
  }
  #Create annotation values
  #df_ann <-prep_annotations(df, md)
  
  haRow <- create_row_annotations(df, md)
  
  #Create the heatmap
  Heatmap(matrix = dplyr::select(df, -feature) %>%
            as.matrix() %>%
            t %>%
            scale(scale = FALSE) %>%
            t %>%
            cor(use = "complete.obs",method = "spearman"),
          name = "correlation",
          column_title = paste0("Correlation of ",assay_name, title_suffix),
          top_annotation = create_top_annotations(df, md),
          left_annotation = haRow,
          show_row_names = FALSE,
          show_column_names = FALSE,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
          na_col = "grey")
}
