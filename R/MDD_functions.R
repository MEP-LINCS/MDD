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
rrscaleValues <- function(x, zeros = .001, ncores = 4){
  x_rr <- as.matrix(x) %>%
    rrscale(zeros = zeros, ncores = ncores)
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
    #   drop_na()
    
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
    #   drop_na()
    
    df_as_matrix <- df_sp %>%
      dplyr::select(-feature) %>%
      as.matrix()
    rownames(df_as_matrix) <- df_sp$feature
    return(df_as_matrix)
  } else if(columns == "time") {
    #feature_ligand x time code
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
    #feature_ligand x time code
  } else stop("received columns value of ",columns, " which must be condition, ligand or time")
  
}

get_iheatmap <- function(df_as_matrix, ...) {
  #Create the heatmap
  hm <- main_heatmap(data = df_as_matrix,
                     name = "score")
}

format_hm <- function(hm, k = 6, cluster_method = "hclust",...){
  hm_mod<- hm %>%
    add_row_labels(font = list(size = 4),
                   side = "right") %>%
    add_col_labels() %>%
    add_row_annotation(ann_df,
                       side = "left",
                       size = 0.05)  %>%
    modify_layout(list(margin = list(r = 120)))
  if(!k==0){
    hm_mod <- hm_mod %>%
      add_row_clustering(name = paste0("Cluster(",cluster_method,")"),
                         k = k,
                         method = cluster_method,
      ) %>%
      add_col_summary(size = .1,
                      groups = paste0("Cluster(",cluster_method,")"))
  }
  return(hm_mod)
}

#apply rrscale on all values of an assay
rrscale_assay <- function(df_long, zscore_cutoff = Inf){
  df_sp <- df_long %>%
    spread(feature, value)
  if("experimentalCondition" %in% colnames(df_sp)) {
    rr_objects <- df_sp %>%
      select(-Type, -experimentalCondition) %>%
      as.matrix() %>%
      rrscale(zeros = 0.01,z = zscore_cutoff)
    
    rr_mat <- rr_objects[["RR"]] %>%
      as_tibble()
    
    df_rr <- df_sp %>%
      select(experimentalCondition, Type) %>%
      bind_cols(rr_mat) %>%
      gather("feature", "value", -Type, -experimentalCondition)
  } else {
    rr_objects <- df_sp %>%
      select(-Type, -ligand, -experimentalTimePoint) %>%
      as.matrix() %>%
      rrscale(zeros = 0.01,z = zscore_cutoff)
    
    rr_mat <- rr_objects[["RR"]] %>%
      as_tibble()
    
    df_rr <- df_sp %>%
      select(ligand, experimentalTimePoint, Type) %>%
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

hafun <- function(x, k){
  hm <- get_iheatmap(x, assay_name = "Prior knowledge, rrscale-transformed") %>%
    format_hm(k = k, cluster_method = cluster_method)
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
