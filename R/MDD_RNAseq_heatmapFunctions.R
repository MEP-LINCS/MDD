# These functions are for plotting level3, level3-centered, and level4-centered
# RNAseq data.

HeatmapL3 <- function(mat, ca = ha.RNA.L3, featureName = "log2(fpkm + 1)", bks = c(0, 5, 150), showRow = FALSE, clust = FALSE, ...) {
  Heatmap(mat, 
          top_annotation = ca,
          cluster_columns = clust,
          show_row_names = showRow,
          cluster_row_slices = FALSE,
          row_names_gp = gpar(fontsize = 7),
          show_column_names = FALSE,
          name = featureName,
          col = colorRamp2(bks, c("blue", "white", "red")),
          ...)
}

HeatmapL3C <- function(mat, ca = ha.RNA.L3, featureName = "Z-score", bks = c(-2, 0, 2), showRow = FALSE, clust = FALSE, ...) {
  Heatmap(mat, 
          top_annotation = ca,
          cluster_columns = clust,
          show_row_names = showRow,
          cluster_row_slices = FALSE,
          row_names_gp = gpar(fontsize = 7),
          show_column_names = FALSE,
          name = featureName,
          col = colorRamp2(bks, c("blue", "white", "red")),
          ...)
}

HeatmapL4C <- function(mat, ca = ha.RNA.L4, featureName = "Z-score", bks = c(-3, 0, 3), showRow = FALSE, clust = FALSE, ...) {
  Heatmap(mat, 
          top_annotation = ca,
          cluster_columns = clust,
          show_row_names = showRow,
          cluster_row_slices = FALSE,
          row_names_gp = gpar(fontsize = 7),
          show_column_names = FALSE,
          name = featureName,
          col = colorRamp2(bks, c("blue", "white", "red")),
          ...)
}