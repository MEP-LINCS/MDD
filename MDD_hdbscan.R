#explore hdbscan 
library(tidyverse)
library(dbscan)

#Load the data matrix
dm <- read_csv("MDD_full_dm.csv") %>%
  select(-X1)


#run the clustering algorithm
full_cl <- hdbscan(dm, minPts = 100)
full_cl
plot(full_cl)

# run hdbscan on umap embedding
load("integrated_analysis/MDD_pam_UMAP.rda")

umap_df <- df_UMAP %>%
  select(UMAP_1, UMAP_2)

umap_cl <- hdbscan(umap_df, minPts = 35)
umap_cl
#plot(umap_cl)

plot(umap_df, col=umap_cl$cluster+1, 
     pch=ifelse(umap_cl$cluster == 0, 8, 1), # Mark noise as star
     cex=ifelse(umap_cl$cluster == 0, 0.5, 0.75), # Decrease size of noise
     xlab=NA, ylab=NA)
colors <- sapply(1:length(umap_cl$cluster), 
                 function(i) adjustcolor(palette()[(umap_cl$cluster+1)[i]], alpha.f = umap_cl$membership_prob[i]))
points(umap_df, col=colors, pch=20)

# #Save the gap analysis object to disk
save(gss, file = paste0("gss_PAM_shared_B",b,".rda"))
