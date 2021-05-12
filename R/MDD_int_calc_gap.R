library(cluster)
library(tidyverse)
library(factoextra)

b <- 100

# #Load the data matrix
# dm <- read_csv("../MDD_combined_dm.csv") %>%
#   select(-X1)
# 
# # #run the gap analysis usin PAM
# gss <- clusGap(dm, FUN = pam, K.max = 25, B = b, pamonce = 5, verbose = TRUE)
# # 
# # #Save the gap analysis object to disk
# save(gss, file = paste0("gss_PAM_combined_B",b,".rda"))

###shared dataset
#Load the data matrix
dm <- read_csv("../MDD_shared_dm.csv")

# #run the gap analysis usin PAM
gss <- clusGap(dm, FUN = pam, K.max = 25, B = b, pamonce = 5, verbose = TRUE)
# 
# #Save the gap analysis object to disk
save(gss, file = paste0("gss_PAM_shared_B",b,".rda"))

# # Silhouette method
# fviz_nbclust(dm, kmeans, method = "silhouette",k.max = 25)+
#   labs(subtitle = "Silhouette method, kmeans")
# 
# # Silhouette method
# fviz_nbclust(dm, pam, method = "silhouette",k.max = 25)+
#   labs(subtitle = "Silhouette method, PAM")
# 
# # Elbow method
# fviz_nbclust(dm, kmeans, method = "wss", k.max = 25) +
#   #geom_vline(xintercept = 4, linetype = 2)+
#   labs(subtitle = "Elbow method")

## the same with cluster-wise colours:
# pdf("MDD_silhouette_PAM.pdf")
# for(k in 2:25){
#   plot(silhouette(pam(dm, k=k, pamonce = 5)), main = paste("k = ",k), do.n.k=FALSE)
# }
# res <- dev.off()
# par(op)

