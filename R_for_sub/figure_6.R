library(tidyverse)
library(factoextra)

# Analysis of counts per GO term
to_cluster <- read.table("../data/GO_PA_counts/goa_counts_norm_CFP.tsv", sep="\t", header=TRUE)
row.names(to_cluster) <- to_cluster$pa_group
to_cluster <- to_cluster %>% dplyr::select(-pa_group)
to_cluster_filt <- to_cluster %>% select_if(colSums(.) != 0)

pca <- prcomp(x = to_cluster_filt, scale. = TRUE, center = TRUE)

# how many PCs?
fviz_eig(pca)

# Results for individuals
res.ind <- get_pca_ind(pca)
View(res.ind$contrib)      # Contributions to the PCs of each site

# kmeans clustering
# how many cluster centers?
set.seed(123)
fviz_nbclust(res.ind$contrib[,1:4], kmeans, method = "wss")

# cluster with 5 centers
set.seed(123)
cl_pca <- kmeans(x=res.ind$contrib[,1:4], centers=5, nstart=25, iter.max=10)
View(cl_pca$cluster)

# Visualizations
fviz_pca_ind(pca, axes=c(1, 2), repel = FALSE, col.ind = as.factor(cl_pca$cluster))

names = row.names(pca$x)
as.data.frame(pca$x[,1:2]) %>%
  ggplot(aes(x=PC1, y=PC2, colour=as.factor(cl_pca$cluster), labels=names)) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size=4) +
  # geom_text(aes(label=names),hjust=0, vjust=-1) +
  # theme_classic()
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab("PC 1 (11.2%)") +
  ylab("PC 2 (9.2%)") +
  labs(colour="Cluster")
