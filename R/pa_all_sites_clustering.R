library(klaR)
library(tidyverse)

# too_large = c("kidney_pa_fav", "kidney_pa_unfav", "liver_pa_unfav")
# too_large_two = c("kidney_pa_fav", "kidney_pa_unfav", "liver_pa_unfav", "endometrium_pa_fav", "endometrium_pa_unfav", "pancreas_pa_fav", "pancreas_pa_unfav")
# to_cluster <- to_cluster %>% select(-kidney_pa_fav, -kidney_pa_unfav, -liver_pa_unfav, -endometrium_pa_fav, -endometrium_pa_unfav, -pancreas_pa_fav, -pancreas_pa_unfav)

to_cluster <- read.table("../data/pathology_atlas/pa_str_proteins_wide.tsv", sep="\t", header=TRUE)
ids <- to_cluster$ID
to_cluster <- as.data.frame(t(to_cluster))[-1,]
colnames(to_cluster) <- ids

# kmodes clustering
# set.seed(123)
# cl <- klaR::kmodes(data=to_cluster, modes=6, iter.max=10)
# clusters <- data.frame(names=row.names(to_cluster), cluster=cl$cluster)

# Multiple Correspondence Analysis
mca <- FactoMineR::MCA(to_cluster, ncp=10, graph = FALSE)
ind_contributions <- as.data.frame(mca$ind$contrib)
var_contributions <- as.data.frame(mca$var$contrib)
mca_plot_df <- as.data.frame(mca$ind$coord)
as.data.frame(mca$eig) %>% View()
as.data.frame(mca$eig)[1:5,] %>% ggplot(aes(x=seq(1,5), y=`cumulative percentage of variance`)) + geom_line()

# kmeans clustering of first three PCs from MCA
set.seed(123)
cl_mca <- kmeans(x=ind_contributions[,1:3], centers=3, nstart=25, iter.max=10)
cl_mca$tot.withinss
as.data.frame(cl_mca$cluster) %>% View()

fviz_cluster(cl_mca, data=ind_contributions[,1:3], geom="point")

mca_plot_df$cluster <- cl_mca$cluster
mca_plot_df %>% 
  ggplot(aes(x=`Dim 1`, y=`Dim 2`, colour=as.factor(cluster), label=row.names(.))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size=2) +
  theme_bw() +
  xlim(-1, 2.5) +
  ylim(-1, 2.5) +
  xlab("Dim 1 (20.85 %)") +
  ylab("Dim 2 (18.46 %)") +
  geom_text(aes(label=names), hjust=0, vjust=2)
  

# Analysis of counts per GO term
to_cluster <- read.table("../data/GO_PA_counts/goa_counts_norm_CFP.tsv", sep="\t", header=TRUE)
# to_cluster <- to_cluster %>% filter(!pa_group %in% c("liver_pa_fav", "liver_pa_unfav", "kidney_pa_fav", "kidney_pa_unfav"))
row.names(to_cluster) <- to_cluster$pa_group
to_cluster <- to_cluster %>% dplyr::select(-pa_group)
to_cluster_filt <- to_cluster %>% select_if(colSums(.) != 0)

# pca, center values around 0 and set st.dev to 1
library(factoextra)
pca <- prcomp(x = to_cluster_filt, scale. = TRUE, center = TRUE)
summary(pca)
View(pca$x)

# how many PCs?
fviz_eig(pca)

eig.val <- get_eigenvalue(pca)
eig.val

# Results for Variables
res.var <- get_pca_var(pca)
res.var$coord          # Coordinates
View(res.var$contrib)        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
View(res.ind$contrib)      # Contributions to the PCs
View(res.ind$cos2)           # Quality of representation

# kmeans clustering
set.seed(123)
fviz_nbclust(res.ind$contrib[,1:4], kmeans, method = "wss")

set.seed(123)
cl_pca <- kmeans(x=res.ind$contrib[,1:4], centers=5, nstart=25, iter.max=10)
View(cl_pca$cluster)

# Plot with clusters
fviz_pca_ind(pca, axes=c(1, 2), repel = FALSE, col.ind = as.factor(cl_pca$cluster))
fviz_pca_var(pca, axes=c(1, 2), repel = TRUE, col.ind = as.factor(cl_pca$cluster))

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


         