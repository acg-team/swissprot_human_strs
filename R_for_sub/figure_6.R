library(tidyverse)
library(umap)
library(factoextra)
cbPalette <- c("#D55E00", "#56B4E9", "#F0E442", "#000000", "#CC79A7", "#009E73", "#0072B2", "#E69F00")

to.cluster.df <- read.table("../data_for_sub/goa_counts_CFP.tsv", sep="\t", header=TRUE) # unnormalized count data
# to.cluster.df <- read.table("../data_for_sub/goa_counts_norm_CFP.tsv", sep="\t", header=TRUE) # normalized count data

filter_df <- function(df, term_threshold=5, group_threshold=20){
  # remove columns and rows that contain very little information
  df.filt.cols <- df[, colSums(df != 0) >= term_threshold]
  df.filt.cols.rows <- df.filt.cols[rowSums(df.filt.cols[, -1] != 0) > group_threshold, ]
  
  return(df.filt.cols.rows)
}

to.cluster.filt.df <- filter_df(to.cluster.df)

labels <- to.cluster.filt.df$pa_group
data.df <- to.cluster.filt.df %>% select(-pa_group)

# uncomment next for log transform
# data.df <- log(data.df + 1)

cluster.embedding <- umap(
  data.df,
  n_neighbors=30, # less fine-grained clusters (focus less on noise)
  n_components=10,
  min_dist=0.00000001, # pack points more closely together
  random_state=123
  )

# plot WSS per number of centers to select optimum number of clusters
set.seed(123)
fviz_nbclust(as.data.frame(cluster.embedding$layout), kmeans, method = "wss")

# cluster (3 centers appears optimal)
set.seed(123)
clusters <- kmeans(x=as.data.frame(cluster.embedding$layout), centers=3, nstart=25, iter.max=10)
clusters$tot.withinss

# make embedding in 2 dimensions (default) for visualization
plot.embedding <- umap(data.df)

# plot 2D embedding, colour points by cluster
as.data.frame(plot.embedding$layout) %>%
  ggplot(aes(V1, V2, colour=as.factor(clusters$cluster))) +
  geom_point(size=4) +
  theme_bw() +
  theme(text = element_text(size=20),
        legend.box.background = element_rect(colour = "black")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  xlab(element_blank()) +
  ylab(element_blank()) +
  scale_colour_manual(values=cbPalette) +
  labs(colour="Cluster")
