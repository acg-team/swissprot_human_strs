library(tidyverse)
master_df <- read.table("../results/overlap/local_disorder/the_good_one.tsv", header=TRUE, sep="\t")

strs <- data.frame(data_set = c(rep("non_pa", 2), rep("pa_fav", 2), rep("pa_unfav", 2)),
                     feature=rep("strs", 6),
                     containing=rep(c("without", "with"), 3),
                     count=c(17226, 2575, 316, 35, 194, 48),
                     percentage=c(87.0, 13.0, 90.0, 10.0, 80.2, 19.8))

idrs <- data.frame(data_set = c(rep("non_pa", 2), rep("pa_fav", 2), rep("pa_unfav", 2)),
                   feature=rep("idrs", 6),
                   containing=rep(c("without", "with"), 3),
                   count=c(8868, 10933, 143, 208, 41, 201),
                   percentage=c(44.8, 55.2, 40.8, 59.2, 16.4, 83.6))
bar_df <- rbind(strs, idrs)


pdf("../results/figures/pa_str_occurence.pdf")
bar_df %>% filter(containing == "with") %>% ggplot(aes(x=data_set, y=percentage, fill=feature)) +
  geom_bar(stat="identity", position = position_dodge(width=0.55), width=0.5) +
  theme_classic() +
  theme(axis.title.x = element_blank())
dev.off()

# pa stuff
fav_ids <- read.table("../data/thyroid/thyroid_fav_ids.tsv", header=FALSE, sep="\t")
fav_ids$pa_crc <- "pa_fav"

unfav_ids <- read.table("../data/thyroid/thyroid_unfav_ids.tsv", header=FALSE, sep="\t")
unfav_ids$pa_crc <- "pa_unfav"

all_pa_ids <- rbind(fav_ids, unfav_ids)
colnames(all_pa_ids) <- c("ID", "pa")
master_df <- left_join(master_df, all_pa_ids, by="ID")
master_df %>% select(ID, pa) %>% distinct() %>%  group_by(pa) %>% summarise(count=n())

# Fisher exact tests to see if occurence of STRs is different between groups for PA CRC
uni_ids <- read.table("../data/pathology_atlas/crc/test_uniprot_only.csv", header=TRUE, sep=",")
uni_ids <- uni_ids %>% select(ID, pa)
uni_ids %>% select(ID, pa) %>% distinct() %>% group_by(pa) %>% summarise(count=n())

master_df <- master_df %>% left_join(uni_ids, by="ID")
master_df %>% select(ID, pa) %>% distinct() %>% group_by(pa) %>% summarise(count=n())

source("src/pa_str_counts.R")
pa_str_counts(master_df, all_fav=351, all_unfav=242)
