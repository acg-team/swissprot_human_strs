library(tidyverse)

# loading STR data
master_df <- read.table("../results/overlap/local_disorder/the_good_one.tsv", header=TRUE, sep="\t")

# loading signal peptide data
signal_peptides <- read.table("../results/misc/signal_peptides.tsv", header=FALSE, sep="\t")
colnames(signal_peptides) <- c("ID", "signal_begin", "signal_end")

# merging data frames
master_df <- master_df %>% left_join(signal_peptides, by="ID")

# for proteins without signal peptides, set the signal begin and end locations to 0
master_df$signal_begin[is.na(master_df$signal_begin)] <- 0 
master_df$signal_end[is.na(master_df$signal_end)] <- 0 

# how many proteins left after filtering out signal peptide STRs?
master_df %>% filter(begin > signal_end) %>% select(ID, promoting) %>% distinct() %>% summarise(count=n())

# how many STRs left after filtering out signal peptide STRs? (grouped by intrinsic disorder propensity)
master_df %>% select(ID, promoting) %>% distinct() %>% group_by(promoting) %>% summarise(count=n())

# calculate normalized STR centers
str_location_df <- master_df %>%
  transmute(tr_type,
            center=round(center),
            ID = ID,
            begin = begin,
            signal_end = signal_end,
            l_effective=l_effective,
            promoting=promoting,
            Length=round(Length),
            repeat_region_length=round(repeat_region_length)) %>%
  mutate(usableLen = Length - repeat_region_length,
         pos = (center-ceiling(repeat_region_length/2))/usableLen) %>%
  filter(usableLen > 0, pos <= 1)

# distribution of ordered STRs different from that of disordered and mixed?
t.test(str_location_df %>% filter(promoting=="order") %>% select(pos), str_location_df %>% filter(promoting=="disorder") %>% select(pos))$p.value
t.test(str_location_df %>% filter(promoting=="order") %>% select(pos), str_location_df %>% filter(promoting=="mixed") %>% select(pos))$p.value
t.test(str_location_df %>% filter(promoting=="mixed") %>% select(pos), str_location_df %>% filter(promoting=="disorder") %>% select(pos))$p.value

# Locations of STR in proteins, coloured on disorder promoting propensity 
pdf("../results/figures/strs_by_promoting.pdf")
str_location_df %>% 
  ggplot(aes(x = pos, colour=promoting, fill=promoting)) +
  geom_density(alpha=.7) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank())
dev.off()

# Pie chart of disorder propensity
pie_df <- master_df %>% group_by(promoting) %>% summarise(count=n())
pie_df$data_set <- "all"

pdf("../results/figures/nosig_strs_by_promoting_pie.pdf")
pie_df %>% ggplot(aes(x=data_set, y=count, fill=promoting)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank())
dev.off()

# distribution of ordered STRs different from that of disordered and mixed?
# now with STRs in signal peptides removed
t.test(str_location_df %>% filter(promoting=="order", begin > signal_end) %>% select(pos), str_location_df %>% filter(promoting=="disorder") %>% select(pos))$p.value
t.test(str_location_df %>% filter(promoting=="order", begin > signal_end) %>% select(pos), str_location_df %>% filter(promoting=="mixed") %>% select(pos))$p.value
t.test(str_location_df %>% filter(promoting=="mixed", begin > signal_end) %>% select(pos), str_location_df %>% filter(promoting=="disorder") %>% select(pos))$p.value

# Locations of STR in proteins, coloured on disorder promoting propensity but 
# STRs in signal peptides are removed
pdf("../results/figures/nosig_strs_by_promoting.pdf")
str_location_df %>% filter(begin > signal_end) %>% 
  ggplot(aes(x = pos, colour=promoting, fill=promoting)) +
  geom_density(alpha=.7) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank())
dev.off()

# Pie chart of disorder propensity but STRs in signal peptides are removed
pie_nosig_df <- master_df %>% filter(begin > signal_end) %>% group_by(promoting) %>% summarise(count=n())
pie_nosig_df$data_set <- "all"

pdf("../results/figures/nosig_strs_by_promoting_pie.pdf")
pie_nosig_df %>% ggplot(aes(x=data_set, y=count, fill=promoting)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank())
dev.off()
