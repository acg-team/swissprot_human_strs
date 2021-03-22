library(tidyverse)

# loading STR data
master_df <- read.table("../data_for_sub/str_sp_final.tsv", header=TRUE, sep="\t")

# loading signal peptide data
signal_peptides <- read.table("../data_for_sub/signal_peptides_all_proteins_with_ptm.tsv", header=FALSE, sep="\t")
colnames(signal_peptides) <- c("ID", "signal_begin", "signal_end", "PTM")

# merging data frames
master_df <- master_df %>% left_join(signal_peptides, by="ID")

# for proteins without signal peptides, set the signal begin and end locations to 0
master_df$signal_begin[is.na(master_df$signal_begin)] <- 0 
master_df$signal_end[is.na(master_df$signal_end)] <- 0 


# calculate normalized STR centers
str_location_df <- master_df %>%
  transmute(tr_type,
            center=round(center),
            ID = ID,
            begin = begin,
            signal_end = signal_end,
            l_effective=l_effective,
            promoting=promoting,
            protein_length=round(protein_length),
            repeat_region_length=round(repeat_region_length)) %>%
  mutate(usableLen = protein_length - repeat_region_length,
         pos = (center-ceiling(repeat_region_length/2))/usableLen) %>%
  filter(usableLen > 0, pos <= 1)


# Locations of STR centers in proteins, coloured on disorder promoting propensity 
str_location_df %>% 
  ggplot(aes(x = pos, colour=promoting, fill=promoting)) +
  geom_density(alpha=.7) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size=20)) +
  xlab("Normalized STR center location")  +
  labs(fill="STR promotes", color="STR promotes")

# Pie chart of disorder propensity
pie_df <- master_df %>% group_by(promoting) %>% summarise(count=n())
pie_df$data_set <- "all"

pie_df %>% ggplot(aes(x=data_set, y=count, fill=promoting)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank())

# Locations of STR in proteins, colored on disorder promoting propensity but 
# STRs in signal peptides are removed
str_location_df %>% filter(begin > signal_end) %>% 
  ggplot(aes(x = pos, colour=promoting, fill=promoting)) +
  geom_density(alpha=.7) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size=20)) +
  xlab("Normalised STR centre location")  +
  labs(fill="STR promotes", color="STR promotes")

# Pie chart of disorder propensity but STRs in signal peptides are removed
pie_nosig_df <- master_df %>% filter(!end <= signal_end) %>% group_by(promoting) %>% summarise(count=n())
pie_nosig_df$data_set <- "all"

pie_nosig_df %>% ggplot(aes(x=data_set, y=count, fill=promoting)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank())
