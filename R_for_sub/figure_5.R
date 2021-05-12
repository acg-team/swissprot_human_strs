library(tidyverse)
cbPalette <- c("#D55E00", "#56B4E9", "#F0E442", "#000000", "#CC79A7", "#009E73", "#0072B2", "#E69F00")

# load overview table of occurrence of STRs and IDRs in pathology atlas proteins
overall_summary_df <- read.table("../data_for_sub/str_idr_pathology_atlas_summary.tsv", header=TRUE, sep="\t")

# For each anatomical site, generate a bar graph of the percentage of STR containing
# proteins that do or do not correlate to patient survival
overall_summary_df %>% filter(feature == "str") %>% 
  ggplot(aes(x=feature, y=perc_with_feature, fill=pa)) +
  geom_bar(stat="identity", position = position_dodge(width=0.55), width = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  facet_wrap(~site, ncol = 6, scale = "free_x") +
  ylab("Percentage with STR") +
  ylim(c(0, 30)) +
  scale_fill_manual(values=c(cbPalette[1], cbPalette[3], cbPalette[2]))

#############################
### SUPPLEMENTAL ANALYSES ###
#############################

# Make plot for percentage of disorder-promoting STR-containing proteins
dispro_summary_df <- read.table("../data_for_sub/dispro_str_pathology_atlas_summary.tsv", header=TRUE, sep="\t")

dispro_summary_df %>% filter(feature == "str", !pa == 0) %>% 
  ggplot(aes(x=feature, y=perc, fill=pa)) +
  geom_bar(stat="identity", position = position_dodge(width=0.55), width = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values=c(cbPalette[1], cbPalette[3], cbPalette[2])) +
  ylim(0, 25) +
  facet_wrap(~site, ncol = 6, scale = "free_x")
