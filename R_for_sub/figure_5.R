library(tidyverse)

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
  facet_wrap(~site, ncol = 6, scale = "free_x")
