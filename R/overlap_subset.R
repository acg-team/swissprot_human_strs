#!/usr/bin/env R

# Setup of packages and functions
library(tidyverse)

# Important paths
base_path = "/Users/maxverbiest/PhD/projects/SP_CRC_Pathway_TRs/"

wnt_overlap_path <- paste0(base_path, "results/overlap/wnt_TR_IDR_overlap.tsv")
wnt_overlap <- read.table(wnt_overlap_path, sep="\t", header=TRUE)
subset <- wnt_overlap %>% filter(l_type == "homo", total_overlap == 0)
