#!/usr/bin/env R

library(dplyr)
source("src/reformat_ids.R")

get_SV_with_overlap <- function(i){
  return_list = list()
  meta_path <- paste0("../data/proteins/meta/merged_", i, "_meta.tsv")
  meta_info <- read.table(meta_path, sep="\t", header=TRUE)
  meta_info <- reformat_ids(meta_info)
  SV_ids <- meta_info %>% filter(SV == "yes") %>% select(ID)
  
  overlap_path <- paste0("../results/overlap/", i, "_TR_IDR_overlap.tsv" )
  overlap <- read.table(overlap_path, sep="\t", header=TRUE)
  overlap <- overlap %>% filter(total_overlap > 0 )
  SV_with_overlap <- inner_join(SV_ids, overlap, by="ID")
  return(list(overlap, SV_with_overlap))
}

pathways <- c(
  "wnt",
  "p53",
  "pi3k",
  "tgf-beta",
  "ras-mapk"
)

data <- get_SV_with_overlap("ras-mapk")
all_overlap <- data[[1]]
SV_overlap <- data[[2]]
