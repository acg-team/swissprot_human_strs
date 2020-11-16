#!/usr/bin/env R

# Setup of packages and functions
library(tidyverse)

source("matteo_functions/load_tr_annotations.R")
source("src/source_all.R")
source_all("src/")

# Important paths
base_path = "/Users/maxverbiest/PhD/projects/SP_CRC_Pathway_TRs/"

pathways <- c(
  "wnt",
  "p53",
  "pi3k",
  "tgf-beta",
  "ras-mapk",
  "mmr",
  "pa_unfav",
  "pa_fav"
)

for(i in pathways){
  data <- load_pathway_data(base_path, i, meta=FALSE, idr=TRUE)
  trs <- data[["trs"]]
  trs_filt <- trs %>% 
    filter(pvalue <= 0.05 & divergence <= 0.1)
  idrs <- data[["idrs"]]  
  
  overlap <- calculate_tr_idr_overlap(trs_filt, idrs)
  
  trs_with_overlap <- overlap %>% 
    filter(total_overlap > 0) %>% 
    summarize(count = n())
  
  tr_in_idr <- overlap %>% 
    filter(TRinDS_overlap > 0) %>% 
    summarize(count = n())
  
  print(paste(i, "TRs with overlap:", trs_with_overlap[1,1], "TRs in IDR:", tr_in_idr[1,1]))
  out_name <- paste0(base_path, "results/overlap/", i, "TR_IDR_overlap.tsv")
  write.table(overlap, file = out_name, row.names = FALSE)
  print(paste("Wrote table to:", out_name))
}