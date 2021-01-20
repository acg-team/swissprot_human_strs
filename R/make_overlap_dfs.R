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

uniprot_rest <- as.character(seq(1, 18))

for(i in uniprot_rest){
  tr_path <- paste0(base_path, "results/TRs/uniprot_rest/part", i, "_merged.tsv")
  tr_df <- load_tr_annotations(tr_path)
  trs_filt_df <- tr_df %>% 
    filter(pvalue <= 0.05 & divergence <= 0.1)
  
  idr_path <- paste0(base_path, "results/disorder/local_annotations/uniprot_rest/part", i, ".tsv")
  idr_df <- read.table(idr_path, header = FALSE, sep = "\t")
  colnames(idr_df) <- c("ID", "begin", "end")
  idr_df <- reformat_ids(idr_df)
  idr_df <- as.data.frame(idr_df)
  
  overlap <- calculate_tr_idr_overlap(trs_filt_df, idr_df)
  
  trs_with_overlap <- overlap %>% 
    filter(total_overlap > 0) %>% 
    summarize(count = n())
  
  tr_in_idr <- overlap %>% 
    filter(TRinDS_overlap > 0) %>% 
    summarize(count = n())
  
  summary <- paste0(i, " TRs:",  dim(overlap[1]), ", with overlap:", trs_with_overlap[1,1], ", TRs in IDR:", tr_in_idr[1,1])
  print(summary)
  out_name <- paste0(base_path, "results/overlap/local_disorder/uniprot_rest/part", i, "_TR_IDR_overlap.tsv")
  write.table(overlap, file = out_name, row.names = FALSE, sep="\t")
  print(paste("Wrote table to:", out_name))
}

# for(i in uniprot){
#   data <- load_pathway_data(base_path, i, meta=FALSE, idr=TRUE)
#   trs <- data[["trs"]]
#   trs_filt <- trs %>% 
#     filter(pvalue <= 0.05 & divergence <= 0.1)
#   idrs <- data[["idrs"]]  
#   
#   overlap <- calculate_tr_idr_overlap(trs_filt, idrs)
#   
#   trs_with_overlap <- overlap %>% 
#     filter(total_overlap > 0) %>% 
#     summarize(count = n())
#   
#   tr_in_idr <- overlap %>% 
#     filter(TRinDS_overlap > 0) %>% 
#     summarize(count = n())
#   
#   summary <- paste0(i, " TRs:",  dim(overlap[1]), ", with overlap:", trs_with_overlap[1,1], ", TRs in IDR:", tr_in_idr[1,1])
#   print(summary)
#   out_name <- paste0(base_path, "results/overlap/local_disorder/", i, "_TR_IDR_overlap.tsv")
#   write.table(overlap, file = out_name, row.names = FALSE, sep="\t")
#   print(paste("Wrote table to:", out_name))
# }