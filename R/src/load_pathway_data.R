#!/usr/bin/env R
# Simple function to load TR, meta (idr) data for CRC pathways
# Very inflexible, highly dependent on directory structure
# Requires functions load_tr_annotations(), reformat_ids()

load_pathway_data <- function(base_path, pathway, meta=FALSE, idr=FALSE){
  options = c("pa_unfav",
              "pa_fav",
              "mmr",
              "wnt",
              "p53",
              "pi3k",
              "tgf-beta",
              "ras-mapk")
  # Check if specified path is valid (error message is really ugly...)
  if(!pathway %in% options){stop("Please supply one of the supported pathways: ", paste(options))}
  
  tr_path <- paste0(base_path, "results/TRs/", pathway, "_TRs.tsv")
  tr_df <- load_tr_annotations(tr_path)
  data_list <- list("trs" = tr_df)
  
  if(meta){
    meta_path <- paste0(base_path, "data/proteins/meta/merged_", pathway, "_meta.tsv")
    meta_df <- read.table(meta_path, header = TRUE, sep = "\t")
    meta_df <- reformat_ids(meta_df)  
    data_list[["meta"]]= meta_df
  }
  
  if(idr){
    idr_path <- paste0(base_path, "results/disorder/", pathway, "_IDRs.tsv")
    idr_df <- read.table(idr_path, header = TRUE, sep = "\t")
    data_list[["idrs"]] = idr_df
  }
  return(data_list)
}
