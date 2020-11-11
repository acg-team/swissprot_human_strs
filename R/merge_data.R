#!/usr/bin/env R

library("dplyr")
library("ggplot2")


test_data <- read.csv("../data/swissrepeats/tr_annotations.csv", header=TRUE)
wnt_ids <- read.table("../data/merged_prot_sets/merged_wnt_meta.tsv", sep = "\t", header = TRUE)
wnt_ids$ID <- as.character(wnt_ids$ID)
mat_wnt <- read.table("../TRAL-Result-Analysis/inst/extdata/TRs_Wnt_proteins_CRC_sp_l100.tsv", sep = "\t", header = TRUE)
mat_wnt_short <- read.table("../TRAL-Result-Analysis/inst/extdata/TRs_Wnt_proteins_CRC_sp.tsv", sep = "\t", header = TRUE)

reformat_ids <- function(df) {
  df <- df %>% 
    rowwise() %>% 
    mutate(
      ID = strsplit(ID, "|", fixed=TRUE)[[1]][2]
    )
  return(df)
}

wnt_ids <- reformat_ids(wnt_ids)
head(wnt_ids, 10)

wnt_data <- inner_join(wnt_ids, test_data, by="ID")
print(length(unique(wnt_data$ID)))
wnt_data_all <- left_join(wnt_ids, test_data, by="ID")
print(length(unique(wnt_data_all$ID)))
head(wnt_data, 10)
wnt_data_all_filter <- wnt_data_all %>% 
  filter(
    !is.na(n_effective) & 
    n_effective >= 2.5 & 
    pvalue <= 0.05 & 
    l_effective <= 3
)
