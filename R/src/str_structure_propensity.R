#!/usr/bin/env R
library(tidyverse)

master_df <- read.table("../results/overlap/local_disorder/master_overlap.tsv", header=TRUE, sep="\t")

aa_df <- data.frame(
  aa = c("W", "F", "Y",	"I", "M", "L", "V", "N", "C", "T", "A", "G", "R", "D", "H", "Q", "K", "S", "E", "P"), 
  propensity = c(rep("order", 10), rep("disorder", 10)),
  score = c(-0.884, -0.697, -0.510, -0.486, -0.397, -0.326, -0.121, 0.007, 0.02, 0.059, 0.06, 0.166, 0.180, 0.192, 0.303, 0.318, 0.586, 0.341, 0.736, 0.987)
)

disorder_propensity <- function(msa) {
  msa <- str_remove_all(msa, ",")
  msa <- str_remove_all(msa, "-")
  msa <- str_remove_all(msa, " ") # No idea why spaces show up sometimes
  ordered <- FALSE
  disordered <- FALSE
  for (i in strsplit(msa, "")[[1]]) {
    propensity <- as.character(aa_df$propensity[which(aa_df$aa == i)])[[1]]
    if (propensity == "order") {
      ordered = TRUE
    } else if (propensity == "disorder") {
      disordered = TRUE
    }
  }
  if (disordered == TRUE){
    if (ordered == TRUE) {
      return("mixed")
    }
    return("disorder")
  }
  return("order")
}

# master_df <- master_df %>% mutate(promoting = sapply(msa_original, disorder_propensity))
