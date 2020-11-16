#!/usr/bin/env R
# Function to reformat SP standard headers into UniProt ID only
# Example
# In:   >sp|Q13315|ATM_HUMAN Serine-protein kinase ATM OS=Homo sapiens OX=9606 GN=ATM PE=1 SV=4
# Out:  Q13315


reformat_ids <- function(df) {
  df <- df %>% 
    rowwise() %>% 
    mutate(
      ID = strsplit(as.character(ID), "|", fixed=TRUE)[[1]][2]
    )
  return(df)
}