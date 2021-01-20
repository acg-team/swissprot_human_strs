library(tidyverse)

str_purity <- function(msa, return_val){
  msa <- as.character(msa)
  frequencies <- as.data.frame(table(strsplit(msa, ",", fixed=TRUE)))
  frequencies <- frequencies %>% mutate(fraction = Freq/sum(Freq)) %>% arrange(desc(fraction))
  dom <- as.character(frequencies[1, "Var1"])
  # return fraction of the TR that is made up of the dominant unit
  if(return_val == "frac"){
    return(frequencies[1, "fraction"])
  } 
  # return the dominant unit
  else if(return_val == "dom"){
    return(dom)
  } 
  # return the longest uninterrupted stretch of the dominant unit
  else if(return_val == "longest"){
    max_count <- 0
    cur_count <- 0
    for(i in strsplit(msa, ",")[[1]]){
      if(i == dom){
        cur_count = cur_count + 1
      } else {
        if(cur_count > max_count){
          max_count = cur_count
        }
        cur_count = 0
      }
      if(cur_count > max_count){
        max_count = cur_count
      }
    }
  return(max_count)
  } 
  # wrong return_val specified
  else{
    stop(sprintf("in str_purity(msa, return_val=\"%s\"): unsupported return_val specified", return_val), call. = FALSE)
  }
}

master_df <- read.table("../results/overlap/local_disorder/the_good_one.tsv", header=TRUE, sep="\t")

master_df <- master_df %>% mutate(
  dominant_unit = sapply(msa_original, str_purity, return_val="dom"),
  frac_dominant = sapply(msa_original, str_purity, return_val="frac"),
  longest_stretch = sapply(msa_original, str_purity, return_val="longest")
)
