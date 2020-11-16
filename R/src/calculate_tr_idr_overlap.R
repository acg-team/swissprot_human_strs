#!/usr/bin/env R
# Function is heavily based on:
# https://github.com/acg-team/swissrepeats/blob/master/results/compute_tr_all_overlap.R

calculate_tr_idr_overlap <- function(tr_df, idr_df) {
  # unique id for each IDR
  idr_df <- idr_df %>% 
    mutate(IDR_id = row_number())
  
  # unique id for each TR, add end of TR
  tr_df <- tr_df %>% 
    mutate(end = begin + repeat_region_length - 1,
           TR_id = row_number())

  # initialize empty columns for overlap values
  tr_df <- tr_df %>% 
    mutate(DSinTR_overlap = 0,
           TRinDS_overlap = 0,
           tail_overlap = 0,
           head_overlap = 0,
           total_overlap = 0)
  
  for (prot_id in unique(tr_df$ID)){ # for each protein that has a TR
    #print(c("Protein", prot_id))
    for (current_tr in tr_df$TR_id[which(tr_df$ID == prot_id)]){ # for each TR for this protein
      tr_begin <- tr_df$begin[which(tr_df$TR_id == current_tr)]
      tr_end <- tr_df$end[which(tr_df$TR_id == current_tr)]
      #print(c("TR", tr_begin, tr_end))
      for (current_idr in idr_df$IDR_id[which(idr_df$ID == prot_id)]){ # for all IDRs for this protein
        idr_begin <- idr_df$begin[which(idr_df$IDR_id == current_idr)]
        idr_end <- idr_df$end[which(idr_df$IDR_id == current_idr)]
        #print(c("IDR", idr_begin, idr_end))
                                  
        # check DS in TR -> There can be multiple (?) 
        # so overlap should be added instead of set for this one
        if(tr_begin <= idr_begin & tr_end >= idr_end){
          tr_df$DSinTR_overlap[which(tr_df$TR_id == current_tr)] = tr_df$DSinTR_overlap[which(tr_df$TR_id == current_tr)] + (idr_end - idr_begin + 1)
        }
        # check TR in DS
        else if(tr_begin >= idr_begin & tr_end <= idr_end){
          tr_df$TRinDS_overlap[which(tr_df$TR_id == current_tr)] = tr_end - tr_begin + 1
        }
        # check head overlap
        else if(tr_begin < idr_begin & tr_end > idr_begin){
          tr_df$head_overlap[which(tr_df$TR_id == current_tr)] = tr_end - idr_begin + 1
        }
        # check tail overlap
        else if(tr_begin > idr_begin & tr_begin < idr_end){
          tr_df$tail_overlap[which(tr_df$TR_id == current_tr)] = idr_end - tr_begin + 1
        }
      }
    }
  }
  tr_df$total_overlap <- rowSums(tr_df[,c("DSinTR_overlap", "TRinDS_overlap", "head_overlap", "tail_overlap")])
  return(tr_df)
}
