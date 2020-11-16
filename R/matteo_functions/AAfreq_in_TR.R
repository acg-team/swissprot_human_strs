#' Calculate the Amino Acid Frequency in a given TR-MSA
#'
#' @param tr_all_sp a data.frame with TR results imported through \code{\link{load_tr_annotations}}. Requires to have the MSA included as column.
#'
#' @return A data.frame with AA, it's absolute frequency and it's ratios.
#' @export
AAfreq_in_TR <- function(tr_all_sp, includeSpaces=FALSE, aa_ignore = c("U", "O", "B", "J", "Z", "X", "*","-",".","+", "other")){
  tr_all_sp <- tr_all_sp[, grepl("msa", colnames(tr_all_sp))] # Look for the column, which contains the MSA.
  aa_freq <- alphabetFrequency(AAStringSet(tr_all_sp)) # in the 'other' column, are the spaces in between the TRs

  if (includeSpaces){
    aa_freq <- colSums(aa_freq)
  } else {
    aa_freq <- colSums(aa_freq[,-ncol(aa_freq)])
  }
  # aa_freq <-  Reduce(`+`, aa_freq)

  df_aafreq  <- as.data.frame(aa_freq)
  df_aafreq$aa <- rownames(df_aafreq)
  # Remove NAs (*,-,.,+)
  df_aafreq <- subset(df_aafreq, !(df_aafreq$aa %in% aa_ignore)) #TODO add
  # Calculate aa composition ratio in TR
  for (i in 1:nrow(df_aafreq)){
    df_aafreq[i,3] <- round(df_aafreq[i,1]/sum(df_aafreq[1]), 4)
  }
  colnames(df_aafreq) <- c("aa_freq", "aa", "aa_ratio")
  return(df_aafreq)
}
