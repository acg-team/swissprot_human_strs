#' Calculate the Amino Acid Frequency in a given AA-Sequence
#'
#' @description  For AA-frequency in a TR-MSA use \code{\link{AAfreq_in_TR}}.
#' @param sp a XStringSet objects
#'
#' @return A data.frame with AA, it's absolute frequency and it's ratios.
#' @export
AAfreq_in_prot <- function(sp){
  aa_freq <- alphabetFrequency(sp) # in the other column, are the spaces in between the TRs
  aa_freq <- colSums(aa_freq)
  # aa_freq <-  Reduce(`+`, aa_freq)

  df_aafreq  <- as.data.frame(aa_freq)
  df_aafreq$aa <- rownames(df_aafreq)
  # Remove NAs (*,-,.,+)
  # df_aafreq <- subset(df_aafreq, !(df_aafreq$aa %in% aa_ignore)) #TODO add
  # Calculate aa composition ratio in TR
  for (i in 1:nrow(df_aafreq)){
    df_aafreq[i,3] <- round(df_aafreq[i,1]/sum(df_aafreq[1]), 4)
  }
  colnames(df_aafreq) <- c("aa_freq", "aa", "aa_ratio")
  return(df_aafreq)
}
