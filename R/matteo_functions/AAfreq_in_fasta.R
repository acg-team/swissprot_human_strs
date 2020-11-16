#' AA frequency from a fasta file
#'
#' @param path str. Path to the fasta file.
#' @param aa_ignore character vector of symbols to ignore in the count.
#'
#' @return A data.frame with AA, it's absolute frequency and it's ratios.
#' @export
AAfreq_in_fasta <- function(path = "/home/delt/ZHAW/CRC_TRs/data/nfkappab_proteins_CRC_sp.fasta", aa_ignore = c("U", "O", "B", "J", "Z", "X", "*","-",".","+", "other")){
  fasta_file <- readAAStringSet(filepath = path, format = "fasta")

  aa_freq <- alphabetFrequency(fasta_file) # in the other column, are the spaces in between the TRs
  aa_freq <- colSums(aa_freq)
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
