#' Get the Protein ID from the Number of TR
#'
#' @param tr_dataframe a data.frame with TR results imported through \code{\link{function}}
#' @param no_TR number of TRs per protein.
#'
#' @return the Protein ID toghether with the number of detected TRs
#' @export
protein_id_by_number_of_TR <- function(tr_dataframe, no_TR){
  conttab <- table(tr_dataframe$ID)
  # get ID of the one with == no_TR
  df <- as.data.frame(which(conttab == no_TR))
  df <- cbind(rownames(df))
  colnames(df) <- c("prot_id")
  return(df)
}
