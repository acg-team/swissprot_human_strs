#' Download Protein Sequence
#'
#' @param prot_id str. with the Protein identifier.
#'
#' @return a XStringSet objects
#' @export
download_prot_sequence <- function(prot_id){
  #################
  ### Download all protein sequences from swissprot and store them as fasta
  #################
  sp_url <- paste0("https://www.uniprot.org/uniprot/", prot_id, ".fasta")
  seq_sp <- getURL(sp_url)
  temp <- tempfile() # store it in a temporary file
  write.table(seq_sp, file = temp, sep = "", row.names = FALSE, col.names = FALSE, quote = FALSE)
  seq_sp <- readAAStringSet(temp, format = "fasta")
  unlink(temp)
  return(seq_sp)
}
