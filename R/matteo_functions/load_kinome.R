#' Load Human Proteome Kinome
#'
#' @param url string of the query to the UniProt REST API to download the protein sequences as fasta file.
#' @param path string with the filepath to the .fasta file of the protein sequences
#' @param OnlyIDs if TRUE, only the the ProteinIDs are return. If FALSE a XString class object of the AAStringSet is returned.
#'
#' @export
#'
#' @examples
#' url_kin <- "https://www.uniprot.org/uniprot/?query=ec:2.7.10.-%20OR%20ec:2.7.11.-%20OR%20ec:2.7.12.-%20OR%20ec:2.7.13.-%20OR%20ec:2.7.14.-%20OR%20ec:2.7.99.-&format=fasta&sort=score&fil=proteome:UP000005640%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22"
#' load_kinome(url_kin, "/myProject/url_kin.fasta")
load_kinome <- function(url, path, OnlyIDs = TRUE){
  # Download the swissprot file only if it doesn't already exist.
  # Uncomment this, if you want to use the most recent available data!
  if(!file.exists(path)){
    download.file(url, destfile = path)
  }

  # read the .fasta file
  sp_kinome <- Biostrings::readAAStringSet(path)

  if(!OnlyIDs){
    return(sp_kinome)
  } else if(OnlyIDs){
    # split fasta header to separate the prot IDs from the rest and convert it to a dataframe
    fastaheader <- base::strsplit(names(sp_kinome), split = "[|]")
    fastaheader <- data.frame(matrix(unlist(fastaheader), nrow=length(fastaheader), byrow = TRUE), stringsAsFactors = FALSE)
    # keep only the protIDs
    kinIDs <- fastaheader[,2]
    return(kinIDs)
  }
}
