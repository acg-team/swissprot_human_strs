#' Filter for Tandem Repeats by Gene Name
#'
#' @param tr_sp a data.frame with TR results imported through \code{\link{load_tr_annotations}} and left-joined meta data from SwissProt through \code{\link{load_swissprot}}.
#' @param genes \code{str} of gene names. For details see \url{https://www.uniprot.org/help/gene_name}
#'
#' @export
#'
getTRbyGene <- function(tr_sp, genes){
  for (gene in genes){
    tr_sp_subset <- tr_sp[grepl(gene, tr_sp$gene_names),]
  }

  if (nrow(tr_sp_subset) != 0){
    return(tr_sp_subset)
  } else if ((nrow(tr_sp_subset) == 0)) {
    return(print("No protein with TRs found in this gene(s)."))
  } else {
    return(print("There is a bug. Somewhere... Find it!"))
  }
}
