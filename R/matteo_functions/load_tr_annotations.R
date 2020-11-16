#' Load Tandem Repeat Annotation Data
#'
#' @param path
#'
#' @return Dataframe object with the results from TRAL and a characterisation for the TR-unit size.
#' @export
load_tr_annotations <- function(path) {
  # path = paste(local_base_path, path, sep=local_path_separator)
  tr_all = read.csv(path, sep="\t", header=TRUE, quote="", stringsAsFactors = FALSE)
  #tr_all = subset(tr_all, pvalue < 0.01)
  tr_all$total_repeat_length = (tr_all$n_effective * tr_all$l_effective)
  tr_all$center = tr_all$begin + (tr_all$total_repeat_length - 1)/2
  tr_all$l_type = ifelse(tr_all$l_effective ==1, "homo",
                         ifelse(tr_all$l_effective >1 & tr_all$l_effective <= 3, "micro",
                                ifelse(tr_all$l_effective < 15, "small",
                                       "domain")))
  # tr_all$fraction_disordered_chars = tr_all$disordered_overlap / (tr_all$l_effective * tr_all$n_effective)
  return(as.data.frame(tr_all))
}
