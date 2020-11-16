#' Tandem Repeat Distribution in Proteome
#'
#' @param tr_all_sp a data.frame with TR results imported through \code{\link{load_tr_annotations}}
#' @param paperfigure boolean. If TRUE, adapts figure for publication (larger font size, ...)
#' @param max_pvalue upper threshold of p-value. \code{tr_all_sp} is subsetted for pvalue <= max_pvalue.
#' @param max_l_effective upper threshold of l_effective \code{tr_all_sp} is subsetted for l_effective <= max_l_effective
#' @param max_n_effective upper threshold of n_effective \code{tr_all_sp} is subsetted for n_effective <= max_n_effective
#'
#' @return Distribution (Heatmap) of tandem repeats (TRs) in Swiss-Prot as a function of their repeat unit length $l_{effective} <= 80$ (x-Axis, x1) and their number of repeat units $n_{effective} <= 40$ (x2, y-Axis). Darker colour indicates a larger number of TRs with a specific length and number of repeats.
#' @export
TRdist_in_proteome <- function(tr_all_sp, paperfigure = FALSE, max_pvalue = 0.01, max_l_effective = 80, max_n_effective = 40){
  refactor = function(col) {
    factor(col, levels = min(col):max(col))
  }

  d <- subset(tr_all_sp,  pvalue <= max_pvalue & l_effective <= max_l_effective & n_effective <= max_n_effective)
  d_summary <- d %>%
    mutate(n_effective_rounded = round(n_effective)) %>%
    group_by(l_effective, n_effective_rounded) %>%
    summarize(count=n()) %>%
    ungroup() %>%
    transmute(l_effective = refactor(l_effective),
              n_effective_rounded = refactor(n_effective_rounded),
              log10count=log10(count))

  p <- ggplot(d_summary, aes(x=l_effective, y=n_effective_rounded)) +
    geom_point(aes(color= log10count, size = 10))+
    # scale_color_continuous(high = "#AA3939", low = "#2D882D")+
    scale_color_continuous(low = "#F9C73F", high = "#A21212")+
    labs(x="Repeat Unit Length", #expression(l[effective]),
         y="Number of repeats", #expression(n[effective]),
         color=expression(log[10](count))) +
    scale_x_discrete(breaks=c(1,seq(0,80,5),80)) +
    scale_y_discrete(breaks=c(2,seq(0,40,5),40)) +
    guides(size = FALSE) # remove size legend

  p <- beautifier(p)
  if (paperfigure){
    p <- paper.figure(p)
  }
  return(p)
}
