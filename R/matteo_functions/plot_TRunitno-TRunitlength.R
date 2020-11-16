#' Plot TR unit Number vs TR unit Length colored by total Repeat Length
#'
#' @param tr_sp a data.frame with TR results imported through \code{\link{load_tr_annotations}} and left-joined meta data from SwissProt through \code{\link{load_swissprot}}.
#' @param ngeq int. Subsets TRs with greater-equal the given integer on their TR-unit number for labelling.
#' @param lgeq int. Subsets TRs with greater-equal the given integer on their TR-unit length for labelling.
#'
#' @return ggplot object
plot_TRunitnoTRunitlength <- function(tr_sp = tr_nfkappab_sp, ngeq = 15, lgeq = 5){
  # Subset the data for labelling
  plot_data <- subset(tr_sp, xor(n_effective >= ngeq, l_effective >= lgeq))
  plot_data <- subset(tr_sp, subset = n_effective >= ngeq | l_effective >= lgeq)

  p1 <- ggplot(data = tr_sp, aes(x=tr_sp$l_effective, y=tr_sp$n_effective))+
    geom_point(aes(color = tr_sp$repeat_region_length), size = 4)+
    scale_color_gradient2( low="yellow", mid="#2D882D",
                           high="#AA3939", space ="Lab" )+
    labs(x= "TR unit length",
         y= "TR unit number",
         color= "Repeat Region Length")+
    geom_text_repel(data          = plot_data,
                    aes(x=plot_data$l_effective, y=plot_data$n_effective, label = plot_data$prot_name),
                    # nudge_y       = 50 - plot_data$n_effective,
                    size          = 3,
                    box.padding   = 0.5,
                    point.padding = 0.5,
                    force         = 100,
                    segment.size  = 0.2,
                    segment.color = "grey50"
    )
  p1 <- beautifier(p1)+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  return(p1)
}
