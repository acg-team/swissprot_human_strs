#' Amino Acid Ratio Plotted Against Disorderpropensity
#'
#' @param tr_all a data.frame with TR results imported through \code{\link{load_tr_annotations}}. Requires to have the MSA included as column.
#' @param plot_title (str) Title of the plot
#' @return ggplot object
#' @export
AAratio_vs_Disorderpropensity <- function(tr_all, plot_title = NULL, sp_overall = FALSE){
  if (sp_overall){
    p1 <- ggplot(AAfreqSP)+
      geom_point(aes(x = as.factor(disorderpropensity), y = aa_ratio_sp))+
      scale_x_discrete(labels = as.character(AAfreqSP$aa),
                       breaks = AAfreqSP$disorderpropensity)+
      labs(x ="Amino Acid",
           y = "AA ratio",
           title = "All Swiss-Prot Proteins")+
      theme_minimal()
    p1 <- beautifier(p1)+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    return(p1)
  } else {
    if (nrow(tr_all) >=25 && colnames(tr_all) %in% c("aa_freq", "aa", "aa_ratio") ){
      p1 <- ggplot(AAfreq_all)+
        geom_point(aes(x = as.factor(disorderpropensity), y = aa_ratio))+
        scale_x_discrete(labels = as.character(AAfreq_all$aa),
                         breaks = AAfreq_all$disorderpropensity)+
        theme_minimal()

      if(!is.null(plot_title)){
        p1 <- p1 +
          labs(x ="Amino Acid",
               y = "AA ratio",
               title = plot_title)
      } else {
        p1 <- p1 +
          labs(x ="Amino Acid",
               y = "AA ratio")
      }

      p1 <- beautifier(p1)+
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
      # p1 <- paper.figure(p1)+
      #   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
      return(p1)
    } else {
      AAfreq_all <- AAfreq_in_TR(tr_all)
      AAfreq_all <- base::merge(x = AAfreq_all, y=AAfreqSP, by = "aa") # drop rows which have no disorderpropensity

      p1 <- ggplot(AAfreq_all)+
        geom_point(aes(x = as.factor(disorderpropensity), y = aa_ratio))+
        scale_x_discrete(labels = as.character(AAfreq_all$aa),
                         breaks = AAfreq_all$disorderpropensity)+
        theme_minimal()

      if(!is.null(plot_title)){
        p1 <- p1 +
          labs(x ="Amino Acid",
               y = "AA ratio",
               title = plot_title)
      } else {
        p1 <- p1 +
          labs(x ="Amino Acid",
               y = "AA ratio")
      }

      p1 <- beautifier(p1)+
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
      # p1 <- paper.figure(p1)+
      #   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
      return(p1)
    }
  }
}
