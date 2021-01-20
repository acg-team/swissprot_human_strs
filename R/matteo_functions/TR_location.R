#' Display TR center location
#'
#' @param tr_all_sp a data.frame with TR results imported through \code{\link{load_tr_annotations}}
#' @param plot_title (str) Title of the plot
#' @param byTRtype boolean. If TRUE, function returns plot filled by TR-type
#' @param byProt boolean. IF TRUE, function returns plot filled by each TR.
#' @param paperfigure boolean. If TRUE, adapts figure for publication (larger font size, ...). Modifications can be done in \code{\link{paper.figure}}
#' @description The TR center is normalized to reduce boundary effects. \cr
#' The center-bias comes from boundary effects from the simple center/Length metric employed; \cr
#' if the TR covers a significant fraction of the protein, then it's center necessarily falls near the middle.\cr
#'  We can compensate for this by normalizing over only the valid center locations:\cr
#'  \deqn{x = \frac{center - l_{eff}*n_{eff}/2}{N-l_{eff}*n_{eff}}} \cr
#'  \cr TR center locations of 0 and 1 refere to the N- and C-terminus respectively.
#'
#' @return ggplot object
#' @export
TR_location <- function(tr_all_sp, plot_title = NULL, byTRtype = FALSE, byProt = TRUE, paperfigure = FALSE){
  if (byTRtype){
    # Only show one of the plot per time
    # TODO: implement showing both at the same time.
    byProt = FALSE

    # renormalized version accounting for true length
    tr_all_sp$repeat_type_bin <- with(tr_all_sp, cut(l_effective, breaks = c(0,1,3, 15, 2000), dig.lab = 12))
    # tr_all_sp_positions <- tr_all_sp %>%
    #   transmute(repeat_type_bin,
    #             center=round(center),
    #             Length=round(Length),
    #             total_repeat_length=round(total_repeat_length)) %>%
    #   mutate(usableLen = Length - total_repeat_length,
    #          pos = (center-ceiling(total_repeat_length/2))/usableLen) %>%
    #   filter(usableLen > 0) %>% #otherwise causes NaN or Inf
    #   tbl_df
    tr_all_sp_positions <- tr_all_sp %>%
      transmute(repeat_type_bin,
                center=round(center),
                l_effective=l_effective,
                promoting=promoting,
                Length=round(Length),
                repeat_region_length=round(repeat_region_length)) %>%
      mutate(usableLen = Length - repeat_region_length,
             pos = (center-ceiling(repeat_region_length/2))/usableLen) %>%
      filter(usableLen > 0) %>% 
      tbl_df

    # possible bug: some centers lie outside the usable length
    # tr_all_sp_positions %>%
    #   mutate(diff= center - floor(Length+1/2 - total_repeat_length/2)) %>%
    #   arrange(-pos) %>%
    #   top_n(1,diff) %>%
    #   invisible

    cols1.4 <- c("#2D882D", "#AA3939", "#AA7939", "#29506D")
    cols2.4 <- c("#2D882D", "#AA3939", "#AA7939", "#29506D")
    
    p1 <- ggplot(tr_all_sp_positions %>% filter(pos <= 1, l_effective <= 2), aes(x = pos, colour=promoting, fill=promoting)) +
      geom_density(alpha=.7)

    # p1 <- ggplot(tr_all_sp_positions %>% filter(pos <= 1), aes(x = pos, colour=repeat_type_bin, fill=repeat_type_bin)) +
    #   geom_density(alpha=.7) +
    #   scale_color_manual( values = cols2.4, #values = c("grey","violet","turquoise4"),
    #                       name  ="TR type",
    #                       breaks=c("(0,1]","(1,3]","(3,15]", "(15,2000]"),
    #                       labels=c("Homo", "Micro", "Small","Domain")) +
    #   scale_fill_manual( values = cols1.4,
    #                      name  ="TR type",
    #                      breaks=c("(0,1]","(1,3]","(3,15]", "(15,2000]"),
    #                      labels=c("Homo", "Micro", "Small","Domain"))

    if(!is.null(plot_title)){
      p1 <- p1 +
        labs(x ="TR center location",
             y = "Density",
             title = plot_title)
    } else {
      p1 <- p1 +
        labs(x ="TR center location",
             y = "Density")
    }

    p1 <- beautifier(p1, x.axis.text.angle = 0)
    if (paperfigure){
      p1 <- paper.figure(p1, x.axis.text.angle = 0)
    }
    returnplot <- p1
  }

  if (byProt){
    # renormalized version accounting for true length
    # tr_all_sp_positions <- tr_all_sp %>%
    #   transmute(ID,
    #             center=round(center),
    #             Length=round(Length),
    #             total_repeat_length=round(total_repeat_length)) %>%
    #   mutate(usableLen = Length - total_repeat_length,
    #          pos = (center-ceiling(total_repeat_length/2))/usableLen) %>%
    #   filter(usableLen > 0) %>% #otherwise causes NaN or Inf
    #   tbl_df
    tr_all_sp_positions <- tr_all_sp %>%
      transmute(repeat_type_bin,
                center=round(center),
                Length=round(Length),
                repeat_region_length=round(repeat_region_length)) %>%
      mutate(usableLen = Length - repeat_region_length,
             pos = (center-ceiling(repeat_region_length/2))/usableLen) %>%
      filter(usableLen > 0) %>% 
      tbl_df

    # possible bug: some centers lie outside the usable length
    # tr_all_sp_positions %>%
    #   mutate(diff= center - floor(Length+1/2 - total_repeat_length/2)) %>%
    #   arrange(-pos) %>%
    #   top_n(1,diff) %>%
    #   invisible

    p2 = ggplot(tr_all_sp_positions %>% filter(pos <= 1), aes(x = pos, colour=ID, fill=ID)) +
      geom_density(alpha=.7) +
      scale_color_discrete(name = "Protein ID")+
      scale_fill_discrete(name = "Protein ID")+
      labs(x="TR center location", y="Density")

    if(!is.null(plot_title)){
      p2 <- p2 +
        labs(x ="TR center location",
             y = "Density",
             title = plot_title)
    } else {
      p2 <- p2 +
        labs(x ="TR center location",
             y = "Density")
    }

    p2 <- beautifier(p2, x.axis.text.angle = 0)
    if (paperfigure){
      p2 <- paper.figure(p2, x.axis.text.angle = 0)
    }

    if (byTRtype){
      returnplot <- returnplot + p2
    } else {
      returnplot <- p2
    }
  }
  return(returnplot)
}
