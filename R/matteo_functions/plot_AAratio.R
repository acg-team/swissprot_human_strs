#' Plot AA ratio
#'
#' @param df data.frame. AA frequency count obtained through \code{\link{AAfreq_in_TR}}.
#' @param plot_title str. Plot title.
#' @param ylim int. Scale the y-axis maximum value.
#'
#' @return ggplot object.
#' @export
#'
#' @examples plot_AAratio(df = AAfreq_in_TR(tr_all_sp))
plot_AAratio <- function(df = AAfreq, plot_title = NULL, ylim =1){
  p1 <- ggplot(df, aes(x = as.factor(df$aa), y = df$aa_ratio))+
    geom_segment( aes(x=as.factor(df$aa) ,xend=as.factor(df$aa), y=0, yend=df$aa_ratio), color="grey")+
    geom_point()+
    scale_x_discrete(labels = as.character(df$aa))+
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
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
    ylim(0,ylim)
  return(p1)
}
