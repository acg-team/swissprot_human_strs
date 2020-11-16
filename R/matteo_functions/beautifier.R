#' ggplot Beautifier Function
#'
#' @param p ggplot object
#' @param x.axis.text.angle angle of x-axis text
#' @param x.axis.text.hjust hjust for adjusting the height of the x-axis text
#'
#' @return ggplot object
#' @export
beautifier <- function(p, x.axis.text.angle = 90, x.axis.text.hjust = NULL){
  p <- p + theme(panel.background = element_rect(fill = 'transparent', colour = NA),
                 text = element_text(),
                 legend.background = element_rect(colour = "white"),
                 legend.text = element_text(family = "sans", face='italic', hjust=0),
                 legend.key = element_rect(colour = 'white', fill = 'white'),
                 strip.background = element_rect(fill = 'transparent', colour = NA), # colour='red', fill='#CCCCFF'
                 strip.text.x = element_text(family = "sans", angle = 0),
                 strip.text.y = element_text(family = "sans", angle = 270, margin = margin(r=30)),
                 axis.text.x = element_text(family = "sans", angle = x.axis.text.angle, hjust = x.axis.text.hjust, margin=margin(1,1,2,1,"pt")),
                 axis.text.y = element_text(family = "sans", margin=margin(1,1,2,1,"pt")),
                 axis.ticks.length = unit(0.05, "cm"),
                 plot.background = element_rect(fill = 'transparent', colour = NA)
  )
  return(p)
}
