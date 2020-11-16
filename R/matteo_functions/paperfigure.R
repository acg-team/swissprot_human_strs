#' ggplot function to prepare figures for publication
#'
#' @param p ggplot object
#' @param x.axis.text.angle angle of x-axis text
#' @param x.axis.text.hjust hjust for adjusting the height of the x-axis text
#'
#' @return ggplot object
#'
paper.figure <- function(p, x.axis.text.angle = 90, x.axis.text.hjust = NULL){
  p <- p + theme(panel.background = element_rect(fill = 'transparent', colour = NA),
                 text = element_text(size=25),
                 #panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                 legend.background = element_rect(colour = "white"),
                 legend.text = element_text(family = "sans", size=15, face='italic', hjust=0),
                 legend.key = element_rect(colour = 'white', fill = 'white'),
                 strip.background = element_rect(fill = 'transparent', colour = NA), # colour='red', fill='#CCCCFF'
                 strip.text.x = element_text(family = "sans",size=15, angle = 0),
                 strip.text.y = element_text(family = "sans",size=15, angle = 270, margin = margin(r=30)),
                 axis.text.x = element_text(family = "sans",size=15, angle = x.axis.text.angle, hjust = x.axis.text.hjust, margin=margin(1,1,2,1,"pt")),
                 axis.text.y = element_text(family = "sans",size=15, margin=margin(1,1,2,1,"pt")),
                 plot.background = element_rect(fill = 'transparent', colour = NA)
  )
  return(p)
}
