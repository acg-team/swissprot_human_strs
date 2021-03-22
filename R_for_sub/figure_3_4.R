library(tidyverse)

data_df <- read.table("../data_for_sub/gpprofiler_ora.tsv", sep="\t", header=TRUE, quote="")

# Arrange data frame by descending log2(fold change)
plot_df <- data_df %>% arrange(log2fc)
fold_levels <- factor(plot_df$term_name, ordered=TRUE)
plot_df$term_name <- factor(plot_df$term_name, levels=fold_levels)

# For the four different data sources, plot the 25 over-represented terms with the
# lowest p value. Bars will be shaded based on IDR content in the term
plot_df %>% group_by(source) %>% slice_max(n=25, order_by=negative_log10_of_adjusted_p_value) %>% ungroup() %>% 
  ggplot(aes(y=term_name, x=log2fc, fill=idr_fraction)) +
  scale_fill_gradient2(limits=c(0, 1)) +
  geom_bar(stat="identity") +
  theme_classic() +
  labs(fill="Fraction\nwith IDR") +
  xlab("log2(Fold change)") +
  theme(axis.title.y = element_blank()) +
  facet_wrap(~source, scales="free_y", nrow = 2)

data_df %>% 
  ggplot(aes(x=idr_fraction, y=signal_fraction)) +
  geom_point(size = 4, aes(colour=source)) +
  geom_smooth() +
  theme_classic() +
  theme(text=element_text(size=25), 
        legend.background = element_rect(linetype = "solid", colour="black")) +
  coord_fixed() +
  labs(colour = "Data source") +
  xlab("Fraction with IDR") +
  ylab("Fraction with signal peptide")
