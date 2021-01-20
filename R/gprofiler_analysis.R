library(tidyverse)
data_df <- read.table("../results/functional_enrichment/gProfiler/gprofiler_reformat_processed_signals.tsv", sep="\t", header=TRUE, quote="")

plot_df <- data_df %>% arrange(log2fc)
fold_levels <- factor(plot_df$term_name, ordered=TRUE)
plot_df$term_name <- factor(plot_df$term_name, levels=fold_levels)

pdf("../results/figures/18012021_update/gprofiler_ora.pdf")
plot_df %>% group_by(source) %>% slice_max(n=25, order_by=negative_log10_of_adjusted_p_value) %>% ungroup() %>% 
  ggplot(aes(y=term_name, x=log2fc, fill=idr_fraction)) +
  scale_fill_gradient2(limits=c(0, 1)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  facet_wrap(~source, scales="free_y", nrow = 2)
dev.off()

# Basic dot plot
data_df %>% 
  ggplot(aes(x=idr_fraction, y=log2fc, colour=source)) +
  geom_point() +
  theme_classic()

# Plot with signal peptide fraction as smoothed line with second axis
coeff <- 3 # look at plot to determine proper coefficient for transformation
data_df %>% 
  ggplot(aes(x=idr_fraction)) +
  geom_smooth(aes(y=signal_fraction*coeff)) +
  geom_point(aes( y=log2fc, colour=source)) +
  scale_y_continuous(
    name = "log2(Fold change)",
    sec.axis = sec_axis(trans=~./coeff, name="Fraction with signal peptide")
  ) +
  theme_classic()

low_idr <- data_df %>% filter(idr_fraction < 0.5)
low_idr %>% select(term_name, idr_fraction, signal_fraction) %>% View()
low_idr %>% summarise(mean_sig_frac = mean(signal_fraction))
data_df %>% anti_join(low_idr, by="term_name") %>% summarise(mean_sig_frac = mean(signal_fraction))
