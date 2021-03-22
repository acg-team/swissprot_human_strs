library(tidyverse)
data_df <- read.table("../results/functional_enrichment/gProfiler/gprofiler_reformat_processed_signals.tsv", sep="\t", header=TRUE, quote="")
str_in_sig_df <- read.table("../results/functional_enrichment/gProfiler/gprofiler_reformat_processed_str_in_sp.tsv", sep="\t", header=TRUE, quote="")

data_df <- data_df %>% left_join(str_in_sig_df %>% mutate(str_sig_frac=signal_fraction, str_sig_num=intersection_size*signal_fraction) %>% dplyr::select(term_name, str_sig_frac, str_sig_num), by="term_name")

plot_df <- data_df %>% arrange(log2fc)
fold_levels <- factor(plot_df$term_name, ordered=TRUE)
plot_df$term_name <- factor(plot_df$term_name, levels=fold_levels)

pdf("../results/figures/18012021_update/gprofiler_ora.pdf")
plot_df %>% group_by(source) %>% slice_max(n=25, order_by=negative_log10_of_adjusted_p_value) %>% ungroup() %>% 
  ggplot(aes(y=term_name, x=log2fc, fill=idr_fraction)) +
  scale_fill_gradient2(limits=c(0, 1)) +
  geom_bar(stat="identity") +
  theme_classic() +
  labs(fill="Fraction\nwith IDR") +
  xlab("log2(Fold change)") +
  theme(axis.title.y = element_blank()) +
  facet_wrap(~source, scales="free_y", nrow = 2)
dev.off()

# Basic dot plots
# log2(fc) vs idr fraction 
data_df %>% 
  ggplot(aes(x=idr_fraction, y=log2fc, colour=source)) +
  geom_point(size = 2) +
  theme_classic()

# signal_fraction vs idr_fraction
pdf("../results/figures/raw_figures/signal_v_idr_fraction.pdf")
data_df %>% 
  ggplot(aes(x=idr_fraction, y=signal_fraction)) +
  geom_point(size = 4, aes(colour=source)) +
  geom_smooth() +
  theme_classic() +
  theme(text=element_text(size=25), 
        legend.background = element_rect(linetype = "solid", colour="black")) +
  # ylim(0, 1) +
  # facet_wrap(~source) +
  coord_fixed() +
  labs(colour = "Data source") +
  xlab("Fraction with IDR") +
  ylab("Fraction with signal peptide")
dev.off()

# absolute values of idr and signal containing
data_df %>% 
  ggplot(aes(x=idr_fraction * intersection_size, y=signal_fraction * intersection_size)) +
  geom_smooth() +
  geom_point(size = 2, aes(colour=source)) +
  theme_classic()

# Plot with signal peptide fraction as smoothed line with second axis
coeff <- 3 # look at plot to determine proper coefficient for transformation
data_df %>% 
  ggplot(aes(x=idr_fraction)) +
  # geom_line(aes(y=signal_fraction*coeff)) +
  geom_smooth(aes(y=signal_fraction*coeff)) +
  geom_point(aes(y=log2fc, colour=source), size = 2) +
  scale_y_continuous(
    name = "log2(Fold change)",
    sec.axis = sec_axis(trans=~./coeff, name="")
  ) +
  theme_classic()

low_idr <- data_df %>% filter(idr_fraction < 0.5)
low_idr %>% select(source, term_name, idr_fraction, signal_fraction) %>% View()
low_idr %>% select(source, term_name, idr_fraction, signal_fraction, str_sig_frac, str_sig_num) %>% View()
low_idr %>% summarise(mean_sig_frac = mean(signal_fraction))
data_df %>% anti_join(low_idr, by="term_name") %>% summarise(mean_sig_frac = mean(signal_fraction))

high_sig <- data_df %>% filter(signal_fraction >= 0.1)
high_sig %>% select(source, term_name, idr_fraction, signal_fraction) %>% View()

data_df %>% filter(idr_fraction > 0.5, signal_fraction > 0.5) %>% 
  select(source, term_name, idr_fraction, signal_fraction) %>% View()

data_df %>% 
  select(source, term_name, idr_fraction, signal_fraction) %>% View()

