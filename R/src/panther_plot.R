library(tidyverse)

# Read in table of differentially represented GO terms
panther <- read.table("../results/functional_enrichment/PANTHER/str_vs_all_mfcc_complete.tsv", header=TRUE, sep="\t")
panther <- panther %>% mutate(Fold.Enrichment = as.numeric(as.character(Fold.Enrichment)))
panther[is.na(panther)] <- 0.01 # Do not make 0, because otherwise log2(0) = -inf

# Calculate log2-fold enrichments, order the table using these values (low -> high)
panther <- panther %>% mutate(log2_fold = log(Fold.Enrichment, base = 2)) %>% arrange(log2_fold)

GO_levels <- factor(panther$term, ordered=TRUE)
panther$term <- factor(panther$term, levels=GO_levels)

# Make bar chart (ordered on log2-fold, coloured on GO domain) and write to file
pdf(file = "../results/functional_enrichment/PANTHER/plots/ccmf_complete.pdf")
panther %>% 
  ggplot(aes(x=term, y=log2_fold, fill=domain)) + 
  geom_bar(stat="identity") + 
  labs(y="log2-fold change") +
  scale_fill_discrete(name = "GO domain", labels = c("Cellular component", "Molecular function")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size=6),
    plot.margin = unit(c(0.2, 0.2, 0.2, 2), units="cm"),
    legend.position = c(0.85, 0.1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 45)) +
  geom_hline(yintercept=0)
dev.off()

# GO_id_levels <- factor(panther$GO_id, ordered=TRUE)
# panther$GO_id <- factor(panther$GO_id, levels=GO_id_levels)

# Make bar chart (ordered on log2-fold, coloured on GO domain) and write to file
# pdf(file = "../results/figures/panther_ORA_complete_overlap_bpccmf.pdf")
# panther %>% 
#   ggplot(aes(x=GO_id, y=log2_fold, fill=domain)) + 
#   geom_bar(stat="identity") + 
#   labs(y="log2-fold change") +
#   scale_fill_discrete(name = "GO domain", labels = c("Biological process", "Cellular component", "Molecular function")) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1, size=6),
#     plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), units="cm"),
#     legend.position = c(0.85, 0.1)) +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
#   geom_hline(yintercept=0)
# dev.off()
