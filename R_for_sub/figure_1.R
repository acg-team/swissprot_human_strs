library(tidyverse)
library(gridExtra)
cbPalette <- c("#D55E00", "#56B4E9", "#F0E442", "#000000", "#CC79A7", "#009E73", "#0072B2", "#E69F00")

str_purity <- function(msa, return_val){
  msa <- as.character(msa)
  frequencies <- as.data.frame(table(strsplit(msa, ",", fixed=TRUE)))
  frequencies <- frequencies %>% mutate(fraction = Freq/sum(Freq)) %>% arrange(desc(fraction))
  dom <- as.character(frequencies[1, "Var1"])
  # return fraction of the TR that is made up of the dominant unit
  if(return_val == "frac"){
    return(frequencies[1, "fraction"])
  } 
  # return the dominant unit
  else if(return_val == "dom"){
    dom <- str_remove_all(dom, "-")
    return(dom)
  } 
  else if(return_val == "cg"){
    aa = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
    # here, cg refers a set of integers with unique subset sums generated using the Conway-Guy sequence
    cg <- c(132568, 199412, 233119, 250115, 258613, 262936, 265136, 266256, 266826, 267111, 267259, 267336, 267376, 267396, 267407, 267413, 267416, 267418, 267419, 267420)
    cg_values <- setNames(as.list(cg), aa)
    dom <- str_remove_all(dom, "-")
    cg_sum <- 0
    for(i in strsplit(dom, "")[[1]]){
      cg_sum <- cg_sum + cg_values[[i]]
    }
    return(cg_sum)
  }
  # return the longest uninterrupted stretch of the dominant unit
  else if(return_val == "longest"){
    max_count <- 0
    cur_count <- 0
    for(i in strsplit(msa, ",")[[1]]){
      if(i == dom){
        cur_count = cur_count + 1
      } else {
        if(cur_count > max_count){
          max_count = cur_count
        }
        cur_count = 0
      }
      if(cur_count > max_count){
        max_count = cur_count
      }
    }
    return(max_count)
  } 
  # wrong return_val specified
  else{
    stop(sprintf("in str_purity(msa, return_val=\"%s\"): unsupported return_val specified", return_val), call. = FALSE)
  }
}

master_df <- read.table("../data_for_sub/str_sp_final.tsv", header=TRUE, sep="\t")

master_df <- master_df %>% mutate(
  dominant_unit = sapply(msa_original, str_purity, return_val="dom"),
  cg_sum = sapply(msa_original, str_purity, return_val="cg")
  # frac_dominant = sapply(msa_original, str_purity, return_val="frac"),
  # longest_stretch = sapply(msa_original, str_purity, return_val="longest")
)

# Counts of how many STRs were found per region length <= 20, split on unit length
master_df %>% 
  filter(repeat_region_length <= 20) %>% 
  group_by(l_effective, repeat_region_length) %>% count() %>%
  ggplot(aes(x=repeat_region_length, y=n, colour=as.factor(l_effective))) +
  geom_point(size=4) +
  geom_line(size=1.5) +
  theme_classic() +
  theme(text = element_text(size=20)) +
  xlab("STR length (amino acids)") +
  ylab("Count") +
  labs(colour = "Repeat unit length") +
  scale_colour_manual(values=c(cbPalette[1], cbPalette[2])) +
  theme(legend.position = c(0.8, 0.8), legend.background = element_rect(linetype = "solid", colour="black"))

# Investigate the length distributios for homorepeats, keeping in mind the 
# AA disorder propensity

# Mapping disorder propensity to amino acids, source of data: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2676888/
aa_df <- data.frame(
  aa = c("W", "F", "Y",	"I", "M", "L", "V", "N", "C", "T", "A", "G", "R", "D", "H", "Q", "K", "S", "E", "P"), 
  propensity = c(rep("order", 10), rep("disorder", 10)),
  score = c(-0.884, -0.697, -0.510, -0.486, -0.397, -0.326, -0.121, 0.007, 0.02, 0.059, 0.06, 0.166, 0.180, 0.192, 0.303, 0.318, 0.586, 0.341, 0.736, 0.987)
)

homo_lengths <- master_df %>% filter(l_effective == 1) %>% 
  group_by(dominant_unit, repeat_region_length) %>% count()
homo_lengths$dominant_unit <- factor(homo_lengths$dominant_unit, levels = aa_df$aa)

# Plot homorepeat regions lengths per AA, do not include AA's with < 10 occurrences
homo_lengths %>% 
  filter(!dominant_unit %in% c("Y", "N", "M", "V", "C"), repeat_region_length <= 50) %>% 
  ggplot(aes(x=repeat_region_length, y=log(n))) +
  geom_point(aes()) +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab("STR length (amino acids)") +
  ylab("log(Count)") +
  facet_wrap(~dominant_unit, nrow = 5)

#############################
### SUPPLEMENTAL ANALYSES ###
#############################

# cbPalette <- c("#D55E00", "#56B4E9", "#F0E442", "#000000", "#CC79A7", "#009E73", "#0072B2", "#E69F00")

# Counts obtained from counting dipeptide occurrences in non-repeating sequence in SwissProt proteins
background_counts_df <- data.frame("dominant_unit"=c('PS', 'SP', 'ED', 'DE', 'GA', 'AG', 'KE', 'EK', 'GS', 'SG', 'ER', 'RE', 'AP', 'PA', 'GP', 'PG', 'RS', 'SR', 'GR', 'RG'),
                                "bg_count"=c(65013, 67656, 49441, 37708, 52947, 55680, 56844, 62152, 66492, 64524, 45113, 45048, 50250, 57180, 51094, 62618, 47512, 51326, 42550, 41482))

master_df %>% filter(l_effective ==2) %>% 
  group_by(cg_sum, dominant_unit) %>% 
  count() %>% View()
master_df %>% filter(l_effective ==2) %>% 
  group_by(cg_sum) %>% 
  count() %>% View()

# add dipeptide counts from the non-repeat background sequence
counts_comparison_df <- master_df %>% filter(l_effective ==2) %>% 
  group_by(cg_sum, dominant_unit) %>% summarise(str_count=n()) %>%
  left_join(background_counts_df, by="dominant_unit") %>% 
  filter(!is.na(bg_count)) %>% 
  pivot_longer(cols=c(str_count, bg_count), names_to="count")

# calculate percentages to facilitate plotting
get_sums <- counts_comparison_df %>% group_by(cg_sum, count) %>% summarize(sum=sum(value))
counts_comparison_df <- counts_comparison_df %>% 
  left_join(get_sums, by=c("cg_sum", "count")) %>% 
  mutate(percentage=(value/sum)*100)

# plots
# remove combinations observed < 30 times in all STRs
counts_comparison_df <- counts_comparison_df %>% filter(!cg_sum %in% c(516941, 517522))
out2 <- by(data=counts_comparison_df, INDICES=counts_comparison_df$cg_sum, FUN=function(m){
  ggplot(m, aes(count, percentage, fill=dominant_unit)) +
  geom_bar(stat="identity", position="dodge2") + theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.ticks.x=element_blank())+ ylab("Percentage") +
  scale_x_discrete(labels=c("bg_count" = "non-STR", "str_count" = "STRs")) +
  scale_fill_manual(values=cbPalette) +
  labs(fill = "Dominant\nunit") 
})
do.call(grid.arrange, c(out2, ncol=3))

# Fisher's exact tests to determine if distribution over the dipeptides is different
## in the STRs compared to non-repeating sequence
by(data=counts_comparison_df, INDICES=counts_comparison_df$cg_sum, FUN=function(m){
  # select relevant columns
  m <- m %>% ungroup() %>% select(dominant_unit, count, value)
  # spread columns, necessary for subsequent test
  m <- m %>% pivot_wider(id_cols=dominant_unit, names_from=count, values_from=value)
  # prepare df for Fisher's exact test
  m <- data.frame(m[,-1], row.names=m$dominant_unit)
  fisher.test(m) # Fisher's exact test
})
