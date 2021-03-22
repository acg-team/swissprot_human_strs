library(tidyverse)

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
    cg <- c(0, 1, 2, 4, 7, 13, 24, 44, 84, 161, 309, 594, 1164, 2284, 4484, 8807, 17305, 34301, 68008, 134852)
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
  cg_sum = sapply(msa_original, str_purity, return_val="cg"),
  frac_dominant = sapply(msa_original, str_purity, return_val="frac"),
  longest_stretch = sapply(msa_original, str_purity, return_val="longest")
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
  xlab("STR length") +
  ylab("Count") +
  labs(colour = "Repeat unit length") +
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
  xlab("Repeat region length") +
  ylab("log(Count)") +
  facet_wrap(~dominant_unit, nrow = 5)

# For the 4 most abundant combinations of AAs to make up STR with unit length 2,
# plot the frequency of the two possible orders (e.g. AP vs PA)
library(gridExtra)
to_keep <- c(8820, 13291, 1164, 4497) # will select only STRs made up of {A, P}, {G, R}, {G, S} or {R, S}
x <- master_df %>% filter(l_effective == 2, cg_sum %in% to_keep) %>%
  group_by(cg_sum, dominant_unit) %>%
  summarise(count=n())

out <- by(data = x, INDICES = x$cg_sum, FUN = function(m) {
  m <- droplevels(m)
  m <- ggplot(m, aes(as.factor(cg_sum), count, group=1, fill = dominant_unit)) + 
    geom_bar(stat="identity", position="dodge2") + theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ ylab("Count") +
    labs(fill = "Dominant\nunit") +
    theme(text = element_text(size=20))
})
do.call(grid.arrange, c(out2, ncol=1))
