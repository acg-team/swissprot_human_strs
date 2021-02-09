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

master_df <- read.table("../results/overlap/local_disorder/the_good_one.tsv", header=TRUE, sep="\t")

master_df <- master_df %>% mutate(
  dominant_unit = sapply(msa_original, str_purity, return_val="dom"),
  cg_sum = sapply(msa_original, str_purity, return_val="cg"),
  frac_dominant = sapply(msa_original, str_purity, return_val="frac"),
  longest_stretch = sapply(msa_original, str_purity, return_val="longest")
)

#
master_df %>% 
  filter(repeat_region_length <= 20) %>% 
  group_by(l_effective, repeat_region_length) %>% count() %>% 
  ggplot(aes(x=repeat_region_length, y=n, colour=as.factor(l_effective))) +
  # geom_bar(stat="identity") +
  geom_point(size=2) +
  geom_line() +
  theme_classic() +
  xlab("STR length") +
  ylab("Count") +
  labs(colour = "Repeat unit length") +
  theme(legend.position = c(0.8, 0.8), legend.background = element_rect(linetype = "solid", colour="black"))
  
# Number of homo repeats per length
master_df %>% filter(l_effective == 1) %>% group_by(dominant_unit) %>% count() %>% View()

homo_lengths <- master_df %>% filter(l_effective == 1) %>% 
  group_by(dominant_unit, repeat_region_length) %>% count()
homo_lengths$dominant_unit <- factor(homo_lengths$dominant_unit, levels = aa_df$aa)

homo_lengths %>% 
  filter(!dominant_unit %in% c("Y", "N", "M", "V", "C"), repeat_region_length <= 50) %>% 
  ggplot(aes(x=repeat_region_length, y=log(n))) +
  geom_point(aes()) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~dominant_unit, nrow = 5)

# Overview of counts
master_df %>% 
  filter(l_effective == 1) %>% 
  group_by(dominant_unit) %>% 
  summarise(
    total = n(),
    max_len = max(repeat_region_length),
    mean_len = mean(repeat_region_length),
    max_longest = max(longest_stretch),
    mean_longest = mean(longest_stretch),
    max_frac = max(frac_dominant),
    mean_frac = mean(frac_dominant)
  ) %>% 
  View()


# functions of l == 2 repeats
master_df %>% 
  filter(l_effective == 2) %>% 
  group_by(dominant_unit) %>% 
  count(name="STR frequency") %>% View()
  group_by(`STR frequency`) %>% count(name="Number of times observed") %>% 
  ggplot(aes(x=`STR frequency`, y=`Number of times observed`)) +
  geom_bar(stat="identity")


master_df %>% 
  filter(dominant_unit == "RS") %>% 
  select(ID, n_effective, begin) %>% 
  # arrange(n_effective) %>% 
  left_join(load_swissprot("../data/swissrepeats/swissprot_human.tsv"), by="ID") %>% 
  arrange(protein_name) %>% 
  View()


# to_keep <- master_df %>% filter(l_effective == 2) %>% 
#   select(cg_sum, dominant_unit) %>% group_by(cg_sum) %>% distinct() %>% 
#   count() %>% filter(n == 2) %>% ungroup() %>% select(cg_sum)
# to_keep <- as.vector(to_keep)[[1]]

master_df %>% filter(l_effective == 2) %>% select(cg_sum) %>% group_by(cg_sum) %>% count() %>% View()

to_keep <- c(8820, 13291, 1164, 4497)
x <- master_df %>% filter(l_effective == 2, cg_sum %in% to_keep) %>% 
  select(cg_sum, dominant_unit, n_effective) %>% 
  group_by(cg_sum, dominant_unit, n_effective) %>% count()

out <- by(data = x, INDICES = x$cg_sum, FUN = function(m) {
  m <- droplevels(m)
  m <- ggplot(m, aes(n_effective, n, group=1, fill = dominant_unit)) + 
    geom_bar(stat="identity") + theme_classic() + xlab("Number of units") + ylab("Count") +
    labs(fill = "Dominant unit") + theme(legend.position = c(0.8, 0.8)) + scale_x_continuous(seq(3, 11, 0.5))
})

x2 <- master_df %>% filter(l_effective == 2, cg_sum %in% to_keep) %>%
  group_by(cg_sum, dominant_unit) %>%
  summarise(count=n())
out2 <- by(data = x2, INDICES = x2$cg_sum, FUN = function(m) {
  m <- droplevels(m)
  m <- ggplot(m, aes(as.factor(cg_sum), count, group=1, fill = dominant_unit)) + 
    geom_bar(stat="identity", position="dodge2") + theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ ylab("Count") +
    labs(fill = "Dominant unit") 
})
do.call(grid.arrange, c(out2, ncol=1))


master_df %>% filter(l_effective == 2, cg_sum %in% to_keep) %>%
  group_by(cg_sum, dominant_unit) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=as.factor(cg_sum), y=count, fill=dominant_unit)) +
  geom_bar(stat="identity", position = position_dodge())
  


