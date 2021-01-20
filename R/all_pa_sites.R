library(tidyverse)

base_bath <- "/Users/maxverbiest/PhD/projects/SP_CRC_Pathway_TRs/data/pathology_atlas/"
sites <- c("bladder",
           "crc",
           "liver",
           "prostate",
           "thyroid",
           "brain",
           "endometrium",
           "lung",
           "skin",
           "breast",
           "hnn",
           "ovary",
           "stomach",
           "cervix",
           "kidney",
           "pancreas",
           "testis")

all_sites_df <- data.frame()
for(i in sites){
  path <- paste0(base_bath, i, "/", i, "_uniprot_only.csv")
  current_site_df <- read.table(path, header=TRUE, sep=",", stringsAsFactors=FALSE)
  if(dim(all_sites_df)[[1]] == 0){
    all_sites_df <- current_site_df
  } else{
    all_sites_df <- all_sites_df %>% full_join(current_site_df, by = "ID")
  }
}

all_sites_df$generality <- rowSums(!is.na(all_sites_df[2:18]))

all_sites_df %>% select(pa_crc, generality) %>% filter(!is.na(pa_crc)) %>% group_by(pa_crc, generality) %>% count() %>% 
  ggplot(aes(x=generality, y=n, colour=pa_crc))  +
  geom_line()

str_df <- read.table("../results/overlap/local_disorder/the_good_one.tsv", header=TRUE, sep="\t")
all_sites_strs_df <- all_sites_df %>% left_join(str_df %>% select(ID) %>% distinct() %>% mutate(str="yes"), by="ID")
all_sites_strs_df %>% group_by(str) %>% summarise(mean = mean(generality))

idr_df <- read.table("../results/disorder/local_annotations/unique_idr_containing.tsv", sep="\t", header=FALSE)
colnames(idr_df) <- c("ID")
idr_df$idr <- "yes"
all_sites_strs_idrs_df <- all_sites_strs_df %>% left_join(idr_df, by="ID")

all_sites_strs_df$fav <- rowSums(all_sites_strs_df[2:18] == "pa_fav", na.rm = TRUE)
all_sites_strs_df$unfav <- rowSums(all_sites_strs_df[2:18] == "pa_unfav", na.rm = TRUE)
all_sites_strs_df <- all_sites_strs_df %>% mutate(directionality = fav - unfav)


overall_summary_df <- read.table("../results/overlap/local_disorder/str_idr_all_pa_summary.tsv", header=TRUE, sep="\t")
overall_summary_df <- data.frame()
for(i in seq(2, 18)){
  current_ids_df <- as.data.frame(all_sites_strs_idrs_df[i])
  
  current_str_df <- current_ids_df
  current_str_df$str <- all_sites_strs_idrs_df$str
  colnames(current_str_df) <- c("pa", "str")
  current_str_df <- current_str_df %>% mutate(feature = case_when(str == "yes" ~ "str")) %>% 
    select(-str)
  
  current_str_summary_df <- current_str_df %>% filter(!is.na(pa)) %>% group_by(pa, feature) %>% 
    count() %>% group_by(pa) %>% mutate(sum = sum(n), perc = n/sum*100) %>% 
    ungroup() %>% filter(!is.na(feature))
  
  current_str_summary_df <- current_str_summary_df %>% 
    add_row(pa = "no_pa", 
            feature="str", 
            n=2658-sum(current_str_summary_df$n), 
            sum=20394-sum(current_str_summary_df$sum)) %>% 
    mutate(perc = n/sum*100, site=sites[i-1])
  
  
  current_idr_df <- current_ids_df
  current_idr_df$idr <- all_sites_strs_idrs_df$idr
  colnames(current_idr_df) <- c("pa", "idr")
  current_idr_df <- current_idr_df %>% mutate(feature = case_when(idr == "yes" ~ "idr")) %>% 
    select(-idr)
  
  current_idr_summary_df <- current_idr_df %>% filter(!is.na(pa)) %>% group_by(pa, feature) %>% 
    count() %>% group_by(pa) %>% mutate(sum = sum(n)) %>% 
    ungroup() %>% filter(!is.na(feature))
  
  current_idr_summary_df <- current_idr_summary_df %>% 
    add_row(pa = "no_pa", 
          feature="idr", 
          n=11277-sum(current_idr_summary_df$n), 
          sum=20394-sum(current_idr_summary_df$sum)) %>% 
    mutate(perc = n/sum*100, site=sites[i-1])
  
  if(dim(overall_summary_df)[[1]] == 0){
    overall_summary_df <- rbind(current_str_summary_df, current_idr_summary_df)
  } else{
    overall_summary_df <- rbind(overall_summary_df, current_str_summary_df, current_idr_summary_df)
  }
}

overall_summary_df$feature <- factor(overall_summary_df$feature, levels=c("str", "idr"))
pdf("../results/figures/18012021_update/str_idr_per_site.pdf")
overall_summary_df %>% filter(feature == "idr") %>% 
  ggplot(aes(x=feature, y=perc, fill=pa)) +
  geom_bar(stat="identity", position = position_dodge(width=0.55), width = 0.5) +
  theme_classic() +
  facet_wrap(~site, ncol = 6, scale = "free_x")
dev.off()

x <- overall_summary_df %>% filter(feature == "str") %>% select(n)
x <- x[[1]]
y <- overall_summary_df %>% filter(feature == "idr") %>% select(n)
y <- y[[1]]
cor.test(x, y, method="pearson")

# statistical tests
stat_df <- data.frame(site="tmp", feature="tmp", comparison="tmp", pvalue=0.00)
for(current_site in sites){
  site_sum <- overall_summary_df %>% filter(site == current_site)
  for(current_feat in c("str")){ #c("idr", "str")){
    feat_sum <- site_sum %>% filter(feature == current_feat) %>% mutate(sum = sum - n)
    pval <- fisher.test(matrix(c(t(feat_sum[1, 3:4])[1:2], t(feat_sum[2, 3:4])[1:2]), nrow=2))$p.value
    stat_df <- stat_df %>% add_row(site=current_site, feature=current_feat, comparison="f_uf", pvalue=pval)
    
    pval <- fisher.test(matrix(c(t(feat_sum[1, 3:4])[1:2], t(feat_sum[3, 3:4])[1:2]), nrow=2))$p.value
    stat_df <- stat_df %>% add_row(site=current_site, feature=current_feat, comparison="f_n", pvalue=pval)

    pval <- fisher.test(matrix(c(t(feat_sum[2, 3:4])[1:2], t(feat_sum[3, 3:4])[1:2]), nrow=2))$p.value
    stat_df <- stat_df %>% add_row(site=current_site, feature=current_feat, comparison="uf_n", pvalue=pval)
  }
}
stat_df <- stat_df %>% filter(!site == "tmp")
# FDR
stat_df <- stat_df %>% arrange(pvalue) %>% mutate(rank=seq(1, dim(.)[1]), FDR=rank/dim(.)[1]*0.05)
stat_df %>% filter(pvalue <= FDR) %>% View()
stat_df %>% filter(pvalue <= FDR) %>% arrange(site) %>% View()
