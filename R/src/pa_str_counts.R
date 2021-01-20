pa_str_counts <- function(df, all_fav, all_unfav, all_prots=20394, all_with_str=2658){
  df <- df %>% select(ID, pa) %>% distinct()
  # number of pa favourable proteins that contain an str
  str_fav <- df %>% select(ID, pa) %>% filter(pa=="pa_fav") %>% summarise(count=n())
  str_fav <- str_fav[[1]]
  # number of pa unfavourable proteins that contain an str
  str_unfav <- df %>% select(ID, pa) %>% filter(pa=="pa_unfav") %>% summarise(count=n())
  str_unfav <- str_unfav[[1]]
  # number of proteins not in current pa set
  all_prots_no_pa <- all_prots - all_fav - all_unfav
  # number of proteins not in pa, with and without str
  no_pa_with_str <- all_with_str - str_fav - str_unfav
  no_pa_no_str <- all_prots_no_pa - no_pa_with_str
  
  # Fisher exact tests, print both table and p-value
  non_vs_fav <- matrix(c(no_pa_no_str, no_pa_with_str, all_fav-str_fav, str_fav), nrow=2, dimnames=list(c("without", "with"), c("non_pa", "pa_fav")))
  print(non_vs_fav)
  print(fisher.test(non_vs_fav)$p.value)
  
  non_vs_unfav <- matrix(c(no_pa_no_str, no_pa_with_str, all_unfav-str_unfav, str_unfav), nrow=2, dimnames=list(c("without", "with"), c("non_pa", "pa_unfav")))
  print(non_vs_unfav)
  print(fisher.test(non_vs_unfav)$p.value)
  
  fav_vs_unfav <- matrix(c(all_fav-str_fav, str_fav, all_unfav-str_unfav, str_unfav), nrow=2, dimnames=list(c("without", "with"), c("pa_fav", "pa_unfav")))
  print(fav_vs_unfav)
  print(fisher.test(fav_vs_unfav)$p.value)
}
