library(tidyverse)
library(this.path)
library(ggsci)
library(lsa)

setwd(this.dir())
source("../../bin/theme_setup.R")
###supplementary Figure10A
#The distance in mutational spectra (defined as the 1- Cosine similarity) between spectra from 
#the earliest time (the first ten months) and later times (sliding window) for the SARS-CoV-2.
df <- read_table("../../data/sars_coding_mutation.txt") %>% filter(sub_type != "sub_type")
meta_df <- read_csv("../../data/all_meta_corrected_host.csv")
all_data_df <- meta_df %>% select(c("accession","date")) %>% 
  right_join(df,by=c("accession"="accession")) %>% na.omit() 
data_df_date <- all_data_df %>% mutate_at(vars("date"),list(date_obj = ~as.Date(., format = "%Y/%m/%d"))) %>%
  mutate(days = as.numeric(date_obj - as.Date("2019-12-01"))) %>%  mutate(months=floor(days/30)) %>% na.omit()
standard_mutation <- c('CUAA', 'CUAC', 'CUAG', 'CUAU', 'CUCA', 'CUCC', 'CUCG', 'CUCU', 'CUGA', 'CUGC', 'CUGG', 'CUGU', 'CUUA', 'CUUC', 'CUUG', 'CUUU', 'GAAA', 'GAAC', 'GAAG', 'GAAU', 'GACA', 'GACC', 'GACG', 'GACU', 'GAGA', 'GAGC', 'GAGG', 'GAGU', 'GAUA', 'GAUC', 'GAUG', 'GAUU', 'AGAA', 'AGAC', 'AGAG', 'AGAU', 'AGCA', 'AGCC', 'AGCG', 'AGCU', 'AGGA', 'AGGC', 'AGGG', 'AGGU', 'AGUA', 'AGUC', 'AGUG', 'AGUU', 'UCAA', 'UCAC', 'UCAG', 'UCAU', 'UCCA', 'UCCC', 'UCCG', 'UCCU', 'UCGA', 'UCGC', 'UCGG', 'UCGU', 'UCUA', 'UCUC', 'UCUG', 'UCUU', 'GUAA', 'GUAC', 'GUAG', 'GUAU', 'GUCA', 'GUCC', 'GUCG', 'GUCU', 'GUGA', 'GUGC', 'GUGG', 'GUGU', 'GUUA', 'GUUC', 'GUUG', 'GUUU', 'CAAA', 'CAAC', 'CAAG', 'CAAU', 'CACA', 'CACC', 'CACG', 'CACU', 'CAGA', 'CAGC', 'CAGG', 'CAGU', 'CAUA', 'CAUC', 'CAUG', 'CAUU', 'GCAA', 'GCAC', 'GCAG', 'GCAU', 'GCCA', 'GCCC', 'GCCG', 'GCCU', 'GCGA', 'GCGC', 'GCGG', 'GCGU', 'GCUA', 'GCUC', 'GCUG', 'GCUU', 'CGAA', 'CGAC', 'CGAG', 'CGAU', 'CGCA', 'CGCC', 'CGCG', 'CGCU', 'CGGA', 'CGGC', 'CGGG', 'CGGU', 'CGUA', 'CGUC', 'CGUG', 'CGUU', 'ACAA', 'ACAC', 'ACAG', 'ACAU', 'ACCA', 'ACCC', 'ACCG', 'ACCU', 'ACGA', 'ACGC', 'ACGG', 'ACGU', 'ACUA', 'ACUC', 'ACUG', 'ACUU', 'UGAA', 'UGAC', 'UGAG', 'UGAU', 'UGCA', 'UGCC', 'UGCG', 'UGCU', 'UGGA', 'UGGC', 'UGGG', 'UGGU', 'UGUA', 'UGUC', 'UGUG', 'UGUU', 'AUAA', 'AUAC', 'AUAG', 'AUAU', 'AUCA', 'AUCC', 'AUCG', 'AUCU', 'AUGA', 'AUGC', 'AUGG', 'AUGU', 'AUUA', 'AUUC', 'AUUG', 'AUUU', 'UAAA', 'UAAC', 'UAAG', 'UAAU', 'UACA', 'UACC', 'UACG', 'UACU', 'UAGA', 'UAGC', 'UAGG', 'UAGU', 'UAUA', 'UAUC', 'UAUG', 'UAUU')
mutation_df <- data.frame(mutation_type=standard_mutation)
df_early_syn <- data_df_date %>% filter(months <= 10 & sub_type=="syn") %>% group_by(mutation_type) %>%summarise(early= n())
df_early_syn <- mutation_df %>% full_join(df_early_syn) %>% mutate_all(~replace(., is.na(.), 0)) %>% rowwise() %>% mutate(early=early/sum(df_early_syn$early))
df_early_nonsyn <- data_df_date %>% filter(months <= 10 & sub_type=="non_syn") %>% group_by(mutation_type) %>% summarise(early= n())
df_early_nonsyn <- mutation_df %>% full_join(df_early_nonsyn) %>% mutate_all(~replace(., is.na(.), 0)) %>% rowwise() %>% mutate(early=early/sum(df_early_nonsyn$early))
distance_df <- data.frame()
slide_spectrum <- function(start, length, stride, size, distance_df){
  for( start in seq(as.numeric(start), as.numeric(length), by=stride)) {
  window_start <- start 
  window_end <- start + size-1
  if(window_end > length) {
    window_end <- length
  }
  average_time <- (window_start+window_end)/2
  df_late_syn <- data_df_date %>% filter(months >= window_start & months <= window_end  & sub_type=="syn") %>% 
    group_by(mutation_type) %>% summarise(late=n())
  syn_count <- sum(df_late_syn$late)
  df_late_nonsyn <- data_df_date %>% filter(months >= window_start & months <= window_end  & sub_type=="non_syn") %>% 
    group_by(mutation_type) %>% summarise(late= n())
  nonsyn_count <- sum(df_late_nonsyn$late)
  min_count <- min(syn_count, nonsyn_count)
  nonsyn_vector <- rep(df_late_nonsyn$mutation_type, df_late_nonsyn$late)
  syn_vector <- rep(df_late_syn$mutation_type, df_late_syn$late)
  #Downsampling ensures the same number of mutations for synonymous and non-synonymous mutational spectra
  set.seed(123) 
  downsample_nonsyn_vector <- sample(nonsyn_vector, size = min_count, replace = FALSE) %>% table()
  df_late_nonsyn <- df_late_nonsyn %>% rowwise() %>% 
    mutate(late = ifelse(mutation_type %in% downsample_nonsyn_vector, downsample_nonsyn_vector[[mutation_type]], late))
  downsample_syn_vector <- sample(syn_vector, size = min_count, replace = FALSE) %>% table()
  df_late_syn <- df_late_syn %>% rowwise() %>% 
    mutate(late = ifelse(mutation_type %in% downsample_syn_vector, downsample_syn_vector[[mutation_type]], late))
  df_late_syn <- mutation_df %>% full_join(df_late_syn) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    rowwise() %>% mutate(late=late/sum(df_late_syn$late))
  df_late_nonsyn <- mutation_df %>% full_join(df_late_nonsyn) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    rowwise() %>% mutate(late=late/sum(df_late_nonsyn$late))
  cosine_similarity_syn <- cosine(df_early_syn$early, df_late_syn$late)
  cosine_similarity_nonsyn <- cosine(df_early_nonsyn$early, df_late_nonsyn$late)
  syn_distance <- 1 - cosine_similarity_syn
  nonsyn_distance <- 1-cosine_similarity_nonsyn
  small_df <- data.frame(
   average_time,
   nonsyn_distance,
   syn_distance)
  distance_df <- rbind(distance_df, small_df )
  }
  return (distance_df)
}
time_distance_df <- slide_spectrum(start=min(data_df_date$months) , length=max(data_df_date$months),  
                                  stride=5, size = 10 ,distance_df=distance_df) 
origin_df <- data.frame(average_time=0, nonsyn_distance=0, syn_distance=0)
time_distance_df <- rbind(time_distance_df, origin_df)
model1 <- lm(data=time_distance_df, nonsyn_distance ~ average_time + 0)
summary(model1)
model2 <- lm(data=time_distance_df, syn_distance ~ average_time + 0)
summary(model2)
#compare two slopes
slope1 <- coef(summary(model1))["average_time", "Estimate"]
slope2 <- coef(summary(model2))["average_time", "Estimate"]
se1 <- coef(summary(model1))["average_time", "Std. Error"]
se2 <- coef(summary(model2))["average_time", "Std. Error"]
t_value <- (slope1 - slope2) / sqrt(se1^2 + se2^2)
df <- min(df.residual(model1), df.residual(model2))
p_value <- 1 - pt(t_value, df)
cat("t-value:", t_value, "\n")
cat("p-value:", p_value, "\n")
plot_time <- time_distance_df %>% pivot_longer(cols = starts_with(c("nonsyn","syn")),
                                               names_to = "mutation_type",
                                               values_to = "distance")
p <- ggplot(plot_time, aes(x = average_time, y = distance, color = mutation_type)) +
    geom_jitter(width = 0.1, height = 0.0008) + 
    geom_smooth(method = "lm", se = FALSE,formula = y ~ x + 0, aes(group = mutation_type), show.legend = FALSE) +
    scale_color_npg(labels = c("Nonsynonymous", "Synonymous")) +
    labs(
      y = "Pairwise distance between spectrums",
      x = "Month (since 2019/12)",
      title = "SARS-CoV-2",
      color = "Mutation type"
    ) +
    my_theme()+
    theme(
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.margin = margin(t = -0.3 , b= 0.1, unit = "cm")) +  
    guides(color = guide_legend(nrow = 1)) 
ggsave("suppl_figure10A.pdf",unit="cm", width=6, height=6,dpi = 300)

###supplementary figure10B
#The distance in mutational spectra (defined as the 1- Cosine similarity) between spectra from 
#the earliest time (the first ten months) and later times (sliding window) for H3N2
df <- read_table("../../data/h3n2_coding_mutation.txt")
#convert mutation spectra in genomic strand to positive strand 
standard_mutation <- c('CUAA', 'CUAC', 'CUAG', 'CUAU', 'CUCA', 'CUCC', 'CUCG', 'CUCU', 'CUGA', 'CUGC', 'CUGG', 'CUGU', 'CUUA', 'CUUC', 'CUUG', 'CUUU', 'GAAA', 'GAAC', 'GAAG', 'GAAU', 'GACA', 'GACC', 'GACG', 'GACU', 'GAGA', 'GAGC', 'GAGG', 'GAGU', 'GAUA', 'GAUC', 'GAUG', 'GAUU', 'AGAA', 'AGAC', 'AGAG', 'AGAU', 'AGCA', 'AGCC', 'AGCG', 'AGCU', 'AGGA', 'AGGC', 'AGGG', 'AGGU', 'AGUA', 'AGUC', 'AGUG', 'AGUU', 'UCAA', 'UCAC', 'UCAG', 'UCAU', 'UCCA', 'UCCC', 'UCCG', 'UCCU', 'UCGA', 'UCGC', 'UCGG', 'UCGU', 'UCUA', 'UCUC', 'UCUG', 'UCUU', 'GUAA', 'GUAC', 'GUAG', 'GUAU', 'GUCA', 'GUCC', 'GUCG', 'GUCU', 'GUGA', 'GUGC', 'GUGG', 'GUGU', 'GUUA', 'GUUC', 'GUUG', 'GUUU', 'CAAA', 'CAAC', 'CAAG', 'CAAU', 'CACA', 'CACC', 'CACG', 'CACU', 'CAGA', 'CAGC', 'CAGG', 'CAGU', 'CAUA', 'CAUC', 'CAUG', 'CAUU', 'GCAA', 'GCAC', 'GCAG', 'GCAU', 'GCCA', 'GCCC', 'GCCG', 'GCCU', 'GCGA', 'GCGC', 'GCGG', 'GCGU', 'GCUA', 'GCUC', 'GCUG', 'GCUU', 'CGAA', 'CGAC', 'CGAG', 'CGAU', 'CGCA', 'CGCC', 'CGCG', 'CGCU', 'CGGA', 'CGGC', 'CGGG', 'CGGU', 'CGUA', 'CGUC', 'CGUG', 'CGUU', 'ACAA', 'ACAC', 'ACAG', 'ACAU', 'ACCA', 'ACCC', 'ACCG', 'ACCU', 'ACGA', 'ACGC', 'ACGG', 'ACGU', 'ACUA', 'ACUC', 'ACUG', 'ACUU', 'UGAA', 'UGAC', 'UGAG', 'UGAU', 'UGCA', 'UGCC', 'UGCG', 'UGCU', 'UGGA', 'UGGC', 'UGGG', 'UGGU', 'UGUA', 'UGUC', 'UGUG', 'UGUU', 'AUAA', 'AUAC', 'AUAG', 'AUAU', 'AUCA', 'AUCC', 'AUCG', 'AUCU', 'AUGA', 'AUGC', 'AUGG', 'AUGU', 'AUUA', 'AUUC', 'AUUG', 'AUUU', 'UAAA', 'UAAC', 'UAAG', 'UAAU', 'UACA', 'UACC', 'UACG', 'UACU', 'UAGA', 'UAGC', 'UAGG', 'UAGU', 'UAUA', 'UAUC', 'UAUG', 'UAUU')
complement_mutation <- c('GAUU', 'GAGU', 'GACU', 'GAAU', 'GAUG', 'GAGG', 'GACG', 'GAAG', 'GAUC', 'GAGC', 'GACC', 'GAAC', 'GAUA', 'GAGA', 'GACA', 'GAAA', 'CUUU', 'CUGU', 'CUCU', 'CUAU', 'CUUG', 'CUGG', 'CUCG', 'CUAG', 'CUUC', 'CUGC', 'CUCC', 'CUAC', 'CUUA', 'CUGA', 'CUCA', 'CUAA', 'UCUU', 'UCGU', 'UCCU', 'UCAU', 'UCUG', 'UCGG', 'UCCG', 'UCAG', 'UCUC', 'UCGC', 'UCCC', 'UCAC', 'UCUA', 'UCGA', 'UCCA', 'UCAA', 'AGUU', 'AGGU', 'AGCU', 'AGAU', 'AGUG', 'AGGG', 'AGCG', 'AGAG', 'AGUC', 'AGGC', 'AGCC', 'AGAC', 'AGUA', 'AGGA', 'AGCA', 'AGAA', 'CAUU', 'CAGU', 'CACU', 'CAAU', 'CAUG', 'CAGG', 'CACG', 'CAAG', 'CAUC', 'CAGC', 'CACC', 'CAAC', 'CAUA', 'CAGA', 'CACA', 'CAAA', 'GUUU', 'GUGU', 'GUCU', 'GUAU', 'GUUG', 'GUGG', 'GUCG', 'GUAG', 'GUUC', 'GUGC', 'GUCC', 'GUAC', 'GUUA', 'GUGA', 'GUCA', 'GUAA', 'CGUU', 'CGGU', 'CGCU', 'CGAU', 'CGUG', 'CGGG', 'CGCG', 'CGAG', 'CGUC', 'CGGC', 'CGCC', 'CGAC', 'CGUA', 'CGGA', 'CGCA', 'CGAA', 'GCUU', 'GCGU', 'GCCU', 'GCAU', 'GCUG', 'GCGG', 'GCCG', 'GCAG', 'GCUC', 'GCGC', 'GCCC', 'GCAC', 'GCUA', 'GCGA', 'GCCA', 'GCAA', 'UGUU', 'UGGU', 'UGCU', 'UGAU', 'UGUG', 'UGGG', 'UGCG', 'UGAG', 'UGUC', 'UGGC', 'UGCC', 'UGAC', 'UGUA', 'UGGA', 'UGCA', 'UGAA', 'ACUU', 'ACGU', 'ACCU', 'ACAU', 'ACUG', 'ACGG', 'ACCG', 'ACAG', 'ACUC', 'ACGC', 'ACCC', 'ACAC', 'ACUA', 'ACGA', 'ACCA', 'ACAA', 'UAUU', 'UAGU', 'UACU', 'UAAU', 'UAUG', 'UAGG', 'UACG', 'UAAG', 'UAUC', 'UAGC', 'UACC', 'UAAC', 'UAUA', 'UAGA', 'UACA', 'UAAA', 'AUUU', 'AUGU', 'AUCU', 'AUAU', 'AUUG', 'AUGG', 'AUCG', 'AUAG', 'AUUC', 'AUGC', 'AUCC', 'AUAC', 'AUUA', 'AUGA', 'AUCA', 'AUAA')
mutation_dict <- setNames(complement_mutation, standard_mutation)
df <- df %>%  rowwise() %>% mutate(mutation_type=mutation_dict[[mutation_type]])
sig_df <- df
meta_df <- read_csv("../../data/all_meta_corrected_host.csv") %>% select(accession, year) 
data_df_date <- sig_df %>% left_join(meta_df,by=c("accession"="accession"))  %>% filter(!is.na(year))
mutation_df <- data.frame(mutation_type=standard_mutation)
df_early_syn <- data_df_date %>% filter(year <= 2004 & year >=2000 & sub_type=="syn") %>% group_by(mutation_type) %>% summarise(early= n())
df_early_syn <- mutation_df %>% full_join(df_early_syn) %>% mutate_all(~replace(., is.na(.), 0)) %>% 
  rowwise() %>% mutate(early=early/sum(df_early_syn$early))
df_early_nonsyn <- data_df_date %>% filter(year <= 2004 & year >=2000 & sub_type=="non_syn") %>% group_by(mutation_type) %>% summarise(early= n())
df_early_nonsyn <- mutation_df %>% full_join(df_early_nonsyn) %>% mutate_all(~replace(., is.na(.), 0)) %>%
  rowwise() %>% mutate(early=early/sum(df_early_nonsyn$early))
distance_df <- data.frame()
slide_spectrum <- function(start, length, stride, size, distance_df){
  for( start in seq(as.numeric(start), as.numeric(length), by=stride)) {
    window_start <- start 
    window_end <- start + size-1
    if(window_end > length) {
      window_end <- length
    }
    average_time <- (window_start+window_end)/2
    df_late_syn <- data_df_date %>% filter(year >= window_start & year <= window_end  & sub_type=="syn") %>% 
      group_by(mutation_type) %>% summarise(late=n())
    syn_count <- sum(df_late_syn$late)
    df_late_nonsyn <- data_df_date %>% filter(year >= window_start & year <= window_end  & sub_type=="non_syn") %>% 
      group_by(mutation_type) %>% summarise(late= n())
    nonsyn_count <- sum(df_late_nonsyn$late)
    min_count <- min(syn_count, nonsyn_count)
    nonsyn_vector <- rep(df_late_nonsyn$mutation_type, df_late_nonsyn$late)
    syn_vector <- rep(df_late_syn$mutation_type, df_late_syn$late)
    set.seed(123) 
    downsample_nonsyn_vector <- sample(nonsyn_vector, size = min_count, replace = FALSE) %>% table()
    df_late_nonsyn <- df_late_nonsyn %>% rowwise() %>% 
      mutate(late = ifelse(mutation_type %in% downsample_nonsyn_vector, downsample_nonsyn_vector[[mutation_type]], late))
    downsample_syn_vector <- sample(syn_vector, size = min_count, replace = FALSE) %>% table()
    df_late_syn <- df_late_syn %>% rowwise() %>% 
      mutate(late = ifelse(mutation_type %in% downsample_syn_vector, downsample_syn_vector[[mutation_type]], late))
    df_late_syn <- mutation_df %>% full_join(df_late_syn) %>% mutate_all(~replace(., is.na(.), 0)) %>%
      rowwise() %>% mutate(late=late/sum(df_late_syn$late))
    df_late_nonsyn <- mutation_df %>% full_join(df_late_nonsyn) %>% mutate_all(~replace(., is.na(.), 0)) %>%
      rowwise() %>% mutate(late=late/sum(df_late_nonsyn$late))
    cosine_similarity_syn <- cosine(df_early_syn$early, df_late_syn$late)
    cosine_similarity_nonsy <- cosine(df_early_nonsyn$early, df_late_nonsyn$late)
    syn_distance <- 1 - cosine_similarity_syn
    nonsyn_distance <- 1-cosine_similarity_nonsy
    small_df <- data.frame(
      average_time,
      nonsyn_distance,
      syn_distance)
    distance_df <- rbind(distance_df, small_df )
  }
  return (distance_df)
}

time_distance_df <- slide_spectrum(start=2000 , length=max(data_df_date$year),  
                                   stride=2.5, size = 5 ,distance_df=distance_df) %>% mutate(average_time = average_time-2000)
origin_df <- data.frame(average_time=0, nonsyn_distance=0, syn_distance=0)
time_distance_df <- rbind(time_distance_df, origin_df)
model1 <- lm(data=time_distance_df, nonsyn_distance ~ average_time + 0)
summary(model1)
model2 <- lm(data=time_distance_df, syn_distance ~ average_time + 0)
summary(model2)
#compare two slopes
slope1 <- coef(summary(model1))["average_time", "Estimate"]
slope2 <- coef(summary(model2))["average_time", "Estimate"]
se1 <- coef(summary(model1))["average_time", "Std. Error"]
se2 <- coef(summary(model2))["average_time", "Std. Error"]
t_value <- (slope1 - slope2) / sqrt(se1^2 + se2^2)
df <- min(df.residual(model1), df.residual(model2))
p_value <- 1 - pt(t_value, df)
cat("t-value:", t_value, "\n")
cat("p-value:", p_value, "\n")
plot_time <- time_distance_df %>% pivot_longer(cols = starts_with(c("nonsyn","syn")),
                                               names_to = "mutation_type",
                                               values_to = "distance")
p <- ggplot(plot_time, aes(x = average_time, y = distance, color = mutation_type)) +
    geom_jitter(width = 0.1, height = 0.003) + 
    geom_smooth(method = "lm", se = FALSE,formula = y ~ x + 0, aes(group = mutation_type), show.legend = FALSE) +
    scale_color_npg(labels=c("Nonsynonymous","Synonymous"))+
    labs(y="Pairwise distance between spectrums",
         x="Year (since 2000)",
         title="H3N2",
         color="Mutation type")+
    my_theme()+
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.margin = margin(t = -0.3 , b= 0.1, unit = "cm")) +  
    guides(color = guide_legend(nrow = 1)) 
ggsave("suppl_figure10B.pdf",unit="cm", width=6, height=6,dpi = 300)






