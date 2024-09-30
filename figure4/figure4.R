library(tidyverse)
library(this.path)
library(ggsci)
library(lsa)
library(ggpubr)
setwd(this.dir())
source("../bin/theme_setup.R")

###Figure4A
#The distance in mutational spectra (defined as the 1- Cosine similarity) between spectra from 
#the earliest time (the first ten months) and later times (sliding window) for the SARS-CoV-2.
df <- read_table("../data/sars_coding_mutation.txt") %>% filter(sub_type != "sub_type")
meta_df <- read_csv("../data/all_meta_corrected_host.csv")
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
      group_by(mutation_type) %>%summarise(late=n())
    df_late_syn <- mutation_df %>% full_join(df_late_syn) %>% mutate_all(~replace(., is.na(.), 0)) %>%
      rowwise() %>% mutate(late=late/sum(df_late_syn$late))
    df_late_nonsyn <- data_df_date %>% filter(months >= window_start & months <= window_end  & sub_type=="non_syn") %>% 
      group_by(mutation_type) %>% summarise(late= n())
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
model2 <- lm(data=time_distance_df, syn_distance ~ average_time + 0)
# compare slopes using t-test
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
ggsave("Figure4A.pdf",unit="cm", width=6, height=6,dpi = 300)


###figure4D
#The distance in mutational spectra (defined as the 1- Cosine similarity) between spectra from 
#the earliest time (the first ten months) and later times (sliding window) for the H3N2.
df <- read_table("../data/h3n2_coding_mutation.txt")
standard_mutation <- c('CUAA', 'CUAC', 'CUAG', 'CUAU', 'CUCA', 'CUCC', 'CUCG', 'CUCU', 'CUGA', 'CUGC', 'CUGG', 'CUGU', 'CUUA', 'CUUC', 'CUUG', 'CUUU', 'GAAA', 'GAAC', 'GAAG', 'GAAU', 'GACA', 'GACC', 'GACG', 'GACU', 'GAGA', 'GAGC', 'GAGG', 'GAGU', 'GAUA', 'GAUC', 'GAUG', 'GAUU', 'AGAA', 'AGAC', 'AGAG', 'AGAU', 'AGCA', 'AGCC', 'AGCG', 'AGCU', 'AGGA', 'AGGC', 'AGGG', 'AGGU', 'AGUA', 'AGUC', 'AGUG', 'AGUU', 'UCAA', 'UCAC', 'UCAG', 'UCAU', 'UCCA', 'UCCC', 'UCCG', 'UCCU', 'UCGA', 'UCGC', 'UCGG', 'UCGU', 'UCUA', 'UCUC', 'UCUG', 'UCUU', 'GUAA', 'GUAC', 'GUAG', 'GUAU', 'GUCA', 'GUCC', 'GUCG', 'GUCU', 'GUGA', 'GUGC', 'GUGG', 'GUGU', 'GUUA', 'GUUC', 'GUUG', 'GUUU', 'CAAA', 'CAAC', 'CAAG', 'CAAU', 'CACA', 'CACC', 'CACG', 'CACU', 'CAGA', 'CAGC', 'CAGG', 'CAGU', 'CAUA', 'CAUC', 'CAUG', 'CAUU', 'GCAA', 'GCAC', 'GCAG', 'GCAU', 'GCCA', 'GCCC', 'GCCG', 'GCCU', 'GCGA', 'GCGC', 'GCGG', 'GCGU', 'GCUA', 'GCUC', 'GCUG', 'GCUU', 'CGAA', 'CGAC', 'CGAG', 'CGAU', 'CGCA', 'CGCC', 'CGCG', 'CGCU', 'CGGA', 'CGGC', 'CGGG', 'CGGU', 'CGUA', 'CGUC', 'CGUG', 'CGUU', 'ACAA', 'ACAC', 'ACAG', 'ACAU', 'ACCA', 'ACCC', 'ACCG', 'ACCU', 'ACGA', 'ACGC', 'ACGG', 'ACGU', 'ACUA', 'ACUC', 'ACUG', 'ACUU', 'UGAA', 'UGAC', 'UGAG', 'UGAU', 'UGCA', 'UGCC', 'UGCG', 'UGCU', 'UGGA', 'UGGC', 'UGGG', 'UGGU', 'UGUA', 'UGUC', 'UGUG', 'UGUU', 'AUAA', 'AUAC', 'AUAG', 'AUAU', 'AUCA', 'AUCC', 'AUCG', 'AUCU', 'AUGA', 'AUGC', 'AUGG', 'AUGU', 'AUUA', 'AUUC', 'AUUG', 'AUUU', 'UAAA', 'UAAC', 'UAAG', 'UAAU', 'UACA', 'UACC', 'UACG', 'UACU', 'UAGA', 'UAGC', 'UAGG', 'UAGU', 'UAUA', 'UAUC', 'UAUG', 'UAUU')
complement_mutation <- c('GAUU', 'GAGU', 'GACU', 'GAAU', 'GAUG', 'GAGG', 'GACG', 'GAAG', 'GAUC', 'GAGC', 'GACC', 'GAAC', 'GAUA', 'GAGA', 'GACA', 'GAAA', 'CUUU', 'CUGU', 'CUCU', 'CUAU', 'CUUG', 'CUGG', 'CUCG', 'CUAG', 'CUUC', 'CUGC', 'CUCC', 'CUAC', 'CUUA', 'CUGA', 'CUCA', 'CUAA', 'UCUU', 'UCGU', 'UCCU', 'UCAU', 'UCUG', 'UCGG', 'UCCG', 'UCAG', 'UCUC', 'UCGC', 'UCCC', 'UCAC', 'UCUA', 'UCGA', 'UCCA', 'UCAA', 'AGUU', 'AGGU', 'AGCU', 'AGAU', 'AGUG', 'AGGG', 'AGCG', 'AGAG', 'AGUC', 'AGGC', 'AGCC', 'AGAC', 'AGUA', 'AGGA', 'AGCA', 'AGAA', 'CAUU', 'CAGU', 'CACU', 'CAAU', 'CAUG', 'CAGG', 'CACG', 'CAAG', 'CAUC', 'CAGC', 'CACC', 'CAAC', 'CAUA', 'CAGA', 'CACA', 'CAAA', 'GUUU', 'GUGU', 'GUCU', 'GUAU', 'GUUG', 'GUGG', 'GUCG', 'GUAG', 'GUUC', 'GUGC', 'GUCC', 'GUAC', 'GUUA', 'GUGA', 'GUCA', 'GUAA', 'CGUU', 'CGGU', 'CGCU', 'CGAU', 'CGUG', 'CGGG', 'CGCG', 'CGAG', 'CGUC', 'CGGC', 'CGCC', 'CGAC', 'CGUA', 'CGGA', 'CGCA', 'CGAA', 'GCUU', 'GCGU', 'GCCU', 'GCAU', 'GCUG', 'GCGG', 'GCCG', 'GCAG', 'GCUC', 'GCGC', 'GCCC', 'GCAC', 'GCUA', 'GCGA', 'GCCA', 'GCAA', 'UGUU', 'UGGU', 'UGCU', 'UGAU', 'UGUG', 'UGGG', 'UGCG', 'UGAG', 'UGUC', 'UGGC', 'UGCC', 'UGAC', 'UGUA', 'UGGA', 'UGCA', 'UGAA', 'ACUU', 'ACGU', 'ACCU', 'ACAU', 'ACUG', 'ACGG', 'ACCG', 'ACAG', 'ACUC', 'ACGC', 'ACCC', 'ACAC', 'ACUA', 'ACGA', 'ACCA', 'ACAA', 'UAUU', 'UAGU', 'UACU', 'UAAU', 'UAUG', 'UAGG', 'UACG', 'UAAG', 'UAUC', 'UAGC', 'UACC', 'UAAC', 'UAUA', 'UAGA', 'UACA', 'UAAA', 'AUUU', 'AUGU', 'AUCU', 'AUAU', 'AUUG', 'AUGG', 'AUCG', 'AUAG', 'AUUC', 'AUGC', 'AUCC', 'AUAC', 'AUUA', 'AUGA', 'AUCA', 'AUAA')
mutation_dict <- setNames(complement_mutation, standard_mutation)
df <- df %>%  rowwise() %>% mutate(mutation_type=mutation_dict[[mutation_type]])
sig_df <- df
meta_df <- read_csv("../data/all_meta_corrected_host.csv") %>% select(accession, year) 
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
                    group_by(mutation_type) %>%summarise(late=n())
    df_late_syn <- mutation_df %>% full_join(df_late_syn) %>% mutate_all(~replace(., is.na(.), 0)) %>%
                    rowwise() %>% mutate(late=late/sum(df_late_syn$late))
    df_late_nonsyn <- data_df_date %>% filter(year >= window_start & year <= window_end  & sub_type=="non_syn") %>% 
                        group_by(mutation_type) %>% summarise(late= n())
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
time_distance_df <- slide_spectrum(start=2000 , length=max(data_df_date$year),  
                                   stride=2.5, size = 5 ,distance_df=distance_df) %>% mutate(average_time = average_time-2000)
origin_df <- data.frame(average_time=0, nonsyn_distance=0, syn_distance=0)
time_distance_df <- rbind(time_distance_df, origin_df)
#compare slopes using t-test
model1 <- lm(data=time_distance_df, nonsyn_distance ~ average_time + 0)
model2 <- lm(data=time_distance_df, syn_distance ~ average_time + 0)
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
    geom_jitter(width = 0.1, height = 0.003) +  # 添加抖动效果
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
ggsave("Figure4D.pdf",unit="cm", width=6, height=6,dpi = 300)

###Figure4B 
#Numbers of nonsynonymous and synonymous mutations explained by two mutational signatures in SARS-CoV-2
df <- read_table("../data/sars_coding_mutation.txt") %>% filter(sub_type != "sub_type")
sig_df <- read_table("../data/CH192_De-Novo_Signatures_sars.txt")
#assign each mutation type to a signature
max_sig_df <- sig_df %>%
              rowwise() %>%
              mutate(max_column = names(select(., -MutationsType))[which.max(c_across(-MutationsType))]) %>%
              select(MutationsType,max_column)
max_df <- df %>% left_join(max_sig_df,by=c("mutation_type"= "MutationsType") ) 
group_df <- max_df %>% group_by(max_column,sub_type) %>% summarise(count=n())
all_group_meta_max_df <- max_df %>% group_by(max_column,sub_type) %>% summarise(count=n()) 
wider_all_group_meta_max_df <- all_group_meta_max_df %>% 
                                pivot_wider(
                                names_from = sub_type,
                                values_from = count) 
longer_sars_df <- wider_all_group_meta_max_df %>% pivot_longer(!max_column, names_to = "category", values_to = "count")
sars_exp_site <- read_tsv("../data/sars_exp_site.txt", col_names = F)[5,2][[1]]
sars_dnds <- wider_all_group_meta_max_df  %>% mutate(dnds = (`non_syn`/ syn)/sars_exp_site)
plot <- ggplot(longer_sars_df, aes(x = max_column, y = count, fill = category)) + 
  geom_col(position = "dodge") +
  scale_fill_npg(labels=c("Nonsynonymous","Synonymous"))+
  scale_x_discrete(labels=c("SARS-CoV-2\nsignatureA\n(early)","SARS-CoV-2\nsignatureB\n(late)"))+
  my_theme()+
  labs(y = "Count",
       fill="Mutation type")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.margin = margin(t = -0.5 , l= -1.2, b= 0.1, unit = "cm")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("Figure4B.pdf",unit="cm", width=4.5, height=6,dpi = 300)


####Figure4C
#dn/ds for signatures of SARS-CoV-2
plot <- ggplot(sars_dnds, aes(x = max_column, y = dnds)) + 
  geom_col(fill="lightblue", width = 0.7) +
  scale_x_discrete(labels=c("SARS-CoV-2\nsignatureA\n(early)","SARS-CoV-2\nsignatureB\n(late)"))+
  my_theme()+
  labs(y = "Dn/ds")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  annotate("text", x = 1.5, y = 0.95, label = "Neutral", size = 7/.pt)
ggsave("Figure4C.pdf",unit="cm", width=4, height=6,dpi = 300)

###Figure4E
#Numbers of nonsynonymous and synonymous mutations explained by two mutational signatures in H3N2
df <- read_table("../data/h3n2_coding_mutation.txt")
sig_df <- read_table("../data/CH192_De-Novo_Signatures_h3n2.txt")
standard_mutation <- c('CUAA', 'CUAC', 'CUAG', 'CUAU', 'CUCA', 'CUCC', 'CUCG', 'CUCU', 'CUGA', 'CUGC', 'CUGG', 'CUGU', 'CUUA', 'CUUC', 'CUUG', 'CUUU', 'GAAA', 'GAAC', 'GAAG', 'GAAU', 'GACA', 'GACC', 'GACG', 'GACU', 'GAGA', 'GAGC', 'GAGG', 'GAGU', 'GAUA', 'GAUC', 'GAUG', 'GAUU', 'AGAA', 'AGAC', 'AGAG', 'AGAU', 'AGCA', 'AGCC', 'AGCG', 'AGCU', 'AGGA', 'AGGC', 'AGGG', 'AGGU', 'AGUA', 'AGUC', 'AGUG', 'AGUU', 'UCAA', 'UCAC', 'UCAG', 'UCAU', 'UCCA', 'UCCC', 'UCCG', 'UCCU', 'UCGA', 'UCGC', 'UCGG', 'UCGU', 'UCUA', 'UCUC', 'UCUG', 'UCUU', 'GUAA', 'GUAC', 'GUAG', 'GUAU', 'GUCA', 'GUCC', 'GUCG', 'GUCU', 'GUGA', 'GUGC', 'GUGG', 'GUGU', 'GUUA', 'GUUC', 'GUUG', 'GUUU', 'CAAA', 'CAAC', 'CAAG', 'CAAU', 'CACA', 'CACC', 'CACG', 'CACU', 'CAGA', 'CAGC', 'CAGG', 'CAGU', 'CAUA', 'CAUC', 'CAUG', 'CAUU', 'GCAA', 'GCAC', 'GCAG', 'GCAU', 'GCCA', 'GCCC', 'GCCG', 'GCCU', 'GCGA', 'GCGC', 'GCGG', 'GCGU', 'GCUA', 'GCUC', 'GCUG', 'GCUU', 'CGAA', 'CGAC', 'CGAG', 'CGAU', 'CGCA', 'CGCC', 'CGCG', 'CGCU', 'CGGA', 'CGGC', 'CGGG', 'CGGU', 'CGUA', 'CGUC', 'CGUG', 'CGUU', 'ACAA', 'ACAC', 'ACAG', 'ACAU', 'ACCA', 'ACCC', 'ACCG', 'ACCU', 'ACGA', 'ACGC', 'ACGG', 'ACGU', 'ACUA', 'ACUC', 'ACUG', 'ACUU', 'UGAA', 'UGAC', 'UGAG', 'UGAU', 'UGCA', 'UGCC', 'UGCG', 'UGCU', 'UGGA', 'UGGC', 'UGGG', 'UGGU', 'UGUA', 'UGUC', 'UGUG', 'UGUU', 'AUAA', 'AUAC', 'AUAG', 'AUAU', 'AUCA', 'AUCC', 'AUCG', 'AUCU', 'AUGA', 'AUGC', 'AUGG', 'AUGU', 'AUUA', 'AUUC', 'AUUG', 'AUUU', 'UAAA', 'UAAC', 'UAAG', 'UAAU', 'UACA', 'UACC', 'UACG', 'UACU', 'UAGA', 'UAGC', 'UAGG', 'UAGU', 'UAUA', 'UAUC', 'UAUG', 'UAUU')
complement_mutation <- c('GAUU', 'GAGU', 'GACU', 'GAAU', 'GAUG', 'GAGG', 'GACG', 'GAAG', 'GAUC', 'GAGC', 'GACC', 'GAAC', 'GAUA', 'GAGA', 'GACA', 'GAAA', 'CUUU', 'CUGU', 'CUCU', 'CUAU', 'CUUG', 'CUGG', 'CUCG', 'CUAG', 'CUUC', 'CUGC', 'CUCC', 'CUAC', 'CUUA', 'CUGA', 'CUCA', 'CUAA', 'UCUU', 'UCGU', 'UCCU', 'UCAU', 'UCUG', 'UCGG', 'UCCG', 'UCAG', 'UCUC', 'UCGC', 'UCCC', 'UCAC', 'UCUA', 'UCGA', 'UCCA', 'UCAA', 'AGUU', 'AGGU', 'AGCU', 'AGAU', 'AGUG', 'AGGG', 'AGCG', 'AGAG', 'AGUC', 'AGGC', 'AGCC', 'AGAC', 'AGUA', 'AGGA', 'AGCA', 'AGAA', 'CAUU', 'CAGU', 'CACU', 'CAAU', 'CAUG', 'CAGG', 'CACG', 'CAAG', 'CAUC', 'CAGC', 'CACC', 'CAAC', 'CAUA', 'CAGA', 'CACA', 'CAAA', 'GUUU', 'GUGU', 'GUCU', 'GUAU', 'GUUG', 'GUGG', 'GUCG', 'GUAG', 'GUUC', 'GUGC', 'GUCC', 'GUAC', 'GUUA', 'GUGA', 'GUCA', 'GUAA', 'CGUU', 'CGGU', 'CGCU', 'CGAU', 'CGUG', 'CGGG', 'CGCG', 'CGAG', 'CGUC', 'CGGC', 'CGCC', 'CGAC', 'CGUA', 'CGGA', 'CGCA', 'CGAA', 'GCUU', 'GCGU', 'GCCU', 'GCAU', 'GCUG', 'GCGG', 'GCCG', 'GCAG', 'GCUC', 'GCGC', 'GCCC', 'GCAC', 'GCUA', 'GCGA', 'GCCA', 'GCAA', 'UGUU', 'UGGU', 'UGCU', 'UGAU', 'UGUG', 'UGGG', 'UGCG', 'UGAG', 'UGUC', 'UGGC', 'UGCC', 'UGAC', 'UGUA', 'UGGA', 'UGCA', 'UGAA', 'ACUU', 'ACGU', 'ACCU', 'ACAU', 'ACUG', 'ACGG', 'ACCG', 'ACAG', 'ACUC', 'ACGC', 'ACCC', 'ACAC', 'ACUA', 'ACGA', 'ACCA', 'ACAA', 'UAUU', 'UAGU', 'UACU', 'UAAU', 'UAUG', 'UAGG', 'UACG', 'UAAG', 'UAUC', 'UAGC', 'UACC', 'UAAC', 'UAUA', 'UAGA', 'UACA', 'UAAA', 'AUUU', 'AUGU', 'AUCU', 'AUAU', 'AUUG', 'AUGG', 'AUCG', 'AUAG', 'AUUC', 'AUGC', 'AUCC', 'AUAC', 'AUUA', 'AUGA', 'AUCA', 'AUAA')
#convert mutation type in negative strand to mutation type in positive strand
mutation_dict <- setNames(complement_mutation, standard_mutation)
df <- df %>%  rowwise() %>% mutate(mutation_type=mutation_dict[[mutation_type]])
#assign each mutation type to a signature
max_sig_df <- sig_df %>%
              rowwise() %>%
              mutate(max_column = names(select(., -MutationsType))[which.max(c_across(-MutationsType))]) %>%
              select(MutationsType,max_column)
max_df <- df %>% left_join(max_sig_df,by=c("mutation_type"= "MutationsType") ) 
group_df <- max_df %>% group_by(max_column,sub_type) %>% summarise(count=n())
wider_all_group_meta_max_df <- group_df %>% pivot_wider(
                                names_from = sub_type,
                                values_from = count)
longer_h3n2_df <- wider_all_group_meta_max_df %>% pivot_longer(!max_column, names_to = "category", values_to = "count")
h3n2_exp_site <- read_tsv("../data/h3n2_exp_site.txt", col_names = F)[5,2][[1]]
h3n2_dnds <- wider_all_group_meta_max_df %>% mutate(dnds = (`non_syn`/ syn)/h3n2_exp_site)
plot <- ggplot(longer_h3n2_df, aes(x = max_column, y = count, fill = category)) + 
      geom_col(position = "dodge") +
      scale_fill_npg(labels=c("Nonsynonymous","Synonymous"))+
      scale_x_discrete(labels=c("H3N2   \nsignatureA","H3N2   \nsignatureB\n(late)","H3N2   \nsignatureC\n(early)"))+
      my_theme()+
      labs(y = "Count",
           fill="Mutation type")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            axis.title.x = element_blank(),
            legend.title = element_blank(),
            legend.margin = margin(t = -0.5 , l= -1.2, b= 0.1, unit = "cm")) +  
      guides(fill = guide_legend(nrow = 1)) 
ggsave("Figure4E.pdf",unit="cm", width=5, height=6,dpi = 300)

###Figure4F
#dn/ds for signatures in H3N2
plot <- ggplot(h3n2_dnds, aes(x = max_column, y = dnds)) + 
        geom_col(fill="lightblue",width = 0.7) +
        scale_fill_npg(labels=c("Nonsynonymous","Synonymous"))+
        scale_x_discrete(labels=c("H3N2   \nsignatureA","H3N2   \nsignatureB\n(late)","H3N2   \nsignatureC\n(early)"))+
        my_theme()+
        labs(y = "Dn/ds")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank())+
        geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
        annotate("text", x = 2, y = 0.95, label = "Neutral", size = 7/.pt)
ggsave("Figure4F.pdf",unit="cm", width=5, height=6 ,dpi = 300)


###Figure4G
#Cosine similarity between the mutation spectra of iSNVs at different frequencies and the mutation spectrum derived from polymorphism data for SARS-CoV-2
low_frequency_isnv <- c("sars_isnv_0.2_0-count", "sars_isnv_0.4_0.2-count")
high_frequency_isnv <- c("sars_isnv_0.6_0.4-count", "sars_isnv_0.8_0.6-count", "sars_isnv_0.95_0.8-count")
data_df <- read_csv("../data/all_pairwise_result1.csv") %>% rowwise() %>% 
            filter(!(spectrumA %in% low_frequency_isnv && spectrumB %in% low_frequency_isnv)) %>%
            filter(!(spectrumA %in% high_frequency_isnv && spectrumB %in% high_frequency_isnv)) %>% 
            ungroup()
data_df[data_df$correct_p<=0.001,"p.signif"]<-"***"
data_df[data_df$correct_p>0.001 & data_df$correct_p<=0.01,"p.signif"]<-"**"
data_df[data_df$correct_p>0.01 & data_df$correct_p<=0.05,"p.signif"]<-"*"
data_df[data_df$correct_p>0.05,"p.signif"]<-"NS"
snp_df1 <- data_df %>% filter(spectrumB=="sars_polymorphism_count")
snv_df1 <- data_df %>% filter(spectrumB!="sars_polymorphism_count") %>% mutate(y.position=seq(0.95, by = 0.07/2, length.out = 6))
data_df <- read_csv("../data/syn_pairwise_result2.csv") %>% rowwise() %>% 
          filter(!(spectrumA %in% low_frequency_isnv && spectrumB %in% low_frequency_isnv)) %>%
          filter(!(spectrumA %in% high_frequency_isnv && spectrumB %in% high_frequency_isnv)) %>% 
          ungroup()
data_df[data_df$correct_p<=0.001,"p.signif"]<-"***"
data_df[data_df$correct_p>0.001 & data_df$correct_p<=0.01,"p.signif"]<-"**"
data_df[data_df$correct_p>0.01 & data_df$correct_p<=0.05,"p.signif"]<-"*"
data_df[data_df$correct_p>0.05,"p.signif"]<-"NS"
snp_df2 <- data_df %>% filter(spectrumB=="sars_polymorphism_syn-count")
snv_df2 <- data_df %>% filter(!spectrumB=="sars_polymorphism_syn-count") %>% mutate(y.position=seq(0.95, by = 0.07/2, length.out = 6))
data_df <- read_csv("../data/nonsyn_pairwise_result3.csv") %>% rowwise() %>% 
            filter(!(spectrumA %in% low_frequency_isnv && spectrumB %in% low_frequency_isnv)) %>%
            filter(!(spectrumA %in% high_frequency_isnv && spectrumB %in% high_frequency_isnv)) %>% 
            ungroup()
data_df[data_df$correct_p<=0.001,"p.signif"]<-"***"
data_df[data_df$correct_p>0.001 & data_df$correct_p<=0.01,"p.signif"]<-"**"
data_df[data_df$correct_p>0.01 & data_df$correct_p<=0.05,"p.signif"]<-"*"
data_df[data_df$correct_p>0.05,"p.signif"]<-"NS"
snp_df3 <- data_df %>% filter(spectrumB=="sars_polymorphism_non_syn-count")
snv_df3 <- data_df %>% filter(!spectrumB=="sars_polymorphism_non_syn-count") %>% mutate(y.position=seq(0.95, by = 0.07/2, length.out = 6))
original_signature <- c("sars_isnv_0.2_0-count","sars_isnv_0.4_0.2-count","sars_isnv_0.6_0.4-count",
                        "sars_isnv_0.8_0.6-count","sars_isnv_0.95_0.8-count","sars_polymorphism_count",
                        "sars_polymorphism_non_syn-count","sars_polymorphism_syn-count") 
correct_signature <- c("iSNV(0.05-0.2)","iSNV(0.2-0.4)","iSNV(0.4-0.6)",
                       "iSNV(0.6-0.8)","iSNV(0.8-0.95)","polymorphism","polymorphism","polymorphism")
signature_dict <- setNames(correct_signature,original_signature)
all_snp_df <- rbind(snp_df1,snp_df2,snp_df3) %>% dplyr::rename("group1"="spectrumA","group2"="spectrumB")  %>% rowwise() %>%
  mutate(group1=signature_dict[[group1]])  %>% mutate(group2=signature_dict[[group2]]) %>% ungroup()
all_snv_df <- rbind(snv_df1,snv_df2,snv_df3) %>% dplyr::rename("group1"="spectrumA","group2"="spectrumB")  %>% rowwise() %>%
  mutate(group1=signature_dict[[group1]]) %>%  mutate(group2=signature_dict[[group2]]) %>% ungroup()
isnv_color <- pal_npg()(10)
isnv_color <- c("#8491B4FF","#3C5488FF","#F39B7FFF","#FF7F50","#E41A1C")
all_snv_df1 <- all_snv_df %>% select(group,group1,group2,p.signif,y.position) %>% 
  filter(!group2 %in% c("sars_polymorphism_count","sars_polymorphism_non_syn-count",
                        "sars_polymorphism_syn-count"))
p1 <- ggplot(all_snp_df, aes(x = group1, y = cosine_similarity, fill=group1)) +
      geom_col(color="black") +
      facet_wrap(.~group, scales = "free_y") +
      scale_fill_manual(values=isnv_color)+
      stat_pvalue_manual(all_snv_df1,
                         size=3,
                         label= "p.signif",
                         color = "midnightblue")+
      coord_cartesian(ylim = c(0.7, 1.15 )) +
      scale_y_continuous( breaks = seq(0.7, 1, by = 0.1))+
      labs(
        x = "Mutation type",
        y = "Pairwise cosine similarity\nbetween polymorhism and iSNV")+
      my_theme()+
      theme(
        # axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      )+
      labs(fill="iSNV frequency")
ggsave("Figure4G.pdf",width = 18, height = 6, units = "cm",dpi = 300)

###Figure4H-K
#prepare data
isnv_pair <- read_csv("../data/isnv_pair_ruan_andreas.csv")
cutoff = 0.08
high_nonsyn <- isnv_pair %>% filter(donor_alt_frequency >= cutoff, recipient_alt_frequency >=0.03) %>%
              filter(func == "A") %>% 
              mutate(
                all_mean_frequency = mean(donor_alt_frequency),
                all_mean_recepient_frequency = mean(recipient_alt_frequency),
                mutation_type = "High-frequency nonsynonymous")
high_syn <- isnv_pair %>% filter(donor_alt_frequency >= cutoff, recipient_alt_frequency >=0.03) %>%
            filter(func == "S") %>%  
            mutate(
              all_mean_frequency = mean(donor_alt_frequency),
              all_mean_recepient_frequency = mean(recipient_alt_frequency),
              mutation_type = "High-frequency synonymous")
low_nonsyn <- isnv_pair %>% filter(donor_alt_frequency < cutoff  & donor_alt_frequency >= 0.03  & recipient_alt_frequency >=0.03) %>% 
              filter(func == "A") %>% 
              mutate(
                all_mean_frequency = mean(donor_alt_frequency),
                all_mean_recepient_frequency = mean(recipient_alt_frequency),
                mutation_type = "low-frequency nonsynonymous")
low_syn <- isnv_pair %>% filter(donor_alt_frequency < cutoff  & donor_alt_frequency >= 0.03 & recipient_alt_frequency >=0.03 ) %>% 
            filter(func == "S") %>% 
            mutate(
              all_mean_frequency = mean(donor_alt_frequency),
              all_mean_recepient_frequency = mean(recipient_alt_frequency),
              mutation_type = "low-frequency synonymous")
df <- rbind(high_nonsyn, low_nonsyn, high_syn, low_syn)
df_long <- data.frame(
            category = rep(c("donor", "recipient"), each = nrow(df)),
            frequency = c(df$donor_alt_frequency, df$recipient_alt_frequency),
            mutatio_type = c(df$mutation_type, df$mutation_type),
            group = rep(1:nrow(df), times = 2))

###Figure4H
#Frequency changes of high-frequency (>=0.08 in donors) non-synonymous iSNVs for SARS-CoV-2 before and after transmission
df <- high_nonsyn %>% filter(donor_alt_frequency < 0.9)
df_long <- data.frame(
          category = rep(c("donor", "recipient"), each = nrow(df)),
          frequency = c(df$donor_alt_frequency, df$recipient_alt_frequency),
          group = rep(1:nrow(df), times = 2))
wilcox_test_result  <- wilcox.test(df$donor_alt_frequency, df$recipient_alt_frequency, paired = TRUE, alternative = "greater")
p <- ggplot(df_long, aes(x = category, y = frequency, group = group)) +
    geom_boxplot(aes(group = category, color=category), linewidth = 0.1, outlier.shape = NA)+
    scale_color_npg()+
    geom_point( aes(group = category, color=category), size = 0.2) +
    geom_line(aes(group = group), linewidth=0.05) +
    annotate("text",
             x = 0.5, y = max(df_long$frequency) + 0.01, 
             label = paste0("p-value : ", format(wilcox_test_result$p.value, digits = 3)), 
             hjust = 0, vjust = 1, size = 7 /.pt)+
    labs(
      y = "iSNV frequency",
      title = "High-frequency\nnonsynonymous"
    ) +
    my_theme()+
    theme(  
      axis.title.x = element_blank(),
      plot.title = element_text(size=7, margin=margin(t= 0.1, l =-1, unit="cm")),
      axis.text.x = element_text(margin=margin(b= 0.1, unit="cm")),
      legend.position = "none"
    )
ggsave("Figure4H.pdf",width = 4.5, height = 5, units = "cm",dpi = 300)

###Figure4I
#Frequency changes of low-frequency (0.03-0.08 in donors) non-synonymous iSNVs for SARS-CoV-2 before and after transmission
df <- low_nonsyn
df_long <- data.frame(
          category = rep(c("donor", "recipient"), each = nrow(df)),
          frequency = c(df$donor_alt_frequency, df$recipient_alt_frequency),
          group = rep(1:nrow(df), times = 2))
wilcox_test_result  <- wilcox.test(df$donor_alt_frequency, df$recipient_alt_frequency, paired = TRUE, alternative = "greater")
p <- ggplot(df_long, aes(x = category, y = frequency, group = group)) +
            geom_boxplot(aes(group = category, color=category), linewidth = 0.1, outlier.shape = NA)+
            scale_color_npg()+
            geom_point( aes(group = category, color=category), size = 0.2) +
            geom_line(aes(group = group), linewidth=0.05) +
            annotate("text",
                     x = 0.5, y = max(df_long$frequency) + 0.01, 
                     label = paste0("p-value : ", format(wilcox_test_result$p.value, digits = 3)), 
                     hjust = 0, vjust = 1, size = 7 /.pt)+
            labs(
              y = "iSNV frequency",
              title = "Low-frequency\nnonsynonymous"
            ) +
            my_theme()+
            theme(  
              axis.title.x = element_blank(),
              plot.title = element_text(margin=margin(t= 0.1, unit="cm")),
              axis.text.x = element_text(margin=margin(b= 0.1, unit="cm")),
              legend.position = "none")

ggsave("Figure4I.pdf",width = 4.5, height = 5, units = "cm",dpi = 300)


###Figure4J
#Frequency changes of High-frequency (>=0.08 in donors) synonymous iSNVs for SARS-CoV-2 before and after transmission
df <- high_syn
df_long <- data.frame(
          category = rep(c("donor", "recipient"), each = nrow(df)),
          frequency = c(df$donor_alt_frequency, df$recipient_alt_frequency),
          group = rep(1:nrow(df), times = 2))
wilcox_test_result  <- wilcox.test(df$donor_alt_frequency, df$recipient_alt_frequency, paired = TRUE, alternative = "greater")
p <- ggplot(df_long, aes(x = category, y = frequency, group = group)) +
      geom_boxplot(aes(group = category, color=category), linewidth = 0.1, outlier.shape = NA)+
      scale_color_npg()+
      geom_point( aes(group = category, color=category), size = 0.2) +
      geom_line(aes(group = group), linewidth=0.05) +
      annotate("text",
               x = 0.5, y = max(df_long$frequency) + 0.08, 
               label = paste0("p-value : ", format(wilcox_test_result$p.value, digits = 3)), 
               hjust = 0, vjust = 1, size = 7 /.pt)+
      labs(
        y = "iSNV frequency",
        title = "High-frequency\nsynonymous"
      ) +
      my_theme()+
      theme(  
        axis.title.x = element_blank(),
        plot.title = element_text(margin=margin(t= 0.1, unit="cm")),
        axis.text.x = element_text(margin=margin(b= 0.1, unit="cm")),
        legend.position = "none")
ggsave("Figure4J.pdf",width = 4.5, height = 5, units = "cm",dpi = 300)



###Figure4K
#Frequency changes of low-frequency (0.03-0.08 in donors) synonymous iSNVs for SARS-CoV-2 before and after transmission
df <- low_syn
df_long <- data.frame(
  category = rep(c("donor", "recipient"), each = nrow(df)),
  frequency = c(df$donor_alt_frequency, df$recipient_alt_frequency),
  group = rep(1:nrow(df), times = 2)
)
wilcox_test_result  <- wilcox.test(df$donor_alt_frequency, df$recipient_alt_frequency, paired = TRUE, alternative = "greater")
p <- ggplot(df_long, aes(x = category, y = frequency, group = group)) +
    geom_boxplot(aes(group = category, color=category), linewidth = 0.1, outlier.shape = NA)+
    scale_color_npg()+
    geom_point( aes(group = category, color=category), size = 0.2) +
    geom_line(aes(group = group), linewidth=0.05) +
    annotate("text",
             x = 0.5, y = max(df_long$frequency) + 0.01, 
             label = paste0("p-value : ", format(wilcox_test_result$p.value, digits = 3)), 
             hjust = 0, vjust = 1, size = 7 /.pt)+
    labs(
      y = "iSNV frequency",
      title = "Low-frequency\nsynonymous"
    ) +
    my_theme()+
    theme(  
      axis.title.x = element_blank(),
      plot.title = element_text(margin=margin(t= 0.1, unit="cm")),
      axis.text.x = element_text(margin=margin(b= 0.1, unit="cm")),
      legend.position = "none")
ggsave("Figure4K.pdf",width = 4.5, height = 5, units = "cm",dpi = 300)




