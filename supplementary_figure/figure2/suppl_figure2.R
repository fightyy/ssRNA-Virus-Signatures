library(tidyverse)
library(this.path)
library(ggsci)
library(scales)
setwd(this.dir())
source("../../bin/theme_setup.R")

###supplementary 2A
###bar plot of the log ratio between early and late spectra for SARS-COV2
df <- read_tsv("../../data/all_signature_sample_csv.txt") %>% filter(! order %in% c("ssRNA_+_", "ssRNA_-_")) %>% 
      filter(host != "internal") %>% rowwise()%>% 
      mutate(start_year = suppressWarnings(as.numeric(strsplit(year, "-")[[1]][1])),
             end_year = suppressWarnings(as.numeric(strsplit(year, "-")[[1]][2]))) %>% na.omit() %>% 
      mutate(mean_between_years = mean(c(start_year, end_year))) %>% 
      select(c(start_year, end_year, mean_between_years), everything())
df_sars <- df %>% filter(species=="sars-cov-2") %>% arrange(mean_between_years) 
df_sars_early <- df_sars %>% filter(mean_between_years <= 10) %>%
              select((ncol(df_sars)-191):ncol(df_sars)) %>% colSums() %>% as.data.frame() 
df_sars_early <- df_sars %>% filter(mean_between_years <= 10) %>%
                select((ncol(df_sars)-191):ncol(df_sars)) %>% colSums() %>% as.data.frame()  %>% 
                mutate(across(starts_with("."), ~ . / sum(`.`))) %>%
                mutate(MutationType=rownames(df_sars_early)) %>% dplyr::rename("SARS-COV-2\n(early)"=".")
df_sars_late <- df_sars %>% filter(mean_between_years >= 38) %>%
                select((ncol(df_sars)-191):ncol(df_sars)) %>% colSums() %>% as.data.frame()  %>% 
                mutate(across(starts_with("."), ~ . / sum(`.`))) %>%
                mutate(MutationType=rownames(df_sars_early)) %>% dplyr::rename("SARS-COV-2\n(late)"=".")
df_sars_ratio <- df_sars_late %>% right_join(df_sars_early)  %>% 
                  mutate(log_ratio=log(`SARS-COV-2\n(late)`/`SARS-COV-2\n(early)`)) %>%
                  mutate(log_ratio=ifelse((`SARS-COV-2\n(late)`+`SARS-COV-2\n(early)`) <= 0.01,0,log_ratio)) %>%
                  filter(abs(log_ratio) >0 )
df_plot_sars <- df_sars_ratio %>% select("MutationType","log_ratio")
df_plot_sars$trinuc <- apply(df_plot_sars ,1,function(x){
  return(paste0(substring(x["MutationType"],3,3),
                "[",
                substring(x["MutationType"],1,1),
                ">",
                substring(x["MutationType"],2,2),
                "]",
                substring(x["MutationType"],4,4)))
})

trinuc_level <- c('A[C>U]A','A[C>U]C','A[C>U]G','A[C>U]U','C[C>U]A','C[C>U]C','C[C>U]G','C[C>U]U','G[C>U]A','G[C>U]C','G[C>U]G','G[C>U]U','U[C>U]A','U[C>U]C','U[C>U]G','U[C>U]U','A[G>A]A','A[G>A]C','A[G>A]G','A[G>A]U','C[G>A]A','C[G>A]C','C[G>A]G','C[G>A]U','G[G>A]A','G[G>A]C','G[G>A]G','G[G>A]U','U[G>A]A','U[G>A]C','U[G>A]G','U[G>A]U','A[A>G]A','A[A>G]C','A[A>G]G','A[A>G]U','C[A>G]A','C[A>G]C','C[A>G]G','C[A>G]U','G[A>G]A','G[A>G]C','G[A>G]G','G[A>G]U','U[A>G]A','U[A>G]C','U[A>G]G','U[A>G]U','A[U>C]A','A[U>C]C','A[U>C]G','A[U>C]U','C[U>C]A','C[U>C]C','C[U>C]G','C[U>C]U','G[U>C]A','G[U>C]C','G[U>C]G','G[U>C]U','U[U>C]A','U[U>C]C','U[U>C]G','U[U>C]U','A[G>U]A','A[G>U]C','A[G>U]G','A[G>U]U','C[G>U]A','C[G>U]C','C[G>U]G','C[G>U]U','G[G>U]A','G[G>U]C','G[G>U]G','G[G>U]U','U[G>U]A','U[G>U]C','U[G>U]G','U[G>U]U','A[C>A]A','A[C>A]C','A[C>A]G','A[C>A]U','C[C>A]A','C[C>A]C','C[C>A]G','C[C>A]U','G[C>A]A','G[C>A]C','G[C>A]G','G[C>A]U','U[C>A]A','U[C>A]C','U[C>A]G','U[C>A]U','A[G>C]A','A[G>C]C','A[G>C]G','A[G>C]U','C[G>C]A','C[G>C]C','C[G>C]G','C[G>C]U','G[G>C]A','G[G>C]C','G[G>C]G','G[G>C]U','U[G>C]A','U[G>C]C','U[G>C]G','U[G>C]U','A[C>G]A','A[C>G]C','A[C>G]G','A[C>G]U','C[C>G]A','C[C>G]C','C[C>G]G','C[C>G]U','G[C>G]A','G[C>G]C','G[C>G]G','G[C>G]U','U[C>G]A','U[C>G]C','U[C>G]G','U[C>G]U','A[A>C]A','A[A>C]C','A[A>C]G','A[A>C]U','C[A>C]A','C[A>C]C','C[A>C]G','C[A>C]U','G[A>C]A','G[A>C]C','G[A>C]G','G[A>C]U','U[A>C]A','U[A>C]C','U[A>C]G','U[A>C]U','A[U>G]A','A[U>G]C','A[U>G]G','A[U>G]U','C[U>G]A','C[U>G]C','C[U>G]G','C[U>G]U','G[U>G]A','G[U>G]C','G[U>G]G','G[U>G]U','U[U>G]A','U[U>G]C','U[U>G]G','U[U>G]U','A[A>U]A','A[A>U]C','A[A>U]G','A[A>U]U','C[A>U]A','C[A>U]C','C[A>U]G','C[A>U]U','G[A>U]A','G[A>U]C','G[A>U]G','G[A>U]U','U[A>U]A','U[A>U]C','U[A>U]G','U[A>U]U','A[U>A]A','A[U>A]C','A[U>A]G','A[U>A]U','C[U>A]A','C[U>A]C','C[U>A]G','C[U>A]U','G[U>A]A','G[U>A]C','G[U>A]G','G[U>A]U','U[U>A]A','U[U>A]C','U[U>A]G','U[U>A]U')
df_plot_sars$muttype <- apply(df_plot_sars,1,function(x){
return(paste0(substring(x["MutationType"],1,1),
              ">",
              substring(x["MutationType"],2,2)))
})

df_plot_sars$trinuc <- factor(df_plot_sars$trinuc,
                              levels=trinuc_level)
df_plot_sars$muttype <- factor(df_plot_sars$muttype,
                                levels=c("C>U","G>A","A>G","U>C",
                                         "G>U","C>A","G>C","C>G",
                                         "A>C","U>G","A>U","U>A"))
mypal <- pal_npg("nrc")(10)
mypal <- c(mypal, c("#4DAF4A" , "#999999" , "#A65628"))
plot_sars <- ggplot(data =df_plot_sars , aes(x = trinuc, y = log_ratio, fill = muttype)) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  ylab("Log-ratio of mutation frequencies between\nlate and early mutation spectrums") +
  labs(title="SARS-COV-2",
       fill = "Mutation type")+
  scale_fill_manual(values = mypal) +
  my_theme()+
  theme(  
    panel.spacing.x = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),  
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.text.x =element_text(size = 7,angle = 90, hjust = 1),
    axis.text.y =element_text(size = 7),
    axis.title.y = element_text(size = 7),
    plot.margin = margin(0, 0, 0, 0, "null"),
    plot.title = element_text(hjust = 0.5, margin = margin(t = 0.1, unit = "cm"))
  )
plot_sars
ggsave("suppl_figure2A.pdf",width =16, height = 6, units = "cm",dpi = 300)

###supplementary 2B
###bar plot of the log ratio between early and late spectrums for H3N2
df_h3n2 <- df %>% filter(species=="InfluenzaAvirus_h3n2_") %>% arrange(mean_between_years) 
sorted_df_h3n2 <- df_h3n2 %>% arrange(desc(`CUAU`)) %>% select(CUAU, everything())
df_h3n2_early <- df_h3n2 %>% filter(mean_between_years <= 2004 && mean_between_years >= 2000) %>%
                select((ncol(df_sars)-191):ncol(df_sars)) %>% colSums() %>% as.data.frame()  %>% 
                mutate(across(starts_with("."), ~ . / sum(`.`))) %>%
                mutate(MutationType=rownames(df_sars_early)) %>% dplyr::rename("H3N2\n(early)"=".")
df_h3n2_late <- df_h3n2 %>% filter(mean_between_years >= 2019 && mean_between_years <= 2023) %>%
                select((ncol(df_sars)-191):ncol(df_sars)) %>% colSums() %>% as.data.frame()  %>% 
                mutate(across(starts_with("."), ~ . / sum(`.`))) %>%
                mutate(MutationType=rownames(df_sars_late)) %>% dplyr::rename("H3N2\n(late)"=".")
df_h3n2_ratio <- df_h3n2_late %>% right_join(df_h3n2_early)  %>% 
                mutate(log_ratio=log(`H3N2\n(late)`/`H3N2\n(early)`)) %>%
                mutate(log_ratio=ifelse((`H3N2\n(early)`+`H3N2\n(late)`)<= 0.01,0,log_ratio)) %>%
                filter(abs(log_ratio) >0 )
df_plot_h3n2 <- df_h3n2_ratio %>% select("MutationType","log_ratio")
df_plot_h3n2$muttype <- apply(df_plot_h3n2 ,1,function(x){
  return(paste0(substring(x["MutationType"],1,1),
                ">",
                substring(x["MutationType"],2,2)))
})
df_plot_h3n2$trinuc <- apply(df_plot_h3n2,1,function(x){
  return(paste0(substring(x["MutationType"],3,3),
                "[",
                substring(x["MutationType"],1,1),
                ">",
                substring(x["MutationType"],2,2),
                "]",
                substring(x["MutationType"],4,4)))
})
df_plot_h3n2$trinuc <- factor(df_plot_h3n2$trinuc,
                              levels=trinuc_level)
df_plot_h3n2$muttype <- factor(df_plot_h3n2$muttype,
                                levels=c("C>U","G>A","A>G","U>C",
                                          "G>U","C>A","G>C","C>G",
                                          "A>C","U>G","A>U","U>A"))
mypal <- pal_npg("nrc")(10)
mypal <- c(mypal, c("#4DAF4A" , "#999999" , "#A65628"))
plot_h3n2 <- ggplot(data = df_plot_h3n2 , aes(x = trinuc, y = log_ratio, fill = muttype)) +
              geom_bar(stat = "identity") + 
              theme_bw() +
              ylab("Log-ratio of mutation frequencies between\nlate and early mutation spectrums") +
              labs(title="H3N2",
                   fill = "Mutation type")+
              scale_fill_manual(values = mypal) +
              coord_cartesian(ylim = c(-1,1))+
              my_theme()+
              theme(  
                panel.spacing.x = unit(0, "lines"),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_blank(), 
                axis.ticks.x = element_blank(), 
                axis.title.x = element_blank(),
                axis.text.x =element_text(size = 7,angle = 90, hjust = 1),
                axis.text.y =element_text(size = 7),
                axis.title.y = element_text(size = 7),
                plot.margin = margin(0, 0, 0, 0, "null"),
                plot.title = element_text(hjust = 0.5, margin = margin(t = 0.1, unit = "cm")))
ggsave("suppl_figure2B.pdf",width =16, height = 6, units = "cm",dpi = 300)

