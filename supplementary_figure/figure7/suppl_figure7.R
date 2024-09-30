library(tidyverse)
library(this.path)
library(rstatix)
library(ggpubr)
library(cowplot)
library(ggpmisc)
setwd(this.dir())
source("../../bin/theme_setup.R")
###supplementary figure7A
###The probability of the mutations in the spike protein predicted to be bounded by MHC molecules for the SARS-COV2
mhc_count <- read_table("../data/mhc_table_sars.txt",col_names=c("strong","weak","all","accession_pos","protein")) %>% 
             mutate("weak_strong"=weak+strong) %>% rowwise() %>% 
             mutate(accession = str_split(accession_pos, "-")[[1]][1],
                    pos = str_split(accession_pos, "-")[[1]][2])

mhc_proportion <- mhc_count %>% group_by(accession,pos,protein) %>% 
                  summarise_at(vars(weak, strong, weak_strong, all), funs(sum)) %>% 
                  mutate_at(vars(weak, strong, weak_strong), funs(./all))
meta_df <- read_csv("../../data/all_meta_corrected_host.csv") %>% select("accession","year","date")
data_df <- mhc_proportion %>% left_join(meta_df) %>% na.omit
mhc_data_df_date <- data_df %>%
                    mutate_at(vars("date"),list(date_obj = ~as.Date(., format = "%Y/%m/%d"))) %>%
                    mutate(days = as.numeric(date_obj - as.Date("2019-12-01"))) %>%  mutate(months=4*floor(days/120))
result1 <- lm(data=mhc_data_df_date, strong ~ months)
summary(result1)
result2 <- lm(data=mhc_data_df_date, weak_strong ~ months)
summary(result2)
p <- ggplot(mhc_data_df_date, aes(x=months, y=weak_strong)) +
      geom_violin(aes(fill= months,group=factor(months)),color=NA) +
      scale_fill_gradient(low = "#4DBBD5FF", high = "#E64B35FF")+
      stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1),
                   geom="pointrange", color="black",
                   shape = 18, size = 0.75)+
      geom_smooth(method="lm", se=FALSE,color="red",linewidth=0.6)+
      stat_poly_eq(aes(label = paste(..eq.label..)),
                   label.x = 0.01 , label.y = 0.9, size = 7 / .pt,
                   formula = y ~ x, parse = TRUE) +
      stat_cor(aes(label = paste(..p.label..)), 
               label.x = -2, label.y = 0.07, size = 7 /.pt) +
      coord_cartesian(ylim = c(0, 0.07))+
      labs(title="SARS-CoV-2",
           fill="Month since\n2019/12",
           x = "Month (since 2019)",
           y = "Proportion of peptides\nbound by MHC class I")+
      my_theme()+
      theme(plot.title =element_text(margin=margin(t= 0.1 ,unit="cm")))
ggsave("suppl_figure7A.pdf",width = 18, height = 6, units = "cm",dpi = 300)

###supplementary figure7B
###The probability of the mutations in the HA protein predicted to be bounded by MHC molecules for the H3N2 influenza virus
mhc_count <- read_table("../data/mhc_table_h3n2.txt",col_names=c("strong","weak","all","accession_pos"))%>% 
             mutate("weak_strong"=weak+strong) %>% rowwise() %>% 
             mutate(accession = str_split(accession_pos, "-")[[1]][1],
                       pos = str_split(accession_pos, "-")[[1]][2])
mhc_proportion <- mhc_count %>% group_by(accession,pos) %>% 
                  summarise_at(vars(weak, strong, weak_strong, all), funs(sum)) %>% 
                  mutate_at(vars(weak, strong, weak_strong), funs(./all))
meta_df <- read_csv("../../data/all_meta_corrected_host.csv") %>% select("accession","year")
data_df <- mhc_proportion %>% left_join(meta_df) %>% na.omit %>% filter(year >=2000) 
data_df$year <- as.numeric(data_df$year)
result1 <- lm(data=data_df, strong ~ year)
summary(result1)
result2 <- lm(data=data_df, weak_strong ~ year)
summary(result2)
p <- ggplot(data_df, aes(x=year, y=weak_strong)) +
    geom_violin(aes(fill= year,group=factor(year)), color=NA) +
    scale_fill_gradient(low = "#4DBBD5FF", high = "#E64B35FF") +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "pointrange", color = "black",
                 shape = 18, size = 0.75) +
    geom_smooth(method="lm", se=FALSE,color="red",linewidth=0.6)+
    stat_poly_eq(aes(label = paste(..eq.label..)),
                 label.x = 0.05 , label.y = 0.9, size = 7 / .pt,
                 formula = y ~ x, parse = TRUE) +
    stat_cor(aes(label = paste(..p.label..)),
             label.x = 2000, label.y =0.08, size = 7 /.pt) +
    labs(title = "H3N2",
         fill = "Year",
         x = "Year",
         y = "Proportion of peptides\nbound by MHC class I") +
    my_theme()+
    theme(plot.title =element_text(margin=margin(t= 0.5 ,unit="cm")))
ggsave("suppl_figure7B.pdf",width = 18, height = 6, units = "cm",dpi = 300)
