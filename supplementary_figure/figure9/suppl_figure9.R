library(tidyverse)
library(this.path)
library(ggsci)
setwd(this.dir())
source("../../bin/theme_setup.R")

###supplementary figure9
###Frequency changes of high-frequency (>=0.08 in donors) non-synonymous iSNVs for SARS-CoV-2 before and after transmission
isnv_pair <- read_csv("../../data/isnv_pair_ruan_andreas.csv")
cutoff = 0.08
high_nonsyn <- isnv_pair %>% filter(donor_alt_frequency >= cutoff, recipient_alt_frequency >=0.03) %>%
              filter(func == "A") %>% 
              mutate(
                all_mean_frequency = mean(donor_alt_frequency),
                all_mean_recepient_frequency = mean(recipient_alt_frequency),
                mutation_type = "High-frequency nonsynonymous")
df <- high_nonsyn
df_long <- data.frame(category = rep(c("donor", "recipient"), each = nrow(df)),
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
        legend.position = "none")
ggsave("suppl_figure9.pdf",width = 4.5, height = 5, units = "cm",dpi = 300)


