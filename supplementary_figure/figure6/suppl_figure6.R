library(tidyverse)
library(this.path)
library(ggplot2)
setwd(this.dir())
source("../../bin/theme_setup.R")
###supplementary figure6
df <- read_tsv("../data/COSMIC_SBS96_Activities.txt")
sbs_df <- df %>% select(-c("Samples"))
col_vector <- colSums(sbs_df)
col_frequency_vector <- col_vector/sum(col_vector)
col_frequency_df <- data.frame(t(col_frequency_vector))
longer_col_frequency_df <- col_frequency_df %>% 
                            pivot_longer(everything(),names_to = "signature", values_to = "proportion") 
longer_col_frequency_df$signature <- factor(longer_col_frequency_df$signature,
                                            levels=c("SBS1","SBS2","SBS3","SBS5",
                                                     "SBS10b","SBS12","SBS16","SBS17a",
                                                     "SBS17b","SBS32","SBS48","SBS85",
                                                     "SBS96D","SBS96H"))
plot <- ggplot(data=longer_col_frequency_df,aes(x=signature,y=proportion)) +
  geom_col(color="black",fill="#F39B7FFF")+
  my_theme()+
  labs(
    x = "Signature",
    y = "Relative contribution of the signature to \n3867 viral mutation spectrums")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  )
  
ggsave("suppl_figure6.pdf",width =8, height = 9, units = "cm",dpi = 300)
