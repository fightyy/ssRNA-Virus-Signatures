library(tidyverse)
library(this.path)
library(ggplot2)
setwd(this.dir())
source("../bin/theme_setup.R")
###Figure5D
#distribution of different mutational signatures across various viral families
meta_df <- read_tsv("../data/context_all_signature_sample_csv.txt") %>% select(c(sample, year, host, species, order,family))
signature_df <- read_tsv("../data/SBS96_De-Novo_Activities_refit.txt")
data_df <- signature_df %>% left_join(meta_df, by=c("Samples"="sample")) 
group_df <- data_df %>% rowwise() %>% 
            mutate(across(matches("SBS"), ~./sum(c_across(matches("SBS")))))
data_df <- group_df %>%  na.omit() %>% 
           mutate(family = gsub("viridae", "", family))
sbs_columns <- grep("SBS", names(group_df), value = TRUE)
sbs_df <- data_df %>% select(sbs_columns) 
sbs_columns <- subset(sbs_columns, apply(sbs_df, 2, function(x) max(x) > 0.2 | mean(x) > 0.05))
all_group_df <- data.frame()
for (familys in unique(data_df$family) ) {
  specie_df <- data_df %>% filter(family == familys)
  sample_count <- nrow(specie_df)
  sbs_vector <- data.frame()
  for (sbs_name in sbs_columns) {
    sbs_df <- specie_df %>% filter(eval(parse(text = sbs_name)) > 0.05 )
    median_proportion <- median(specie_df %>% select(sbs_name) %>% as.vector() %>% unlist())
    proportion <- round(nrow(sbs_df) /nrow(specie_df), digits = 1) 
    row_vector <- c(sbs_name, familys, proportion, sample_count, paste0(familys), median_proportion)
    all_group_df<- rbind(all_group_df,row_vector)
  }
}
colnames(all_group_df) <- c("signature", "family", "proportion", "count", "family_count", "median_proportion")
all_group_df$family_count <- factor(all_group_df$family_count, levels = c("Arena", "Filo", "Hanta", "Orthomyxo","Paramyxo", "Peribunya", 
                                    "Phenui", "Rhabdo", "Calici", "Corona", "Flavi","Picorna", "Toga","Pox"))
plot <- ggplot(all_group_df, aes(x = factor(signature), y = factor(family_count), fill = as.numeric(median_proportion))) +
        geom_tile(fill = NA, colour = 'grey') +
        geom_point( aes( size = proportion), shape = 21) +
        scale_size_discrete(range = c(0, 5),breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
        scale_fill_gradient(low = "#4DBBD5FF", high = "red") +
        labs(x = "Mutational signature", y = "Virus family ") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        my_theme()+
        labs(fill="Median signature exposure \n among the viral families",
             size="Proportion of samples \n with the signature")
ggsave("Figure5D.pdf",width =16, height = 8, units = "cm",dpi = 300)

