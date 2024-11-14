library(tidyverse)
library(this.path)
library(ggsci)
library(ggplot2)
library(proxy)
library(parallel)
library(factoextra)
library(car)
setwd(this.dir())
source("../../bin/theme_setup.R")

###prepare pca data
df <- read_tsv("../../data/context_all_signature_sample_csv.txt") 
data_df <- df  %>% mutate(family = gsub("viridae", "", family)) %>% filter(host != "Others")
num_cores=5
df_normalized <- mclapply(seq_len(nrow(data_df)), function(i) {
  row <- data_df[i, ]
  row[-(1:6)] <- row[-(1:6)] / sum(row[-(1:6)])
  return(row)
}, mc.cores = num_cores)

df_normalized <- do.call(rbind, df_normalized) 
df_normalized <- df_normalized %>% filter(host != "internal")
df_normalized <- df_normalized %>% filter(! species %in% c("Monkeypox_virus","Human_orthorubulavirus_2"))

pca_df <- df_normalized[,(ncol(df_normalized)-191):ncol(df_normalized)] 
species_vector <-  c("InfluenzaAvirus_h5n1_", "InfluenzaAvirus_h3n2_", "Lyssavirus_rabies", 
                     "Mumps_orthorubulavirus", "sars-cov-2" )
balti_dict <- setNames(c(rep("-ssRNA", times = 4), rep("+ssRNA", times = 1)), species_vector)
df_normalized_good <- df_normalized %>% filter(! order %in% c("ssRNA_+_", "ssRNA_-_"))  %>% 
                      filter(species %in% species_vector) %>% rowwise() %>% 
                      mutate(baltimore_class= balti_dict[[species]]) %>% select(baltimore_class,everything())
pca_df <- df_normalized_good %>% select((ncol(df_normalized_good)-191):ncol(df_normalized_good))
pca_result <- prcomp(pca_df)
pc1 <- pca_result$sdev[1]^2 / sum(pca_result$sdev^2)
pc2 <- pca_result$sdev[2]^2 / sum(pca_result$sdev^2)
pca_plot_data <- as_tibble(get_pca_ind(pca_result)$coord) %>% 
  select("Dim.1", "Dim.2") %>% 
  cbind(df_normalized_good[,1:(ncol(df_normalized_good)-192)]) 


###Figure1D
#PCA plot of mutational spectra for eight viral species(colored by species)
species_colors <- c("#E41A1C","#FF4500", "#F781BF", "#D2691E", "#8DA0CB")
species_dict <- setNames(species_colors, species_vector)


p3 <- ggplot(pca_plot_data, aes(x = Dim.1, y = Dim.2,shape=species)) +
  geom_point(aes(color = species,shape=baltimore_class),size=1) +
  scale_color_manual(values = species_dict,
                     limits= c("InfluenzaAvirus_h5n1_", "InfluenzaAvirus_h3n2_",  
                               "Lyssavirus_rabies",
                               "sars-cov-2",
                               "Mumps_orthorubulavirus"),
                     labels=c("H5N1", "H3N2",
                              "Rabies lyssavirus", "SARS-CoV-2",
                              "Mumps\northorubulavirus")) +
  scale_shape_manual(values = c("-ssRNA" = 16, 
                                "+ssRNA" = 17))+
  labs(x = paste("PC1","(",round(pc1*100,1),"%)",sep="") , 
       y = paste("PC2","(",round(pc2*100,1),"%)",sep="")) +
  my_theme()+
  theme(legend.key.height = unit(0.3, "cm"))+
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))+
  labs(shape="Baltimore class",
       color="Viral species")
ggsave("suppl_figure2A.pdf",width = 8, height = 6, units = "cm",dpi = 300)


###Figure1E
#PCA plot of mutational spectra for eight viral species(colored by host)
host_colors <- c( "#E7B800", "#CF4E9CFF", "#2F509EFF", "#8C57A2FF")
p2 <- ggplot(pca_plot_data, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = host,shape=baltimore_class),size=1) +
  scale_color_manual(values=host_colors,
                     labels=c("Aves", "Homo sapiens", "Mammalia"))+
  scale_shape_manual(values = c("-ssRNA" = 16, 
                                "+ssRNA" = 17))+
  labs(x = paste("PC1","(",round(pc1*100,1),"%)",sep="") , y = paste("PC2","(",round(pc2*100,1),"%)",sep="")) +
  my_theme()+
  theme(legend.key.height = unit(0.3, "cm"))+
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))+
  labs(color="Host",
       shape="Baltimore class")
ggsave("suppl_figure2B.pdf",width = 8, height = 6, units = "cm",dpi = 300)


###Figure 1F
#percentage of variance explained by vial species and hosts for different PC
data_df <- df_normalized 
short_good_df <- data_df %>% filter(! order %in% c("ssRNA_+_", "ssRNA_-_")) %>% filter(host != "internal") %>% 
                 filter(! species %in% c("Monkeypox_virus","Human_orthorubulavirus_2","Human_orthorubulavirus_4",
                        "Zaire_ebolavirus","West_Nile_virus", "Zika_virus"))
pca_df <- short_good_df %>% select((ncol(short_good_df)-191):ncol(short_good_df))
pca_result <- prcomp(pca_df)
pc1 <- pca_result$sdev[1]^2 / sum(pca_result$sdev^2)
pc2 <- pca_result$sdev[2]^2 / sum(pca_result$sdev^2)
pca_plot_data <- as_tibble(get_pca_ind(pca_result)$coord) %>% 
                 select("Dim.1", "Dim.2","Dim.3") %>% 
                 cbind(short_good_df[,1:(ncol(short_good_df)-192)])
#calculate variance for unbalanced data
model3 <- lm(Dim.1 ~ species + host, data = pca_plot_data)
pc1_variance <- suppressWarnings(Anova(model3, type = "III")[[1]]/sum(Anova(model3, type = "III")[[1]])[1])
pc1_variance <- pc1_variance[2:3]

model4 <- lm(Dim.2 ~ species + host, data = pca_plot_data)
pc2_variance <- suppressWarnings(Anova(model4, type = "III")[[1]]/sum(Anova(model4, type = "III")[[1]])[1])
pc2_variance <- pc2_variance[2:3]

model5 <- lm(Dim.3 ~ species + host, data = pca_plot_data)
pc3_variance <- suppressWarnings(Anova(model5, type = "III")[[1]]/sum(Anova(model5, type = "III")[[1]])[1])
pc3_variance <- pc3_variance[2:3]

variance_df <- data.frame(variable=c("viral species","host"),
                          PC1=pc1_variance,
                          PC2=pc2_variance,
                          PC3=pc3_variance)
longer_variance_df <- variance_df %>% pivot_longer(!variable, names_to="pc_axis",values_to="proportion")

plot2 <- ggplot(longer_variance_df, aes(x = pc_axis , y = proportion, fill=variable)) + 
  geom_col(position = "dodge") +
  scale_fill_npg(labels=c("Host", "Viral species"))+
  labs(x="PC dimensions",
       y="Proportion of variance\n explained by a variable",
       fill="Variable")+
  my_theme()
ggsave("suppl_figure2C.pdf", plot2, width = 8, height = 6, units = "cm",dpi = 300)

