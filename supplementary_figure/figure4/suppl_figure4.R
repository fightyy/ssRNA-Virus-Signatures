library(tidyverse)
library(this.path)
library(MutationalPatterns)
library(mSigTools)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(factoextra)
library(ggtree)
library(ggsci)
library(ggtreeExtra)
library(ggnewscale)
library(cols4all)
setwd(this.dir())
source("../../bin/theme_setup.R")

###supplementary 4A
###Clustering results of mutational signatures across different viral species
data_df <- read_table("../data/all_signature_species.txt")
data_df <- data_df %>% select(-contains("Human"))
data <- data_df[, -1]
hclust_signatures <- cluster_signatures(data) 
species_color <- pal_npg("nrc")(10)
species_color <- c(species_color, c("#4DAF4A" , "#999999" , "#A65628"))
group_color <- c4a("palette36", 36)
group_color <- c(group_color, c4a("brewer.paired", 12))
host_color <- c( "#E7B800", "#CF4E9CFF", "#97A1A7FF", "#2F509EFF", "#8C57A2FF")
labels <- c("H3N2-Homo sapiens-A",
            "H3N2-Homo sapiens-B",
            "H3N2-Mammalia-A",
            "H5N1-Aves-A",
            "Rabies lyssavirus-Homo sapiens-A",
            "Rabies lyssavirus-Mammalia-A",
            "Mumps orthorubulavirus-Homo sapiens-A",
            "West Nile virus-Aves-A",
            "West Nile virus-Homo sapiens-A",
            "West Nile virus-Insecta-A",
            "West Nile virus-Insecta-B",
            "Zika virus-Homo sapiens-A",
            "Sars cov2-Homo sapiens-A",
            "Sars cov2-Homo sapiens-B")
species_vector <- c("H3N2","H5N1","Rabies lyssavirus","Mumps orthorubulavirus",
                    "West Nile virus","Zika virus","Sars cov2")
balti_dict <- setNames(c(rep("-ssRNA", times = 4), rep("+ssRNA", times = 3)),species_vector)
hclust_signatures$labels <- labels 
labels(hclust_signatures)[hclust_signatures$order]
labels(hclust_signatures)
infor_df <- data.frame(labels=labels)
infor_df <- infor_df %>% rowwise() %>%
            mutate(species= strsplit(labels, "-")[[1]][1]) %>% 
            mutate(host= strsplit(labels, "-")[[1]][2]) %>% 
            mutate(baltimore_class= balti_dict[[species]])
clus <- cutree(hclust_signatures, h=0.15)
g <- split(names(clus), clus)
p <- ggtree(hclust_signatures,branch.length="none")
clades <- sapply(g, function(n) MRCA(p, n))
tree <- groupClade(p, clades, group_name='subtree') +aes(color=subtree)+
        scale_color_manual(values=group_color)+
        guides(color = "none")
species_colors <- c("#E41A1C","#FF4500", "#F781BF", "#FF7F50", "#D2691E", "#6BAED6", "#9ECAE1", "#8DA0CB")
species <-  c("InfluenzaAvirus_h5n1_", "InfluenzaAvirus_h3n2_", "Lyssavirus_rabies", 
              "Mumps_orthorubulavirus", "Zaire_ebolavirus", "West_Nile_virus", "Zika_virus","sars-cov-2" )
species <-  c("H3N2", "H5N1", "Rabies lyssavirus", 
              "Mumps orthorubulavirus", "Zaire_ebolavirus", "West Nile virus", "Zika virus","Sars cov2" )
species_dict <- setNames(species_colors, species)
plot <- tree + geom_fruit(
        data=infor_df,
        geom=geom_tile,
        width=0.5,
        offset = 0.1,
        mapping=aes(y=labels,fill=baltimore_class))+
        scale_fill_manual(values=species_color,
                          guide=guide_legend(title="Baltimore class",keywidth=0.5, keyheight=0.5, order=3),
                          na.translate=FALSE)+
        new_scale_fill() 
plot1 <- plot + geom_fruit(
        data=infor_df,
        geom=geom_tile,
        width=0.5,
        offset = 0.1,
        mapping=aes(y=labels,fill=species))+
        scale_fill_manual(values=species_dict,
                          guide=guide_legend(title="Viral specie",keywidth=0.5, keyheight=0.5, order=3),
                          na.translate=FALSE)+
        new_scale_fill() 

plot2 <- plot1 + geom_fruit(
        data=infor_df,
        geom=geom_tile,
        offset = 0.1,
        width=0.5,
        mapping=aes(y=labels,fill=host))+
        scale_fill_manual(values=host_color,
                          guide=guide_legend(title="Host",keywidth=0.5, keyheight=0.5, order=3),
                          na.translate=FALSE)
ggsave("suppl_figure4A.pdf",unit="cm", width=8, height=8,dpi = 300)

###supplementary 4B
###Clustering results of mutational signatures across different viral families
data_df <- read_table("../data/all_signature_family.txt")
data_df <- data_df %>% select(-contains("Human"))
data <- data_df[, -1]
hclust_signatures <- cluster_signatures(data) 
species_color <- pal_npg("nrc")(10)
species_color <- c(species_color, c("#4DAF4A" , "#999999" , "#A65628"))
group_color <- c4a("palette36", 36)
group_color <- c(group_color, c4a("brewer.paired", 12))
host_color <- c( "#E7B800", "#CF4E9CFF", "#97A1A7FF", "#2F509EFF", "#8C57A2FF")
host_dict <- setNames(c("Homo sapiens", "Reptilia", "Mammalia", "Insecta","Aves"),
                      c("Homo", "Reptilia", "Mammalia", "Insecta","Aves"))
species_vector <- c("Arena", "Filo", "Hanta", "Orthomyxo","Paramyxo", "Peribunya", 
                    "Phenui", "Rhabdo", "Calici", "Corona", "Flavi","Picorna", "Toga")
balti_dict <- setNames(c(rep("-ssRNA", times = 8), rep("+ssRNA", times = 5)),species_vector)
infor_df <- data.frame(labels=hclust_signatures$labels)
infor_df <- infor_df %>% rowwise() %>%
            mutate(species= strsplit(labels, "_")[[1]][1]) %>% 
            mutate(host= strsplit(labels, "_")[[1]][2]) %>%
            mutate(host= host_dict[[host]]) %>%
            mutate(labels = factor(labels, levels =labels(hclust_signatures))) %>% 
            mutate(baltimore_class= balti_dict[[species]])
clus <- cutree(hclust_signatures, h=0.15)
g <- split(names(clus), clus)
p <- ggtree(hclust_signatures,branch.length="none")
clades <- sapply(g, function(n) MRCA(p, n))
tree <- groupClade(p, clades, group_name='subtree') +aes(color=subtree)+
        scale_color_manual(values=group_color)+
        guides(color = "none")
plot <- tree + geom_fruit(
        data=infor_df,
        geom=geom_tile,
        width=1,
        offset = 0.01,
        mapping=aes(y=labels,fill=baltimore_class))+
        scale_fill_manual(values=species_color,
                          guide=guide_legend(title="Baltimore class",keywidth=0.5, keyheight=0.5, order=3),
                          na.translate=FALSE)+
        new_scale_fill() 
plot1 <- plot + geom_fruit(
        data=infor_df,
        geom=geom_tile,
        width=1,
        offset = 0.08,
        mapping=aes(y=labels,fill=species))+
        scale_fill_manual(values=species_color,
                          guide=guide_legend(title="Viral family",keywidth=0.5, keyheight=0.5, order=3),
                          na.translate=FALSE)+
        new_scale_fill() 
plot2 <- plot1 + geom_fruit(
          data=infor_df,
          geom=geom_tile,
          offset = 0.08,
          width=1,
          mapping=aes(y=labels,fill=host))+
          scale_fill_manual(values=host_color,
                            guide=guide_legend(title="Host",keywidth=0.5, keyheight=0.5, order=3),
                            na.translate=FALSE)
ggsave("suppl_figure4B.pdf",unit="cm", width=8, height=16,dpi = 300)

