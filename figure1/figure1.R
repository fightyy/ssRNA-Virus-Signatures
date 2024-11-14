library(tidyverse)
library(this.path)
library(ggsci)
library(ggplot2)
library(proxy)
library(parallel)
library(factoextra)
library(car)
setwd(this.dir())
source("../bin/theme_setup.R")

###prepare spectrum data
standard_mutation <- c('CUAA', 'CUAC', 'CUAG', 'CUAU', 'CUCA', 'CUCC', 'CUCG', 'CUCU', 'CUGA', 'CUGC', 'CUGG', 'CUGU', 'CUUA', 'CUUC', 'CUUG', 'CUUU', 'GAAA', 'GAAC', 'GAAG', 'GAAU', 'GACA', 'GACC', 'GACG', 'GACU', 'GAGA', 'GAGC', 'GAGG', 'GAGU', 'GAUA', 'GAUC', 'GAUG', 'GAUU', 'AGAA', 'AGAC', 'AGAG', 'AGAU', 'AGCA', 'AGCC', 'AGCG', 'AGCU', 'AGGA', 'AGGC', 'AGGG', 'AGGU', 'AGUA', 'AGUC', 'AGUG', 'AGUU', 'UCAA', 'UCAC', 'UCAG', 'UCAU', 'UCCA', 'UCCC', 'UCCG', 'UCCU', 'UCGA', 'UCGC', 'UCGG', 'UCGU', 'UCUA', 'UCUC', 'UCUG', 'UCUU', 'GUAA', 'GUAC', 'GUAG', 'GUAU', 'GUCA', 'GUCC', 'GUCG', 'GUCU', 'GUGA', 'GUGC', 'GUGG', 'GUGU', 'GUUA', 'GUUC', 'GUUG', 'GUUU', 'CAAA', 'CAAC', 'CAAG', 'CAAU', 'CACA', 'CACC', 'CACG', 'CACU', 'CAGA', 'CAGC', 'CAGG', 'CAGU', 'CAUA', 'CAUC', 'CAUG', 'CAUU', 'GCAA', 'GCAC', 'GCAG', 'GCAU', 'GCCA', 'GCCC', 'GCCG', 'GCCU', 'GCGA', 'GCGC', 'GCGG', 'GCGU', 'GCUA', 'GCUC', 'GCUG', 'GCUU', 'CGAA', 'CGAC', 'CGAG', 'CGAU', 'CGCA', 'CGCC', 'CGCG', 'CGCU', 'CGGA', 'CGGC', 'CGGG', 'CGGU', 'CGUA', 'CGUC', 'CGUG', 'CGUU', 'ACAA', 'ACAC', 'ACAG', 'ACAU', 'ACCA', 'ACCC', 'ACCG', 'ACCU', 'ACGA', 'ACGC', 'ACGG', 'ACGU', 'ACUA', 'ACUC', 'ACUG', 'ACUU', 'UGAA', 'UGAC', 'UGAG', 'UGAU', 'UGCA', 'UGCC', 'UGCG', 'UGCU', 'UGGA', 'UGGC', 'UGGG', 'UGGU', 'UGUA', 'UGUC', 'UGUG', 'UGUU', 'AUAA', 'AUAC', 'AUAG', 'AUAU', 'AUCA', 'AUCC', 'AUCG', 'AUCU', 'AUGA', 'AUGC', 'AUGG', 'AUGU', 'AUUA', 'AUUC', 'AUUG', 'AUUU', 'UAAA', 'UAAC', 'UAAG', 'UAAU', 'UACA', 'UACC', 'UACG', 'UACU', 'UAGA', 'UAGC', 'UAGG', 'UAGU', 'UAUA', 'UAUC', 'UAUG', 'UAUU')
complement_mutation <- c('GAUU', 'GAGU', 'GACU', 'GAAU', 'GAUG', 'GAGG', 'GACG', 'GAAG', 'GAUC', 'GAGC', 'GACC', 'GAAC', 'GAUA', 'GAGA', 'GACA', 'GAAA', 'CUUU', 'CUGU', 'CUCU', 'CUAU', 'CUUG', 'CUGG', 'CUCG', 'CUAG', 'CUUC', 'CUGC', 'CUCC', 'CUAC', 'CUUA', 'CUGA', 'CUCA', 'CUAA', 'UCUU', 'UCGU', 'UCCU', 'UCAU', 'UCUG', 'UCGG', 'UCCG', 'UCAG', 'UCUC', 'UCGC', 'UCCC', 'UCAC', 'UCUA', 'UCGA', 'UCCA', 'UCAA', 'AGUU', 'AGGU', 'AGCU', 'AGAU', 'AGUG', 'AGGG', 'AGCG', 'AGAG', 'AGUC', 'AGGC', 'AGCC', 'AGAC', 'AGUA', 'AGGA', 'AGCA', 'AGAA', 'CAUU', 'CAGU', 'CACU', 'CAAU', 'CAUG', 'CAGG', 'CACG', 'CAAG', 'CAUC', 'CAGC', 'CACC', 'CAAC', 'CAUA', 'CAGA', 'CACA', 'CAAA', 'GUUU', 'GUGU', 'GUCU', 'GUAU', 'GUUG', 'GUGG', 'GUCG', 'GUAG', 'GUUC', 'GUGC', 'GUCC', 'GUAC', 'GUUA', 'GUGA', 'GUCA', 'GUAA', 'CGUU', 'CGGU', 'CGCU', 'CGAU', 'CGUG', 'CGGG', 'CGCG', 'CGAG', 'CGUC', 'CGGC', 'CGCC', 'CGAC', 'CGUA', 'CGGA', 'CGCA', 'CGAA', 'GCUU', 'GCGU', 'GCCU', 'GCAU', 'GCUG', 'GCGG', 'GCCG', 'GCAG', 'GCUC', 'GCGC', 'GCCC', 'GCAC', 'GCUA', 'GCGA', 'GCCA', 'GCAA', 'UGUU', 'UGGU', 'UGCU', 'UGAU', 'UGUG', 'UGGG', 'UGCG', 'UGAG', 'UGUC', 'UGGC', 'UGCC', 'UGAC', 'UGUA', 'UGGA', 'UGCA', 'UGAA', 'ACUU', 'ACGU', 'ACCU', 'ACAU', 'ACUG', 'ACGG', 'ACCG', 'ACAG', 'ACUC', 'ACGC', 'ACCC', 'ACAC', 'ACUA', 'ACGA', 'ACCA', 'ACAA', 'UAUU', 'UAGU', 'UACU', 'UAAU', 'UAUG', 'UAGG', 'UACG', 'UAAG', 'UAUC', 'UAGC', 'UACC', 'UAAC', 'UAUA', 'UAGA', 'UACA', 'UAAA', 'AUUU', 'AUGU', 'AUCU', 'AUAU', 'AUUG', 'AUGG', 'AUCG', 'AUAG', 'AUUC', 'AUGC', 'AUCC', 'AUAC', 'AUUA', 'AUGA', 'AUCA', 'AUAA')
spectrum_df <- read_tsv("../data/all_mpxv_sars_phyvirus_sample.txt") %>% select(c(name, brLen,standard_mutation))
history_df <- read_csv("../data/all_meta_corrected_host.csv")
data_df <- history_df %>% select(accession, year, corrected_host, species, order) %>% 
           right_join(spectrum_df, c("accession" = "name"), relationship = "many-to-many")  %>% 
           na.omit() %>% mutate_at( vars(-c("accession", "corrected_host", "species", "order")),as.numeric)
#convert raw spectrum of -ssRNA virus to the genome strand
corrected_data_df <- data.frame()
for (order_name in (data_df %>% .$order %>% unique()))  {
  small_df1 <- data_df %>% filter(order==order_name) 
  if (order_name %in% c("Mononegavirales", "Articulavirales", "ssRNA_-_") ){
    small_df2 <- small_df1
    colnames(small_df2) <- c("accession", "year","corrected_host", "species" ,"order" ,"brLen", complement_mutation)
    small_df3 <- small_df2[c("accession", "year","corrected_host", "species" ,"order" ,"brLen", standard_mutation)]
  }
  else{
    small_df3 <- small_df1 
  }
  corrected_data_df <- rbind(small_df3, corrected_data_df)
}
species_vector <-  c(  "InfluenzaAvirus_h5n1_", 
                       "InfluenzaAvirus_h3n2_",  
                       "Lyssavirus_rabies", 
                       "Mumps_orthorubulavirus",
                       "Zaire_ebolavirus",
                       "West_Nile_virus", 
                       "Zika_virus",
                       "sars-cov-2" )
corrected_data_df <- corrected_data_df %>% filter(species %in% species_vector)
grouped_df <- corrected_data_df %>%
              filter(corrected_host != "Others") %>% 
              group_by(corrected_host, species, order) %>% 
              summarise_at(colnames(data_df)[7:198],sum) 
enough_grouped_df <- grouped_df %>% ungroup() %>% filter(rowSums(select(., 4:195))>100)

#adjust viral spectra to human composition
df <- read_tsv("../data/all_phyvirus_human_compos.txt") %>% filter(name != "name") %>% mutate(across(-1, ~as.numeric(.)))
normalise_df  <- df %>% rowwise() %>% mutate(across(-1, ~ ./ sum(c_across(-1))))
Human_row <- normalise_df[normalise_df$name == "Human", -1] %>% as.vector()
df_result <- normalise_df %>% select(-name)
df_result2 <- Human_row / df_result  
df_result_normalized <- cbind(name = normalise_df$name, df_result2)

correct_spectrum <- function(factor_df, spectrum_df, species_name) {
  factor_row <- factor_df[factor_df$name == species_name, -1] %>% mutate_all(as.numeric) %>% as.vector() %>% unlist()
  species_spectrum_df <- spectrum_df %>% filter(species == species_name) 
  spectrum_df <- species_spectrum_df %>% select(-c(corrected_host,species, order)) %>% mutate_all(as.numeric)
  corrected_spectrum_df  <- apply(spectrum_df, 1, function(row) row * factor_row)
  corrected_spectrum_df  <- as.data.frame(t(corrected_spectrum_df)) %>% mutate_all(round)
  name_corrected_spectrum_df <- species_spectrum_df %>% select(corrected_host,species, order)
  all_corrected_spectrum_df <- cbind(name_corrected_spectrum_df, corrected_spectrum_df)
  return(all_corrected_spectrum_df)
} 
all_corrected_df <- data.frame()
for (specie in species_vector){
  corrected_df <- correct_spectrum(df_result_normalized, enough_grouped_df, specie)
  all_corrected_df <- rbind(all_corrected_df , corrected_df)
}


###Figure1A
#plot mutational spectra for H3N2 in humans, H3N2 in mammalian and SARS-CoV-2 in humans
small_spectrum_df <- grouped_df %>% ungroup() %>% filter(species %in% c("InfluenzaAvirus_h3n2_", "sars-cov-2")) %>%
                      filter(corrected_host %in% c("Mammalia", "Homo_sapiens")) 
normalised_df <- small_spectrum_df [-4,] %>% select(-c("species","corrected_host", "order")) %>% 
                  rowwise() %>% mutate(across(everything(), ~./ sum(c_across(everything())))) %>% t() %>% as.data.frame()
colnames(normalised_df) <- c("H3N2\n(Homo sapiens)", "SARS-CoV-2\n(Homo sapiens)", "H3N2\n(Mammalia)")
normalised_df$MutationType <- rownames(normalised_df)
SBS192.input.catalog <- normalised_df
SBS192.input.catalog$muttype <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["MutationType"],1,1),
                ">",
                substring(x["MutationType"],2,2)))
})
SBS192.input.catalog$trinuc <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["MutationType"],3,3),
                substring(x["MutationType"],1,1),
                substring(x["MutationType"],4,4)))
})
long_rename_SBS192.input.catalog <- SBS192.input.catalog  %>% 
                                    pivot_longer( cols = c("H3N2\n(Homo sapiens)", "H3N2\n(Mammalia)", "SARS-CoV-2\n(Homo sapiens)"),
                                                  names_to = "signature",
                                                  values_to = "exposure")
long_rename_SBS192.input.catalog$muttype <- factor(long_rename_SBS192.input.catalog$muttype,
                                                   levels=c("C>U","G>A","A>G","U>C",
                                                            "G>U","C>A","G>C","C>G",
                                                            "A>C","U>G","A>U","U>A"))
mypal <- pal_npg("nrc")(10)
mypal <- c(mypal, c("#4DAF4A" , "#999999" , "#A65628"))
plot <- ggplot(data =long_rename_SBS192.input.catalog , aes(x = trinuc, y = exposure, fill = muttype)) +
        geom_bar(stat = "identity") +
        facet_grid(signature ~ muttype, scales = "free_x", space = "free_x") +
        theme_bw() +
        ylab("Proportion of mutation") +
        scale_fill_manual(values = mypal) +
        scale_y_continuous(expand = c(0, 0, 0.05, 0))+
        my_theme()+
        theme(
          panel.spacing.x = unit(0, "lines"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 0),
          axis.text.y =element_text(size = 7),
          axis.title.y = element_text(size = 9,margin = margin(r=3)),
          strip.text = element_text(size = 9),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0, "null"))+
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("Figure1A.pdf",width =18, height = 9, units = "cm",dpi = 300)


###Figure1B
#calculate pairwise cosine similarity of viral mutational spectra between different host form the same viral species
calculate_similarity_host <- function (new_grouped_df) {
order_vector <- new_grouped_df$species %>% unique()
order_test_df <- c()
for (order_name in order_vector ) {
  test_df <- new_grouped_df  %>% filter(species == order_name) 
  host_vector <- unique(test_df$corrected_host)
  if (length(host_vector) <=1) {
    next 
  }
  host_combine <- combn(host_vector, 2)
  for (index in rep(1:ncol(host_combine),1)) {
    host_name <-  host_combine[, index]
    matrix <- test_df %>% filter(corrected_host %in% host_name) %>% ungroup() %>% select(4:195)
    cosine_similarity <- proxy::simil(x = matrix[1,]/sum(matrix[1,]), y = matrix[2,]/sum(matrix[2,]), method = "cosine")
    result_df <- c(order_name, host_name[1], host_name[2], cosine_similarity)
    order_test_df <- rbind(order_test_df, result_df)
  }
}
colnames(order_test_df) <- c("order", "hostA", "hostB", "cosine_similarity") 
order_test_df <- as.data.frame(order_test_df) 
order_test_df <- rownames_to_column(order_test_df, var = "rowname") %>% select(-rowname) %>% mutate_at("cosine_similarity", as.numeric)
return(order_test_df )
}
df_100 <- calculate_similarity_host(all_corrected_df) 
p1 <- ggplot(df_100 , aes(x = factor(order), y = cosine_similarity)) +
      geom_boxplot(aes(group = factor(order))) +
      geom_point(
        position=position_jitter(width = 0.2), 
        pch=21, 
        size=2,
        fill = "#F39B7FFF")+
      scale_x_discrete(labels=c( "H3N2", "H5N1 ",
                                 "Rabies lyssavirus", "SARS-CoV-2", "West Nile virus","Zika virus"))+
      labs(
        x = "Viral species",
        y = "Pairwise cosine similarity between \n viral hosts' mutation spectrums")+
      my_theme()+
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("Figure1B.pdf",width = 5.5, height = 6, units = "cm",dpi = 300)




###Figure1C
#calculate pairwise cosine similarity of viral mutational specta between different viral species from the same host
calculate_similarity_species <- function (new_grouped_df) {
  order_vector <- new_grouped_df$corrected_host %>% unique()
  order_test_df <- c()
  for (host_name in order_vector ) {
    test_df <- new_grouped_df  %>% filter(corrected_host == host_name) 
    species_vector <- unique(test_df$species)
    if (length(species_vector) <=1) {
      next 
    }
    species_combine <- combn(species_vector, 2)
    for (index in rep(1:ncol(species_combine),1)) {
      species_name <-  species_combine[, index]
      matrix <- test_df %>% filter(species %in% species_name ) %>% ungroup() %>% select(4:195)
      result <- fisher.test(matrix,simulate.p.value=TRUE) 
      cosine_similarity <- proxy::simil(x = matrix[1,]/sum(matrix[1,]), y = matrix[2,]/sum(matrix[2,]), method = "cosine")
      result_df <- c(host_name, species_name[1], species_name[2], result$p.value, cosine_similarity)
      order_test_df <- rbind(order_test_df, result_df)
    }
  }
  colnames(order_test_df) <- c("host", "speciesA", "speciesB", "p_value", "cosine_similarity") 
  order_test_df <- as.data.frame(order_test_df)
  order_test_df <- rownames_to_column(order_test_df, var = "rowname") %>% select(-rowname) %>% 
    mutate_at("cosine_similarity", as.numeric)
  return(order_test_df )
}
host_df_100 <- calculate_similarity_species(all_corrected_df)
p2 <- ggplot(host_df_100, aes(x = factor(host), y = cosine_similarity)) +
      geom_boxplot(aes(group = factor(host))) +
      geom_point(fill = "#F39B7FFF",
                 position=position_jitter(width = 0.2), 
                 pch=21, 
                 size=2)+
      scale_x_discrete(labels=c("Aves", "Homo sapiens",
                                "Insecta", "Mammalia"))+
      labs(
        x = "Host",
        y = "Pairwise cosine similarity between \n species' mutation spectrums")+
      my_theme()+
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("Figure1C.pdf",width = 5.5, height = 6, units = "cm",dpi = 300)


###prepare pca data
df <- read_tsv("../data/context_all_signature_sample_csv.txt") 
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
                     "Zaire_ebolavirus", "Mumps_orthorubulavirus", "West_Nile_virus", "Zika_virus","sars-cov-2" )
balti_dict <- setNames(c(rep("-ssRNA", times = 5), rep("+ssRNA", times = 3)), species_vector)
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
species_colors <- c("#E41A1C","#FF4500", "#F781BF", "#D2691E", "#FF7F50",  "#6BAED6", "#9ECAE1", "#8DA0CB")
species_dict <- setNames(species_colors, species_vector)
balti_dict <- setNames(c(rep("-ssRNA", times = 5), rep("+ssRNA", times = 3)),species_vector)

p3 <- ggplot(pca_plot_data, aes(x = Dim.1, y = Dim.2,shape=species)) +
  geom_point(aes(color = species,shape=baltimore_class),size=1) +
  scale_color_manual(values = species_dict,
                     limits= c("InfluenzaAvirus_h5n1_", "InfluenzaAvirus_h3n2_",  
                               "Lyssavirus_rabies",
                               "Zaire_ebolavirus",
                               "West_Nile_virus", "Zika_virus",
                               "sars-cov-2",
                               "Mumps_orthorubulavirus"),
                     labels=c("H5N1", "H3N2",
                              "Rabies lyssavirus", "Zaire ebolavirus",
                              "West Nile virus","Zika virus","SARS-CoV-2",
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
ggsave("Figure1D.pdf",width = 8, height = 6, units = "cm",dpi = 300)


###Figure1E
#PCA plot of mutational spectra for eight viral species(colored by host)
host_colors <- c( "#E7B800", "#CF4E9CFF", "#97A1A7FF", "#2F509EFF", "#8C57A2FF")
p2 <- ggplot(pca_plot_data, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = host,shape=baltimore_class),size=1) +
  scale_color_manual(values=host_colors,
                     labels=c("Aves", "Homo sapiens",
                              "Insecta", "Mammalia"))+
  scale_shape_manual(values = c("-ssRNA" = 16, 
                                "+ssRNA" = 17))+
  labs(x = paste("PC1","(",round(pc1*100,1),"%)",sep="") , y = paste("PC2","(",round(pc2*100,1),"%)",sep="")) +
  my_theme()+
  theme(legend.key.height = unit(0.3, "cm"))+
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))+
  labs(color="Host",
       shape="Baltimore class")
ggsave("Figure1E.pdf",width = 8, height = 6, units = "cm",dpi = 300)


###Figure 1F
#percentage of variance explained by vial species and hosts for different PC
data_df <- df_normalized 
short_good_df <- data_df %>% filter(! order %in% c("ssRNA_+_", "ssRNA_-_")) %>% filter(host != "internal") %>% 
                filter(! species %in% c("Monkeypox_virus","Human_orthorubulavirus_2","Human_orthorubulavirus_4"))
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
ggsave("Figure1F.pdf", plot2, width = 8, height = 6, units = "cm",dpi = 300)

