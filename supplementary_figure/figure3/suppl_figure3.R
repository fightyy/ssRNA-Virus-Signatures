library(tidyverse)
library(this.path)
library(proxy)
library(parallel)
library(factoextra)
library(car)
library(ggsci)
setwd(this.dir())
source("../../bin/theme_setup.R")

#prepare data
standard_mutation <- c('CUAA', 'CUAC', 'CUAG', 'CUAU', 'CUCA', 'CUCC', 'CUCG', 'CUCU', 'CUGA', 'CUGC', 'CUGG', 'CUGU', 'CUUA', 'CUUC', 'CUUG', 'CUUU', 'GAAA', 'GAAC', 'GAAG', 'GAAU', 'GACA', 'GACC', 'GACG', 'GACU', 'GAGA', 'GAGC', 'GAGG', 'GAGU', 'GAUA', 'GAUC', 'GAUG', 'GAUU', 'AGAA', 'AGAC', 'AGAG', 'AGAU', 'AGCA', 'AGCC', 'AGCG', 'AGCU', 'AGGA', 'AGGC', 'AGGG', 'AGGU', 'AGUA', 'AGUC', 'AGUG', 'AGUU', 'UCAA', 'UCAC', 'UCAG', 'UCAU', 'UCCA', 'UCCC', 'UCCG', 'UCCU', 'UCGA', 'UCGC', 'UCGG', 'UCGU', 'UCUA', 'UCUC', 'UCUG', 'UCUU', 'GUAA', 'GUAC', 'GUAG', 'GUAU', 'GUCA', 'GUCC', 'GUCG', 'GUCU', 'GUGA', 'GUGC', 'GUGG', 'GUGU', 'GUUA', 'GUUC', 'GUUG', 'GUUU', 'CAAA', 'CAAC', 'CAAG', 'CAAU', 'CACA', 'CACC', 'CACG', 'CACU', 'CAGA', 'CAGC', 'CAGG', 'CAGU', 'CAUA', 'CAUC', 'CAUG', 'CAUU', 'GCAA', 'GCAC', 'GCAG', 'GCAU', 'GCCA', 'GCCC', 'GCCG', 'GCCU', 'GCGA', 'GCGC', 'GCGG', 'GCGU', 'GCUA', 'GCUC', 'GCUG', 'GCUU', 'CGAA', 'CGAC', 'CGAG', 'CGAU', 'CGCA', 'CGCC', 'CGCG', 'CGCU', 'CGGA', 'CGGC', 'CGGG', 'CGGU', 'CGUA', 'CGUC', 'CGUG', 'CGUU', 'ACAA', 'ACAC', 'ACAG', 'ACAU', 'ACCA', 'ACCC', 'ACCG', 'ACCU', 'ACGA', 'ACGC', 'ACGG', 'ACGU', 'ACUA', 'ACUC', 'ACUG', 'ACUU', 'UGAA', 'UGAC', 'UGAG', 'UGAU', 'UGCA', 'UGCC', 'UGCG', 'UGCU', 'UGGA', 'UGGC', 'UGGG', 'UGGU', 'UGUA', 'UGUC', 'UGUG', 'UGUU', 'AUAA', 'AUAC', 'AUAG', 'AUAU', 'AUCA', 'AUCC', 'AUCG', 'AUCU', 'AUGA', 'AUGC', 'AUGG', 'AUGU', 'AUUA', 'AUUC', 'AUUG', 'AUUU', 'UAAA', 'UAAC', 'UAAG', 'UAAU', 'UACA', 'UACC', 'UACG', 'UACU', 'UAGA', 'UAGC', 'UAGG', 'UAGU', 'UAUA', 'UAUC', 'UAUG', 'UAUU')
complement_mutation <- c('GAUU', 'GAGU', 'GACU', 'GAAU', 'GAUG', 'GAGG', 'GACG', 'GAAG', 'GAUC', 'GAGC', 'GACC', 'GAAC', 'GAUA', 'GAGA', 'GACA', 'GAAA', 'CUUU', 'CUGU', 'CUCU', 'CUAU', 'CUUG', 'CUGG', 'CUCG', 'CUAG', 'CUUC', 'CUGC', 'CUCC', 'CUAC', 'CUUA', 'CUGA', 'CUCA', 'CUAA', 'UCUU', 'UCGU', 'UCCU', 'UCAU', 'UCUG', 'UCGG', 'UCCG', 'UCAG', 'UCUC', 'UCGC', 'UCCC', 'UCAC', 'UCUA', 'UCGA', 'UCCA', 'UCAA', 'AGUU', 'AGGU', 'AGCU', 'AGAU', 'AGUG', 'AGGG', 'AGCG', 'AGAG', 'AGUC', 'AGGC', 'AGCC', 'AGAC', 'AGUA', 'AGGA', 'AGCA', 'AGAA', 'CAUU', 'CAGU', 'CACU', 'CAAU', 'CAUG', 'CAGG', 'CACG', 'CAAG', 'CAUC', 'CAGC', 'CACC', 'CAAC', 'CAUA', 'CAGA', 'CACA', 'CAAA', 'GUUU', 'GUGU', 'GUCU', 'GUAU', 'GUUG', 'GUGG', 'GUCG', 'GUAG', 'GUUC', 'GUGC', 'GUCC', 'GUAC', 'GUUA', 'GUGA', 'GUCA', 'GUAA', 'CGUU', 'CGGU', 'CGCU', 'CGAU', 'CGUG', 'CGGG', 'CGCG', 'CGAG', 'CGUC', 'CGGC', 'CGCC', 'CGAC', 'CGUA', 'CGGA', 'CGCA', 'CGAA', 'GCUU', 'GCGU', 'GCCU', 'GCAU', 'GCUG', 'GCGG', 'GCCG', 'GCAG', 'GCUC', 'GCGC', 'GCCC', 'GCAC', 'GCUA', 'GCGA', 'GCCA', 'GCAA', 'UGUU', 'UGGU', 'UGCU', 'UGAU', 'UGUG', 'UGGG', 'UGCG', 'UGAG', 'UGUC', 'UGGC', 'UGCC', 'UGAC', 'UGUA', 'UGGA', 'UGCA', 'UGAA', 'ACUU', 'ACGU', 'ACCU', 'ACAU', 'ACUG', 'ACGG', 'ACCG', 'ACAG', 'ACUC', 'ACGC', 'ACCC', 'ACAC', 'ACUA', 'ACGA', 'ACCA', 'ACAA', 'UAUU', 'UAGU', 'UACU', 'UAAU', 'UAUG', 'UAGG', 'UACG', 'UAAG', 'UAUC', 'UAGC', 'UACC', 'UAAC', 'UAUA', 'UAGA', 'UACA', 'UAAA', 'AUUU', 'AUGU', 'AUCU', 'AUAU', 'AUUG', 'AUGG', 'AUCG', 'AUAG', 'AUUC', 'AUGC', 'AUCC', 'AUAC', 'AUUA', 'AUGA', 'AUCA', 'AUAA')
spectrum_df <- read_tsv("../../data/all_mpxv_sars_phyvirus_sample.txt") %>% select(c(name, brLen,standard_mutation))
colnames(spectrum_df) <- c("name", "brLen", standard_mutation)
history_df <- read_csv("../../data/all_meta_corrected_host.csv")
data_df <- history_df %>% select(accession, year, corrected_host, species, family, order) %>% 
            right_join(spectrum_df, c("accession" = "name"), relationship = "many-to-many")  %>% na.omit() %>% 
            mutate_at(vars(-c("accession","year","corrected_host", "species","family","order")),as.numeric) %>% 
            filter(brLen <= 0.5) %>% left_join(history_df %>% select(family, species) %>% distinct(family, species)) %>% 
            select(family, everything()) %>% mutate(family = gsub("viridae", "", family)) 
#for -ssRNA virus, we convert the raw spectrum to the positive strand
corrected_data_df <- data.frame()
for (order_name in (data_df %>% .$order %>% unique()))  {
  small_df1 <- data_df %>% filter(order==order_name) 
  if (order_name %in% c("Mononegavirales", "Articulavirales", "ssRNA_-_") ){
    small_df2 <- small_df1
    colnames(small_df2) <- c("family","accession","year","corrected_host","species","order","brLen", complement_mutation)
    small_df3 <- small_df2[c("family","accession","year","corrected_host","species","order","brLen", standard_mutation)]
  }
  else{
    small_df3 <- small_df1 
  }
  corrected_data_df <- rbind(small_df3, corrected_data_df)
}
species_vector <-  c("InfluenzaAvirus_h5n1_",
                       "InfluenzaAvirus_h3n2_",
                       "Lyssavirus_rabies",
                       "Mumps_orthorubulavirus",
                       "Zaire_ebolavirus",
                       "West_Nile_virus",
                       "Zika_virus",
                       "sars-cov-2")
phyvirus_species_vector <- data_df  %>% filter(order %in% c("ssRNA_-_","ssRNA_+_")) %>%
                        select(species) %>% unique() %>% as.vector() %>% unlist()
species_vector <- species_vector  %>% append(phyvirus_species_vector) %>% unique()
data_corrected_data_df <- corrected_data_df %>% filter(species %in% species_vector)
grouped_df <- data_corrected_data_df %>%filter(corrected_host != "Others") %>% 
              group_by(corrected_host, family, species) %>% 
              summarise_at(colnames(data_df)[8:199],sum) %>% ungroup()
#adjust viral spectrum to human genome composition
df <- read_tsv("../../data/all_phyvirus_human_compos.txt") %>% filter(name != "name") %>% mutate(across(-1, ~as.numeric(.)))
normalise_df  <- df %>% rowwise() %>% mutate(across(-1, ~ ./ sum(c_across(-1))))
Human_row <- normalise_df[normalise_df$name == "Human", -1] %>% as.vector()
df_result <- normalise_df %>% select(-name)
df_result2 <- Human_row / df_result  
df_result_normalized <- cbind(name = normalise_df$name, df_result2)
correct_spectrum <- function(factor_df, spectrum_df, species_name) {
  factor_row <- factor_df[factor_df$name == species_name, -1] %>% mutate_all(as.numeric) %>% as.vector() %>% unlist()
  species_spectrum_df <- spectrum_df %>% filter(species == species_name) 
  if (nrow(species_spectrum_df) == 0) {
    return()
  }
  spectrum_df <- species_spectrum_df %>% select(-c(corrected_host,species,family)) %>% 
    mutate_all(as.numeric)
  corrected_spectrum_df  <- apply(spectrum_df, 1, function(row) row * factor_row)
  corrected_spectrum_df  <- as.data.frame(t(corrected_spectrum_df)) %>% mutate_all(round)
  name_corrected_spectrum_df <- species_spectrum_df %>% select(corrected_host,species,family)
  all_corrected_spectrum_df <- cbind(name_corrected_spectrum_df, corrected_spectrum_df)
  return(all_corrected_spectrum_df)
} 
all_corrected_df <- data.frame()
for (specie in species_vector){
  # print(specie)
  corrected_df <- correct_spectrum(df_result_normalized, grouped_df, specie)
  all_corrected_df <- rbind(all_corrected_df , corrected_df)
}

###supplementary figure1A
###Pairwise cosine similarity of viral mutational spectrums between different hosts of the same viral family
calculate_similarity_family <- function (new_grouped_df) {
  order_vector <- new_grouped_df$family %>% unique()
  order_test_df <- c()
  for (order_name in order_vector ) {
    test_df <- new_grouped_df  %>% filter(family == order_name) 
    host_vector <- unique(test_df$corrected_host)
    if (length(host_vector) <=1) {
      next 
    }
    host_combine <- combn(host_vector, 2)
    for (index in rep(1:ncol(host_combine),1)) {
      host_name <-  host_combine[, index]
      matrix <- test_df %>% filter(corrected_host %in% host_name) %>% ungroup() %>% select(3:194)
      cosine_similarity <- proxy::simil(x = matrix[1,]/sum(matrix[1,]), y = matrix[2,]/sum(matrix[2,]), method = "cosine")
      result_df <- c(order_name, host_name[1], host_name[2], cosine_similarity)
      order_test_df <- rbind(order_test_df, result_df)
    }
  }
  colnames(order_test_df) <- c("order", "hostA", "hostB",  "cosine_similarity") 
  order_test_df <- as.data.frame(order_test_df) 
  order_test_df <- rownames_to_column(order_test_df, var = "rowname") %>% select(-rowname) %>% 
                    mutate_at("cosine_similarity", as.numeric)
  return(order_test_df )
}
family_grouped_df <- all_corrected_df %>%
                    filter(corrected_host != "Others") %>% 
                    group_by(corrected_host, family) %>% 
                    summarise_at(colnames(data_df)[8:199],sum) 
family_grouped_df_100 <- family_grouped_df %>% ungroup() %>% filter(rowSums(select(., 3:194))>100) 
family_df_100 <- calculate_similarity_family(family_grouped_df_100)
p1 <- ggplot(family_df_100, aes(x = factor(order), y = cosine_similarity)) +
      geom_boxplot(outlier.color=NA) +
      geom_point(
        position=position_jitter(width = 0.2), 
        pch=21, #圆形
        size=1,
        fill = "#F39B7FFF")+
      labs(
        x = "Viral family",
        y = "Pairwise cosine similarity \n between hosts' mutation spectra")+
      my_theme()+
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("suppl_figure1A.pdf",width = 5.5, height = 6, units = "cm",dpi = 300)


###supplementary figure1B
###Pairwise cosine similarity of viral mutational spectrums between different viral families of the same host
calculate_similarity_host <- function (new_grouped_df) {
  order_vector <- new_grouped_df$corrected_host %>% unique()
  order_test_df <- c()
  for (host_name in order_vector ) {
    test_df <- new_grouped_df  %>% filter(corrected_host == host_name) 
    family_vector <- unique(test_df$family)
    if (length(family_vector) <=1) {
      next 
    }
    family_combine <- combn(family_vector, 2)
    for (index in rep(1:ncol(family_combine),1)) {
      family_name <- family_combine[, index]
      matrix <- test_df %>% filter(family %in% family_name ) %>% ungroup() %>% select(3:194)
      cosine_similarity <- proxy::simil(x = matrix[1,]/sum(matrix[1,]), y = matrix[2,]/sum(matrix[2,]), method = "cosine")
      result_df <- c(host_name, family_name[1], family_name[2], cosine_similarity)
      order_test_df <- rbind(order_test_df, result_df)
    }
  }
  colnames(order_test_df) <- c("host", "familyA", "familyB", "cosine_similarity") 
  order_test_df <- as.data.frame(order_test_df)
  order_test_df <- rownames_to_column(order_test_df, var = "rowname") %>% select(-rowname) %>%
                   mutate_at("cosine_similarity", as.numeric)
  return(order_test_df )
}
host_df_100  <- calculate_similarity_host(family_grouped_df_100)
p2 <- ggplot(host_df_100 , aes(x = factor(host), y = cosine_similarity)) +
  geom_boxplot(outlier.color=NA) +
  geom_point(
    position=position_jitter(width = 0.2), 
    pch=21, #圆形
    size=1,
    fill = "#F39B7FFF")+
  scale_x_discrete(labels=c("Aves", "Homo sapiens",
                            "Insecta", "Mammalia", "Reptilia"))+
  labs(
    x = "Host",
    y = "Pairwise cosine similarity \n between families' mutation spectra")+
  my_theme()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("suppl_figure1B.pdf",width = 5.5, height = 6, units = "cm",dpi = 300)


###supplementary figure1C
#PCA plot of mutational spectrums for 13 viral families (colored by viral families)
df <- read_tsv("../../data/context_all_signature_sample_csv.txt") 
data_df <- df  %>% mutate(family = gsub("viridae", "", family)) %>% filter(host != "Others")
num_cores=5
df_normalized <- mclapply(seq_len(nrow(data_df)), function(i) {
  row <- data_df[i, ]
  row[-(1:6)] <- row[-(1:6)] / sum(row[-(1:6)])
  return(row)
}, mc.cores = num_cores)
df_normalized <- do.call(rbind, df_normalized) 
family_vector <- c("Arena", "Filo", "Hanta", "Orthomyxo","Paramyxo", "Peribunya", 
                   "Phenui", "Rhabdo", "Calici", "Corona", "Flavi","Picorna", "Toga")
balti_dict <- setNames(c(rep("-ssRNA", times = 8), rep("+ssRNA", times = 5)),family_vector)
df_normalized <- df_normalized %>% filter(host != "internal") 
df_normalized <- df_normalized %>% filter(family != "Pox") %>% rowwise() %>% 
                mutate(baltimore_class= balti_dict[[family]]) %>% select(baltimore_class,everything())
pca_df <- df_normalized[,(ncol(df_normalized)-191):ncol(df_normalized)] 
pca_result <- prcomp(pca_df)
pc1 <- pca_result$sdev[1]^2 / sum(pca_result$sdev^2)
pc2 <- pca_result$sdev[2]^2 / sum(pca_result$sdev^2)
pca_plot_data <- as_tibble(get_pca_ind(pca_result)$coord) %>% 
                select("Dim.1", "Dim.2", "Dim.3") %>% 
                cbind(df_normalized[,1:(ncol(df_normalized)-192)])
mypal <- pal_npg("nrc")(10)
mypal <- c(mypal, c("#4DAF4A" , "#999999" , "#A65628"))
p3 <- ggplot(pca_plot_data, aes(x = Dim.1, y = Dim.2)) +
      geom_point(aes(color = family,shape=baltimore_class),alpha=0.7,size=0.5) +
      scale_shape_manual(values = c("-ssRNA" = 16, 
                                    "+ssRNA" = 17))+
      scale_color_manual(values = mypal)+
      labs(x = paste("PC1","(",round(pc1*100,1),"%)",sep="") , y = paste("PC2","(",round(pc2*100,1),"%)",sep="")) +
      my_theme()+
      labs(color="Viral family",
           shape="Baltimore class")
ggsave("suppl_figure1C.pdf",width = 8, height = 6, units = "cm",dpi = 300)

###supplementary figure1D
#PCA plot of mutational spectrums for viral species across 13 viral families (colored by hosts)
p2 <- ggplot(pca_plot_data, aes(x = Dim.1, y = Dim.2)) +
      geom_point(aes(color = host,shape=baltimore_class),alpha=0.7,size=0.5) +
      scale_shape_manual(values = c("-ssRNA" = 16, 
                                    "+ssRNA" = 17))+
      scale_color_npg(labels=c("Aves", "Homo sapiens",
                               "Insecta", "Mammalia", "Reptilia"))+
      labs(x = paste("PC1","(",round(pc1*100,1),"%)",sep="") , y = paste("PC2","(",round(pc2*100,1),"%)",sep="")) +
      my_theme()+
      labs(color="Host",
           shape="Baltimore class")
ggsave("suppl_figure1D.pdf",width = 8, height = 6, units = "cm",dpi = 300)

###supplementary figure1E
#Percentage of variance explained by viral family and hosts for different principal components
model3 <- lm(Dim.1 ~ family + host, data = pca_plot_data)
pc1_variance <- Anova(model3, type = "III")[[1]]/sum(Anova(model3, type = "III")[[1]])[1]
pc1_variance <- pc1_variance[2:3]

model4 <- lm(Dim.2 ~ family + host, data = pca_plot_data)
pc2_variance <- Anova(model4, type = "III")[[1]]/sum(Anova(model4, type = "III")[[1]])[1]
pc2_variance <- pc2_variance[2:3]

model5 <- lm(Dim.3 ~ family + host, data = pca_plot_data)
pc3_variance <- Anova(model5, type = "III")[[1]]/sum(Anova(model5, type = "III")[[1]])[1]
pc3_variance <- pc3_variance[2:3]
variance_df <- data.frame(variable=c("viral family","host"),
                          PC1=pc1_variance,
                          PC2=pc2_variance,
                          PC3=pc3_variance)
longer_variance_df <- variance_df %>% pivot_longer(!variable, names_to="pc_axis",values_to="proportion")
plot2 <- ggplot(longer_variance_df, aes(x = pc_axis , y = proportion, fill=variable)) + 
          geom_col(position = "dodge") +
          scale_fill_npg(labels=c("Host", "Viral family"))+
          labs(x="PC dimensions",
               y="Proportion of variance\n explained by a variable",
               fill="Variable")+
          my_theme()
ggsave("suppl_figure1E.pdf", plot2, width = 8, height = 6, units = "cm",dpi = 300)







