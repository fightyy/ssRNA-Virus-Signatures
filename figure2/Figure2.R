library(tidyverse)
library(this.path)
library(parallel)
library(ggsci)
library(factoextra)
setwd(this.dir())
source("../bin/theme_setup.R")

###Figure2A
#spectra of SARS-CoV-2 at the early and late period
standard_mutation <- c('CUAA', 'CUAC', 'CUAG', 'CUAU', 'CUCA', 'CUCC', 'CUCG', 'CUCU', 'CUGA', 'CUGC', 'CUGG', 'CUGU', 'CUUA', 'CUUC', 'CUUG', 'CUUU', 'GAAA', 'GAAC', 'GAAG', 'GAAU', 'GACA', 'GACC', 'GACG', 'GACU', 'GAGA', 'GAGC', 'GAGG', 'GAGU', 'GAUA', 'GAUC', 'GAUG', 'GAUU', 'AGAA', 'AGAC', 'AGAG', 'AGAU', 'AGCA', 'AGCC', 'AGCG', 'AGCU', 'AGGA', 'AGGC', 'AGGG', 'AGGU', 'AGUA', 'AGUC', 'AGUG', 'AGUU', 'UCAA', 'UCAC', 'UCAG', 'UCAU', 'UCCA', 'UCCC', 'UCCG', 'UCCU', 'UCGA', 'UCGC', 'UCGG', 'UCGU', 'UCUA', 'UCUC', 'UCUG', 'UCUU', 'GUAA', 'GUAC', 'GUAG', 'GUAU', 'GUCA', 'GUCC', 'GUCG', 'GUCU', 'GUGA', 'GUGC', 'GUGG', 'GUGU', 'GUUA', 'GUUC', 'GUUG', 'GUUU', 'CAAA', 'CAAC', 'CAAG', 'CAAU', 'CACA', 'CACC', 'CACG', 'CACU', 'CAGA', 'CAGC', 'CAGG', 'CAGU', 'CAUA', 'CAUC', 'CAUG', 'CAUU', 'GCAA', 'GCAC', 'GCAG', 'GCAU', 'GCCA', 'GCCC', 'GCCG', 'GCCU', 'GCGA', 'GCGC', 'GCGG', 'GCGU', 'GCUA', 'GCUC', 'GCUG', 'GCUU', 'CGAA', 'CGAC', 'CGAG', 'CGAU', 'CGCA', 'CGCC', 'CGCG', 'CGCU', 'CGGA', 'CGGC', 'CGGG', 'CGGU', 'CGUA', 'CGUC', 'CGUG', 'CGUU', 'ACAA', 'ACAC', 'ACAG', 'ACAU', 'ACCA', 'ACCC', 'ACCG', 'ACCU', 'ACGA', 'ACGC', 'ACGG', 'ACGU', 'ACUA', 'ACUC', 'ACUG', 'ACUU', 'UGAA', 'UGAC', 'UGAG', 'UGAU', 'UGCA', 'UGCC', 'UGCG', 'UGCU', 'UGGA', 'UGGC', 'UGGG', 'UGGU', 'UGUA', 'UGUC', 'UGUG', 'UGUU', 'AUAA', 'AUAC', 'AUAG', 'AUAU', 'AUCA', 'AUCC', 'AUCG', 'AUCU', 'AUGA', 'AUGC', 'AUGG', 'AUGU', 'AUUA', 'AUUC', 'AUUG', 'AUUU', 'UAAA', 'UAAC', 'UAAG', 'UAAU', 'UACA', 'UACC', 'UACG', 'UACU', 'UAGA', 'UAGC', 'UAGG', 'UAGU', 'UAUA', 'UAUC', 'UAUG', 'UAUU')
df <- read_tsv("../data/all_signature_sample_csv.txt") %>% filter(! order %in% c("ssRNA_+_", "ssRNA_-_")) %>% filter(host != "internal") %>% 
      rowwise()%>% mutate(start_year = suppressWarnings(as.numeric(strsplit(year, "-")[[1]][1])),
                          end_year = suppressWarnings(as.numeric(strsplit(year, "-")[[1]][2]))) %>%
      na.omit() %>% mutate(mean_between_years = mean(c(start_year, end_year))) %>% 
      select(c(start_year, end_year, mean_between_years), everything())
df_sars <- df %>% filter(species=="sars-cov-2") %>% arrange(mean_between_years) 
df_sars_early <- df_sars %>% filter(mean_between_years <= 10) %>%
                select((ncol(df_sars)-191):ncol(df_sars)) %>% colSums() %>% as.data.frame()  %>% 
                mutate(across(starts_with("."), ~ . / sum(`.`))) %>%
                mutate(MutationType=standard_mutation) %>% dplyr::rename("SARS-CoV-2\n(early)"=".")
df_sars_late <- df_sars %>% filter(mean_between_years >= 38) %>%
                select((ncol(df_sars)-191):ncol(df_sars)) %>% colSums() %>% as.data.frame()  %>% 
                mutate(across(starts_with("."), ~ . / sum(`.`))) %>%
                mutate(MutationType=standard_mutation) %>% dplyr::rename("SARS-CoV-2\n(late)"=".")
df_plot_sars <- df_sars_early %>% right_join(df_sars_late, by="MutationType")
df_plot_sars$muttype <- apply(df_plot_sars ,1,function(x){
  return(paste0(substring(x["MutationType"],1,1),
                ">",
                substring(x["MutationType"],2,2)))
})
df_plot_sars$trinuc <- apply(df_plot_sars,1,function(x){
  return(paste0(substring(x["MutationType"],3,3),
                substring(x["MutationType"],1,1),
                substring(x["MutationType"],4,4)))
})
long_df_plot_sars <- df_plot_sars %>% 
                    pivot_longer( cols = contains("Sars"),
                                  names_to = "signature",
                                  values_to = "exposure")
long_df_plot_sars$muttype <- factor(long_df_plot_sars$muttype,
                                    levels=c("C>U","G>A","A>G","U>C",
                                             "G>U","C>A","G>C","C>G",
                                             "A>C","U>G","A>U","U>A"))
mypal <- pal_npg("nrc")(10)
mypal <- c(mypal, c("#4DAF4A" , "#999999" , "#A65628"))
plot_sars <- ggplot(data =long_df_plot_sars , aes(x = trinuc, y = exposure, fill = muttype)) +
              geom_bar(stat = "identity") + 
              facet_grid(signature ~ muttype, scales = "free_x", space = "free_x") +
              theme_bw() +
              ylab("Proportion of mutation") +
              scale_fill_manual(values = mypal) +
              scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
              my_theme()+
              theme(  
                panel.spacing.x = unit(0, "lines"),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                # axis.title.x = element_blank(),
                axis.ticks.x = element_blank(), 
                axis.title.x = element_blank(),
                axis.text.x =element_blank(),
                axis.text.y =element_text(size = 7),
                axis.title.y = element_text(size = 9),
                strip.text = element_text(size = 8),
                legend.position = "none",
                plot.margin = margin(0, 0, 0, 0, "null"))
ggsave("Figure2A.pdf",width =16, height = 5, units = "cm",dpi = 300)


###Figure2D
#spectra of SARS-CoV-2 at the early and late period
df_h3n2 <- df %>% filter(species=="InfluenzaAvirus_h3n2_")
df_h3n2_early <- df_h3n2 %>% filter(mean_between_years <= 2004 && mean_between_years >= 2000) %>%
                select((ncol(df_h3n2)-191):ncol(df_h3n2)) %>% colSums() %>% as.data.frame()  %>% 
                mutate(across(starts_with("."), ~ . / sum(`.`))) %>%
                mutate(MutationType=standard_mutation) %>% dplyr::rename("H3N2\n(early)"=".")
df_h3n2_late <- df_h3n2 %>% filter(mean_between_years >= 2019 && mean_between_years <= 2023) %>%
                select((ncol(df_h3n2)-191):ncol(df_h3n2)) %>% colSums() %>% as.data.frame()  %>% 
                mutate(across(starts_with("."), ~ . / sum(`.`))) %>%
                mutate(MutationType=standard_mutation) %>% dplyr::rename("H3N2\n(late)"=".")
df_plot_h3n2 <- df_h3n2_early %>% right_join(df_h3n2_late,by="MutationType")
df_plot_h3n2$muttype <- apply(df_plot_h3n2 ,1,function(x){
  return(paste0(substring(x["MutationType"],1,1),
                ">",
                substring(x["MutationType"],2,2)))
})
df_plot_h3n2$trinuc <- apply(df_plot_h3n2,1,function(x){
  return(paste0(substring(x["MutationType"],3,3),
                substring(x["MutationType"],1,1),
                substring(x["MutationType"],4,4)))
})
long_df_plot_h3n2 <- df_plot_h3n2 %>% 
                    pivot_longer( cols = contains("H3N2"),
                                  names_to = "signature",
                                  values_to = "exposure")
long_df_plot_h3n2$muttype <- factor(long_df_plot_h3n2$muttype,
                                    levels=c("C>U","G>A","A>G","U>C",
                                             "G>U","C>A","G>C","C>G",
                                             "A>C","U>G","A>U","U>A"))
mypal <- pal_npg("nrc")(10)
mypal <- c(mypal, c("#4DAF4A" , "#999999" , "#A65628"))
plot_h3n2 <- ggplot(data =long_df_plot_h3n2 , aes(x = trinuc, y = exposure, fill = muttype)) +
            geom_bar(stat = "identity") + 
            facet_grid(signature ~ muttype, scales = "free_x", space = "free_x") +
            theme_bw() +
            ylab("Proportion of mutation") +
            scale_fill_manual(values = mypal) +
            scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
            my_theme()+
            theme(  
              panel.spacing.x = unit(0, "lines"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              # axis.title.x = element_blank(),
              axis.ticks.x = element_blank(), 
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y =element_text(size = 7),
              axis.title.y = element_text(size = 9),
              strip.text = element_text(size = 8),
              legend.position = "none",
              plot.margin = margin(0, 0, 0, 0, "null"))
ggsave("Figure2D.pdf",width =16, height = 4.5, units = "cm",dpi = 300)

###Figure2BCEF
#prepare pca data
df <- read_tsv("../data/all_signature_sample_csv.txt") 
data_df <- df  %>% mutate(family = gsub("viridae", "", family)) %>% filter(host != "Others")
num_cores=5
df_normalized <- mclapply(seq_len(nrow(data_df)), function(i) {
  row <- data_df[i, ]
  row[-(1:6)] <- row[-(1:6)] / sum(row[-(1:6)])
  return(row)
}, mc.cores = num_cores)
df_normalized <- do.call(rbind, df_normalized) 
df <- df_normalized  %>% filter(! order %in% c("ssRNA_+_", "ssRNA_-_")) %>% filter(host != "internal") %>% 
  filter(species %in% c("sars-cov-2","InfluenzaAvirus_h3n2_"))
for (specie in unique(df$species)) {
  #figurr2B and Figure2E
  #pca plot of sars-cov-2 and h3n2
  df_normalized <- df %>% filter(species==specie)
  short_good_df <- df_normalized 
  pca_df <- short_good_df %>% select((ncol(short_good_df)-191):ncol(short_good_df))
  pca_result <- prcomp(pca_df)
  pc1 <- pca_result$sdev[1]^2 / sum(pca_result$sdev^2)
  pc2 <- pca_result$sdev[2]^2 / sum(pca_result$sdev^2)
  pca_plot_data <- as_tibble(get_pca_ind(pca_result)$coord) %>%
                   select("Dim.1", "Dim.2") %>%
                   cbind(short_good_df[,1:(ncol(short_good_df)-192)])
  
  new_pca_plot_data <- suppressWarnings(pca_plot_data %>% rowwise()%>% mutate(start_year = as.numeric(strsplit(year, "-")[[1]][1]),
                                                                                end_year = as.numeric(strsplit(year, "-")[[1]][2])) %>%
                                            na.omit() %>% mutate(mean_between_years = mean(c(start_year, end_year)))) 
  if(specie == "sars-cov-2") {
    file_name <- "Figure2B.pdf"
    standard_name <- "SARS-CoV-2"
    standard_time <- "Month since\n  2019/12"
  }
  else{
    file_name <- "Figure2E.pdf"
    standard_name <- "H3N2"
    standard_time <- "Year"
    new_pca_plot_data <- new_pca_plot_data %>% filter(mean_between_years> 1995)
  }
  
  plot1 <- ggplot(new_pca_plot_data, aes(x=Dim.1, y=Dim.2))+
    geom_point(aes(color = mean_between_years)) +
    scale_color_gradient(low = "blue" , high = "red") +
    labs(title=standard_name) +
    labs(x = paste("PC1","(",round(pc1*100,1),"%)",sep="") , y = paste("PC2","(",round(pc2*100,1),"%)",sep="")) +
    labs(color=standard_time)+
    my_theme()
  ggsave(file_name, plot1, width = 8, height = 6, units = "cm",dpi = 300)
  
  #figurr2C and Figure2F
  #pca loading of sars-cov-2 and h3n2
  pc1_loading <- pca_result$rotation[, 1]
  sorted_pc1_loading <- pc1_loading[order(-abs(pc1_loading))] %>% as.data.frame() %>% dplyr::rename("loading"=".")
  sorted_pc1_loading$Mutationtype <- rownames(sorted_pc1_loading)
  sorted_pc1_loading <- sorted_pc1_loading %>% rowwise() %>% 
                        mutate(trinuc=paste0(substring(Mutationtype, 3, 3),
                                             substring(Mutationtype, 1, 1),
                                             substring(Mutationtype, 4, 4))) %>% 
                        mutate(muttype=paste0(substring(Mutationtype,1,1),
                                              ">",
                                              substring(Mutationtype,2,2))) %>%
                        mutate(Mutationtype=paste0(substring(Mutationtype, 3, 3),
                                                   "[",
                                                   substring(Mutationtype, 1, 1),
                                                   ">",
                                                   substring(Mutationtype, 2, 2),
                                                   "]",
                                                   substring(Mutationtype, 4, 4)))
                      
  sorted_pc1_loading$Mutationtype <- factor(sorted_pc1_loading$Mutationtype,
                                            levels=sorted_pc1_loading$Mutationtype)
  
  mypal <- pal_npg("nrc")(10)
  mypal <- c(mypal, c("#4DAF4A" , "#999999" , "#A65628"))
  muttype_vector <- c("C>U","G>A","A>G","U>C",
                      "G>U","C>A","G>C","C>G",
                      "A>C","U>G","A>U","U>A")
  color_dict <- setNames(mypal,muttype_vector)
  if(specie == "sars-cov-2") {
    file_name2 <- "Figure2C.pdf"
  }
  else{
    file_name2 <- "Figure2F.pdf"
  }
  plot2 <- ggplot(sorted_pc1_loading[1:10,], aes(x = factor(Mutationtype), y = loading)) + 
            geom_col(aes(fill=muttype)) +
            scale_fill_manual(values=color_dict)+
            labs(x="Mutation type",
                 y="PC1 loading",
                 title=standard_name,
                 fill="Mutation type")+
            my_theme()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  axis.title.x = element_text(size = 7, margin = margin(r=0.1,unit="cm")),
                  axis.title.y = element_text(size = 7, margin = margin(t=0.2,unit="cm"))
            )
  ggsave(file_name2, plot2, width = 8, height = 6, units = "cm",dpi = 300)
}



