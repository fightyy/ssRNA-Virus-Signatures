library(tidyverse)
library(this.path)
library(ggsci)
setwd(this.dir())
source("../../bin/theme_setup.R")

###supplementary figure3A
###nonsynonymous mutational spectra for different iSNVs and SNPs of the SARS-CoV-2
SBS192.input.catalog <- read_csv("../data/sars_nosnyn_spectrum.csv") 
matrix <- SBS192.input.catalog %>% select(-1)
normalized_matrix <- apply(matrix,1, function(x) x/colSums(matrix)) %>% t() %>% as.data.frame()
normalized_matrix$Mutationtype <-  SBS192.input.catalog$Mutationtype
SBS192.input.catalog <- normalized_matrix
SBS192.input.catalog$muttype <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["Mutationtype"],1,1),
                ">",
                substring(x["Mutationtype"],2,2)))
})
SBS192.input.catalog$trinuc <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["Mutationtype"],3,3),
                substring(x["Mutationtype"],1,1),
                substring(x["Mutationtype"],4,4)))
})

long_rename_SBS192.input.catalog <- SBS192.input.catalog %>% 
                                    pivot_longer( cols = contains("count"),
                                                  names_to = "signature",
                                                  values_to = "exposure")

long_rename_SBS192.input.catalog$muttype <- factor(long_rename_SBS192.input.catalog$muttype,
                                                   levels=c("C>U","G>A","A>G","U>C",
                                                            "G>U","C>A","G>C","C>G",
                                                            "A>C","U>G","A>U","U>A"))
original_signature <- c("sars_isnv_0.2_0-count","sars_isnv_0.4_0.2-count","sars_isnv_0.6_0.4-count",
                        "sars_isnv_0.8_0.6-count","sars_isnv_0.95_0.8-count","sars_polymorphism_count",
                        "sars_polymorphism_non_syn-count","sars_polymorphism_syn-count") 
correct_signature <- c("iSNV\n(0.05-0.2)","iSNV\n(0.2-0.4)","iSNV\n(0.4-0.6)",
                       "iSNV\n(0.6-0.8)","iSNV\n(0.8-0.95)","polymorphism","polymorphism","polymorphism")
signature_dict <- setNames(correct_signature,original_signature)
long_rename_SBS192.input.catalog <- long_rename_SBS192.input.catalog %>% rowwise() %>%
                                    mutate(signature=signature_dict[[signature]])
mypal <- pal_npg("nrc")(10)
mypal <- c(mypal, c("#4DAF4A" , "#999999" , "#A65628"))
plot <- ggplot(data =long_rename_SBS192.input.catalog , aes(x = trinuc, y = exposure, fill = muttype)) +
  geom_bar(stat = "identity") + 
  facet_grid(signature ~ muttype, scales = "free_x", space = "free_x") +
  theme_bw() +
  ylab("Proportion of mutation") +
  labs(title="Mutation spectra of nonsynonymous mutation")+
  scale_fill_manual(values = mypal) +
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
  my_theme()+
  theme(  
    panel.spacing.x = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y =element_text(size = 7),
    axis.title.y = element_text(size = 9),
    strip.text = element_text(size = 5),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0, "null"))
ggsave("suppl_figure3A.pdf",width =20, height = 9, units = "cm",dpi = 300)


###supplementary figure3B
###synonymous mutational spectra for different iSNVs and SNPs of the SARS-CoV-2
SBS192.input.catalog <- read_csv("../data/sars_syn_spectrum.csv") 
matrix <- SBS192.input.catalog %>% select(-1)
normalized_matrix <- apply(matrix,1, function(x) x/colSums(matrix)) %>% t() %>% as.data.frame()
normalized_matrix$Mutationtype <-  SBS192.input.catalog$Mutationtype
SBS192.input.catalog <- normalized_matrix
SBS192.input.catalog$muttype <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["Mutationtype"],1,1),
                ">",
                substring(x["Mutationtype"],2,2)))
})
SBS192.input.catalog$trinuc <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["Mutationtype"],3,3),
                substring(x["Mutationtype"],1,1),
                substring(x["Mutationtype"],4,4)))
})

long_rename_SBS192.input.catalog <- SBS192.input.catalog %>% 
                                    pivot_longer( cols = contains("count"),
                                                  names_to = "signature",
                                                  values_to = "exposure")

long_rename_SBS192.input.catalog$muttype <- factor(long_rename_SBS192.input.catalog$muttype,
                                                   levels=c("C>U","G>A","A>G","U>C",
                                                            "G>U","C>A","G>C","C>G",
                                                            "A>C","U>G","A>U","U>A"))
original_signature <- c("sars_isnv_0.2_0-count","sars_isnv_0.4_0.2-count","sars_isnv_0.6_0.4-count",
                        "sars_isnv_0.8_0.6-count","sars_isnv_0.95_0.8-count","sars_polymorphism_count",
                        "sars_polymorphism_non_syn-count","sars_polymorphism_syn-count") 
correct_signature <- c("iSNV\n(0.05-0.2)","iSNV\n(0.2-0.4)","iSNV\n(0.4-0.6)",
                       "iSNV\n(0.6-0.8)","iSNV\n(0.8-0.95)","polymorphism","polymorphism","polymorphism")
signature_dict <- setNames(correct_signature,original_signature)
long_rename_SBS192.input.catalog <- long_rename_SBS192.input.catalog %>% rowwise() %>%
                                    mutate(signature=signature_dict[[signature]])
mypal <- pal_npg("nrc")(10)
mypal <- c(mypal, c("#4DAF4A" , "#999999" , "#A65628"))
plot <- ggplot(data =long_rename_SBS192.input.catalog , aes(x = trinuc, y = exposure, fill = muttype)) +
        geom_bar(stat = "identity") + 
        facet_grid(signature ~ muttype, scales = "free_x", space = "free_x") +
        theme_bw() +
        ylab("Proportion of mutation") +
        labs(title="Mutation spectra of synonymous mutation")+
        scale_fill_manual(values = mypal) +
        scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
        my_theme()+
        theme(  
          panel.spacing.x = unit(0, "lines"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y =element_text(size = 7),
          axis.title.y = element_text(size = 9),
          strip.text = element_text(size = 5),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0, "null"))
ggsave("suppl_figure3B.pdf",width =20, height = 9, units = "cm",dpi = 300)

