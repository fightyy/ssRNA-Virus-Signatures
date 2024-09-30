library(tidyverse)
library(ggsci)
library(this.path)
setwd(this.dir())
source("../bin/theme_setup.R")

###Figure3A
#signatures of SARS-CoV-2
SBS192.input.catalog <- read_tsv("../data/CH192_De-Novo_Signatures_sars.txt") 
SBS192.input.catalog$muttype <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["MutationsType"],1,1),
                ">",
                substring(x["MutationsType"],2,2)))
})
SBS192.input.catalog$trinuc <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["MutationsType"],3,3),
                substring(x["MutationsType"],1,1),
                substring(x["MutationsType"],4,4)))
})
rename_SBS192.input.catalog <- SBS192.input.catalog %>% 
                                dplyr::rename( "SARS-COV-2\nsignatureA"="SBS192A",
                                                 "SARS-COV-2\nsignatureB"="SBS192B")
long_rename_SBS192.input.catalog <- rename_SBS192.input.catalog %>% 
                                    pivot_longer( cols = c("SARS-COV-2\nsignatureA","SARS-COV-2\nsignatureB"),
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
        ylab("Probability of mutation") +
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
          strip.text = element_text(size = 9),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0, "null"))
ggsave("Figure3A.pdf",width =16, height = 6, units = "cm",dpi = 300)

###figure3C
#signatures of H3N2
SBS192.input.catalog <- read_tsv("../data/CH192_De-Novo_Signatures_h3n2.txt") 
SBS192.input.catalog$muttype <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["MutationsType"],1,1),
                ">",
                substring(x["MutationsType"],2,2)))
})

SBS192.input.catalog$trinuc <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["MutationsType"],3,3),
                substring(x["MutationsType"],1,1),
                substring(x["MutationsType"],4,4)))
})

rename_SBS192.input.catalog <- SBS192.input.catalog %>% dplyr::rename("H3N2\nsignatureA"="SBS192A",
                                                               "H3N2\nsignatureB"="SBS192B",
                                                               "H3N2\nsignatureC"="SBS192C"
)

long_rename_SBS192.input.catalog <- rename_SBS192.input.catalog %>% 
                                    pivot_longer( cols = c("H3N2\nsignatureA", "H3N2\nsignatureB","H3N2\nsignatureC"),
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
        ylab("Probability of mutation") +
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
          strip.text = element_text(size = 9),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0, "null"))
ggsave("Figure3C.pdf",width =16, height = 9, units = "cm",dpi = 300)

###figure3B
#number of mutations explained by signatures in SARS-CoV-2
meta_df <- read_tsv("../data/all_signature_sample_csv.txt") %>% rowwise() %>% mutate() %>% select(c(sample, year, host, species, order))
signature_df <- read_tsv("../data/CH192_De-Novo_Activities_refit_sars.txt") %>% rename_all(~gsub("-", "_", .))
data_df <- signature_df %>% left_join(meta_df, by=c("Samples"="sample")) %>% filter(host != "internal") %>% filter(! order %in% c("ssRNA_+_", "ssRNA_-_"))
time_group_df <- suppressWarnings(data_df %>% rowwise() %>% 
                                    mutate(start_year = as.numeric(unlist(strsplit(year, "-"))[[1]]),
                                           end_year = as.numeric(unlist(strsplit(year, "-"))[[2]])) %>%
                                    mutate(mean_between_years = floor(mean(c(start_year, end_year)))))
signatrue_vector <- colnames(time_group_df)[grep("192", colnames(time_group_df))]
group_time_group_df <- time_group_df %>% select(signatrue_vector, "mean_between_years") %>% group_by(mean_between_years) %>%
                       summarise(across(everything(), ~mean(.x, na.rm = TRUE))) %>% na.omit()
group_data_df <- group_time_group_df %>% pivot_longer(cols= signatrue_vector ,
                                                      values_to =  "count" , 
                                                      names_to = "signature")
signatrue_color <- c("#E64B35FF","#3C5488FF","#FFA500")
plot <- ggplot(group_data_df, aes(x = mean_between_years, y=count , fill = signature)) + 
        geom_col(position = "stack") +
        scale_fill_manual(values=signatrue_color,
                          labels=c("SARS-COV-2\nsignatureA", "SARS-COV-2\nsignatureB","SARS-COV-2\nsignatureC"))+
        labs(x = "Month", y = "Count")+
        labs(fill="Mutation signature")+
        my_theme()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.key.height=unit(1.2,"cm"))
ggsave("Figure3B.pdf",unit="cm", width=8, height=6,dpi = 300)




###figure3D
#number of mutations explained by signatures in H3N2
meta_df <- read_tsv("../data/all_signature_sample_csv.txt") %>% rowwise() %>% mutate() %>% select(c(sample, year, host, species, order))
signature_df <- read_tsv("../data/CH192_De-Novo_Activities_refit_h3n2.txt") %>% rename_all(~gsub("-", "_", .))
data_df <- signature_df %>% left_join(meta_df, by=c("Samples"="sample")) %>% filter(host != "internal") %>% filter(! order %in% c("ssRNA_+_", "ssRNA_-_"))
time_group_df <- suppressWarnings(data_df %>% rowwise() %>% 
                                    mutate(start_year = as.numeric(unlist(strsplit(year, "-"))[[1]]),
                                    end_year = as.numeric(unlist(strsplit(year, "-"))[[2]])) %>%
                                    mutate(mean_between_years = floor(mean(c(start_year, end_year)))))
signatrue_vector <- colnames(data_df)[grep("192", colnames(data_df))]
group_time_group_df <- time_group_df %>% select(signatrue_vector, "mean_between_years") %>% group_by(mean_between_years) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE))) %>% filter(mean_between_years >2000)
group_data_df <- group_time_group_df %>% pivot_longer(cols= signatrue_vector ,
                                                      values_to =  "count" , 
                                                      names_to = "signature")
signatrue_color <- c("#E64B35FF","#3C5488FF","#FFA500")
plot <- ggplot(group_data_df, aes(x = mean_between_years, y=count , fill = signature)) + 
  geom_col(position = "stack") +
  scale_fill_manual(values=signatrue_color,
                    labels=c("H3N2\nsignatureA", "H3N2\nsignatureB","H3N2\nsignatureC"))+
  labs(x = "Year",y = "Count")+
  labs(fill="Mutation signature")+
  my_theme()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height=unit(1.2,"cm"))
ggsave("Figure3D.pdf",unit="cm", width=8, height=6,dpi = 300)


