library(tidyverse)
library(this.path)
library(ggsci)
setwd(this.dir())
source("../../bin/theme_setup.R")
###supplementary figure8
###signatures extracted by mSigHdp
SBS192.input.catalog <- data.table::fread("signatures_mSigHdp.txt") 
matrix <- SBS192.input.catalog %>% select(-1)
normalized_matrix <- apply(matrix,1, function(x) x/colSums(matrix)) %>% t() %>% as.data.frame()
normalized_matrix$MutationType <-  SBS192.input.catalog$MutationType
SBS192.input.catalog <- normalized_matrix
SBS192.input.catalog$muttype <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["MutationType"],3,3),
                ">",
                substring(x["MutationType"],5,5)))
})

SBS192.input.catalog$trinuc <- apply(SBS192.input.catalog,1,function(x){
  return(paste0(substring(x["MutationType"],1,1),
                substring(x["MutationType"],3,3),
                substring(x["MutationType"],7,7)))
})
long_rename_SBS192.input.catalog <- SBS192.input.catalog %>% 
                                    pivot_longer( cols = contains("SBS"),
                                                  names_to = "signature",
                                                  values_to = "exposure") %>% rowwise() %>% 
                                    mutate(signature=paste0("SBS\n96",substr(signature[[1]], 6, 6)))

long_rename_SBS192.input.catalog$muttype <- factor(long_rename_SBS192.input.catalog$muttype,
                                                   levels=c("C>A","C>G", "C>T",
                                                            "T>A","T>C","T>G"))
mypal <- pal_npg("nrc")(10)
mypal <- c("#4DBBD5FF","black","#E64B35FF","#999999","#97C758", "#E8BFBC")
plot <- ggplot(data =long_rename_SBS192.input.catalog , aes(x = trinuc, y = exposure, fill = muttype)) +
  geom_bar(stat = "identity") + 
  facet_grid(signature ~ muttype, scales = "free") +
  theme_bw() +
  ylab("Probability of mutation") +
  scale_fill_manual(values = mypal) +
  scale_y_continuous(expand = c(0, 0, 0.05, 0))+
  my_theme()+
  theme(  
    panel.spacing.x = unit(0, "lines"),
    panel.spacing.y = unit(0.1, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x =element_text(size = 5),
    axis.text.y =element_text(size = 6),
    axis.title.y = element_text(size = 8),
    strip.text.y = element_text(size = 7),
    strip.text.x = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0, "null"))+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("suppl_figure8.pdf",width =16, height = 18, units = "cm",dpi = 300)

