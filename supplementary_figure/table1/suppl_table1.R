library(tidyverse)
library(this.path)
library(ggsci)
setwd(this.dir())
source("../../bin/theme_setup.R")
###supplementary table1
###numbers of failures and successes for nonsynonymous iSNVs and synonymous iSNVs 
isnv_pair <- read_csv("../../data/isnv_pair_ruan_andreas.csv")
cutoff = 0.08
#high frequency
high_nonsyn <- isnv_pair %>% filter(donor_alt_frequency >= cutoff) %>%
                mutate(transmit=ifelse(recipient_alt_frequency >=0.03, "success", "failure")) %>%
                filter(func == "A") %>% group_by(donor , ref_pos_alt, donor_alt_frequency) %>% 
                summarise(success = sum(transmit == "success"), all = mean(recipient_count), 
                          mean_frequency = mean(donor_alt_frequency)*mean(recipient_count))  %>% 
                mutate(failure = sum(all-success))  %>% ungroup() %>%
                mutate(all_mean_frequency = sum(mean_frequency)/sum(all))
high_syn <- isnv_pair %>% filter(donor_alt_frequency >= cutoff) %>%
              mutate(transmit=ifelse(recipient_alt_frequency >=0.03, "success", "failure")) %>%
              filter(func == "S") %>% group_by(donor , ref_pos_alt, donor_alt_frequency) %>% 
              summarise(success = sum(transmit == "success"), all = mean(recipient_count), 
                        mean_frequency = mean(donor_alt_frequency)*mean(recipient_count))  %>% 
              mutate(failure = sum(all-success))  %>% ungroup() %>%
              mutate(all_mean_frequency = sum(mean_frequency)/sum(all))
high_freq_chi_table <- data.frame(success=c(sum(high_nonsyn$success) , sum(high_syn$success)), 
                                  failure=c(sum(high_nonsyn$failure) , sum(high_syn$failure)))
rownames(high_freq_chi_table) = c("nonsyn", "syn")
#low frequency 
low_nonsyn <- isnv_pair %>% filter(donor_alt_frequency < cutoff  & donor_alt_frequency >= 0.03 ) %>% 
              mutate(transmit=ifelse(recipient_alt_frequency >=0.03, "success", "failure")) %>%
              filter(func == "A") %>% group_by(donor , ref_pos_alt, donor_alt_frequency) %>% 
              summarise(success = sum(transmit == "success"), all = mean(recipient_count), 
                        mean_frequency = mean(donor_alt_frequency)*mean(recipient_count))  %>% 
              mutate(failure = sum(all-success))  %>% ungroup() %>%
              mutate(all_mean_frequency = sum(mean_frequency)/sum(all))
low_syn <- isnv_pair %>% filter(donor_alt_frequency < cutoff  & donor_alt_frequency >= 0.03 ) %>% 
            mutate(transmit=ifelse(recipient_alt_frequency >=0.03, "success", "failure")) %>%
            filter(func == "S") %>% group_by(donor , ref_pos_alt, donor_alt_frequency) %>% 
            summarise(success = sum(transmit == "success"), all = mean(recipient_count), 
                      mean_frequency = mean(donor_alt_frequency)*mean(recipient_count))  %>% 
            mutate(failure = sum(all-success))  %>% ungroup() %>%
            mutate(all_mean_frequency = sum(mean_frequency)/sum(all))
low_high_nonsyn_table <- data.frame(mutation_type = c("nonsyn", "nonsyn") ,
                                    observe_success=c(sum(low_nonsyn$success), sum(high_nonsyn$success)), 
                                    observe_failure=c(sum(low_nonsyn$failure), sum(high_nonsyn$failure)))
rownames(low_high_nonsyn_table) = c(paste0("low(0.03-",cutoff, ")"), paste0("high(>=", cutoff, ")"))
low_high_syn_table <- data.frame(mutation_type = c("syn", "syn") ,
                                 observe_success=c(sum(low_syn$success), sum(high_syn$success)), 
                                 observe_failure=c(sum(low_syn$failure), sum(high_syn$failure)))
rownames(low_high_syn_table) = c(paste0("low(0.03-",cutoff, ")"), paste0("high(>=", cutoff, ")"))
###merge
all_table <- rbind(low_high_nonsyn_table, low_high_syn_table)
all_table$frequency <- rownames(all_table)
all_table <- all_table %>% select(frequency, everything())
frequency_dict <- setNames(c("low(0.03-0.08)", "high(>=0.08)", "low(0.03-0.08)", "high(>=0.08)"), 
                           c("low(0.03-0.08)1", "high(>=0.08)1","low(0.03-0.08)", "high(>=0.08)"))
result_table <- all_table %>% mutate(success_prop = observe_success / (observe_success+observe_failure)) %>% 
                rowwise() %>% mutate(frequency = frequency_dict[[frequency]]) %>% 
                dplyr::rename( "iSNV frequency" = "frequency", "Muation type"= "mutation_type", 
                        "Success" = "observe_success", "Failure" = "observe_failure",
                        "Success proportion" =  "success_prop")
write_csv(result_table, "suppl_table1.csv")



