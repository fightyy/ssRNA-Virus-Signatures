library(tidyverse)
library(this.path)
library(ggsci)
library(segmented)
setwd(this.dir())
source("../../bin/theme_setup.R")

###suppl_figure11B
#The distance in synonymous mutational spectra between the earliest time and later times for SARS-CoV-2
dnds_df <- read_csv("../../data/sars_dnds.csv")
#synonymous spectrum distance regression
lm_df <- read_csv("../data/sars_spectrum_distance.csv")
plot_lm_df <- lm_df %>% filter(mutation_type == "syn_distance") %>% filter(average_time!=0)
dnds_min <- 0.55
dnds_max <- 0.75
distance_range <- max(plot_lm_df$distance)-min(plot_lm_df$distance)
dnds_range <- dnds_max - dnds_min
scale_factor <- distance_range / dnds_range
mypal <- pal_npg("nrc")(10)
lm_fit <- lm(distance ~ average_time, data = plot_lm_df)
#piecewise regression
segment_fit <- segmented(lm_fit, npsi=1)
model <- segment_fit  
breakpoints <- model$psi[, "Est."]
slopes <- slope(model)$average_time[, "Est."]
segments <- c(0, breakpoints, 47) 
mid_points <- (segments[-1] + segments[-length(segments)]) / 2
slope_labels <- data.frame(
                x = mid_points,   
                y = max(plot_lm_df$distance) * 1.2,
                label = round(slopes, 4))

distance_plot <- ggplot() +
                geom_point(data = plot_lm_df, aes(x = average_time, y = distance), color=mypal[2], size=1) +
                geom_line(data = plot_lm_df, aes(x = average_time, y = fitted(segment_fit)),
                          color = mypal[2], size = 0.5) +
                geom_text(data = slope_labels,
                          aes(x = x, y = y, label = label),
                          color = "red", size = 2, fontface = "bold") +
                geom_vline(xintercept = segment_fit[["psi"]][, "Est."], linetype = "dashed", color = "black") +
                geom_point(
                  data = dnds_df,
                  aes(x = start, y = (dnds - dnds_min) * scale_factor + min(plot_lm_df$distance)),
                  color = "black", size = 0.5) +
                geom_line(data = dnds_df,
                          aes(x = start, y = (dnds - dnds_min) * scale_factor + min(plot_lm_df$distance)),
                          color = "black", size = 0.1)+
                scale_y_continuous(
                  name = "Pairwise distance between\nsynonymous spectra",
                  sec.axis = sec_axis(
                    ~ (. - min(plot_lm_df$distance))/ scale_factor + dnds_min,
                    name = "dN/dS")) +
                labs(x = "Month (since 2019/12)",  title="SARS-CoV-2") +
                my_theme()+
                theme(
                  legend.position = "bottom",
                  legend.title = element_blank(),
                  axis.title.y.right = element_text(angle = 90,face="italic"),
                  legend.margin = margin(t = -0.3 , b= 0.1, unit = "cm")) 

ggsave("suppl_figure11B.pdf",unit="cm", width=8, height=6,dpi = 300)

###suppl_figure11A
#The distance in nonsynonymous mutational spectra between the earliest time and later times for SARS-CoV-2
plot_lm_df <- lm_df %>% filter(mutation_type == "nonsyn_distance") %>% filter(average_time!=0)
dnds_min <- 0.55
dnds_max <- 0.75
distance_range <- max(plot_lm_df$distance)-min(plot_lm_df$distance)
dnds_range <- dnds_max - dnds_min
scale_factor <- distance_range / dnds_range
mypal <- pal_npg("nrc")(10)
lm_fit <- lm(distance ~ average_time, data = plot_lm_df)
#piecewise regression
segment_fit <- segmented(lm_fit, npsi=1)
model <- segment_fit 
breakpoints <- model$psi[, "Est."]
slopes <- slope(model)$average_time[, "Est."]
segments <- c(0, breakpoints, 47)
mid_points <- (segments[-1] + segments[-length(segments)]) / 2 
slope_labels <- data.frame(
                x = mid_points, 
                y = max(plot_lm_df$distance) * 1.05+seq(0,0.025,length.out= length(mid_points)),
                label = round(slopes, 4))
distance_plot <- ggplot() +
                geom_point(data = plot_lm_df, aes(x = average_time, y = distance), color=mypal[1], size=1) +
                geom_line(data = plot_lm_df, aes(x = average_time, y = fitted(segment_fit)),
                          color = mypal[1], size = 0.5) +
                geom_text(data = slope_labels,
                          aes(x = x, y = y, label = label),
                          color = "red", size = 2, fontface = "bold") +
                geom_vline(xintercept = segment_fit[["psi"]][, "Est."], linetype = "dashed", color = "black") +
                geom_point(
                  data = dnds_df,
                  aes(x = start, y = (dnds - dnds_min) * scale_factor + min(plot_lm_df$distance)),
                  color = "black", size = 0.5) +
                geom_line(data = dnds_df,
                          aes(x = start, y = (dnds - dnds_min) * scale_factor + min(plot_lm_df$distance)),
                          color = "black", size = 0.1)+
                scale_y_continuous(
                  name = "Pairwise distance between\nnonsynonymous spectra",
                  sec.axis = sec_axis(
                    ~ (. - min(plot_lm_df$distance))/ scale_factor + dnds_min,
                    name = "dN/dS")) +
                labs(x = "Month (since 2019/12)",  title="SARS-CoV-2") +
                my_theme()+
                theme(
                  legend.position = "bottom",
                  legend.title = element_blank(),
                  axis.title.y.right = element_text(angle = 90,face="italic"),
                  legend.margin = margin(t = -0.3 , b= 0.1, unit = "cm")) 

ggsave("suppl_figure11A.pdf",unit="cm", width=8, height=6,dpi = 300)

###suppl_fig11D
#The distance in synonymous mutational spectra between the earliest time and later times for H3N2
dnds_df <- read_csv("../../data/h3n2_dnds.csv")
lm_df <- read_csv("../data/h3n2_spectrum_distance.csv")
plot_lm_df <- lm_df %>% filter(mutation_type == "syn_distance") %>% filter(average_time!=0)
dnds_min <- 0.35
dnds_max <- 0.55
distance_range <- max(plot_lm_df$distance)-min(plot_lm_df$distance)
dnds_range <- dnds_max - dnds_min
scale_factor <- distance_range / dnds_range
mypal <- pal_npg("nrc")(10)
lm_fit <- lm(distance ~ average_time, data = plot_lm_df)
#piecewise regression
segment_fit <- segmented(lm_fit, npsi=1)
model <- segment_fit 
breakpoints <- model$psi[, "Est."]
slopes <- slope(model)$average_time[, "Est."]
segments <- c(0, breakpoints, 24) 
mid_points <- (segments[-1] + segments[-length(segments)]) / 2  
slope_labels <- data.frame(
                  x = mid_points,   
                  y = max(plot_lm_df$distance) * 1.2,
                  label = round(slopes, 4))
distance_plot <- ggplot() +
                geom_point(data = plot_lm_df, aes(x = average_time, y = distance), color=mypal[2], size=1) +
                geom_line(data = plot_lm_df, aes(x = average_time, y = fitted(segment_fit)),
                          color = mypal[2], size = 0.5) +
                geom_text(data = slope_labels,
                          aes(x = x, y = y, label = label),
                          color = "red", size = 2, fontface = "bold") +
                geom_vline(xintercept = segment_fit[["psi"]][, "Est."], linetype = "dashed", color = "black") +
                geom_point(
                  data = dnds_df,
                  aes(x = start, y = (dnds - dnds_min) * scale_factor + min(plot_lm_df$distance)),
                  color = "black", size = 0.5
                ) +
                geom_line(data = dnds_df,
                          aes(x = start, y = (dnds - dnds_min) * scale_factor + min(plot_lm_df$distance)),
                          color = "black", size = 0.1)+
                scale_y_continuous(
                  name = "Pairwise distance between\nsynonymous spectra",
                  sec.axis = sec_axis(
                    ~ (. - min(plot_lm_df$distance))/ scale_factor + dnds_min,
                    name = "dN/dS"
                  )
                ) +
                labs(x = "Year (since 2000)",  title="H3N2") +
                my_theme()+
                theme(
                  legend.position = "bottom",
                  legend.title = element_blank(),
                  axis.title.y.right = element_text(angle = 90,face="italic"),
                  legend.margin = margin(t = -0.3 , b= 0.1, unit = "cm")) 

ggsave("suppl_figure11D.pdf",unit="cm", width=8, height=6,dpi = 300)

###suppl_fig11C
##The distance in nonsynonymous mutational spectra between the earliest time and later times for H3N2
plot_lm_df <- lm_df %>% filter(mutation_type == "nonsyn_distance") %>% filter(average_time!=0)
dnds_min <- 0.35
dnds_max <- 0.55
distance_range <- max(plot_lm_df$distance)-min(plot_lm_df$distance)
dnds_range <- dnds_max - dnds_min
scale_factor <- distance_range / dnds_range
mypal <- pal_npg("nrc")(10)
lm_fit <- lm(distance ~ average_time, data = plot_lm_df)
segment_fit <- segmented(lm_fit, npsi=1)
model <- segment_fit
breakpoints <- model$psi[, "Est."]
slopes <- slope(model)$average_time[, "Est."]
segments <- c(0, breakpoints, 24) 
mid_points <- (segments[-1] + segments[-length(segments)]) / 2
slope_labels <- data.frame(
                x = mid_points,      
                y = max(plot_lm_df$distance) * 1.05+seq(0,0.025,length.out= length(mid_points)),
                label = round(slopes, 4))
distance_plot <- ggplot() +
                geom_point(data = plot_lm_df, aes(x = average_time, y = distance), color=mypal[1], size=1) +
                geom_line(data = plot_lm_df, aes(x = average_time, y = fitted(segment_fit)),
                          color = mypal[1], size = 0.5) +
                geom_text(data = slope_labels,
                          aes(x = x, y = y, label = label),
                          color = "red", size = 2, fontface = "bold") +
                geom_vline(xintercept = segment_fit[["psi"]][, "Est."], linetype = "dashed", color = "black") +
                geom_point(
                  data = dnds_df,
                  aes(x = start, y = (dnds - dnds_min) * scale_factor + min(plot_lm_df$distance)),
                  color = "black", size = 0.5) +
                geom_line(data = dnds_df,
                          aes(x = start, y = (dnds - dnds_min) * scale_factor + min(plot_lm_df$distance)),
                          color = "black", size = 0.1)+
                scale_y_continuous(
                  name = "Pairwise distance between\nnonsynonymous spectra",
                  sec.axis = sec_axis(
                    ~ (. - min(plot_lm_df$distance))/ scale_factor + dnds_min,
                    name = "dN/dS")) +
                labs(x = "Year (since 2000)",  title="H3N2") +
                my_theme()+
                theme(
                  legend.position = "bottom",
                  legend.title = element_blank(),
                  axis.title.y.right = element_text(angle = 90,face="italic"),
                  legend.margin = margin(t = -0.3 , b= 0.1, unit = "cm")) 
ggsave("suppl_figure11C.pdf",unit="cm", width=8, height=6,dpi = 300)



