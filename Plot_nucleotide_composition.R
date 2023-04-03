setwd("C:/Users/Teresa/OneDrive - University Of Cambridge/St Johns Cambridge/1-Part II Biochemistry/Part II Project/Files/Randomer sequencing")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr) 
library(RColorBrewer)


data<-read.csv("Nucleotides_triplets.csv", stringsAsFactors = FALSE)
data <- as.data.frame(data)
data <- subset(data, Pool == 'FN17L' | Pool =='FN73L' | Pool == 'F201N' | Pool == 'F205N')
data <- subset(data, Types == 'All triplets in sequence')

data$Nucleotide[data$Nucleotide == 'T'] <- 'U'

theme_set(theme_bw()+
            theme(axis.text.x = element_text(size=8), 
                  axis.text.y = element_text(size=8), 
                  axis.title.x = element_text(size=9), 
                  axis.title.y = element_text(size=9), 
                  text = element_text(family = "sans"),
                  legend.title = element_blank(),
                  legend.box.margin = margin (-8,-8,-8,-8),
                  strip.text = element_text(size = 9),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

p1 <- ggplot(data, aes(x = factor (Position, levels = c("1st", "2nd", "3rd")), 
                       y = Value ,
                       
                       fill = factor (Nucleotide, levels = c('A','U','G','C'))))+
  
  stat_summary(geom = "bar", 
               fun.min = mean,
               fun.max = mean,
               fun = mean,  
               size = 0,
               position = "stack") +
  
  geom_hline(yintercept=0.1, linetype = "dashed",
             color = "black", linewidth=0.02)+
  geom_hline(yintercept=0.2, linetype = "dashed",
             color = "black", linewidth=0.02)+
  geom_hline(yintercept=0.3, linetype = "dashed",
             color = "black", linewidth=0.02)+
  geom_hline(yintercept=0.4, linetype = "dashed",
             color = "black", linewidth=0.02)+
  geom_hline(yintercept=0.5, linetype = "dashed",
             color = "black", linewidth=0.02)+
  geom_hline(yintercept=0.6, linetype = "dashed",
             color = "black", linewidth=0.02)+
  geom_hline(yintercept=0.7, linetype = "dashed",
             color = "black", linewidth=0.02)+
  geom_hline(yintercept=0.8, linetype = "dashed",
             color = "black", linewidth=0.02)+
  geom_hline(yintercept=0.9, linetype = "dashed",
             color = "black", linewidth=0.02)+
  
  
  ylab("Percentage of Sequences") +
  xlab("Nucleotide Position in Triplet (5' --> 3')")+
  scale_fill_manual(values = c("skyblue3",  "indianred1", "burlywood2","grey70"))+
  
  facet_wrap(vars(Pool), 
             ncol =4
  )+ 
  scale_y_continuous(expand = c(0,0), labels = scales::percent, breaks = scales::pretty_breaks(n = 10) )+
  scale_x_discrete(expand = c(0,0))

print(p1)
ggsave(paste0("Nucleotides_triplets", ".png"), 
       p1, width = 12, height = 5.5, units = "cm", dpi = 400)