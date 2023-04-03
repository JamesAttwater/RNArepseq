setwd("C:/Users/Teresa/OneDrive - University Of Cambridge/St Johns Cambridge/1-Part II Biochemistry/Part II Project/Files/Randomer sequencing")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr) 
library(RColorBrewer)


data4 <- read.csv("Ribozyme_align.csv")


data <- subset(data4, pool== "RZ5N" | pool== "N17L"| pool== "N73L" )

#data$score = log10(data$score)
supp.labs <- c( "Type 5 complement (3' --> 5')", "Type 1 complement (3' --> 5')",
                "Type 5 (5' --> 3')", "Type 1 (5' --> 3')")
names(supp.labs) <- c("3_type5_reversecomp","4_type1_reversecomp","1_type5",  "2_type1" )

labels_seq <- data$position_seq

pool <- c("RZ5N", "F201N", "F205N", "N17L", "N73L", "N73H", "FN17L","FN73L","FN73H")
intensity <- c(10892154, 10000000, 39507018.98, 62852360, 246946355.5, 95373574, 113600458.1, 174104558.3, 48852872) #need to correct F201N

(df <- data.frame(pool, intensity))


data_merge <- merge(data, df, by = "pool")

data_merge$normalized <- data_merge$score*data_merge$intensity/10**5


theme_set(theme_bw()+
            theme(axis.text.x = element_text(size=6), 
                  axis.text.y = element_text(size=6), 
                  axis.title.x = element_text(size=9), 
                  axis.title.y = element_text(size=9), 
                  legend.text = element_text(size=8),
                  legend.title = element_blank(),
                  legend.box.margin = margin (-8,-8,-8,-8),
                  text = element_text(family = "sans"),
                  strip.text = element_text(size = 8),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))
p <- ggplot(data=data_merge, 
            aes(x=position, 
                y=normalized, 
                group= factor(pool, levels=pools),
                color = factor(pool, levels=pools))) +
  geom_line(size=0.8)+
  labs( x= "Position in Sequence",
        y= "Normalized number of sequences")+
  facet_wrap(vars(ribozyme), 
             ncol =2, 
             labeller = labeller(ribozyme=supp.labs))+#,
  #scales = "free")+ 
  scale_y_continuous(expand = c(0,0), labels = scales::comma )+  
  scale_x_continuous(breaks = seq(0,152,10))+
  scale_color_manual("pool", values = c("cyan3",  "darkolivegreen4", "coral2","darkgoldenrod2"))


print(p)


#save plot to file
ggsave(paste0("Ribozyme_align_plot", ".png"), 
       p, width = 16.8, height = 8, units = "cm", dpi = 400)