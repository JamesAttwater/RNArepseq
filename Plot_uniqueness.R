setwd("C:/Users/Teresa/OneDrive - University Of Cambridge/St Johns Cambridge/1-Part II Biochemistry/Part II Project/Files/Randomer sequencing")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr) 
library(RColorBrewer)

data <- read.csv("Uniqueness_counts_splits.csv")

data$split <- as.character(data$split)
data$split[data$split == "1-ribo"] = "Ribozyme aligned"

data$split[data$split == "2-grey"] = "Grey area"

data$split[data$split == "3-non_ribo"] = "Not ribozyme aligned"

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

p = ggplot (data,
            aes(x= factor(pool, levels = c("RZ5N","F201N","F205N","N17L","N73L","N73H","FN17L","FN73L","FN73H")),
                y = value, fill = factor(split, levels = c("Ribozyme aligned","Grey area", "Not ribozyme aligned"))))+
  
  geom_bar(position="dodge", stat="identity")+
  
  
  #axis labels
  labs(x = "Sequence Pool", y = "Percentage of Sequences")+
  
  scale_y_continuous(expand = c(0,0),labels = scales::percent_format(accuracy = 5L), limits = c(0,1) ) + #expand removes white space above and below, without percentage get 0 to 1.0
  
  #adjust colours (it only takes as many as you have categories but you always need to supply enough) - these are just some that I had used. you can find more here: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
  scale_fill_manual("cell_type", values = c(  "skyblue3", "grey","indianred1", "darkslategray1", "aquamarine2")) 

print (p)

#save plot to file --> can change resolution with dpi, file type, dimensions, and file name
ggsave(paste0("Uniqueness_split_barplot", ".png"), 
       p, width = 16, height = 7, units = "cm", dpi = 400)