setwd("C:/Users/Teresa/OneDrive - University Of Cambridge/St Johns Cambridge/1-Part II Biochemistry/Part II Project/Files/Randomer sequencing")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr) 
library(RColorBrewer)

data <- read.csv("Complementarity_counts_splits.csv")


data <- subset(data, category == 'fraction_unique_complementary_seq' | category =='fraction_complementary_seq_with_repeat')

data$category <- as.character(data$category)

data$category[data$category == "fraction_unique_complementary_seq"] <- "unique"
data$category[data$category == "fraction_complementary_seq_with_repeat"] <- "with repeat"

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

#labels for in the plot etc
cat = c("unique","with repeats")
facet_labels <- c( "Ribozyme aligned","Grey area","Not ribozyme aligned")
names(facet_labels) <- c("1-ribo",  "2-grey", "3-non_ribo" ) #this assigns the names respectively to those categories in the data

p = ggplot (data,
            aes(x= factor(pool, levels = c("RZ5N","F201N","F205N","N17L","N73L","N73H","FN17L","FN73L","FN73H")),
                y = value,
                fill = category))+
  
  geom_bar(position="dodge", stat="identity")+
  
  #facet_grid is responsible for those three parts to the plot
  facet_wrap(vars(split), 
             ncol =1,
             labeller = labeller(split=facet_labels))+ 
  
  #axis labels
  labs(x = "Sequence Pool", y = "Percentage of Sequences")+
  
  scale_y_continuous(expand = c(0,0),labels = scales::percent_format(accuracy = 5L), limits = c(0,0.5) ) + 
  scale_fill_manual("cell_type", values = c(  "skyblue3", "indianred1","grey", "darkslategray1", "aquamarine2")) 

print (p)

#save plot to file --> can change resolution with dpi, file type, dimensions, and file name
ggsave(paste0("Complementarity_split_barplot", ".png"), 
       p, width = 13, height = 8, units = "cm", dpi = 400)