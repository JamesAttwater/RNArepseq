setwd("C:/Users/Teresa/OneDrive - University Of Cambridge/St Johns Cambridge/1-Part II Biochemistry/Part II Project/Files/Randomer sequencing")

#libraries --> if not installed need to do install.packages("libraryname") in console
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr) 
library(RColorBrewer)


data <- read.csv ("Length_allselec_unfiltered.csv")
data$Pool <- as.character(data$Pool)

#remove unwanted pools
data <- data[which(data$Pool == c("RZ5N","F201N","F205N","N17L","N73L","N73H", "FN17L","FN73L","FN73H")),]

theme_set(theme_classic()+
            theme(axis.text.x = element_text(size=7), 
                  axis.text.y = element_text(size=7), 
                  axis.title.x = element_text(size=9), 
                  axis.title.y = element_text(size=9), 
                  text = element_text(family = "sans"),
                  legend.title = element_blank(),
                  legend.box.margin = margin (-8,-8,-8,-8),
                  #legend.position = 'bottom',
                  
                  legend.key.size = unit(8, "point"),
                  strip.text = element_text(size = 8),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))
palette <- colorRampPalette(brewer.pal(8, "Paired"))(9)
pools <- c("F201N","F205N","FN17L","FN73L","FN73H", "RZ5N", "N17L","N73L","N73H")
p <- ggplot(data=data, 
            aes(Length, fill= factor(Pool, levels = pools), col = factor(Pool, levels = pools))) +
  geom_density(position = 'identity', alpha=.4)+
  scale_fill_manual(values = rev(palette)) +
  scale_color_manual(values = rev(palette))+
  
  labs( x= "Length of Sequence",
        y= "Percentage of Sequences")+ 
  scale_y_continuous(expand = c(0,0) , labels = scales::percent_format(accuracy = 5L)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,70,5))


print(p)

#save plot to file
ggsave(paste0("Length_distribution", ".png"), 
       p, width = 10.5, height = 4.5, units = "cm", dpi = 400)
