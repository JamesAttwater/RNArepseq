setwd("C:/Users/Teresa/OneDrive - University Of Cambridge/St Johns Cambridge/1-Part II Biochemistry/Part II Project/Files/Randomer sequencing")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr) 
library(RColorBrewer)

### Sequence pool explanation:
#F201N:A, F205N: B, FN17L: C, FN73L: D, FN73H: E, RZ5N: F, N17L: G, N73L: H, N73H: I


pool <- c("RZ5N", "F201N", "F205N", "N17L", "N73L", "N73H", "FN17L","FN73L","FN73H")
intensity <- c(10892154, 10000000, 39507018.98, 62852360, 246946355.5, 95373574, 113600458.1, 174104558.3, 48852872) #need to correct F201N
filtered <- c(0.455,0.362,0.499,0.649,0.689,0.600, 0.673,0.674,0.538)

intense_filt = intensity*filtered

#(df <- data.frame(pool, intense_filt))
(df <- data.frame(pool, intensity))
options(digits = 18)

data <- read.csv("Alignment_count_all.csv")
data4 <- as.data.frame(data)
data <- subset(data4, type== "grey_area" | type== "not_ribozyme"| type== "plus_aligned" | type== "minus_aligned"  | type== "plus_and_minus")
data_merge <- merge(data, df, by = "pool")

pools <- c("F201N","F205N","FN17L","FN73L","FN73H", "RZ5N", "N17L","N73L","N73H")

types <- c("plus_aligned",  "minus_aligned","plus_and_minus", "grey_area", "not_ribozyme")

data_merge$sum_pool = with(data_merge, ave(value, pool, FUN = function(i)i / sum(i)))

data_merge$normalized <- data_merge$sum_pool*data_merge$intensity/10**6


theme_set(theme_classic()+
            theme(axis.text.x = element_text(size=7), 
                  axis.text.y = element_blank(), 
                  axis.ticks.y=element_blank(),
                  axis.title.x = element_text(size=9), 
                  axis.title.y = element_text(size=9), 
                  legend.text = element_text(size=8),
                  legend.title = element_blank(),
                  #  legend.position='bottom',
                  legend.key.size = unit(8, "point"),
                  legend.box.margin = margin (-10,-10,-10,-10),
                  text = element_text(family = "sans"),
                  strip.text = element_text(size = 8),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))


###Plot adding up to 100%
p <- ggplot(data=data_merge, 
            aes(x=factor(pool, levels = pools), 
                y=normalized,
                fill = factor(type, levels = types), 
                group = factor(type, levels = types))) +
  
  geom_bar( position = 'fill',stat="identity")+
  labs( x= "Sequence Pool",
        y= "Percentage of Sequences")+ 
  scale_fill_manual("type", values = c("steelblue4", "skyblue3", "skyblue", "grey", "indianred1"), labels = c("(+)-strand only aligned", "(-)-strand only aligned","(+)- and (-)-strand aligned", "grey area","not ribozyme")) +
  
  scale_y_continuous(expand = c(0,0), labels = scales::comma)  #can insert limits = c(lower, upper)



print(p)

#save plot to file
ggsave(paste0("Alignment_count_percent", ".png"), 
       p, width = 13.7, height = 3.6, units = "cm", dpi = 400)

### Plot stacking bars

p <- ggplot(data=data_merge, 
            aes(x=factor(pool, levels = pools), 
                y=normalized,
                fill = factor(type, levels = types), 
                group = factor(type, levels = types))) +
  
  geom_bar( position = 'stack',stat="identity")+
  labs( x= "Sequence Pool",
        y= "Normalized number\nof sequences")+ 
  scale_fill_manual("type", values = c("steelblue4", "skyblue3", "skyblue", "grey", "indianred1"), labels = c("(+)-strand only aligned", "(-)-strand only aligned","(+)- and (-)-strand aligned", "grey area","not ribozyme")) +
  
  scale_y_continuous(expand = c(0,0), labels = scales::comma)  



print(p)

#save plot to file
ggsave(paste0("Alignment_count_stack", ".png"), 
       p, width = 13.7, height = 3.6, units = "cm", dpi = 400)



### 9-nt comparison in longer sequence pools

pool <- c("RZ5N", "F201N", "F205N", "N17L", "N73L", "N73H", "FN17L","FN73L","FN73H")

intensity <- c(10892154, 10000000, 39507018.98, 62852360, 246946355.5, 95373574, 113600458.1, 174104558.3, 48852872) #need to correct F201N

(df <- data.frame(pool, intensity))
options(digits = 18)

data<- read.csv("Alignment_count_longer_9nt.csv")
data4 <- as.data.frame(data)
data <- subset(data4, type== "grey_area" | type== "not_ribozyme"| type== "plus_aligned" | type== "minus_aligned"  | type== "plus_and_minus")
data_merge <- merge(data, df, by = "pool")

pools <- c("F201N","F205N","FN17L","FN73L","FN73H", "RZ5N", "N17L","N73L","N73H")

types <- c("plus_aligned",  "minus_aligned","plus_and_minus", "grey_area", "not_ribozyme")

data_merge$sum_pool = with(data_merge, ave(value, pool, FUN = function(i)i / sum(i)))

data_merge$normalized <- data_merge$sum_pool*data_merge$intensity/10**6

theme_set(theme_classic()+
            theme(axis.text.x = element_text(size=7), 
                  axis.text.y = element_text(size=7), 
                  axis.title.x = element_text(size=9), 
                  axis.title.y = element_text(size=9), 
                  legend.text = element_text(size=8),
                  legend.title = element_blank(),
                  #  legend.position='bottom',
                  legend.key.size = unit(8, "point"),
                  legend.box.margin = margin (-10,-10,-10,-10),
                  text = element_text(family = "sans"),
                  strip.text = element_text(size = 8),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

p <- ggplot(data=data_merge, 
            aes(x=factor(pool, levels = pools), 
                y=normalized,
                fill = factor(type, levels = types), 
                group = factor(type, levels = types))) +
  
  geom_bar( position = 'fill',stat="identity")+
  labs( x= "Sequence Pool",
        y= "Percentage of Sequences")+ #darkslategray1
  scale_fill_manual("type", values = c("steelblue4", "skyblue3", "skyblue", "grey", "indianred1"), labels = c("(+)-strand only aligned", "(-)-strand only aligned","(+)- and (-)-strand aligned", "grey area","not ribozyme")) +
  
  scale_y_continuous(expand = c(0,0), labels = scales::comma) 



print(p)

#save plot to file
ggsave(paste0("Alignment_count_9nt_percent", ".png"), 
       p, width = 9.3, height = 3.6, units = "cm", dpi = 400)


p <- ggplot(data=data_merge, 
            aes(x=factor(pool, levels = pools), 
                y=normalized,
                fill = factor(type, levels = types), 
                group = factor(type, levels = types))) +
  
  geom_bar( position = 'stack',stat="identity")+
  labs( x= "Sequence Pool",
        y= "Normalised Number of Sequences")+
  scale_fill_manual("type", values = c("steelblue4", "skyblue3", "skyblue", "grey", "indianred1"), labels = c("(+)-strand only aligned", "(-)-strand only aligned","(+)- and (-)-strand aligned", "grey area","not ribozyme")) +
  
  scale_y_continuous(expand = c(0,0), labels = scales::comma) 



print(p)

#save plot to file
ggsave(paste0("Alignment_count_9nt_stack", ".png"), 
       p, width = 9.3, height = 3.6, units = "cm", dpi = 400)


###Random Sequences alignment


pool <- c("RZ5N", "F201N", "F205N", "N17L", "N73L", "N73H", "FN17L","FN73L","FN73H")

options(digits = 18)

data<- read.csv("Alignment_count_9nt_random.csv")
data4 <- as.data.frame(data)
data <- subset(data4, type== "grey_area" | type== "not_ribozyme"| type== "plus_aligned" | type== "minus_aligned"  | type== "plus_and_minus")

pools <- c("F201N","F205N","FN17L","FN73L","FN73H", "RZ5N", "N17L","N73L","N73H")

types <- c("plus_aligned",  "minus_aligned","plus_and_minus", "grey_area", "not_ribozyme")


theme_set(theme_classic()+
            theme(axis.text.x = element_text(size=7), 
                  axis.text.y = element_text(size=7), 
                  axis.title.x = element_text(size=9), 
                  axis.title.y = element_text(size=9), 
                  legend.text = element_text(size=8),
                  legend.title = element_blank(),
                  #  legend.position='bottom',
                  legend.key.size = unit(8, "point"),
                  legend.box.margin = margin (-10,-10,-10,-10),
                  text = element_text(family = "sans"),
                  strip.text = element_text(size = 8),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

p <- ggplot(data=data, 
            aes(x=factor(pool, levels = pools), 
                y=value,
                fill = factor(type, levels = types), 
                group = factor(type, levels = types))) +
  
  geom_bar( position = 'fill',stat="identity")+
  labs( x= "Sequence Pool",
        y= "Percentage of Sequences")+ #darkslategray1
  scale_fill_manual("type", values = c("steelblue4", "skyblue3", "skyblue", "grey", "indianred1"), labels = c("(+)-strand only aligned", "(-)-strand only aligned","(+)- and (-)-strand aligned", "grey area","not ribozyme")) +
  
  scale_y_continuous(expand = c(0,0), labels = scales::percent) #can insert limits = c(lower, upper)



print(p)

#save plot to file
ggsave(paste0("Alignment_count_random_9nt", ".png"), 
       p, width = 9.3, height = 3.6, units = "cm", dpi = 400)


