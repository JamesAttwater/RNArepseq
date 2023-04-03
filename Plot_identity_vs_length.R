setwd("C:/Users/Teresa/OneDrive - University Of Cambridge/St Johns Cambridge/1-Part II Biochemistry/Part II Project/Files/Randomer sequencing")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr) 
library(RColorBrewer)

###HI sequences exluding plus and plusminus aligned sequences

df = read.csv("Length_vs_alignment_score_H_I.csv")

df2 =read.csv("Intensities_HI.csv") #normalising with gel intensities from different lengths

sample = "H_I"

df3 = df[!(df$classification=='plus_aligned'|df$classification=='plus_and_minus'),]

library(plyr)
counts <- ddply(df3, .(df3$len_seq, df3$score, df3$class), nrow)
names(counts) <- c("len_seq", "score", "class", "Freq")

counts$len_seq <- as.factor(counts$len_seq)
merge = merge(df2,counts, by ='len_seq')
merge$Freq_intensity = merge$Freq*merge$Intensity_factor


classes = c("aligned", "grey_area", "not_ribozyme")

theme_set(theme_classic()+
            theme(axis.text.x = element_text(size=20), 
                  axis.text.y = element_text(size=20), 
                  axis.title.x = element_text(size=20), 
                  axis.title.y = element_text(size=20), 
                  text = element_text(family = "sans"),
                  legend.title = element_blank(),
                  legend.text = element_text(size =20),
                  legend.key.size = unit(40, "point"),
                  legend.box.margin = margin (-10,-10,-10,-10),
                  strip.text = element_text(size = 5),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

p=ggplot(merge,
         aes(x = len_seq,
             y = score,
             fill = factor(class, levels = classes), 
             color = factor(class, levels = classes)
         ))+
  geom_point(aes(size = log10(Freq_intensity)))+
  scale_size_continuous(range=c(0.01,5))+
  labs( x= "Sequence Length (nt)",
        y= "Matching Nucleotides with TPR")+
  ylim(0,66)+
  scale_x_continuous(breaks = seq(6,66,3), limits = c(6,66))+
  
  scale_fill_manual("class", values = c( "skyblue3",  "grey", "indianred1"), labels = c("aligned", "grey area","not ribozyme"))+
  scale_color_manual("class", values = c( "skyblue3",  "grey", "indianred1"), labels = c("aligned", "grey area","not ribozyme"))


ggsave(paste0("Length_Homology_",sample, "_noplus_noplusminus.png"), 
       p, width = 40, height = 24, units = "cm", dpi = 600)


###HI plus and plusminus aligned sequences

df = read.csv("Length_vs_alignment_score_H_I.csv")

df2 =read.csv("Intensities_HI.csv") #normalising with gel intensities from different lengths

sample = "H_I"

df3 = df[(df$classification=='plus_aligned'|df$classification=='plus_and_minus'),]

library(plyr)
counts <- ddply(df3, .(df3$len_seq, df3$score, df3$class), nrow)
names(counts) <- c("len_seq", "score", "class", "Freq")

counts$len_seq <- as.factor(counts$len_seq)
merge = merge(df2,counts, by ='len_seq')
merge$Freq_intensity = merge$Freq*merge$Intensity_factor


classes = c("aligned", "grey_area", "not_ribozyme")

theme_set(theme_classic()+
            theme(axis.text.x = element_text(size=20), 
                  axis.text.y = element_text(size=20), 
                  axis.title.x = element_text(size=20), 
                  axis.title.y = element_text(size=20), 
                  text = element_text(family = "sans"),
                  legend.title = element_blank(),
                  legend.text = element_text(size =20),
                  legend.key.size = unit(40, "point"),
                  legend.box.margin = margin (-10,-10,-10,-10),
                  strip.text = element_text(size = 5),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

p=ggplot(merge,
         aes(x = len_seq,
             y = score,
             fill = factor(class, levels = classes), 
             color = factor(class, levels = classes)
         ))+
  geom_point(aes(size = log10(Freq_intensity)))+
  scale_size_continuous(range=c(0.01,5))+
  labs( x= "Sequence Length (nt)",
        y= "Matching Nucleotides with TPR")+
  ylim(0,66)+
  scale_x_continuous(breaks = seq(6,66,3), limits = c(6,66))+
  
  scale_fill_manual("class", values = c( "skyblue3",  "grey", "indianred1"), labels = c("aligned", "grey area","not ribozyme"))+
  scale_color_manual("class", values = c( "skyblue3",  "grey", "indianred1"), labels = c("aligned", "grey area","not ribozyme"))



ggsave(paste0("Length_Homology_",sample, "_plus_plusminus.png"), 
       p, width = 40, height = 24, units = "cm", dpi = 600)


### random sequences all but plus and plus minus aligned sequences

df = read.csv("Length_vs_alignment_score_random_HI.csv")

df2 =read.csv("Intensities_HI.csv") #normalising with gel intensities from different lengths

sample = "random_HI" #still has same length distribution as HI though

df3 = df[!(df$classification=='plus_aligned'|df$classification=='plus_and_minus'),]

library(plyr)
counts <- ddply(df3, .(df3$len_seq, df3$score, df3$class), nrow)
names(counts) <- c("len_seq", "score", "class", "Freq")

counts$len_seq <- as.factor(counts$len_seq)
merge = merge(df2,counts, by ='len_seq')
merge$Freq_intensity = merge$Freq*merge$Intensity_factor


classes = c("aligned", "grey_area", "not_ribozyme")

theme_set(theme_classic()+
            theme(axis.text.x = element_text(size=20), 
                  axis.text.y = element_text(size=20), 
                  axis.title.x = element_text(size=20), 
                  axis.title.y = element_text(size=20), 
                  text = element_text(family = "sans"),
                  legend.title = element_blank(),
                  legend.text = element_text(size =20),
                  legend.key.size = unit(40, "point"),
                  legend.box.margin = margin (-10,-10,-10,-10),
                  strip.text = element_text(size = 5),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

p=ggplot(merge,
         aes(x = len_seq,
             y = score,
             fill = factor(class, levels = classes), 
             color = factor(class, levels = classes)
         ))+
  geom_point(aes(size = log10(Freq_intensity)))+
  scale_size_continuous(range=c(0.01,5))+
  labs( x= "Sequence Length (nt)",
        y= "Matching Nucleotides with TPR")+
  ylim(0,66)+
  scale_x_continuous(breaks = seq(6,66,3), limits = c(6,66))+
  
  scale_fill_manual("class", values = c( "skyblue3",  "grey", "indianred1"), labels = c("aligned", "grey area","not ribozyme"))+
  scale_color_manual("class", values = c( "skyblue3",  "grey", "indianred1"), labels = c("aligned", "grey area","not ribozyme"))



ggsave(paste0("Length_Homology_",sample, "_plus_plusminus.png"), 
       p, width = 40, height = 24, units = "cm", dpi = 600)


#HI random only plotting plus and plus minus aligned sequences
df = read.csv("Length_vs_alignment_score_random_HI.csv")

df2 =read.csv("Intensities_HI.csv") #normalising with gel intensities from different lengths

sample = "random_HI" #still has same length distribution as HI though

df3 = df[(df$classification=='plus_aligned'|df$classification=='plus_and_minus'),]

library(plyr)
counts <- ddply(df3, .(df3$len_seq, df3$score, df3$class), nrow)
names(counts) <- c("len_seq", "score", "class", "Freq")

counts$len_seq <- as.factor(counts$len_seq)
merge = merge(df2,counts, by ='len_seq')
merge$Freq_intensity = merge$Freq*merge$Intensity_factor


classes = c("aligned", "grey_area", "not_ribozyme")

theme_set(theme_classic()+
            theme(axis.text.x = element_text(size=20), 
                  axis.text.y = element_text(size=20), 
                  axis.title.x = element_text(size=20), 
                  axis.title.y = element_text(size=20), 
                  text = element_text(family = "sans"),
                  legend.title = element_blank(),
                  legend.text = element_text(size =20),
                  legend.key.size = unit(40, "point"),
                  legend.box.margin = margin (-10,-10,-10,-10),
                  strip.text = element_text(size = 5),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

p=ggplot(merge,
         aes(x = len_seq,
             y = score,
             fill = factor(class, levels = classes), 
             color = factor(class, levels = classes)
         ))+
  geom_point(aes(size = log10(Freq_intensity)))+
  scale_size_continuous(range=c(0.01,5))+
  labs( x= "Sequence Length (nt)",
        y= "Matching Nucleotides with TPR")+
  ylim(0,66)+
  scale_x_continuous(breaks = seq(6,66,3), limits = c(6,66))+
  
  scale_fill_manual("class", values = c( "skyblue3",  "grey", "indianred1"), labels = c("aligned", "grey area","not ribozyme"))+
  scale_color_manual("class", values = c( "skyblue3",  "grey", "indianred1"), labels = c("aligned", "grey area","not ribozyme"))



ggsave(paste0("Length_Homology_",sample, "_noplus_noplusminus.png"), 
       p, width = 40, height = 24, units = "cm", dpi = 600)