setwd("C:/Users/Teresa/OneDrive - University Of Cambridge/St Johns Cambridge/1-Part II Biochemistry/Part II Project/Files/Randomer sequencing")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr) 
library(RColorBrewer)


df = read.csv("Peptide_length_family_box_length_HI.csv")
df2 = read.csv("Peptide_length_family_box_length_HI_codon_random.csv")
df_merge = merge(df,df2, by ='Length')

names(df_merge)[names(df_merge) == "Frequency.x"] <- "Sequences.HI"
names(df_merge)[names(df_merge) == "Frequency.y"] <- "Random.Sequences"

long = gather(df_merge, Sample, Frequency, Sequences.HI:Random.Sequences, factor_key=TRUE)
long$Length <- as.factor(long$Length)

theme_set(theme_classic()+
            theme(axis.text.x = element_text(size=5), 
                  axis.text.y = element_text(size=5), 
                  axis.title.x = element_text(size=5), 
                  axis.title.y = element_text(size=5), 
                  text = element_text(family = "sans"),
                  legend.title = element_blank(),
                  legend.text = element_text(size =5),
                  legend.key.size = unit(8, "point"),
                  legend.box.margin = margin (-10,-10,-10,-10),
                  strip.text = element_text(size = 5),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

p=ggplot(long,
         aes(x = Length, y = Frequency, fill = Sample))+
  # facet_wrap(vars(Sample))+
  geom_bar(stat = 'identity', position='dodge')+
  labs( x= "Length of Peptides in Amino acids",
        y= "Percentage of Peptides")+
  
  scale_y_continuous(expand = c(0,0), labels = scales::percent)+
  
  scale_fill_manual("Sample", values = c( "indianred1",  "grey"), labels = c("Sequences HI", "Random Sequences"))

print(p)


ggsave(paste0("Family_box_codons", ".png"), 
       p, width = 10, height = 6, units = "cm", dpi = 300)