
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyverse)
library(ggprism)

# load data--------------------------------------------
load("3.combined_ArgMrgIntTnp_allMappedSamples.RData")
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")

# analysis
regGenes <- unique(regGenes.v2_df$ARGsubtype) 
totalARG_df.all_nonReg <- combined_deepargByGene_df %>% filter(!ARG %in% regGenes)

totalARG_number_df <- totalARG_df.all_nonReg %>% group_by(sampleID) %>% dplyr::filter(RPKM!=0) %>% summarise(richness=n())

ARG_number <- transform(totalARG_number_df, Sample_type=rep(c("ALRT", "ALRT", "ULRT", "UT", "ULRT","UT"), each=3))
ARG_number <- transform(ARG_number, remediateYear=rep(c("3", "4", "4", "4", "3","3"), each=3))
Average <- read.csv("ARG_average_number.csv", header = 1, row.names = 1)

# plot
ARG_number$Sample_type <- factor(ARG_number$Sample_type, levels = c("UT", "ULRT", "ALRT"))
p <- ggplot() +
  geom_point(data=ARG_number, aes(x=Sample_type, y=richness, color=Sample_type, shape=factor(remediateYear)),size=9, alpha=0.6) + 
  geom_point(data=Average, aes(x=Sample_type, y=average),size=1, alpha=0.6) + 
  theme_prism() + theme(axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm')))) +
  scale_shape_manual(values = c(21,16)) + 
  ylim(150,500) + xlab("") +ylab("No.of ARG subtypes/Gb") + scale_x_discrete(expand = expansion(0.2,0)) + 
  scale_color_manual(values =  c('#3CA183', '#4DA0CE', '#D7A641')) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 18)) + 
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 15))
p
p <- p + geom_segment(aes(x = "UT",y = 1759.072,xend = "UT" ,yend = 1760))
p
