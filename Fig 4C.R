setwd("G:/5.repaire_ARG/repaire_importance/PathoFact/bins_VF")

library(dplyr)
library(tidyverse)
library(reshape2)

# load data------------------------------------
VFs_bins <- read.csv("VFS_bins.csv")
VFs_bins

# integrate data-------------------------------
VFs_bins <- VFs_bins %>% group_by(sample, virus.factors) %>% summarise(sum=sum(sum))
VFs_bins <- transform(VFs_bins,layer.abb=c("ALRT","ALRT","ALRT","ALRT","ALRT","ALRT","ALRT","ALRT","ALRT","ALRT","ALRT","ALRT","ULRT","ULRT","ULRT",
                                           "ULRT","ULRT","ULRT","UT","UT","UT","UT","UT","UT","ULRT","ULRT","ULRT","ULRT","ULRT","ULRT","UT","UT","UT","UT","UT","UT"))
VFs_bins_all <- VFs_bins %>% group_by(layer.abb,virus.factors) %>% summarise(avg=mean(sum),sd=sd(sum)) #求平均值和偏差

library(ggplot2)
library(ggsci)
library(reshape2)
library(ggprism)

VFs_bins_all$layer.abb <- factor(VFs_bins_all$layer.abb, levels=c("UT","ULRT","ALRT"))
VFs_bins_all$virus.factors <- factor(VFs_bins_all$virus.factors, levels = c("VF1","VF2"))

p <-ggplot(VFs_bins_all, aes(x=layer.abb,y=avg)) + 
  geom_bar(aes(fill=virus.factors), stat= "identity", position = position_dodge(0.65), width = 0.6, alpha=0.7) +
  geom_pointrange( aes(x=layer.abb, y=avg, ymin=avg-sd, ymax=avg+sd), stat= "identity", position = position_dodge2(width = 0.65),colour="darkgray", alpha=1, size=1.3)  +
  #geom_errorbar(aes(ymin = avg-sd, ymax = avg+sd), position = position_dodge2(width = 0.25))+
  scale_fill_manual(values = c('#698B97', '#C2D3CD')) +
  scale_x_discrete(expand = expansion(0.2,0)) + 
  ylim(0,40000) +
  theme_prism(base_size = 14,
              base_line_size = 1,
              base_rect_size = 2) + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm')))) + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        #legend.position = 'right',
        legend.title = element_blank()) +
  xlab('') + ylab('Abundance of VF in AC-MAGs')
p

