setwd("G:/5.repaire_ARG/repaire_importance/ARG_mobility_Gene")

library(dplyr)
library(ggpmisc)
library(ggplot2)
library(tidyverse)
library(ggprism)


# load data
load("3.combined_ArgMrgIntTnp_allMappedSamples.RData")
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")

# analysis
regGenes <- unique(regGenes.v2_df$ARGsubtype) 
totalARG_df.all_nonReg <- combined_deepargByGene_df %>% filter(!ARG %in% regGenes)

totalARG_df <- totalARG_df.all_nonReg %>% group_by(sampleID) %>% summarise(totalARG = sum(DepthPG))
totalARG_number_df <- totalARG_df.all_nonReg %>% group_by(sampleID) %>% dplyr::filter(RPKM!=0) %>% summarise(richness=n()) 

combined.dat <- merge(totalARG_df,
                      combined_IntTnp_df %>% 
                        mutate(MobilityGene = sapply( strsplit(target,"_", fixed = T),"[[",1) , 
                               MobGene.abund = DepthPG) %>% 
                        dplyr::select(sampleID,MobilityGene, MobGene.abund, sampleType), by="sampleID")


combined.dat2 <- merge(totalARG_df,
                       combined_IntTnp_df %>% 
                         dplyr::filter(target == "transposase_Eval1e-10"|target == "resolvase_Eval1e-10"|target == "recombinase_Eval1e-10"|target == "integrase_Eval1e-10") %>%
                         mutate(MobilityGene = sapply( strsplit(target,"_", fixed = T),"[[",1) , 
                                MobGene.abund = DepthPG) %>% 
                         dplyr::select(sampleID,MobilityGene, MobGene.abund, sampleType), by="sampleID") 
combined.dat2
totalMGE_df <- combined.dat2 %>% group_by(sampleID) %>% summarise(totalMGE = sum(MobGene.abund)) 

# plot

Average <- read.csv("Total_average_MGE.csv", header = 1, row.names = 1)

totalMGE_df$Sample_type <- c("ALRT","ALRT","ALRT","ALRT","ALRT","ALRT","ULRT", "ULRT","ULRT","UT","UT","UT","ULRT", "ULRT","ULRT","UT","UT","UT")#新增加一列列名
totalMGE_df$Remediate_year <- c("3","3","3","4","4","4","4","4","4","4","4","4","3","3","3","3","3","3")
totalMGE_df$Sample_type <- factor(totalMGE_df$Sample_type, levels = c("UT", "ULRT", "ALRT")) #指定想要呈现的形式

p <- ggplot() + 
  geom_point(data=totalMGE_df, aes(x=Sample_type, y=totalMGE, color=Sample_type, shape=Remediate_year),size=7, alpha=0.7) + 
  geom_point(data=Average, aes(x=Sample_type, y=average),size=1, alpha=0.7) + 
  theme_prism() + theme(axis.line = element_line(arrow = arrow(length = unit(0.3, 'cm')))) +
  scale_shape_manual(values = c(21,16)) +
  ylim(0,15000) + xlab("") +ylab("MGE (Coverage, ×/Gb)") + scale_x_discrete(expand = expansion(0.2,0)) + 
  scale_color_manual(values =  c('#3CA183', '#4DA0CE', '#D7A641')) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 15)) + #将X轴旋转,以及修改相关字体
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.text = element_text(size = 13))
p
p <- p +theme(axis.title = element_text(size = 18, color = "gray10"))
p <- p + theme(axis.title.x = element_text(size = 13)) + theme(axis.title.y = element_text(size = 14))
p 
