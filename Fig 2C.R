

library(data.table)
library(dplyr)
library(reshape2)
library(scales)
library(ggplot2)
library(tidyverse)
library(ggprism)
library(tidyr)

# load data-----------------------------
load("3.combined_ArgMrgIntTnp_allMappedSamples.RData")
load("2.ARGtypes.RData") 
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")

# analysis
regGenes <- unique(regGenes.v2_df$ARGsubtype) 
totalARG_df.all_nonReg <- combined_deepargByGene_df %>% filter(!ARG %in% regGenes)

avgRelAbund_df <- totalARG_df.all_nonReg %>% select(sampleID, drug_type, RPKM) %>% filter(drug_type != "unclassified") %>%
  group_by(sampleID,drug_type) %>% summarise(rpkm=sum(RPKM)) %>%
  group_by(drug_type) %>% summarise(avgRpkm = mean(rpkm)) %>% mutate(relAbund=avgRpkm/sum(avgRpkm)) %>% arrange(desc(relAbund))


ARGtype_rank <- c(avgRelAbund_df$drug_type[1:10],"other")
otherTypes <- avgRelAbund_df$drug_type[11:nrow(avgRelAbund_df)]

tmpDat <- totalARG_df.all_nonReg %>% select(ARG,sampleID,DepthPG,drug_type)

tmpDat$drug_type[which(tmpDat$drug_type %in% otherTypes)] <- "other" 

ARG_dominated_composition <- tmpDat %>% filter(drug_type!="unclassified") %>% 
  group_by(sampleID, drug_type) %>% summarise(DepthPG=sum(DepthPG)) %>% mutate(relAbund=DepthPG/sum(DepthPG))

sample.info <- fread("sampleInfo.txt") 
sample.info$sampleID <- tolower(sample.info$sampleID) 

ARG_dominated_composition$sampleID[which(ARG_dominated_composition$sampleID =='TAS-2-Q1')] <- 'tas2'
ARG_dominated_composition$sampleID[which(ARG_dominated_composition$sampleID =='TAS-4-Q1')] <- 'tas4'
ARG_dominated_composition$sampleID[which(ARG_dominated_composition$sampleID =='TAS-6-Q1')] <- 'tas6'
ARG_dominated_composition$sampleID[which(ARG_dominated_composition$sampleID =='TAX-2_Q1')] <- 'tax2'
ARG_dominated_composition$sampleID[which(ARG_dominated_composition$sampleID =='TAX-4_Q1')] <- 'tax4'
ARG_dominated_composition$sampleID[which(ARG_dominated_composition$sampleID =='TAX-6_Q1')] <- 'tax6'
ARG_dominated_composition$sampleID[which(ARG_dominated_composition$sampleID =='TUT1_Q1')] <- 'tut1'
ARG_dominated_composition$sampleID[which(ARG_dominated_composition$sampleID =='TUT21_Q1')] <- 'tut21'
ARG_dominated_composition$sampleID[which(ARG_dominated_composition$sampleID =='TUT4_Q1')] <- 'tut4'

ARG_dominated_composition$layer.abb <- sapply(ARG_dominated_composition$sampleID,
                                              function(x) sample.info$lay.abb[which(sample.info$sampleID == x)])

ARG_dominated_composition_arrange <- dplyr::arrange(ARG_dominated_composition, desc(layer.abb), drug_type)
ARG_average_dominated_composition_arrange <- ARG_dominated_composition_arrange %>% group_by(layer.abb, drug_type) %>% summarise(avrage=mean(relAbund))
ARG_average_dominated_composition_arrange %>% group_by(layer.abb) %>% summarise(sum(avrage))

ARG_average_dominated_composition_arrange$drug_type <- factor(ARG_average_dominated_composition_arrange$drug_type, levels = c("other", "mupirocin", "sulfonamide", "phenicol", "MLS", "aminoglycoside","glycopeptide", "bacitracin", "beta-lactam", "tetracycline","multidrug"))
ARG_average_dominated_composition_arrange$layer.abb <- factor(ARG_average_dominated_composition_arrange$layer.abb, levels = c("UT", "ULRT", "ALRT"))

# plot
p <- ggplot(data = ARG_average_dominated_composition_arrange,aes(x = layer.abb,y = avrage)) + 
  geom_bar(aes(fill = drug_type),stat = 'identity',
           color = 'black',size =0.5,
           width = 0.5) + 
  scale_x_discrete(expand = expansion(0.2,0)) + #+ scale_y_discrete(expand = expansion(0,0.2))
  #scale_fill_npg() + 
  #scale_fill_aaas() + 
  scale_fill_manual(values =  c('ForestGreen', 'OliveDrab', 'DarkSeaGreen', 'LightCyan3',
                                'LightBlue2', 'SkyBlue4', '#4DA0CE','Beige',
                                'Wheat', 'LightGoldenrod', '#D7A641'))  + 
  #'#3CA183', '#4DA0CE', '#D7A641'
  theme_prism(base_size = 14,
              base_line_size = 0.5,
              base_rect_size = 2) + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm')))) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.position = 'right',
        legend.text = element_text(size = 15)) +
  xlab('') + ylab('Composition of ARG types') 
#scale_y_continuous(expand = expansion(0,0))
p
