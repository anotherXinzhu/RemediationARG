

library(dplyr)
library(data.table)
library(ggplot2)

# load data----------------------------
load("1.closeConn.dat.RData")
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")
meta <- fread("sampleInfo.txt")

# integrate data-----------------------
regGenes <- unique(regGenes.v2_df$ARGsubtype) 
closeConn_dat.all <- closeConn_dat.all %>% filter(!ARG %in% regGenes)

meta$sampleID <- toupper(meta$sampleID) 

plotDat_all <- closeConn_dat.all %>%
  dplyr::select(sample_name, ARGtype, minDist) %>% 
  dplyr::mutate(layer = sapply(sample_name,function(x) meta$layer[meta$sampleID == x])) 

plotDat_all <- closeConn_dat.all %>%
  dplyr::select(sample_name, minDist) %>%
  mutate(layer = sapply(sample_name,function(x) meta$layer[meta$sampleID == x]))

plotDat_all_average <- plotDat_all %>% group_by(layer) %>% summarise(mean = mean(minDist))

# plot

plotDat_all$minDist <- plotDat_all$minDist/1000 

p <- plotDat_all %>%
  ggplot(aes(x = layer, y = minDist)) +
  stat_summary(fun.y = "mean", geom = "bar", position = "identity",fill= c('#3CA183', '#9BC4DA','#D7A641'), alpha = 0.7, width = 0.5) + 
  geom_point(shape=1,color="darkgray",alpha=0.7, size=1, position = position_jitter(width = 0.1)) + 
  scale_x_discrete(expand = expansion(0.2,0)) +
  coord_cartesian(ylim = c(0,250)) +
  theme_prism(base_size = 14,
              base_line_size = 1,
              base_rect_size = 2) + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm')))) +
  theme(axis.text.x = element_text (hjust = 1 )) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_blank()) +
  xlab("") + ylab("ARG-MGE neareast distance (kb)") 
p
