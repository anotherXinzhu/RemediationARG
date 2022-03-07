
library(ggplot2)
library(tidyverse)
library(ggsci)
library(reshape2)

# load data
composition <- read.csv("16sphylum.csv", sep = ",", header = 1)

# analysis
composition_df <- composition %>% arrange(desc(sum))
percent(composition_df$sum, accuracy = 0.1)
composition_rank <- c(composition_df$X[1:10],"other")
otherTypes <- composition_df$X[11:nrow(composition_df)]

tmpDat <- composition %>% select(-sum)
tmpDat$X[which(tmpDat$X %in% otherTypes)] <- "other" 
tmpDat_width <- melt(tmpDat)

tmpDat <- transform(tmpDat_width, layer.abb=rep(c("ALRT", "ULRT", "UT", "ALRT", "ULRT","UT"), each=183))
tmpDat1 <- tmpDat %>% group_by(X, layer.abb) %>% summarise(sum=sum(value))
tmpDat2 <- tmpDat1 %>% group_by(layer.abb) %>% summarise(sum(sum))
tmpDat2 <- tmpDat1 %>% group_by(layer.abb) %>% dplyr::mutate(freq = sum/sum(sum))
aa <- tmpDat2 %>% group_by(layer.abb) %>% summarise(sum(freq)) 

# plot
tmpDat2$X <- factor(tmpDat2$X, levels = c("other","Chloroflexi","Verrucomicrobia","Firmicutes","Planctomycetes","Bacteroidetes","Nitrospirae","Acidobacteria","Actinobacteria","Euryarchaeota","Proteobacteria"))
tmpDat2$layer.abb <- factor(tmpDat2$layer.abb, levels = c("UT", "ULRT", "ALRT"))

p <- ggplot(data = tmpDat2,aes(x = layer.abb,y = freq)) + 
  geom_bar(aes(fill = tmpDat2$X),stat = 'identity',
           color = 'black',size =0.5,
           width = 0.5) + 
  scale_x_discrete(expand = expansion(0.2,0)) + 
  scale_fill_manual(values =  c('ForestGreen', 'OliveDrab', 'DarkSeaGreen', 'LightCyan3',
                                'LightBlue2', 'SkyBlue4', '#4DA0CE','Beige',
                                'Wheat', 'LightGoldenrod', '#D7A641'))  + 
  theme_prism(base_size = 14,
              base_line_size = 0.5,
              base_rect_size = 2) + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm')))) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.position = 'right',
        legend.text = element_text(size = 15)) +
  xlab('') + ylab('Phyla composition') 
p
