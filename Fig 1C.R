
library(vegan)
library(ggplot2)
library(reshape2)
library(grid)

# load data
fo <- read.csv("diversity_da1.csv", header = 1, row.names = 1)

# analysis
bray_dis <- vegdist(fo, method = 'bray')
pcoa <- cmdscale(bray_dis, k = (nrow(fo) - 1), eig = T)
group <- read.csv("group.csv",header = 1, row.names = 1)
anosim.result <- anosim(bray_dis,group$group,permutations =999)
plot(anosim.result)
two_px_exp <- (pcoa$eig)[1:2]/sum(pcoa$eig)*100

pcoa_points <- data.frame((pcoa$points)[,1:2]) 
colnames(pcoa_points) <- c('PCoA1', 'PCoA2')

# plot
group$group  <- factor(group$group ,levels = c('UT','ULRT','ALRT'))
ggplot(data = pcoa_points, aes(PCoA1 ,PCoA2)) +
  geom_point(aes(fill = group$group) ,color = 'black' ,size =6.5 ,alpha = 0.8, shape = 21) +
  stat_ellipse(aes(fill = group$group ), geom = "polygon", level = 0.95, alpha = 0.2, show.legend = TRUE) +	
  scale_fill_manual(values = c('#3CA183', '#4DA0CE', '#D7A641')) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(hjust = 0.5), legend.position = 'right') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  theme(legend.title =element_blank())+
  labs(x = "PCoA1", y = "PCoA2")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 17)) + 
  theme(axis.text.y = element_text(size = 17)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 13)) +
  theme(title = element_text(size = 17)) 
