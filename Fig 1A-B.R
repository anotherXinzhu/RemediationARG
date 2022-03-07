
library(vegan)
library(dplyr)
library(scales)
library(ggpubr)
library(ggsci)
library(ggprism)

# load data ------------------------------
diversity_da <- read.csv("diversity_da.csv", header = 1, row.names = 1)

# analysis
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)    
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)   
}
richness <- alpha_index(t(diversity_da), method = 'richness', base = exp(1))

shannon <- alpha_index(t(diversity_da), method = 'shannon', base = exp(1))
all_diversity <- data.frame(cbind(richness,shannon))

all_diversity_richeness_average <- all_diversity %>% group_by(layer.abb) %>% summarise(average=mean(richness))
all_diversity_shannon_average <- all_diversity %>% group_by(layer.abb) %>% summarise(average=mean(shannon))

# plot for Richeness
p1 <- ggplot() + 
  geom_point(data=all_diversity, aes(x=layer.abb, y=richness, color=layer.abb, shape=remediation_year),size=9, alpha=0.6) + 
  geom_point(data=all_diversity_richness_average, aes(x=layer.abb, y=average),size=1, alpha=0.7) + 
  theme_prism() + theme(axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm')))) +
  scale_shape_manual(values = c(21,16)) +
  ylim(0,650) + xlab("") +ylab("Richness") + scale_x_discrete(expand = expansion(0.2,0)) + 
  scale_color_manual(values =  c('#3CA183', '#4DA0CE', '#D7A641')) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 18)) + 
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 15))
p1

p1 <- p1 + geom_segment(aes(x = "UT",y = 1759.072,xend = "UT" ,yend = 1760)) 


# plot for Shannon
p2 <- ggplot() + 
  geom_point(data=all_diversity, aes(x=layer.abb, y=shannon, color=layer.abb, shape=remediation_year),size=9, alpha=0.6) + 
  geom_point(data=all_diversity_shannon_average, aes(x=layer.abb, y=average),size=1, alpha=0.7) + 
  theme_prism() + theme(axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm')))) +
  scale_shape_manual(values = c(21,16)) +
  ylim(0,6.5) + xlab("") +ylab("Richness") + scale_x_discrete(expand = expansion(0.2,0)) + 
  scale_color_manual(values =  c('#3CA183', '#4DA0CE', '#D7A641')) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 18)) + 
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 15))
p2

p2 <- p2 + geom_segment(aes(x = "UT",y = 1759.072,xend = "UT" ,yend = 1760)) 
