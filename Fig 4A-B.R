
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

# load data-----------------------------------------
load("1.binAbund_diversity.RData")
load("2.ARG.inAllBins.RData")
load("2.Bin.taxon.df.RData")
sample.info <- fread("sampleInfo.txt") 
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")
sample.info$sampleID <- tolower(sample.info$sampleID) 

# iintegrate data-----------------------------------

regGenes <- unique(regGenes.v2_df$ARGsubtype)
nonRegARG.bins_df <- ARG.bins_df %>% filter(!ARG %in% regGenes) 
length(unique(nonRegARG.bins_df$bin.name))

ACbin_dat <-  merge(bin.dat_combined %>% mutate(bin.id = sub("\\.filtered\\.fa$","", bin)),
                    bin.taxon_df %>% dplyr::select(user_genome, phylum, family, genus, lastTaxon),
                    by.x = "bin.id", by.y = "user_genome") %>% 
  filter(bin.id %in% nonRegARG.bins_df$bin.name) 
ACbin_dat1 <- ACbin_dat %>% group_by(sample) %>% dplyr::select(bin, DepthPG,sample)  

ACbin_dat1 <- spread(ACbin_dat1, bin, DepthPG) 
ACbin_dat1$sumbin <- rowSums(ACbin_dat1[,2:615]) 
ACbin_dat2 <- ACbin_dat1[,c(1,616)]  

ACbin_dat$layer.abb <- sapply(ACbin_dat$sample,
                              function(x) sample.info$lay.abb[which(sample.info$sampleID == x)]) 

TaxLvl <- "phylum" 

dat.tmp <- ACbin_dat
colnames(dat.tmp)[colnames(dat.tmp) == TaxLvl] <- "taxon"  

# correct taxon
corrected_taxon <- vector("character",nrow(dat.tmp))  

for(i in c(1:nrow(dat.tmp))){ 
  #i=3
  taxon <- dat.tmp$taxon[i] 
  
  if(taxon=="") next
  parts <- strsplit(taxon," ",fixed = T)[[1]]
  
  for(pt in parts){
    j=which(parts == pt)
    
    parts[j] <- sub("_[A-Z]+","",pt)
    
  }
  
  corrected_taxon[i] <- paste(parts,collapse = " ")
}
dat.tmp$taxon <- corrected_taxon      

# Figure 4a
plotDat <- dat.tmp %>% 
  group_by(sample) %>%
  summarise(ACbin.depthPG = sum(DepthPG))

plotDat <- transform(plotDat, layer.abb=c("ALRT","ALRT","ALRT","ALRT","ALRT","ALRT","ULRT","ULRT","ULRT","UT","UT","UT","ULRT","ULRT","ULRT","UT","UT","UT"))
plotDat <- plotDat %>% group_by(layer.abb) %>% summarise(avg=mean(ACbin.depthPG),sd=sd(ACbin.depthPG)) 
plotDat$layer.abb <- factor(plotDat$layer.abb, levels = c("UT","ULRT","ALRT"))

p2<-ggplot(plotDat) +
  geom_col(aes(x=layer.abb,y=avg), stat="identity", width = 0.5, size = 0.05, color = "black", fill = c('#D7A641', '#9BC4DA', '#3CA183'), alpha = 0.7) +
  geom_pointrange( aes(x=layer.abb, y=avg, ymin=avg-sd, ymax=avg+sd), stat= "identity", position = position_dodge2(width = 0.65),colour="darkgray", alpha=1, size=1.3)  +
  scale_x_discrete(expand = expansion(0.2,0)) + 
  ylim(0,150) +
  theme_prism(base_size = 14,
              base_line_size = 1,
              base_rect_size = 2) + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm')))) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_blank()) +
  xlab('') + ylab('Average bin (coverage, ×/Gb)')
p2

# Figure 4b
plotDat <- dat.tmp %>% 
  group_by(layer.abb, taxon) %>% 
  summarise(totalDepthPG = sum(DepthPG)) %>%
  mutate(freq = totalDepthPG/sum(totalDepthPG))

plotDat %>% group_by(layer.abb) %>% summarise(t=sum(freq)) 
N=10
taxRank_Df <- plotDat %>% group_by(taxon) %>% summarise(avg.freq = mean(freq), sd = sd (freq)) %>% arrange(desc(avg.freq))
sum(taxRank_Df$avg.freq)
otherTaxon <- taxRank_Df$taxon[(N+1):nrow(taxRank_Df)]

plotDat$taxon[plotDat$taxon %in% otherTaxon] <- "other"
plotDat <- plotDat %>% group_by(layer.abb, taxon) %>% summarise(freq= sum(freq))

# plotting

plotDat$taxon <- factor(plotDat$taxon,
                        levels = c("other", "Acidobacteriota","Actinobacteriota", "Bacteroidota","Firmicutes","Gemmatimonadota",
                                   "Myxococcota", "Nanoarchaeota","Nitrospirota","Proteobacteria","Thermoplasmatota"))

library(ggplot2)
library(ggsci)
library(reshape2)
library(ggprism)
plotDat$layer.abb <- factor(plotDat$layer.abb, levels = c("UT","ULRT","ALRT"))

p <- ggplot(data = plotDat,aes(x=layer.abb, y=freq)) + 
  geom_bar(aes(fill = taxon),stat = 'identity',
           color = 'black',size =0.5,
           width = 0.5) + 
  scale_x_discrete(expand = expansion(0.2,0)) + 
  scale_fill_manual(values =  c('ForestGreen', 'OliveDrab', 'DarkSeaGreen', 'LightCyan3',
                                'LightBlue2', 'SkyBlue4', '#4DA0CE','Beige',
                                'Wheat', 'LightGoldenrod', '#D7A641'))  +
  
  theme_prism(base_size = 14,
              base_line_size = 1,
              base_rect_size = 2) + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm')))) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'right',
        legend.title = element_blank()) +
  xlab('') + ylab('Phyla composition of AC-MAGs') 

p

