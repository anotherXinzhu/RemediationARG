
library(dplyr)
library(reshape2)
library(grid) 
library(data.table)
library(ggsci)
library(scales)
library(ggplot2)
library(reshape2)

# load quantificaiton data -------------------
load("3.combined_ArgMrgIntTnp_allMappedSamples.RData")
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")

regGenes <- unique(regGenes.v2_df$ARGsubtype) 
totalARG_df.all_nonReg <- combined_deepargByGene_df %>% filter(!ARG %in% regGenes)

dat <- totalARG_df.all_nonReg %>% select(sampleID, ARG, DepthPG) %>% dcast(sampleID ~ ARG)

rownames(dat) <- dat$sampleID
dat <- dat[-1]
any(is.na(dat))

if(any(sapply(dat,sum) == 0)) dat_positive <- dat[,-which(sapply(dat,sum) == 0)] else dat_positive <- dat

# grouping information
meta_df <- combined_deepargByGene_df %>% select(sampleID, sampleType, remediateYear) %>% unique()

meta_df <-  meta_df[match(rownames(dat_positive), meta_df$sampleID),]


# analysis
library(vegan)
dat_arg <- dat_positive
groups_df <- meta_df
rownames(dat_arg) == groups_df$sampleID




arg.dist <- vegdist(dat_arg, method = "bray", binary = FALSE)

pcoa <- cmdscale(arg.dist, k = (nrow(dat_arg) - 1), eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa$eig/sum(pcoa$eig) 

# plot

site <- data.frame(pcoa$point)[1:2] 
site$name <- rownames(site)

meta_df$sampleID == site$name 
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')

#g="ecosystem"
site$group <- meta_df[,"sampleType"] #color by position
arg.ano <- with(dat_arg, anosim(arg.dist, meta_df[,"sampleType"]))
#summary(arg.ano)
arg.ano$statistic # R
arg.ano$signif # p

#ggplot2 
grob1 <- grobTree(textGrob(paste("Anosim.R=",round(arg.ano$statistic, 4),
                                 ", pvalue=",arg.ano$signif,
                                 sep = ""), 
                           x=0.05,  y=0.1, hjust=0, gp=gpar(col="black", fontsize=13)))

p <- ggplot(data = site, aes(X1, X2)) +
  geom_point(aes(color = group),size=7) +
  geom_point(size=7,shape=21) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.2, show.legend = TRUE) +    
  scale_color_manual(values = c('#3CA183', '#9BC4DA', '#D7A641',"steelblue3")) +
  scale_fill_manual(values = c('#3CA183', '#9BC4DA', '#D7A641',"steelblue3")) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(hjust = 0.5), legend.position = 'right') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  labs(x = pcoa1, y = pcoa2, title = paste('PCoA, distance: ',distMethod,", Binary: ",ifBinary,sep = ""))+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 17)) + 
  theme(axis.text.y = element_text(size = 17)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 13)) +
  theme(title = element_text(size = 17)) +
  annotation_custom(grob1)

ggsave(p, filename = "Fig2B.pdf", device = "pdf")
