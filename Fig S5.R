

library(reshape2)
library(tibble)
library(dplyr)
library(data.table)
library(ggplot2)
library(vegan)
library(factoextra) 
library(vegan) 
library(grid) 

# load data
load("3.combined_ArgMrgIntTnp_allMappedSamples.RData")
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")
load("1.ARG_dat.RData")
load("1.phychem.sampleinfo.RData") 

regGenes <- unique(regGenes.v2_df$ARGsubtype) 
combined_deepargByGene_df <- combined_deepargByGene_df %>% filter(!ARG %in% regGenes)

# integrate data
dat <- combined_deepargByGene_df %>% dplyr::select(sampleID, ARG, DepthPG) %>% dcast(sampleID ~ ARG) 

rownames(dat) <- dat$sampleID 
dat <- dat[-1]
any(is.na(dat))

rownames(dat) <- tolower(sub( "-", "", sub("([\\-]?_?Q1)",  "", rownames(dat))))

ARG_dat <- dat
readme <- "ARG quantified as coverage/Gb"
save(ARG_dat,readme, file = "1.ARG_dat.RData") 

# read arg data ------------------------------------

all(rownames(phychem_df) %in% samples.all) 


performance <- NULL
arg <- ARG_dat[match(samples.all, rownames(ARG_dat)),]
otu <- phychem_df[match(samples.all, rownames(phychem_df)),] 

if(any(sapply(arg,sum) == 0)) arg <- arg[,-which(sapply(arg,sum) == 0)] else arg <- arg 
if(any(sapply(otu,sum) == 0)) otu <- otu[,-which(sapply(otu,sum) == 0)] else otu <- otu


# PCA --------
otu_pca <- rda(otu, scale = T) 

arg_hel <- decostand(arg, method = 'hellinger')  
arg_pca <- rda(arg_hel, scale = FALSE)

# PCoA  binary ---------
arg.dist <- vegdist(arg, method = "bray", binary = T)
arg.pcoa.bnr <- cmdscale(arg.dist, k = (nrow(arg) - 1), eig = TRUE)

otu.dist <- vegdist(otu, method = "bray", binary = T)
otu.pcoa.bnr <- cmdscale(otu.dist, k = (nrow(arg) - 1), eig = TRUE)

# PCoA numeric ----------
arg.dist <- vegdist(arg, method = "bray", binary = F)
arg.pcoa <- cmdscale(arg.dist, k = (nrow(arg) - 1), eig = TRUE)

otu.dist <- vegdist(otu, method = "bray", binary = F)
otu.pcoa <- cmdscale(otu.dist, k = (nrow(arg) - 1), eig = TRUE)


# Procrustes 
for(arg.methd in c("pca", "pcoa.binary","pcoa")){
  
  if(arg.methd == "pca"){
    site_arg <- summary(arg_pca, scaling = 1)$site
  }else if(arg.methd == "pcoa.binary"){
    site_arg <- data.frame(arg.pcoa.bnr$point)[1:2] 
  }else{
    site_arg <-  data.frame(arg.pcoa$point)[1:2] 
  }
  
  for(otu.methd in c("pca", "pcoa.binary", "pcoa") ) {
    if(otu.methd == "pca"){
      site_otu <- summary(otu_pca, scaling = 1)$site
    }else if(arg.methd == "pcoa.binary"){
      site_otu <- data.frame(otu.pcoa.bnr$point)[1:2] 
    }else{
      site_otu <-  data.frame(otu.pcoa$point)[1:2] 
    }
    
    proc <- procrustes(X = site_otu, Y = site_arg, symmetric = TRUE)
    
    set.seed(123)

    prot <- protest(X = site_otu, Y = site_arg, permutations = how(nperm = 999))
    
    prot$signif  
    prot$ss  
    
    vec <- c(arg.method = arg.methd, phychem.method = otu.methd,
             pvalue = prot$signif, m2 = round(prot$ss,3 ))
    
    performance <- bind_rows(performance, vec) 
  } 
}




## PCA - PCA ------------------------------

arg.methd = "pca"
if(arg.methd == "pca"){
  site_arg <- summary(arg_pca, scaling = 1)$site
}else if(arg.methd == "pcoa.binary"){
  site_arg <- data.frame(arg.pcoa.bnr$point)[1:2] 
}else{
  site_arg <-  data.frame(arg.pcoa$point)[1:2] 
}

otu.methd = "pca"
if(otu.methd == "pca"){
  site_otu <- summary(otu_pca, scaling = 1)$site
}else if(arg.methd == "pcoa.binary"){
  site_otu <- data.frame(otu.pcoa.bnr$point)[1:2] 
}else{
  site_otu <-  data.frame(otu.pcoa$point)[1:2] 
}  

# Procrustes

site_otu <- summary(otu_pca, scaling = 1)$site

proc <- procrustes(X = site_otu, Y = site_arg, symmetric = TRUE)
summary(proc)  
names(proc)

head(proc$Yrot)  
head(proc$X)  
proc$ss  
proc$rotation 

plot(proc, kind = 2)
Resi <- residuals(proc)
Resi[order(Resi)]  

set.seed(123)
prot <- protest(X = site_otu, Y = site_arg, permutations = how(nperm = 999))
prot

names(prot)
prot$signif  
prot$ss  


Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)


ymax = max(c(Y$X2,Y$PC2))
xmin = min(c(Y$X1,Y$PC1))
Y$sample <- rownames(Y)


Y$layer <- sapply(Y$sample, function(x)meta_df$layer[meta_df$sampleID == x])
Y$year <-  as.character(sapply(Y$sample, function(x)meta_df$year[meta_df$sampleID == x]))

ellipse_dat <- cbind.data.frame(
  sample = c(rownames(Y), rownames(Y)),
  X_area = c(Y$X1, Y$PC1),
  Y_area = c(Y$X2, Y$PC2)
) %>% mutate(layer = sapply(sample, function(x)meta_df$layer[meta_df$sampleID == x]),
             year =  as.character(sapply(sample, function(x)meta_df$year[meta_df$sampleID == x])))


# plot-----------------------------------------------------
ellipse_dat$layer  <- factor(ellipse_dat$layer ,levels = c('contaminated','remediate_inner','remediate_top'))#定义因子

library(ggplot2)
p <- ggplot(Y) +
  geom_point(aes(X1, X2, color = year), size = 5, shape = 16) + 
  geom_point(aes(PC1, PC2, color = year), size = 5, shape = 1) + 
  geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.1, 'cm')), 
               size = 0.3) +
  stat_ellipse(data = ellipse_dat , aes(X_area, Y_area, fill=layer), geom = 'polygon',
               level = 0.95, linetype =  "dashed", size = 1, alpha = 0.1, show.legend = TRUE) +  
  scale_color_manual(values = c('#3CA183', '#4DA0CE', '#D7A641')) +
  scale_fill_manual(values = c('#3CA183', '#4DA0CE', '#D7A641')) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.position = 'right') +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') + 
  theme(legend.title =element_blank())+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 17)) +
  theme(axis.text.y = element_text(size = 17)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 13)) +
  theme(title = element_text(size = 17)) +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  annotate('text', label = sprintf(paste('M^2 == ',round(prot$ss,2), sep = "") ), 
           x = xmin + 0.02, y = ymax - 0.01, size = 5, parse = TRUE) +
  annotate('text', label = paste('P < ', round(prot$signif,3), sep = ""), 
           x = xmin + 0.02, y = ymax - 0.03, size = 5, parse = TRUE) 

p

