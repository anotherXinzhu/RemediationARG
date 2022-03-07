library(dplyr)
library(ggpmisc)
library(ggplot2)
library(ggpubr)

# load data
load("3.combined_ArgMrgIntTnp_allMappedSamples.RData")
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")

regGenes <- unique(regGenes.v2_df$ARGsubtype) 
totalARG_df.all_nonReg <- combined_deepargByGene_df %>% filter(!ARG %in% regGenes)

totalARG_df <- totalARG_df.all_nonReg %>% group_by(sampleID) %>% summarise(totalARG = sum(DepthPG))
totalARG_number_df <- totalARG_df.all_nonReg %>% group_by(sampleID) %>% dplyr::filter(RPKM!=0) %>% summarise(richness=n()) 

combined.dat <- merge(totalARG_df,
                      combined_IntTnp_df %>% 
                        mutate(MobilityGene = sapply( strsplit(target,"_", fixed = T),"[[",1) , 
                               MobGene.abund = DepthPG) %>% 
                        dplyr::select(sampleID,MobilityGene, MobGene.abund, sampleType), by="sampleID")




# plot
combined.dat$sampleType <- factor(combined.dat$sampleType, levels = c("contaminated", "remediate_inner", "remediate_top"))
combined.dat$MobilityGene <- factor(combined.dat$MobilityGene, levels = c("transposase", "integrase", "recombinase","resolvase", "class1Int", "clinicalInt1")) 
combined.dat$MobilityGene
p <- ggplot(combined.dat) +
  geom_point(aes(x = MobGene.abund, y=totalARG, fill=sampleType), shape=21, size=2.5) +
  geom_smooth(aes(x = MobGene.abund, y=totalARG),method = 'lm', formula = y~x, se = TRUE, show.legend = FALSE) + 
  stat_poly_eq(aes(x = MobGene.abund, y=totalARG, label = paste(..eq.label.., ..rr.label.., stat(p.value.label), sep = '~`,`~')),
               formula = y~x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7) +
  facet_wrap(.~MobilityGene, scales = "free") +
  theme(panel.grid = element_blank()) + 
  scale_color_manual(values =  c('#3CA183', '#4DA0CE', '#D7A641')) +
  scale_fill_manual(values =  c('#3CA183', '#4DA0CE', '#D7A641')) +
  ylim(0,2000) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 10)) + 
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 11)) + 
  theme(legend.title = element_text(size = 12)) +
  labs(x = "MGE (coverage,×/Gb)", y="Total ARG (coverage,×/Gb)")+theme_pubr(base_size = 14,border = T)
p

geom_smooth()
ggsave(p, filename = "4.corr_ARG.mobilityGene.pdf", device = "pdf", width = 7, height = 3.5)
ggsave(p, filename = "4.corr_ARG.mobilityGene.png", device = "png", width = 7, height = 3.5)
