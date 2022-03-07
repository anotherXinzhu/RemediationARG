

library(dplyr)
library(data.table)

# load data-----------------------------------------------------------
load("1.binAbund_diversity.RData") 
load("2.ARG.inAllBins.RData")
load("2.Bin.taxon.df.RData")
meta <- fread("sampleInfo.txt") %>% mutate(sampleID = tolower(sampleID))
ARGMapping_df <- fread("deeparg_gene_mapping.txt")


# integrate data-----------------------------------------------------
bin980.abund_df.l <- merge(
  bin980.Coverage_df.w %>% 
    dplyr::select(-secondary_cluster,-cluster.rep.bin) %>% 
    reshape2::melt(variable.name="sample",value.name = "coverage"),
  bin980.DepthPG_df.w %>%
    dplyr::select(-secondary_cluster,-cluster.rep.bin) %>% 
    reshape2::melt(variable.name="sample",value.name = "DepthPG"),
  by=c("genome","sample")
)

bin980.abund_df.l$genome <- sub("\\.filtered\\.fa", "", bin980.abund_df.l$genome)

tmp <- merge(ARG.bins_df %>% group_by(bin.name, ARG) %>% summarise(numORF = n()), 
             bin.taxon_df %>% dplyr::select(-classification), 
             by.x = "bin.name", by.y = "user_genome")

dat<-merge(tmp, 
           bin980.abund_df.l %>% dplyr::select(genome, sample, DepthPG),
           by.x =  "bin.name", by.y = "genome")

dat$ARGDepthPG = dat$numORF * dat$DepthPG



dat.tmp <- dat 
colnames(dat.tmp)[colnames(dat.tmp) == "phylum" ] <- "Taxon"
ARGtype.phylumContrib_df <- dat.tmp %>%
  mutate(ARGtype = sapply(ARG, function(x)ARGMapping_df$drug_type[ARGMapping_df$deeparg_genes == x])) %>%
  group_by(sample, ARGtype, Taxon) %>% summarise(ARGDepthPG = sum(ARGDepthPG)) %>% 
  mutate(layer = sapply(sample, function(x)meta$layer[meta$sampleID == x])) %>%
  group_by(layer, ARGtype, Taxon) %>%
  summarise(sum.ARGDepthPG = sum(ARGDepthPG), avg.ARGDepthPG = mean(ARGDepthPG)) %>%
  group_by(layer, ARGtype) %>%
  mutate(taxonContrib = avg.ARGDepthPG/sum(avg.ARGDepthPG))




#  plot ----------------------
library(ggplot2)
load("2.ARGtypes.RData") 

plotDat <- ARGtype.phylumContrib_df %>%
  filter(!ARGtype %in% c(otherTypes,"unclassified")) 

corrected_taxons <- vector("character",nrow(plotDat))
for(i in c(1:nrow(plotDat))){
  #i=3
  Taxon <- plotDat$Taxon[i]
  
  if(Taxon=="") next
  parts <- strsplit(Taxon," ",fixed = T)[[1]]
  
  for(pt in parts){
    j=which(parts == pt)
    
    parts[j] <- sub("_[A-Z]+","",pt)
    
  }
  
  corrected_taxons[i] <- paste(parts, collapse = " ")
}
plotDat$Taxon <- corrected_taxons

# bar plot -------------------
majorPhyla4 <- c("Thermoplasmatota","Proteobacteria","Actinobacteriota","Acidobacteriota","Others" )

plotDat.bar = plotDat
plotDat.bar$Taxon[!plotDat.bar$Taxon %in% majorPhyla4] <- "Others"

plotDat.bar <- 
  plotDat.bar %>%
  group_by(layer,ARGtype,Taxon) %>% summarise(taxonContrib = sum(taxonContrib))

plotDat.bar %>% group_by(layer, ARGtype) %>% summarise(test=sum(taxonContrib))


plotDat.bar$Taxon <- factor(plotDat.bar$Taxon, levels=rev(majorPhyla4))
plotDat.bar$ARGtype <- factor(plotDat.bar$ARGtype,
                              levels = c("mupirocin","bacitracin","sulfonamide","tetracycline","multidrug",
                                         "aminoglycoside","glycopeptide","MLS","phenicol","beta-lactam"))
Colors <- rev(c("#CA9E44","#DFD180","#658633","#27813A","Gray"))

FigS4.tmp <-  ggplot(plotDat.bar) +
  geom_col(aes(x=layer, y=taxonContrib,  fill=Taxon)) +
  facet_wrap(.~ARGtype ,nrow =2 )+
  scale_fill_manual(values = Colors) +
  scale_x_discrete(labels=c("UT","ULRT","ALRT"))+
  theme_bw() +
  theme(strip.background = element_blank(),panel.grid = element_blank())

FigS4.tmp

ggsave(FigS4.tmp, filename = "5.FigS4tmp.majorTaxonContribution.pdf",
       device = "pdf", height = 6, width = 14)
