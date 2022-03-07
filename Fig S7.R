
library(dplyr)
library(data.table)
library(ggpubr)
library(ggsci)

# load data----------------------------------------------------------
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")
regGenes <- regGenes.v2_df$ARGsubtype %>% unique()
load("1.binAbund_diversity.RData")
load("3.combined_ArgMrgIntTnp_allMappedSamples.RData")
load("2.ARG.inAllBins.RData")
load("parsed_Taxon.RData")

sigBins <- 
  fread("8_2.SigBins_numARG_noCtrlVars_v1.txt",header = F)$V1 
sigBins <- paste(sigBins, ".filtered.fa",sep = "")

# Figure S7A---------------------------------------------------------
# integrate data-----------------------------------------------------
remove(combined_deepargByDrug_df,combined_deepargByMech_df,combined_IntTnp_df, combined_metalByCompound_df,
       combined_metalByGene_df,combined_metalByMech_df, combined_metalByType_df)
numARG.dat <-combined_deepargByGene_df %>%
  filter(!ARG %in% regGenes) %>%
  filter(Coverage > 0) %>% 
  group_by(sampleID) %>%
  summarise(numARG = n()) 

numARG.dat$sampleID <- tolower(sub("Q1", "", gsub("[_\\-]","",numARG.dat$sampleID) ))

bin.dat <- nrBin615.DepthPG_df.l %>%
  reshape2::dcast(sample~bin, value.var = "DepthPG") %>%
  tibble::column_to_rownames("sample") 

bin.dat <- bin.dat[match(numARG.dat$sampleID, rownames(bin.dat)),]

Corr.res <- NULL
for(i in c(1:ncol(bin.dat))){
  ct <- cor.test(numARG.dat$numARG, bin.dat[,i])
  
  vec <- c(colnames(bin.dat)[i],ct$estimate,ct$p.value)
  names(vec)<-c("bin","r",'p')
  
  Corr.res <- bind_rows(Corr.res, vec)
}

# compute ARG number ----------------------

numARG.bins.dat <- ARG.bins_df %>% 
  filter(!ARG %in% regGenes) %>%
  dplyr::select(bin.name, ARG) %>% unique() %>%
  group_by(bin.name) %>% summarise(numARGsubtype = n()) %>%
  filter(bin.name %in% sub("\\.filtered\\.fa", "", nrBin615.DepthPG_df.l$bin) )


# plot: the correlation of bin abundance and ARG number  -----------------------------
tmp <- Corr.res %>% filter(bin %in% sigBins) 
sum(tmp$r > 0,na.rm = T); sum(tmp$r < 0 , na.rm = T) 
sum(tmp$r > 0 & tmp$p<0.05,na.rm = T)

sigBin_r.pos <- (tmp %>% filter(r > 0))$bin
sigBin_r.neg <- (tmp %>% filter(r < 0))$bin

# compare
plotDat_corr <- Corr.res %>% filter(bin %in% sigBins)
plotDat_corr$sigBin <- sapply(plotDat_corr$bin,
                              function(x){
                                if( x %in% sigBin_r.pos) "postive.corr" else if(x %in% sigBin_r.neg) "negative.corr"
                              })
plotDat_corr$r <- as.numeric(plotDat_corr$r)
library(ggprism)
ggplot(plotDat_corr) +
  geom_point(aes(x=sigBin,y=r, color=sigBin),size=2) +
  theme_prism()

# # Figure S7A
p1 <- ggplot(NULL, aes(x = r, y=..count../sum(..count..) * 100)) +  
  geom_histogram(data = plotDat_corr %>% filter(sigBin == "postive.corr"),
                 color = "indianred3", fill = "indianred3",
                 position = "identity", bins = 20, alpha = 0.4) +
  geom_histogram(data = plotDat_corr %>% filter(sigBin == "negative.corr"),
                 color = "skyblue3", fill = "skyblue3",
                 position = "identity", bins = 20, alpha = 0.4) +
  ylab("Frequency (%)") + xlab("Correlation between MAG abundance and No. of ARG subtype in sample") +
  theme_prism()

p1

# plot: selected bins ---------------------
# compare
plotDat.numARG <-  
  numARG.bins.dat %>% filter(bin.name %in% sub("\\.filtered\\.fa","",sigBins)) %>% 
  mutate(sigBins = sapply(bin.name,
                          function(x){
                            if(x %in% sub("\\.filtered\\.fa","",sigBin_r.pos))"postive.corr" else "negative.corr"
                          }))
library(ggrepel)
ggplot(plotDat.numARG) +
  geom_boxplot(aes(x=sigBins,y=numARGsubtype),outlier.shape = NA)+
  geom_jitter(aes(x=sigBins,y=numARGsubtype), width = 0.2) +
  geom_text_repel(aes(x=sigBins,y=numARGsubtype, label=bin.name) )


# plot # Figure S7B
p2 <-ggplot(NULL, aes(x = numARGsubtype, y=..count../sum(..count..) * 100)) +  
  geom_histogram(data = plotDat.numARG %>% filter(sigBins == "postive.corr"),
                 color = "indianred3", fill = "indianred3",
                 position = "identity", bins = 20, alpha = 0.4) +
  geom_histogram(data = plotDat.numARG %>% filter(sigBins == "negative.corr"),
                 color = "skyblue3", fill = "skyblue3",
                 position = "identity", bins = 20, alpha = 0.4) +
  ylab("Frequency (%)") + xlab("No. of ARG subtype in each MAG") +
  theme_prism()

p2



taxon_simpDf <- taxon_982bin_df %>%
  dplyr::select(bin, phylum, genus, lastLevel, taxonName) %>%
  filter(bin %in% sub("\\.filtered\\.fa","", nrBin615.DepthPG_df.l$bin)) %>%
  mutate(sigBins = sapply(bin,
                          function(x){
                            if(x %in% sub("\\.filtered\\.fa","",sigBin_r.pos)){
                              "postive.corr"
                            } else if(x %in% sub("\\.filtered\\.fa","",sigBin_r.neg) ){
                              "negative.corr"
                            }else{
                              "otherBin"
                            }
                          }))


for(TaxonL in c("genus","phylum")){
  TaxonL <- "genus"
  taxon_simpDf.tmp <- taxon_simpDf
  colnames(taxon_simpDf.tmp)[colnames(taxon_simpDf.tmp) == TaxonL] <- "Taxon"
  
  table(taxon_simpDf.tmp$Taxon) 
  
  corrected_taxon <- vector("character",nrow(taxon_simpDf.tmp))
  for(i in c(1:nrow(taxon_simpDf.tmp))){
    txn <- taxon_simpDf.tmp$Taxon[i]
    
    if(txn=="") next
    parts <- strsplit(txn," ",fixed = T)[[1]]
    
    for(pt in parts){
      j=which(parts == pt)
      
      parts[j] <- sub("_[A-Z]+","",pt)
      
    }
    
    corrected_taxon[i] <- paste(parts,collapse = " ")
  }
  unique(corrected_taxon)
  taxon_simpDf.tmp$corrected_taxon <- corrected_taxon
  
  
  
  plotDat_taxon <- taxon_simpDf.tmp %>% 
    group_by(sigBins,corrected_taxon) %>%
    summarise(numBin=n()) %>% mutate(freq=numBin/sum(numBin))
  
  assign(paste("plotDat_",TaxonL,sep = ""), plotDat_taxon, envir = .GlobalEnv)
  #plotDat_phylum <- plotDat_taxon 
  plotDat_genus <- plotDat_taxon
  
}


#Figure S7C
# plot for phylum ---------------------------------
color_phylum_df <- 
  cbind.data.frame(phylum = c("Acidobacteriota", "Actinobacteriota","Dormibacterota", "Firmicutes","Nitrospirota",
                              "Planctomycetota","Proteobacteria","Thermoplasmatota","Bacteroidota","Gemmatimonadota",
                              "Verrucomicrobiota"),
                   Colors = c("#27813A","#658633","#476B83","#AFC7C7","#F0F0D9",
                              "#4891B8","#DFD180","#CA9E44","#8AB18A","#ACD4E1",
                              "#EED8AF"),
                   stringsAsFactors=F)

plotDat_taxon <-  plotDat_phylum %>% dplyr::filter(sigBins != "otherBin")


phylum_ABC <- unique(plotDat_taxon$corrected_taxon)[order(unique(plotDat_taxon$corrected_taxon))]
colors <- sapply(phylum_ABC, function(x) color_phylum_df$Colors[which(color_phylum_df$phylum == x)])
plotDat_taxon <- plotDat_taxon 
p3 <- ggplot(plotDat_taxon) +
  geom_col(aes(x=sigBins,y=freq,fill=corrected_taxon),width = 0.6 ) +  theme(legend.position = "none") +
  scale_fill_manual(values = colors)+
  theme_prism()

p3 

# Figure S7D
# plot for genus ---------------------------------
library(ggsci)
library(RColorBrewer)
nb.cols <- unique((plotDat_genus %>% filter(sigBins != "otherBin"))$corrected_taxon) %>% length()
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

p4<- ggplot(plotDat_genus %>% filter(sigBins != "otherBin")) +
  geom_col(aes(x=sigBins,y=freq,fill=corrected_taxon), width = 0.6) + theme(legend.position = "bottom") +
  
  scale_fill_manual(values = mycolors)+ 
  theme_prism()
p4 

library(ggpubr)
P <- ggarrange(ggarrange(p1,p2,ncol=1),
               p3, p4, ncol = 3, widths = c(0.2,0.35,0.43))

ggsave(P, filename = "ARG-diversity-associated-MAGs.pdf",device = "pdf",
       width = 12, height = 4)
