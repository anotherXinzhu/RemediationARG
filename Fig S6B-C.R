

library(data.table)
library(dplyr)
library(eulerr)
library(ggplot2)
library(ggvenn)
library(patchwork)
library(ggsci)
library(RColorBrewer)
library(ggpubr)
library(cowplot)

# load data-----------------------
deepargMapping_df <- fread("deeparg_gene_mapping.txt",data.table = F)
load("3.combined_ArgMrgIntTnp_allMappedSamples.RData") 
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData") 
load("2.ARGtypes.RData")

# integrate data------------------
regGenes <- unique(regGenes.v2_df$ARGsubtype) 
combined_deepargByGene_df <- combined_deepargByGene_df %>% filter(!ARG %in% regGenes)

dat <- combined_deepargByGene_df %>% dplyr::select(sampleID, ARG, DepthPG, drug_type, res_mechanism, sampleType, remediateYear)

totalARGDepthPG_Df <- 
  combined_deepargByGene_df %>% 
  group_by(sampleID) %>%
  summarise(totalARGDepthPG = sum(DepthPG)) %>% as.data.frame()

dat <- merge(dat, totalARGDepthPG_Df, by="sampleID" ) %>% dplyr::mutate(relAbund = DepthPG/totalARGDepthPG)

dat_original <- dat 

tmp <- dat %>% filter(drug_type != "unclassified") %>% select(-totalARGDepthPG, -relAbund) 
dat <- merge(tmp,
             tmp %>% group_by(sampleID) %>% summarise(totalARGDepthPG = sum(DepthPG)),
             by="sampleID") %>%
  mutate(relAbund = DepthPG/totalARGDepthPG)

AvgRelAbund_df <- dat %>% group_by(sampleType, ARG) %>% summarise(avgRelAbund = mean(relAbund)) %>% as.data.frame() %>%
  arrange(desc(avgRelAbund)) %>% arrange(sampleType)

topN <- 50
topARG_contaminated <- (AvgRelAbund_df %>% filter(sampleType == "contaminated") )$ARG[1:topN]
topARG_remediate_inner <- (AvgRelAbund_df %>% filter(sampleType == "remediate_inner") )$ARG[1:topN]
topARG_remediate_top <- (AvgRelAbund_df %>% filter(sampleType == "remediate_top") )$ARG[1:topN]

minRA <- 0.01
topARG_contaminated <- (AvgRelAbund_df %>% filter(sampleType == "contaminated") %>% filter(avgRelAbund > minRA))$ARG 
topARG_remediate_inner <- (AvgRelAbund_df %>% filter(sampleType == "remediate_inner") %>% filter(avgRelAbund > minRA))$ARG
topARG_remediate_top <- (AvgRelAbund_df %>% filter(sampleType == "remediate_top") %>% filter(avgRelAbund > minRA))$ARG


# venn--------------------------------------------------------------
# plot
# Venn_No.abundant ARGs
tmp <- merge(cbind.data.frame(ARG=topARG_contaminated,
                              Contaminated = T,
                              stringsAsFactors=F), 
             cbind.data.frame(ARG=topARG_remediate_inner,
                              Remediate_inner=T,
                              stringsAsFactors=F),
             by="ARG", all=T) 
VennPlotDat <- merge(tmp, 
                     cbind.data.frame(ARG=topARG_remediate_top,
                                      Remediate_top=T,
                                      stringsAsFactors=F),
                     by="ARG", all=T)

VennPlotDat[is.na(VennPlotDat)] <- F
fit <- euler(VennPlotDat[-1], shape = "ellipse")

F1 <- ggvenn::ggvenn(VennPlotDat, fill_color = c("#3CA183", "#9BC4DA","#D7A641","red"), 
                     fill_alpha = 0.3,
                     stroke_color="black",
                     stroke_alpha = 0.5,
                     stroke_size = 0.2,
                     text_size = 3.5) + 
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(plot.subtitle = element_text(size = 20, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 14)) + theme(plot.subtitle = element_text(size = 10, colour = "black"))
F1 


# Venn_No.ARGs
ARG_contaminated <-unique(( dat %>% filter(sampleType == "contaminated") %>% filter(DepthPG >0) )$ARG)
ARG_remediate_inner <- unique((dat %>% filter(sampleType == "remediate_inner")%>% filter(DepthPG >0))$ARG)
ARG_remediate_top <- unique((dat %>% filter(sampleType == "remediate_top") %>% filter(DepthPG >0))$ARG)

tmp <- merge(cbind.data.frame(ARG=ARG_contaminated,
                              Contaminated = T,
                              stringsAsFactors=F), 
             cbind.data.frame(ARG=ARG_remediate_inner,
                              Remediate_inner=T,
                              stringsAsFactors=F),
             by="ARG",all=T)

VennPlotDat_allARG <- merge(tmp, 
                            cbind.data.frame(ARG=ARG_remediate_top,
                                             Remediate_top=T,
                                             stringsAsFactors=F),
                            by="ARG", all=T)
VennPlotDat_allARG[is.na(VennPlotDat_allARG)] <- F


fit <- euler(VennPlotDat_allARG[-1], shape = "ellipse")

F2 <- ggvenn::ggvenn(VennPlotDat_allARG, fill_color = c("#3CA183", "#9BC4DA","#D7A641","red"), 
                     fill_alpha = 0.3,
                     stroke_color="black",
                     stroke_alpha = 0.5,
                     stroke_size = 0.2,
                     text_size = 3.5) + 
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(plot.subtitle = element_text(size = 20, face = "bold"), 
  legend.text = element_text(size = 12, face = "bold"), 
  legend.title = element_text(size = 14)) + 
  theme(plot.subtitle = element_text(size = 10, colour = "black"))
F2 

# ARG Types -----------------------------------------------------

color_ARGtypes_df <- cbind.data.frame(ARGtype_rank,
                                      color = c("#DC0000FF", "#E64B35FF",   "#00A087FF", "#3C5488FF", "#B09C85FF","#F39B7FFF",
                                                "#4DBBD5FF","#8491B4FF","#91D1C2FF","#7E6148FF","gray"),
                                      stringsAsFactors=F)
# abbreviations for ARG subtype in x axis-------------------------------
ARGsubtypes <- unique(c(topARG_contaminated, topARG_remediate_inner, topARG_remediate_top))
abb  <- unname(sapply(ARGsubtypes, function(x) {
  splited <- strsplit(x,"_",fixed = T)[[1]]
  Abb <- splited[which(nchar(splited) == min(nchar(splited)))]
  
  if(length(Abb) > 1) Abb <- Abb[1]
  return(Abb)
}) )

ARGsubtype_abb_df <- cbind.data.frame(ARGsubtypes,abb,stringsAsFactors=F)


for(i in c(1:length(ARGsubtypes))){
  if(ARGsubtype_abb_df$abb[i] == "to"){
    ARGsubtype_abb_df$abb[i] <- "mupA"
  }else if(ARGsubtypes[i] == "EmrB-QacA_family_major_facilitator_transporter"){
    ARGsubtype_abb_df$abb[i] <- "comD"  # reference: https://www.ncbi.nlm.nih.gov/gene/8624036
  }else if(ARGsubtypes[i] == "major_facilitator_superfamily_transporter"){
    ARGsubtype_abb_df$abb[i] <- "MFS"
  }else if(ARGsubtypes[i] == "cob(I)alamin_adenolsyltransferase"){
    ARGsubtype_abb_df$abb[i] <- "cobI_ATR"
  }else if(ARGsubtype_abb_df$abb[i] == "or"){
    ARGsubtype_abb_df$abb[i] <- "rpsD"
  }
}


# UT ------------------------
dat_contaminated <- 
  dat %>% 
  filter(sampleType == "contaminated") %>% 
  filter(ARG %in% topARG_contaminated) 

dat_contaminated_sd <- dat_contaminated %>%
  group_by(ARG) %>% 
  summarise(MRA=mean(relAbund), sd=sd(relAbund), resMech=unique(res_mechanism), drugType=unique(drug_type)) %>%
  arrange(desc(MRA))

dat_contaminated$ARG <- factor(dat_contaminated$ARG, levels = dat_contaminated_sd$ARG)
dat_contaminated_sd$ARG <- factor(dat_contaminated_sd$ARG, levels = dat_contaminated_sd$ARG)

table(dat_contaminated_sd$drugType) 

# plot

p0_1 <- ggplot(dat_contaminated) + 
  geom_point(aes( x = ARG, y = relAbund),shape=1,alpha=0.7,color="darkgray", size = 4) +
  theme_bw() +
  labs(y = "Relative abundance", x = "", title = "") +
  scale_y_continuous(limits = c(-0.01,0.18),oob=rescale_none) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
  scale_x_discrete(labels=sapply(levels(dat_contaminated$ARG), function(x) ARGsubtype_abb_df$abb[which(ARGsubtype_abb_df$ARGsubtypes == x)] )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15, color = "black")) + 
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.text = element_text(size = 15)) 

p0_1

argTypes <- unique(dat_contaminated_sd$drugType)[order(unique(dat_contaminated_sd$drugType))]
argTypes[which(argTypes %in% otherTypes)] <- "other"
Colors= sapply(argTypes,function(x) color_ARGtypes_df$color[color_ARGtypes_df$ARGtype_rank==x])
p0_2 <- ggplot(dat_contaminated_sd) + 
  geom_pointrange(aes(data= ,x = ARG, y = MRA,ymin=MRA-sd, ymax=MRA+sd, color=drugType),size=1.5) +
  scale_color_manual(values = unname(Colors) ) +  
  scale_y_continuous(limits = c(-0.01,0.2),oob=rescale_none) +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("y.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("y.ticks") +
  rremove("y.axis") + rremove("y.text") + rremove("ylab")
p0_2

aligned_plots <- align_plots(p0_1, p0_2+theme(legend.position = "none"), align="hv", axis="tblr") 
p_contaminated <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_contaminated


# ULRT ------------------------
dat_remediate_inner <- 
  dat %>% 
  filter(sampleType == "remediate_inner") %>% 
  filter(ARG %in% topARG_remediate_inner) 

dat_remediate_inner_sd <- dat_remediate_inner %>%
  group_by(ARG) %>% 
  summarise(MRA=mean(relAbund), sd=sd(relAbund), resMech=unique(res_mechanism), drugType=unique(drug_type)) %>%
  arrange(desc(MRA))

dat_remediate_inner$ARG <- factor(dat_remediate_inner$ARG, levels = dat_remediate_inner_sd$ARG)
dat_remediate_inner_sd$ARG <- factor(dat_remediate_inner_sd$ARG, levels = dat_remediate_inner_sd$ARG)

table(dat_remediate_inner_sd$drugType) 
write.csv(dat_remediate_inner, file = "dat_remediate_inner.csv")

p0_1 <- ggplot(dat_remediate_inner) + 
  geom_point(aes( x = ARG, y = relAbund),shape=1,alpha=0.7,color="darkgray", size = 4) +
  theme_bw() +
  labs(y = "Relative abundance", x = "", title = "") +
  scale_y_continuous(limits = c(-0.01,0.15),oob=rescale_none) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) + 
  scale_x_discrete(labels=sapply(levels(dat_contaminated$ARG), function(x) ARGsubtype_abb_df$abb[which(ARGsubtype_abb_df$ARGsubtypes == x)] )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15, color = "black")) + 
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.text = element_text(size = 15)) 


p0_1

argTypes <- unique(dat_remediate_inner_sd$drugType)[order(unique(dat_remediate_inner_sd$drugType))]
argTypes[which(argTypes %in% otherTypes)] <- "other"
Colors= sapply(argTypes,function(x) color_ARGtypes_df$color[color_ARGtypes_df$ARGtype_rank==x])
p0_2 <- ggplot(dat_remediate_inner_sd) + 
  geom_pointrange(aes(data= ,x = ARG, y = MRA,ymin=MRA-sd, ymax=MRA+sd, color=drugType),size=1.5) +
  scale_color_manual(values = unname(Colors) ) + 
  scale_y_continuous(limits = c(-0.01,0.2),oob=rescale_none) +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("y.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("y.ticks") +
  rremove("y.axis") + rremove("y.text") + rremove("ylab")
p0_2
p0_2

ggsave(p0_2, filename = "1.plots_AbundantARG/legend.pdf", device = "pdf")

aligned_plots <- align_plots(p0_1, p0_2+theme(legend.position = "none"), align="hv", axis="tblr")
p_remediate_inner <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_remediate_inner

# ALRT ------------------------
dat_remediate_top <- 
  dat %>% 
  filter(sampleType == "remediate_top") %>% 
  filter(ARG %in% topARG_remediate_top) 

dat_remediate_top_sd <- dat_remediate_top %>%
  group_by(ARG) %>% 
  summarise(MRA=mean(relAbund), sd=sd(relAbund), resMech=unique(res_mechanism), drugType=unique(drug_type)) %>%
  arrange(desc(MRA))

dat_remediate_top$ARG <- factor(dat_remediate_top$ARG, levels = dat_remediate_top_sd$ARG)
dat_remediate_top_sd$ARG <- factor(dat_remediate_top_sd$ARG, levels = dat_remediate_top_sd$ARG)

table(dat_remediate_top_sd$drugType) 



p0_1 <- ggplot(dat_remediate_top) + 
  geom_point(aes( x = ARG, y = relAbund),shape=1,alpha=0.7,color="darkgray", size = 4) +
  theme_bw() +
  labs(y = "Relative abundance", x = "", title = "") +
  scale_y_continuous(limits = c(-0.01,0.15),oob=rescale_none) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) + 
  scale_x_discrete(labels=sapply(levels(dat_contaminated$ARG), function(x) ARGsubtype_abb_df$abb[which(ARGsubtype_abb_df$ARGsubtypes == x)] )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15, color = "black")) + 
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.text = element_text(size = 15)) 

p0_1

argTypes <- unique(dat_remediate_top_sd$drugType)[order(unique(dat_remediate_top_sd$drugType))]
argTypes[which(argTypes %in% otherTypes)] <- "other"
Colors= sapply(argTypes,function(x) color_ARGtypes_df$color[color_ARGtypes_df$ARGtype_rank==x])
p0_2 <- ggplot(dat_remediate_top_sd) + 
  geom_pointrange(aes(data= ,x = ARG, y = MRA,ymin=MRA-sd, ymax=MRA+sd, color=drugType),size=1) +
  scale_color_manual(values = unname(Colors) ) + 
  scale_y_continuous(limits = c(-0.01,0.2),oob=rescale_none) +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("y.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("y.ticks") +
  rremove("y.axis") + rremove("y.text") + rremove("ylab")
p0_2


aligned_plots <- align_plots(p0_1, p0_2+theme(legend.position = "none"), align="hv", axis="tblr")
p_remediate_top <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_remediate_top

ggsave(ggarrange(p_contaminated,p_remediate_inner,p_remediate_top,nrow = 1,ncol = 3),
       device = "pdf", width = 15,height = 4, filename = "plots_AbundantARG_Fig.pdf")
