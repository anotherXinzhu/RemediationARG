

library(dplyr)
library(reshape2)
library(caret)
library(mediation)
library(vegan)
library(ggplot2)
library(tidyverse)

# load data---------------------------------------
load("3.combined_ArgMrgIntTnp_allMappedSamples.RData")
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")

phychem_dat <- read.csv("phychem_model2_dat.csv",row.names =1)
pollutionIndex_df <- read.csv("pollutionIndex_model1_df.csv", row.names =1)

diversity <- read.csv("16S_all_diversity.csv",header = T, sep = ",")
richness <- diversity[,c(1,2)]
# integrate data
regGenes <- unique(regGenes.v2_df$ARGsubtype) 

numARG.dat <-combined_deepargByGene_df %>%  
  filter(!ARG %in% regGenes) %>%
  filter(Coverage > 0) %>% 
  group_by(sampleID) %>%
  summarise(numARG = n()) 

numARG.dat$sampleID <- toupper(sub("Q1", "",gsub( "[_\\-]", "", numARG.dat$sampleID) ) )

# selcted ARG-rich MAGs with abundance

if(F){
  
  MAGs.dat <- fread("corelated_bin.csv")
  MAGs.dat <- MAGs.dat %>% dplyr::mutate(bin = sub("\\.filtered\\.fa","",sample)) %>% 
    dplyr::select(-sample) %>%
    reshape2::melt(variable.name="sample", value.name = "DepthPG") %>%
    group_by(sample) %>% summarise(ArgRichMAG.abund=sum(DepthPG))

}


if(T){
  
  tmp <- fread("corelated_bin.csv") %>% column_to_rownames("sample") %>% t() %>% data.frame()
  
  pca.mag <- prcomp(tmp, scale. = T) 
  MAGs.dat <- cbind.data.frame(sample= toupper(rownames(tmp)) , ArgRichMAG.abund = pca.mag$x[,1])
  
}


if(T){
  tmp1 <- merge(merge(numARG.dat, richness, by = "sampleID"), MAGs.dat, by.x = "sampleID", by.y = "sample") 
  tmp2 <- merge(phychem_dat,  pollutionIndex_df,by=0)
  dat <- merge(tmp1, tmp2 , 
               by.x = "sampleID", by.y = "Row.names")
  
}

# 2. train and test data set------------------------------------------------
#define training and testing sets
train.dat <- dat[ c(seq(1,16,3),seq(2,17,3)), -1] 
test.respon <- dat[seq(3,18,3), c("numARG")]
test.pred <- dat[seq(3,18,3), -(1:2)] 


#
set.seed(1)
model3 <- train(
  numARG~., 
  data = train.dat,
  method = 'pls',
  preProcess = c("center","scale") 
)
model3 


predictions = predict(model3, newdata = test.pred)

# RMSE
sqrt(mean((test.respon - predictions)^2))  

# R2
cor(test.respon, predictions) ^ 2  


# plot：variable importance --------------

model <- model3
# importance
VarImportance <- varImp(model)$importance  

# direction: negative / positive influence 
directions <- sapply(rownames(VarImportance), 
                     function(x){
                       ct = cor.test(dat[,"numARG"], dat[,x]) 
                       r = ct$estimate
                       direction = if(r > 0) "positive" else "negative"
                     })

# co-linear group 
library(Hmisc)
load("3_0_3.co-linear-envs_w_pollutIndex_16S.MAGs.RData")

# group
Colinear.Grps <- sapply(rownames(VarImportance),
                        function(x){
                          grp = names(which(sapply(Co.linear.Var_list,function(L) x %in% L)))
                          if(length(grp) == 0) NA else grp
                        })

# directions, VarImportance, Colinear.Grps
plotDat <- cbind.data.frame(VarImportance,
                            directions,
                            Colinear.Grps) 
plotDat$Colinear.Grps[is.na(plotDat$Colinear.Grps)] <- "NoGroup"
plotDat$Colinear.Grps[1] <- "CoLinear.Grp1"

# compute average importance
Colinear.Grps_rank <- 
  plotDat %>% group_by(Colinear.Grps) %>% 
  summarise(avg.importance = mean(Overall)) %>% 
  arrange(desc(avg.importance))

plotDat$Colinear.Grps <- 
  factor(plotDat$Colinear.Grps, levels = Colinear.Grps_rank$Colinear.Grps) 

plotDat <- plotDat %>% 
  arrange(desc(Overall )) %>%  arrange(Colinear.Grps) %>% 
  tibble::rownames_to_column("Vars")

plotDat$Vars <- factor(plotDat$Vars, levels = plotDat$Vars)


#ggplot

ggplot(plotDat) +
  geom_col(aes(x=Vars, y=Overall, fill=directions), width = 0.75)+
  facet_grid(.~Colinear.Grps,scales = "free_x",space = "free_x") +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = c("#7179AB","#E44B4B")) +
  theme_bw(base_size = 15) + 
  theme(plot.subtitle = element_text(size = 20,
                                     colour = "black"), plot.title = element_text(hjust = 0.5)) +labs(x = "Factors", y = "Importance") + 
  theme(panel.grid.major = element_line(colour = "ghostwhite"),
        panel.grid.minor = element_line(colour = "ghostwhite"),
        axis.text.x = element_text(colour = "black",vjust = 1,hjust = 1,angle = 90),
        axis.text.y = element_text(colour = "black"),
        panel.background = element_rect(fill = "white",
                                        colour = "white", size = 0, linetype = "solid"),
        plot.background = element_rect(colour = NA),legend.position = 'top') +
  theme(strip.background = element_blank())


#  mediation ---------------------------------------------------------------------------------------------------------------------------
library(vegan)
dat.st <- decostand(dat[,-1], method = "standardize", MARGIN = 2) 
sapply(dat.st, mean); sapply(dat.st,sd) 
# pick most important variables
VarImportance <- VarImportance %>% arrange(desc(Overall)) 

Mediation_res <- NULL 

tmp <- rownames(VarImportance)[1:12]
topN.phych <- tmp[!tmp %in% c("richness","ArgRichMAG.abund")]
pca.ml <- prcomp(dat %>% dplyr::select(all_of(topN.phych)), scale. = T) #

dat.mediation <- cbind.data.frame(dat.st[,c("numARG","richness","ArgRichMAG.abund")],
                                  phychem.PC1 = (pca.ml$x[,1]),
                                  stringsAsFactors=F)


Y = "numARG"
trtr = "phychem.PC1"


for(mdtr in c("richness","ArgRichMAG.abund")){
  id = paste(trtr, mdtr, Y, sep = "_")
  dat.tmp <- dat.mediation
  colnames(dat.tmp)[colnames(dat.tmp) == Y] <- "Y"
  colnames(dat.tmp)[colnames(dat.tmp) == trtr] <- "Treat"
  colnames(dat.tmp)[colnames(dat.tmp) == mdtr] <- "Mediator"
  
  
  model.m <- try(lm(Mediator ~ Treat , data = dat.tmp )) 
  model.y <- try(lm(Y ~ Treat + Mediator, data = dat.tmp))
  
  if(class(model.m) != 'try-error' & class(model.y) != 'try-error' ){
    testModel  <- try(mediation::mediate(model.m,model.y,treat="Treat",mediator="Mediator",boot=F,sims=1000), silent = T)
    
    if(class(testModel) != 'try-error'  ){
      summary = summary(testModel)
      
      res <- capture.output(summary,append=FALSE)
      
      tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
      tmp <- tmp[tmp != "" & tmp!="."]
      tmp <- tmp[!grepl("*",tmp,fixed = T)] 
      ACME.p <- tmp[length(tmp)]
      
      tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
      tmp <- tmp[tmp != "" & tmp!="."]
      tmp <- tmp[!grepl("*",tmp,fixed = T) ]
      ADE.p <- tmp[length(tmp)]
      
      tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
      tmp <- tmp[tmp != "" ]
      i_str = which(grepl("Mediated", tmp))
      prop.mediated <- tmp[(i_str + 1)]
      
      spearman.r = cor(dat.tmp$Treat, dat.tmp$Mediator, method = "spearman")
      
      
      
      vec = c(id, N, ACME.p, ADE.p, prop.mediated, spearman.r)
      names(vec) <- c("Treat_Mediator_Y", "topN.phychem","ACME.p", "ADE.p", "prop.mediated", "spearman.r")
    }else{
      vec = c(id, N,NA,NA,NA,NA)
      names(vec) <- c("Treat_Mediator_Y", "topN.phychem", "ACME.p", "ADE.p", "prop.mediated", "spearman.r")
    }
  }else{
    vec = c(id, N,NA,NA,NA,NA)
    names(vec) <- c("Treat_Mediator_Y", "topN.phychem", "ACME.p", "ADE.p", "prop.mediated", "spearman.r")
  }
  
  Mediation_res <- bind_rows(Mediation_res, vec)
}


names(Mediation_res) <- seq(1,length(Mediation_res),1)
Mediation_res <- data.frame(Mediation_res) %>% t() %>% data.frame()

write.csv(Mediation_res, file = "Mediation.res_richness_MAG_ARG number.csv")

