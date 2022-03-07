
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)

# load data
load("repaire_mapping_df1.Rdata")
load("6_1.regGenes.nonRegGenes_DeepARG.v1.v2.RData")

repaire_mapping_df2 <- repaire_mapping_df1 
read_id <- rownames(repaire_mapping_df1)
rownames(repaire_mapping_df2) <- NULL 
repaire_mapping_df2 <- cbind(read_id, repaire_mapping_df2)

combined <- repaire_mapping_df2 %>% mutate(read_id = sapply(strsplit(read_id, " ", fixed = T), "[[", 1))
duplicated(combined$read_id) 

which(duplicated(combined$read_id))
combined[which(duplicated(combined$read_id)),]

nodupicated <- combined %>%
  dplyr::filter(!grepl('k',read_id)) 

duplicated <- combined %>%
  dplyr::filter(grepl('k',read_id)) 

#TUT1
duplicated1 <- head(combined %>%
                      dplyr::filter(grepl('k',read_id)),37115) 
tail(duplicated1)

duplicated1$read_id1 <- c('_tut1') 
duplicated1$read_id <- str_c(duplicated1$read_id, duplicated1$read_id1)
duplicated1 <- duplicated1[,-20] 

#TUT21
duplicated21 <- (duplicated[37116:82175,])
duplicated21$read_id1 <- c('_tut21')
duplicated21$read_id <- str_c(duplicated21$read_id, duplicated21$read_id1)
duplicated21 <- duplicated21[,-20]

#TUT4
duplicated4 <- tail(combined %>%  
                      dplyr::filter(grepl('k',read_id)),36186)
tail(duplicated4)

duplicated4$read_id1 <- c('_tut4') 
duplicated4$read_id <- str_c(duplicated4$read_id, duplicated4$read_id1)
duplicated4 <- duplicated4[,-20] 



# recombined
all_repaire_mapping <-rbind(nodupicated, duplicated1, duplicated21, duplicated4)

duplicated(all_repaire_mapping$read_id) 


regGenes <- unique(regGenes.v2_df$ARGsubtype)
total_ALRT2_plasmid_ARG_nonReg <- ALRT2_plasmid_ARG %>% filter(!ARG %in% regGenes)
# 
#ALRT2
ALRT2 <- read.csv("ALRT2.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
ALRT2$ARG <- sapply(strsplit(ALRT2$'best.hit', "|", fixed = T),"[", 3) 
ALRT2 <- ALRT2 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
ALRT2 <- unite(ALRT2, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
ALRT2_df <- ALRT2 %>% dplyr::select(read_id, ARG)
ALRT2_plasmid_ARG <- merge(all_repaire_mapping[,c(1,2)], ALRT2_df, by = "read_id")


names(total_ALRT2_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_ALRT2_plasmid_ARG_nonReg$samplename <- rep(c("ALRT2"), each = 174) 
total_ALRT2_plasmid_ARG_nonReg$sampletype <- rep(c("ALRT"), each = 174)

#ALRT8
ALRT8 <- read.csv("ALRT8.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
ALRT8$ARG <- sapply(strsplit(ALRT8$'best.hit', "|", fixed = T),"[", 3) 
ALRT8 <- ALRT8 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
ALRT8 <- unite(ALRT8, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
ALRT8_df <- ALRT8 %>% dplyr::select(read_id, ARG)
ALRT8_plasmid_ARG <- merge(all_repaire_mapping[,c(1,3)], ALRT8_df, by = "read_id")

total_ALRT8_plasmid_ARG_nonReg <- ALRT8_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_ALRT8_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_ALRT8_plasmid_ARG_nonReg$samplename <- rep(c("ALRT8"), each = 212)
total_ALRT8_plasmid_ARG_nonReg$sampletype <- rep(c("ALRT"), each = 212)

#ALRT9
ALRT9 <- read.csv("ALRT9.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
ALRT9$ARG <- sapply(strsplit(ALRT9$'best.hit', "|", fixed = T),"[", 3) 
ALRT9 <- ALRT9 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
ALRT9 <- unite(ALRT9, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
ALRT9_df <- ALRT9 %>% dplyr::select(read_id, ARG)
ALRT9_plasmid_ARG <- merge(all_repaire_mapping[,c(1,4)], ALRT9_df, by = "read_id")

total_ALRT9_plasmid_ARG_nonReg <- ALRT9_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_ALRT9_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_ALRT9_plasmid_ARG_nonReg$samplename <- rep(c("ALRT9"), each = 600)
total_ALRT9_plasmid_ARG_nonReg$sampletype <- rep(c("ALRT"), each = 600)

#TAS2
TAS2 <- read.csv("TAS2.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
TAS2$ARG <- sapply(strsplit(TAS2$'best.hit', "|", fixed = T),"[", 3) 
TAS2 <- TAS2 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
TAS2 <- unite(TAS2, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
TAS2_df <- TAS2 %>% dplyr::select(read_id, ARG)
TAS2_plasmid_ARG <- merge(all_repaire_mapping[,c(1,5)], TAS2_df, by = "read_id")

total_TAS2_plasmid_ARG_nonReg <- TAS2_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_TAS2_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_TAS2_plasmid_ARG_nonReg$samplename <- rep(c("TAS2"), each = 253)
total_TAS2_plasmid_ARG_nonReg$sampletype <- rep(c("ALRT"), each = 253)

#TAS4
TAS4 <- read.csv("TAS4.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
TAS4$ARG <- sapply(strsplit(TAS4$'best.hit', "|", fixed = T),"[", 3) 
TAS4 <- TAS4 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
TAS4 <- unite(TAS4, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
TAS4_df <- TAS4 %>% dplyr::select(read_id, ARG)
TAS4_plasmid_ARG <- merge(all_repaire_mapping[,c(1,6)], TAS4_df, by = "read_id")

total_TAS4_plasmid_ARG_nonReg <- TAS4_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_TAS4_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_TAS4_plasmid_ARG_nonReg$samplename <- rep(c("TAS4"), each = 175)
total_TAS4_plasmid_ARG_nonReg$sampletype <- rep(c("ALRT"), each = 175)


#TAS6
TAS6 <- read.csv("TAS6.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
TAS6$ARG <- sapply(strsplit(TAS6$'best.hit', "|", fixed = T),"[", 3) 
TAS6 <- TAS6 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
TAS6 <- unite(TAS6, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
TAS6_df <- TAS6 %>% dplyr::select(read_id, ARG)
TAS6_plasmid_ARG <- merge(all_repaire_mapping[,c(1,7)], TAS6_df, by = "read_id")

total_TAS6_plasmid_ARG_nonReg <- TAS6_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_TAS6_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_TAS6_plasmid_ARG_nonReg$samplename <- rep(c("TAS6"), each = 208)
total_TAS6_plasmid_ARG_nonReg$sampletype <- rep(c("ALRT"), each = 208)


#TAX2
TAX2 <- read.csv("TAX2.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
TAX2$ARG <- sapply(strsplit(TAX2$'best.hit', "|", fixed = T),"[", 3) 
TAX2 <- TAX2 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
TAX2 <- unite(TAX2, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
TAX2_df <- TAX2 %>% dplyr::select(read_id, ARG)
TAX2_plasmid_ARG <- merge(all_repaire_mapping[,c(1,8)], TAX2_df, by = "read_id")

total_TAX2_plasmid_ARG_nonReg <- TAX2_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_TAX2_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_TAX2_plasmid_ARG_nonReg$samplename <- rep(c("TAX2"), each = 133)
total_TAX2_plasmid_ARG_nonReg$sampletype <- rep(c("ULRT"), each = 133)

#TAX4
TAX4 <- read.csv("TAX4.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
TAX4$ARG <- sapply(strsplit(TAX4$'best.hit', "|", fixed = T),"[", 3) 
TAX4 <- TAX4 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
TAX4 <- unite(TAX4, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
TAX4_df <- TAX4 %>% dplyr::select(read_id, ARG)
TAX4_plasmid_ARG <- merge(all_repaire_mapping[,c(1,9)], TAX4_df, by = "read_id")

total_TAX4_plasmid_ARG_nonReg <- TAX4_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_TAX4_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_TAX4_plasmid_ARG_nonReg$samplename <- rep(c("TAX4"), each = 367)
total_TAX4_plasmid_ARG_nonReg$sampletype <- rep(c("ULRT"), each = 367)


#TAX6
TAX6 <- read.csv("TAX6.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
TAX6$ARG <- sapply(strsplit(TAX6$'best.hit', "|", fixed = T),"[", 3) 
TAX6 <- TAX6 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
TAX6 <- unite(TAX6, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
TAX6_df <- TAX6 %>% dplyr::select(read_id, ARG)
TAX6_plasmid_ARG <- merge(all_repaire_mapping[,c(1,10)], TAX6_df, by = "read_id")

total_TAX6_plasmid_ARG_nonReg <- TAX6_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_TAX6_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_TAX6_plasmid_ARG_nonReg$samplename <- rep(c("TAX6"), each = 378)
total_TAX6_plasmid_ARG_nonReg$sampletype <- rep(c("ULRT"), each = 378)

#ULRT2
ULRT2 <- read.csv("ULRT2.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
ULRT2$ARG <- sapply(strsplit(ULRT2$'best.hit', "|", fixed = T),"[", 3)
ULRT2 <- ULRT2 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
ULRT2 <- unite(ULRT2, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
ULRT2_df <- ULRT2 %>% dplyr::select(read_id, ARG)
ULRT2_plasmid_ARG <- merge(all_repaire_mapping[,c(1,14)], ULRT2_df, by = "read_id")

total_ULRT2_plasmid_ARG_nonReg <- ULRT2_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_ULRT2_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_ULRT2_plasmid_ARG_nonReg$samplename <- rep(c("ULRT2"), each = 390)
total_ULRT2_plasmid_ARG_nonReg$sampletype <- rep(c("ULRT"), each = 390)

#ULRT8
ULRT8 <- read.csv("ULRT8.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
ULRT8$ARG <- sapply(strsplit(ULRT8$'best.hit', "|", fixed = T),"[", 3) 
ULRT8 <- ULRT8 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
ULRT8 <- unite(ULRT8, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
ULRT8_df <- ULRT8 %>% dplyr::select(read_id, ARG)
ULRT8_plasmid_ARG <- merge(all_repaire_mapping[,c(1,15)], ULRT8_df, by = "read_id")

total_ULRT8_plasmid_ARG_nonReg <- ULRT8_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_ULRT8_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_ULRT8_plasmid_ARG_nonReg$samplename <- rep(c("ULRT8"), each = 363)
total_ULRT8_plasmid_ARG_nonReg$sampletype <- rep(c("ULRT"), each = 363)


#ULRT9
ULRT9 <- read.csv("ULRT9.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
ULRT9$ARG <- sapply(strsplit(ULRT9$'best.hit', "|", fixed = T),"[", 3) 
ULRT9 <- ULRT9 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
ULRT9 <- unite(ULRT9, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
ULRT9_df <- ULRT9 %>% dplyr::select(read_id, ARG)
ULRT9_plasmid_ARG <- merge(all_repaire_mapping[,c(1,16)], ULRT9_df, by = "read_id")

total_ULRT9_plasmid_ARG_nonReg <- ULRT9_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_ULRT9_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_ULRT9_plasmid_ARG_nonReg$samplename <- rep(c("ULRT9"), each = 347)
total_ULRT9_plasmid_ARG_nonReg$sampletype <- rep(c("ULRT"), each = 347)


#UT7
UT7 <- read.csv("UT7.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
UT7$ARG <- sapply(strsplit(UT7$'best.hit', "|", fixed = T),"[", 3) 
UT7 <- UT7 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
UT7 <- unite(UT7, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
UT7_df <- UT7 %>% dplyr::select(read_id, ARG)
UT7_plasmid_ARG <- merge(all_repaire_mapping[,c(1,17)], UT7_df, by = "read_id")

total_UT7_plasmid_ARG_nonReg <- UT7_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_UT7_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_UT7_plasmid_ARG_nonReg$samplename <- rep(c("UT7"), each = 68)
total_UT7_plasmid_ARG_nonReg$sampletype <- rep(c("UT"), each = 68)

#UT8
UT8 <- read.csv("UT8.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
UT8$ARG <- sapply(strsplit(UT8$'best.hit', "|", fixed = T),"[", 3)
UT8 <- UT8 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
UT8 <- unite(UT8, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
UT8_df <- UT8 %>% dplyr::select(read_id, ARG)
UT8_plasmid_ARG <- merge(all_repaire_mapping[,c(1,18)], UT8_df, by = "read_id")

total_UT8_plasmid_ARG_nonReg <- UT8_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_UT8_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_UT8_plasmid_ARG_nonReg$samplename <- rep(c("UT8"), each = 71)
total_UT8_plasmid_ARG_nonReg$sampletype <- rep(c("UT"), each = 71)

#UT9
UT9 <- read.csv("UT9.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t")
UT9$ARG <- sapply(strsplit(UT9$'best.hit', "|", fixed = T),"[", 3) 
UT9 <- UT9 %>% mutate(read_id1 = sapply(strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_id11 = sapply(strsplit(read_id,"_", fixed = T),"[[",2)) %>% mutate(read_id111 = sapply(strsplit(read_id,"_", fixed = T),"[[",3)) %>% mutate(read_id1111 = sapply(strsplit(read_id,"_", fixed = T),"[[",4)) %>% mutate(read_id11111 = sapply(strsplit(read_id,"_", fixed = T),"[[",5)) %>% mutate(read_id111111 = sapply(strsplit(read_id,"_", fixed = T),"[[",6)) 
UT9 <- unite(UT9, "read_id", c("read_id1", "read_id11", "read_id111", "read_id1111", "read_id11111", "read_id111111"), sep = "_", remove = T)
UT9_df <- UT9 %>% dplyr::select(read_id, ARG)
UT9_plasmid_ARG <- merge(all_repaire_mapping[,c(1,19)], UT9_df, by = "read_id")

total_UT9_plasmid_ARG_nonReg <- UT9_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_UT9_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_UT9_plasmid_ARG_nonReg$samplename <- rep(c("UT9"), each = 72)
total_UT9_plasmid_ARG_nonReg$sampletype <- rep(c("UT"), each = 72)

#TUT1
TUT1 <- read.csv("TUT1.contig.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t") 
TUT1$ARG <- sapply(strsplit(TUT1$'best.hit', "|", fixed = T),"[", 3)
TUT1_ <- TUT1 %>% mutate(read_idtut1 = sapply( strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_idtut11 = sapply( strsplit(read_id,"_", fixed = T),"[[",2))
TUT1_$read_idtut111 <- c('tut1')
TUT1_ <- TUT1_[,-4] 
TUT1_ <- unite(TUT1_, "read_id", c("read_idtut1", "read_idtut11", "read_idtut111"), sep = "_", remove = T) 
TUT1_df <- TUT1_ %>% dplyr::select(read_id, ARG)
TUT1_plasmid_ARG <- merge(all_repaire_mapping[,c(1,11)], TUT1_df, by = "read_id")

total_TUT1_plasmid_ARG_nonReg <- TUT1_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_TUT1_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_TUT1_plasmid_ARG_nonReg$samplename <- rep(c("TUT1"), each = 164)
total_TUT1_plasmid_ARG_nonReg$sampletype <- rep(c("UT"), each = 164)

#TUT21
TUT21 <- read.csv("TUT21.contig.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t") 
TUT21$ARG <- sapply(strsplit(TUT21$'best.hit', "|", fixed = T),"[", 3)
TUT21_ <- TUT21 %>% mutate(read_idtut21 = sapply( strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_idtut211 = sapply( strsplit(read_id,"_", fixed = T),"[[",2))
TUT21_$read_idtut2111 <- c('tut21')
TUT21_ <- TUT21_[,-4]
TUT21_ <- unite(TUT21_, "read_id", c("read_idtut21", "read_idtut211", "read_idtut2111"), sep = "_", remove = T)
TUT21_df <- TUT21_ %>% dplyr::select(read_id, ARG)
TUT21_plasmid_ARG <- merge(all_repaire_mapping[,c(1,12)], TUT21_df, by = "read_id")

total_TUT21_plasmid_ARG_nonReg <- TUT21_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_TUT21_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_TUT21_plasmid_ARG_nonReg$samplename <- rep(c("TUT21"), each = 197)
total_TUT21_plasmid_ARG_nonReg$sampletype <- rep(c("UT"), each = 197)


#TUT4
TUT4 <- read.csv("TUT4.contig.1kb.gene.deepargout.mapping.ARG",header=1, sep = "\t") 
TUT4$ARG <- sapply(strsplit(TUT4$'best.hit', "|", fixed = T),"[", 3)
TUT4_ <- TUT4 %>% mutate(read_idtut4 = sapply( strsplit(read_id,"_", fixed = T),"[[",1)) %>% mutate(read_idtut44 = sapply( strsplit(read_id,"_", fixed = T),"[[",2))
TUT4_$read_idtut444 <- c('tut4')
TUT4_ <- TUT4_[,-4]
TUT4_ <- unite(TUT4_, "read_id", c("read_idtut4", "read_idtut44", "read_idtut444"), sep = "_", remove = T)
TUT4_df <- TUT4_ %>% dplyr::select(read_id, ARG)
TUT4_plasmid_ARG <- merge(all_repaire_mapping[,c(1,13)], TUT4_df, by = "read_id")

total_TUT4_plasmid_ARG_nonReg <- TUT4_plasmid_ARG %>% filter(!ARG %in% regGenes)
names(total_TUT4_plasmid_ARG_nonReg)[2] <- "DepthPG"
total_TUT4_plasmid_ARG_nonReg$samplename <- rep(c("TUT4"), each = 162)
total_TUT4_plasmid_ARG_nonReg$sampletype <- rep(c("UT"), each = 162)

total_plasmid_ARG <- dplyr::bind_rows(total_ALRT2_plasmid_ARG_nonReg, total_ALRT8_plasmid_ARG_nonReg,
                                      total_ALRT9_plasmid_ARG_nonReg, total_TAS2_plasmid_ARG_nonReg,
                                      total_TAS4_plasmid_ARG_nonReg, total_TAS6_plasmid_ARG_nonReg,
                                      total_TAX2_plasmid_ARG_nonReg, total_TAX4_plasmid_ARG_nonReg,
                                      total_TAX6_plasmid_ARG_nonReg, total_TUT1_plasmid_ARG_nonReg,
                                      total_TUT21_plasmid_ARG_nonReg, total_TUT4_plasmid_ARG_nonReg,
                                      total_ULRT2_plasmid_ARG_nonReg, total_ULRT8_plasmid_ARG_nonReg,
                                      total_ULRT9_plasmid_ARG_nonReg, total_UT7_plasmid_ARG_nonReg,
                                      total_UT8_plasmid_ARG_nonReg, total_UT9_plasmid_ARG_nonReg)
total_plasmid_ARG

total_plasmid_ARG_df <- total_plasmid_ARG %>% group_by(samplename) %>% summarise(DepthPG = sum(Abundance)) 
total_plasmid_ARG_df$Sample_type <- c("ALRT","ALRT","ALRT","ALRT","ALRT","ALRT","ULRT", "ULRT","ULRT","UT","UT","UT","ULRT", "ULRT","ULRT","UT","UT","UT")
total_plasmid_ARG_df$Remediate_year <- c("3","3","3","4","4","4","4","4","4","4","4","4","3","3","3","3","3","3")


# plot

total_plasmid_ARG_df <- read.csv("total_plasmid_ARG_df.csv", header = 1, row.names = 1) 
Average <- read.csv("Total_average_plasmid_ARG.csv", header = 1, row.names = 1)

total_plasmid_ARG_df$Sample_type <- factor(total_plasmid_ARG_df$Sample_type, levels = c("UT", "ULRT", "ALRT")) 

p <- ggplot() +
  geom_point(data=total_plasmid_ARG_df, aes(x=Sample_type, y=DepthPG, color=Sample_type, shape=factor(Remediate_year)),size=7, alpha=0.7) + 
  geom_point(data=Average, aes(x=Sample_type, y=average),size=1, alpha=0.7) + 
  theme_prism() + theme(axis.line = element_line(arrow = arrow(length = unit(0.3, 'cm')))) +
  scale_shape_manual(values = c(21,16)) +
  ylim(0,8e-06) + xlab("") +ylab("Plsmids with ARGs (coverage, ×/Gb)") + scale_x_discrete(expand = expansion(0.2,0)) + 
  scale_color_manual(values =  c('#3CA183', '#4DA0CE', '#D7A641')) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 15)) + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.text = element_text(size = 13))
p
p <- p +theme(axis.title = element_text(size = 18, color = "gray10"))
p <- p + theme(axis.title.x = element_text(size = 13)) + theme(axis.title.y = element_text(size = 14))
p 

