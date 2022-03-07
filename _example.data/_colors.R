library(ggsci)
library(scales)
library(RColorBrewer) #自己根据所需要的颜色数量来制备调色板




# ARG　Types ---------------------------


color_ARGtypes_df <- cbind.data.frame(ARGtype_rank,
                                      color = c("#DC0000FF", "#E64B35FF",   "#00A087FF", "#3C5488FF", "#B09C85FF","#F39B7FFF",
                                                "#4DBBD5FF","#8491B4FF","#91D1C2FF","#7E6148FF","gray"),
                                      stringsAsFactors=F)


# samples -------------------

color_layer_df <- cbind.data.frame(layer = c("contaminated", "remediate_top", "remediate_inner"),
                                   color = c( "#EB6363" ,"#619CFF","#00BA38"),
                                   stringsAsFactors=F)

color_year.layer_df <- cbind.data.frame(layer = c("3_contaminated", "4_contaminated","3_remediate_top","4_remediate_top", "3_remediate_inner","4_remediate_inner"),
                                   color = c("#46698366","#466983FF","#00828066","#008280FF","#99660066","#996600FF"),
                                   stringsAsFactors=F)
