#Supplemental Figures

library(ggplot2)
library(Matrix)
library(cowplot)
library(RColorBrewer)
library(grid)
library(Seurat)
library(plyr)
library(dplyr)
library(hdf5r)
library(data.table)
library(readbitmap)
library(rjson)
library(dataVisEasy)
library(glmGamPoi)

#Figure S1A
setwd("/Data/10Xsc")
data<-readRDS("10xsc_SeuratObject.rds")
setwd("/Annotations")
sample_annot<-read.table("10xsc-Sample_Annotations.txt", sep ="\t", header= T, row.names = 1)
sample_annot<-sample_annot[which(rownames(sample_annot) %in% gsub("\\.", "-",colnames(data))),]
data$Cell_Type<-sample_annot$Cell_Type
data$Cell_Type<-factor(data$Cell_Type, levels=c("Astrocyte", "Endothelial", "Ependymal", "Erythroblast", "Immune", "Microglia", 
                                                "Neuroepithelial", "Neuron", "Oligodendrocyte", "Platelets", "Reticulocyte", "Stem cell"))
Idents(data) <- data$Cell_Type
VlnPlot(data, c("Slc12a5", "Syt1", "Thy1", "Rbfox3", "Nefl", "Tubb3"))

#Figure S1B
setwd("/Data/10Xsc")
Neuron<-readRDS("10xsc_NeuronalSubset_SO.rds")
Neuron <- SCTransform(Neuron, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

anchors <- FindTransferAnchors(reference = Neuron, query = brain.merge, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Neuron$seurat_clusters, prediction.assay = TRUE,
                                  weight.reduction = brain.merge[["pca"]], dims = 1:30)
brain.merge[["predictions"]] <- predictions.assay
DefaultAssay(brain.merge) <- "predictions"
SpatialFeaturePlot(brain.merge, features = c("1"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

#Figure S2
setwd("/Data/IF_quantification/Figure S2")

#Male
FB<-read.csv("Male_BS-127_111_Left_Lateral_DMV_PACAP_CHAT_THColocalizedObjects.csv", sep =",", header= T)
Non<-read.csv("Male_BS-127_111_Left_Lateral_DMV_PACAP_CHAT_THRemovedObjects.csv", sep =",", header= T)
All<-rbind(FB,Non)

data<-data.frame(CHAT = All$Intensity_MeanIntensity_CHAT, PACAP = All$Intensity_MeanIntensity_PACAP, TH = All$Intensity_MeanIntensity_TH)
data<-t(data)
colnames(data)<-c(1:48)

annotations<-data.frame(Connectivity = c(rep("Cardiac-projecting", nrow(FB)), rep("Other", nrow(Non))))
rownames(annotations)<-c(1:48)

initiate_params()
set_annotations(annotations)
set_annot_samps(c("Connectivity"))
annot_cols <- list('Connectivity'=c('Cardiac-projecting'='turquoise3', 'Other'='grey69'))
set_annot_cols(annot_cols)

scatterGenes(data, "PACAP","CHAT", na.fix =  T,  color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)
scatterGenes(data, "TH","CHAT", na.fix =  T,  color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)

table(annotations$Connectivity)

#Female
FB<-read.csv("Female_BS-124_132_Left_Lateral_DMV_PACAP_CHAT_THColocalizedObjects.csv", sep =",", header= T)
Non<-read.csv("Female_BS-124_132_Left_Lateral_DMV_PACAP_CHAT_THRemovedObjects.csv", sep =",", header= T)
All<-rbind(FB,Non)

data<-data.frame(CHAT = All$Intensity_MeanIntensity_CHAT, PACAP = All$Intensity_MeanIntensity_PACAP, TH = All$Intensity_MeanIntensity_TH)
data<-t(data)
colnames(data)<-c(1:32)

annotations<-data.frame(Connectivity = c(rep("Cardiac-projecting", nrow(FB)), rep("Other", nrow(Non))))
rownames(annotations)<-c(1:32)

initiate_params()
set_annotations(annotations)
set_annot_samps(c("Connectivity"))
annot_cols <- list('Connectivity'=c('Cardiac-projecting'='turquoise3', 'Other'='grey69'))
set_annot_cols(annot_cols)

scatterGenes(data, "PACAP","CHAT", na.fix =  T,  color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)
scatterGenes(data, "TH","CHAT", na.fix =  T,  color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)

table(annotations$Connectivity)

#Figure S3E
setwd("/Data/IF_quantification/Figure S3")

Pos<-read.csv("BS-127_137-Rostral-Left-DMV-CARTColocalizedObjects.csv", sep = ",")
length(Pos$Intensity_MeanIntensity_CART)

Neg<-read.csv("BS-127_137-Rostral-Left-DMV-CARTRemovedObjects.csv", sep = ",")
length(Neg$Intensity_MeanIntensity_CART)

df<-data.frame(Intensity = c(Pos$Intensity_MeanIntensity_CART, Neg$Intensity_MeanIntensity_CART), 
               FB = c(rep("Cardiac", nrow(Pos)), rep("Other", nrow(Neg))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = FB)) + geom_density(size =2) +
  scale_x_continuous(breaks=seq(0.02,0.07,by = 0.005),limits = c(0.02,0.07))+
  scale_color_manual(values=c("turquoise3", "grey69"), name = "FB") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure S3F
DMV<-read.csv("BS-127_137-Rostral-Left-DMV-CARTColocalizedObjects.csv", sep = ",")
length(DMV$Intensity_MeanIntensity_CART)

NAm<-read.csv("BS-127_137-Rostral-Left-NA-CARTColocalizedObjects.csv", sep = ",")
length(NAm$Intensity_MeanIntensity_CART)


df<-data.frame(Intensity = c(DMV$Intensity_MeanIntensity_CART, NAm$Intensity_MeanIntensity_CART), 
               Region = c(rep("DMV", nrow(DMV)), rep("NA", nrow(NAm))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Region)) + geom_density(size =2) +
  scale_x_continuous(breaks=seq(0.01,0.07,by = 0.005),limits = c(0.01,0.07))+
  scale_color_manual(values=c("green","orange"), name = "Region") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure S3K
Pos<-read.csv("2NDBS-124_129-rostral-Right-DMV-20x-gCART488ColocalizedObjects.csv", sep = ",")
length(Pos$Intensity_MeanIntensity_CART)

Neg<-read.csv("2NDBS-124_129-rostral-Right-DMV-20x-gCART488RemovedObjects.csv", sep = ",")
length(Neg$Intensity_MeanIntensity_CART)

df<-data.frame(Intensity = c(Pos$Intensity_MeanIntensity_CART, Neg$Intensity_MeanIntensity_CART), 
               FB = c(rep("Cardiac", nrow(Pos)), rep("Other", nrow(Neg))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = FB)) + geom_density(size =2) +
  scale_x_continuous(breaks=seq(0.02,0.05,by = 0.005),limits = c(0.02,0.05))+
  scale_color_manual(values=c("turquoise3", "grey69"), name = "FB") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure S3L
DMV<-read.csv("2NDBS-124_129-rostral-Right-DMV-20x-gCART488ColocalizedObjects.csv", sep = ",")
length(DMV$Intensity_MeanIntensity_CART)

NAm<-read.csv("BS-124_129-rostral-Right-NA-20x-gCART488ColocalizedObjects.csv", sep = ",")
length(NAm$Intensity_MeanIntensity_CART)

df<-data.frame(Intensity = c(DMV$Intensity_MeanIntensity_CART, NAm$Intensity_MeanIntensity_CART), 
               Region = c(rep("DMV", nrow(DMV)), rep("NA", nrow(NAm))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Region)) + geom_density(size =2) +
  scale_x_continuous(breaks=seq(0.02,0.07,by = 0.005),limits = c(0.02,0.07))+
  scale_color_manual(values=c("green","orange"), name = "Region") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure S4
setwd("Data/LCM-RNA-seq")
data<-readRDS("LCMseq_SeuratObject.rds")

zero<-rownames(FindMarkers(data, ident.1 = 0, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE))
one<-rownames(FindMarkers(data, ident.1 = 1, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE))
two<-rownames(FindMarkers(data, ident.1 = 2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE))
three<-rownames(FindMarkers(data, ident.1 = 3, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE))
four<-rownames(FindMarkers(data, ident.1 = 4, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE))
five<-rownames(FindMarkers(data, ident.1 = 5, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE))

FeaturePlot(data, features =c(zero[1:2], one[2:3], two[1:2], three[1:2], four[1:2], five[1:2]))

#Figure S5
setwd("Data/IF_quantification/Figure S5/Left_Right")

#Left
Leftbp<-read.table("BS-127_114-lateral-Left-DMV-SNCAImage.csv", sep = ",", header = T)

#Total FB
L_TFB<-Leftbp$Count_IdentifyPrimaryObjectsFB
#Total SNCA
L_TP<-Leftbp$Count_IdentifyPrimaryObjectsSNCA
#Within FB+, how many are SNCA+?
L_FB_P<-Leftbp$Count_ColocalizedObjects / L_TFB
#Within FB+, how many are SNCA-?
L_Pneg<-Leftbp$Count_IdentifyPrimaryObjectsFB - Leftbp$Count_ColocalizedObjects
L_FB_Pneg<-L_Pneg/L_TFB

Rightbp<-read.table("BS-127_114-lateral-Right-DMV-SNCAImage.csv", sep = ",", header = T)
#Right
#Total FB
R_TFB<-Rightbp$Count_IdentifyPrimaryObjectsFB
#Total SNCA
R_TP<-Rightbp$Count_IdentifyPrimaryObjectsSNCA
#Within FB+, how many are SNCA+?
R_FB_P<-Rightbp$Count_ColocalizedObjects / R_TFB
#Within FB+, how many are SNCA-?
R_Pneg<-Rightbp$Count_IdentifyPrimaryObjectsFB - Rightbp$Count_ColocalizedObjects
R_FB_Pneg<-R_Pneg/R_TFB

dfbp<-data.frame(Percentage = c(L_FB_P, L_FB_Pneg, R_FB_P, R_FB_Pneg), Side = c("Left", "Left", "Right", "Right"), Detection = c("SNCA+", "SNCA-", "SNCA+", "SNCA-") )
dfbp$Detection<-factor(dfbp$Detection, levels = c("SNCA+", "SNCA-"))
ggplot(dfbp, aes(x=Detection, y=Percentage, fill = Side)) + 
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  #scale_x_discrete(limits = c("0", "1","2","3","4","5")) +
  scale_y_continuous(breaks=seq(0,1,by =0.5),limits = c(0,1.1), expand = c(0,0))+
  scale_fill_manual(values=c("tan4","violet"), name = "Side") + 
  labs(y = "Percentage")+
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

###################
###Density Plots###
###################

Left<-read.csv("BS-127_114-lateral-Left-DMV-SNCAColocalizedObjects.csv", sep = ",")
length(Left$Intensity_MeanIntensity_SNCA)
Right<-read.csv("BS-127_114-lateral-Right-DMV-SNCAColocalizedObjects.csv", sep = ",")
length(Right$Intensity_MeanIntensity_SNCA)

df<-data.frame(Intensity = c(Left$Intensity_MeanIntensity_SNCA, Right$Intensity_MeanIntensity_SNCA), 
               Side = c(rep("Left", nrow(Left)), rep("Right", nrow(Right))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Side)) + geom_density(size =2) +
  #+ xlim(2,200) 
  scale_x_continuous(breaks=seq(0,0.08,by = 0.01),limits = c(0,0.08))+
  scale_color_manual(values=c("tan4","violet")) + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure S5F
setwd("Data/IF_quantification/Figure S5/Rostral_Int_Caudal")

Caudalbp<-read.table("BS-127_84-Caudal_Left-DMV-20x-cSNCA555Image.csv", sep = ",", header = T)
#Total FB
mC_TFB<-Caudalbp$Count_IdentifyPrimaryObjectsFB
#Total SNCA
mC_TP<-Caudalbp$Count_IdentifyPrimaryObjectsSNCA
#Within FB+, how many are SNCA+?
mC_FB_P<-Caudalbp$Count_ColocalizedObjects / mC_TFB

#Within FB+, how many are SNCA-?
mC_Pneg<-Caudalbp$Count_IdentifyPrimaryObjectsFB - Caudalbp$Count_ColocalizedObjects
mC_FB_Pneg<-mC_Pneg/mC_TFB


Intermediatebp<-read.table("BS-127_114-lateral-Left-DMV-20x-cSNCA555Image.csv", sep = ",", header = T)
#Total FB
mI_TFB<-Intermediatebp$Count_IdentifyPrimaryObjectsFB
#Total SNCA
mI_TP<-Intermediatebp$Count_IdentifyPrimaryObjectsSNCA
#Within FB+, how many are SNCA+?
mI_FB_P<-Intermediatebp$Count_ColocalizedObjects / mI_TFB
#Within FB+, how many are SNCA-?
mI_Pneg<-Intermediatebp$Count_IdentifyPrimaryObjectsFB - Intermediatebp$Count_ColocalizedObjects
mI_FB_Pneg<-mI_Pneg/mI_TFB

Rostralbp<-read.table("BS-127-144-Rostral_Left_DMV_DMV_SNCAImage.csv", sep = ",", header = T)
#Right
#Total FB
mR_TFB<-Rostralbp$Count_IdentifyPrimaryObjectsFB
#Total SNCA
mR_TP<-Rostralbp$Count_IdentifyPrimaryObjectsSNCA
#Within FB+, how many are SNCA+?
mR_FB_P<-Rostralbp$Count_ColocalizedObjects / R_TFB
#Within FB+, how many are SNCA-?
mR_Pneg<-Rostralbp$Count_IdentifyPrimaryObjectsFB - Rostralbp$Count_ColocalizedObjects
mR_FB_Pneg<-mR_Pneg/mR_TFB

dfbp<-data.frame(Percentage = c(mC_FB_P, mC_FB_Pneg, mI_FB_P, mI_FB_Pneg, mR_FB_P, mR_FB_Pneg), Region = c("Caudal", "Caudal", "Intermediate", "Intermediate", "Rostral", "Rostral"), Detection = c("SNCA+", "SNCA-","SNCA+", "SNCA-", "SNCA+", "SNCA-") )
dfbp$Detection<-factor(dfbp$Detection, levels = c("SNCA+", "SNCA-"))
dfbp$Region<-factor(dfbp$Region, levels = c("Rostral", "Intermediate", "Caudal"))

ggplot(dfbp, aes(x=Detection, y=Percentage, fill = Region)) + 
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  #scale_x_discrete(limits = c("0", "1","2","3","4","5")) +
  scale_y_continuous(breaks=seq(0,1,by =0.5),limits = c(0,1.1), expand = c(0,0))+
  scale_fill_manual(values=c("green","orange","purple"), name = "Region") + 
  labs(y = "Percentage")+
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

###############
#Density Plots#
###############

Caudal<-read.csv("BS-127_84-Caudal_Left-DMV-20x-cSNCA555ColocalizedObjects.csv", sep = ",")
length(Caudal$Intensity_MeanIntensity_SNCA)


Intermediate<-read.csv("BS-127_114-lateral-Left-DMV-20x-cSNCA555ColocalizedObjects.csv", sep = ",")
length(Intermediate$Intensity_MeanIntensity_SNCA)

Rostral<-read.csv("BS-127-144-Rostral_left_DMV_SNCARemovedObjects.csv", sep = ",")
length(Rostral$Intensity_MeanIntensity_SNCA)

df<-data.frame(Intensity = c(Caudal$Intensity_MeanIntensity_SNCA, Intermediate$Intensity_MeanIntensity_SNCA, Rostral$Intensity_MeanIntensity_SNCA), 
               Region = c(rep("Caudal", nrow(Caudal)), rep("Intermediate", nrow(Intermediate)), rep("Rostral", nrow(Rostral))))

range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Region)) + geom_density(size = 2) +
  #+ xlim(2,200) 
  scale_x_continuous(breaks=seq(0,0.08,by = 0.01),limits = c(0,0.08))+
  scale_color_manual(values=c("purple","orange","green"), name = "Region") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))


#FigureS6
setwd("Data/IF_quantification/Figure S6/Left_Right")

#Left
Leftbp<-read.table("BS-124-79-_Left_DMV_SNCAImage.csv", sep = ",", header = T)

#Total FB
L_TFB<-Leftbp$Count_IdentifyPrimaryObjectsFB
#Total SNCA
L_TP<-Leftbp$Count_IdentifyPrimaryObjectsSNCA
#Within FB+, how many are SNCA+?
L_FB_P<-Leftbp$Count_ColocalizedObjects / L_TFB
#Within FB+, how many are SNCA-?
L_Pneg<-Leftbp$Count_IdentifyPrimaryObjectsFB - Leftbp$Count_ColocalizedObjects
L_FB_Pneg<-L_Pneg/L_TFB

Rightbp<-read.table("BS-124-79-_Right_DMV_SNCAImage.csv", sep = ",", header = T)

#Right
#Total FB
R_TFB<-Rightbp$Count_IdentifyPrimaryObjectsFB
#Total SNCA
R_TP<-Rightbp$Count_IdentifyPrimaryObjectsSNCA
#Within FB+, how many are SNCA+?
R_FB_P<-Rightbp$Count_ColocalizedObjects / R_TFB
#Within FB+, how many are SNCA-?
R_Pneg<-Rightbp$Count_IdentifyPrimaryObjectsFB - Rightbp$Count_ColocalizedObjects
R_FB_Pneg<-R_Pneg/R_TFB

dfbp<-data.frame(Percentage = c(L_FB_P, L_FB_Pneg, R_FB_P, R_FB_Pneg), Side = c("Left", "Left", "Right", "Right"), Detection = c("SNCA+", "SNCA-", "SNCA+", "SNCA-") )
dfbp$Detection<-factor(dfbp$Detection, levels = c("SNCA+", "SNCA-"))
ggplot(dfbp, aes(x=Detection, y=Percentage, fill = Side)) + 
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  #scale_x_discrete(limits = c("0", "1","2","3","4","5")) +
  scale_y_continuous(breaks=seq(0,1,by =0.5),limits = c(0,1), expand = c(0,0))+
  scale_fill_manual(values=c("tan4","violet"), name = "Side") + 
  labs(y = "Percentage")+
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

###################
###Density Plots###
###################

Left<-read.csv("BS-124-79-_Left_DMV_SNCAColocalizedObjects.csv", sep = ",")
length(Left$Intensity_MeanIntensity_SNCA)
Right<-read.csv("BS-124-79-_Right_DMV_SNCAColocalizedObjects.csv", sep = ",")
length(Right$Intensity_MeanIntensity_SNCA)

df<-data.frame(Intensity = c(Left$Intensity_MeanIntensity_SNCA, Right$Intensity_MeanIntensity_SNCA), 
               Side = c(rep("Left", nrow(Left)), rep("Right", nrow(Right))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Side)) + geom_density(size =2) +
  #+ xlim(2,200) 
  scale_x_continuous(breaks=seq(0.03,0.09,by = 0.01),limits = c(0.03,0.09))+
  scale_color_manual(values=c("tan4","violet")) + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure S6F
setwd("Data/IF_quantification/Figure S6/Rostral_Int_Caudal")

Caudalbp<-read.table("BS-124-23-_Caudal_DMV_SNCAImage.csv", sep = ",", header = T)
#Total FB
mC_TFB<-Caudalbp$Count_IdentifyPrimaryObjectsFB
#Total SNCA
mC_TP<-Caudalbp$Count_IdentifyPrimaryObjectsSNCA
#Within FB+, how many are SNCA+?
mC_FB_P<-Caudalbp$Count_ColocalizedRegion / mC_TFB

#Within FB+, how many are SNCA-?
mC_Pneg<-Caudalbp$Count_IdentifyPrimaryObjectsFB - Caudalbp$Count_ColocalizedRegion
mC_FB_Pneg<-mC_Pneg/mC_TFB


Intermediatebp<-read.table("BS-124-79-_Left_DMV_SNCAImage.csv", sep = ",", header = T)
#Total FB
mI_TFB<-Intermediatebp$Count_IdentifyPrimaryObjectsFB
#Total SNCA
mI_TP<-Intermediatebp$Count_IdentifyPrimaryObjectsSNCA
#Within FB+, how many are SNCA+?
mI_FB_P<-Intermediatebp$Count_ColocalizedRegion / mI_TFB
#Within FB+, how many are SNCA-?
mI_Pneg<-Intermediatebp$Count_IdentifyPrimaryObjectsFB - Intermediatebp$Count_ColocalizedRegion
mI_FB_Pneg<-mI_Pneg/mI_TFB


Rostralbp<-read.table("BS-124-129-_Rostral_Left_SNCAImage.csv", sep = ",", header = T)
#Right
#Total FB
mR_TFB<-Rostralbp$Count_IdentifyPrimaryObjectsFB
#Total SNCA
mR_TP<-Rostralbp$Count_IdentifyPrimaryObjectsSNCA
#Within FB+, how many are SNCA+?
mR_FB_P<-Rostralbp$Count_ColocalizedRegion / mR_TFB
#Within FB+, how many are SNCA-?
mR_Pneg<-Rostralbp$Count_IdentifyPrimaryObjectsFB - Rostralbp$Count_ColocalizedRegion
mR_FB_Pneg<-mR_Pneg/mR_TFB

dfbp<-data.frame(Percentage = c(mC_FB_P, mC_FB_Pneg, mI_FB_P, mI_FB_Pneg, mR_FB_P, mR_FB_Pneg), Region = c("Caudal", "Caudal", "Intermediate", "Intermediate", "Rostral", "Rostral"), Detection = c("SNCA+", "SNCA-","SNCA+", "SNCA-", "SNCA+", "SNCA-") )
dfbp$Detection<-factor(dfbp$Detection, levels = c("SNCA+", "SNCA-"))
dfbp$Region<-factor(dfbp$Region, levels = c("Rostral", "Intermediate", "Caudal"))

ggplot(dfbp, aes(x=Detection, y=Percentage, fill = Region)) + 
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  #scale_x_discrete(limits = c("0", "1","2","3","4","5")) +
  scale_y_continuous(breaks=seq(0,1,by =0.5),limits = c(0,1.1), expand = c(0,0))+
  scale_fill_manual(values=c("green","orange","purple"), name = "Region") + 
  labs(y = "Percentage")+
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))
###################
###Density plots###
###################

Caudal<-read.csv("BS-124-23-_CaudalL_DMV_SNCAColocalizedRegion.csv", sep = ",")
length(Caudal$Intensity_MeanIntensity_SNCA)


Intermediate<-read.csv("BS-124-79-_Left_DMV_SNCAColocalizedRegion.csv", sep = ",")
length(Intermediate$Intensity_MeanIntensity_SNCA)

Rostral<-read.csv("BS-124-129-_Rostral_Left_DMV_SNCAColocalizedRegion.csv", sep = ",")
length(Rostral$Intensity_MeanIntensity_SNCA)

df<-data.frame(Intensity = c(Caudal$Intensity_MeanIntensity_SNCA, Intermediate$Intensity_MeanIntensity_SNCA, Rostral$Intensity_MeanIntensity_SNCA), 
               Region = c(rep("Caudal", nrow(Caudal)), rep("Intermediate", nrow(Intermediate)), rep("Rostral", nrow(Rostral))))

range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Region)) + geom_density(size = 2) +
  #+ xlim(2,200) 
  scale_x_continuous(breaks=seq(0.01,0.1,by = 0.01),limits = c(0.01,0.1))+
  scale_color_manual(values=c("purple","orange","green"), name = "Region") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))
