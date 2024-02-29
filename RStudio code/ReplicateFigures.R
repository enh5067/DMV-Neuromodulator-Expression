#Code to generate figures in DMV manuscript

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
library(readbitmap)

#Figure 1
brain.merge<-readRDS("Visium_SeuratObject.rds")

#Figure 1A
SpatialDimPlot(brain.merge)

#Figure 1B
DimPlot(brain.merge, reduction = "umap", label = T)

#Figure 1C
SpatialDimPlot(brain.merge, images = "Slice1.5")

#Figure 2
brain.merge<-PrepSCTFindMarkers(brain.merge, assay = "SCT", verbose = TRUE)
DE<-FindMarkers(brain.merge, ident.1 = "10")
SS<-DE[which(DE$p_val < 0.05),]
SS_LFC<-SS[order(SS$avg_log2FC, decreasing = T),]
#SpatialFeaturePlot(brain.merge, features = "Adcyap1")

##10x genomics pipeline to visualize a gene from 10x genomics website (February 2022):
geom_spatial <-  function(mapping = NULL,
                          data = NULL,
                          stat = "identity",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = FALSE,
                          ...) {
  
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

sample_names <- c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8")

image_paths <- c("S1/outs/spatial/tissue_lowres_image.png",
                 "S2/outs/spatial/tissue_lowres_image.png",
                 "S3/outs/spatial/tissue_lowres_image.png",
                 "S4/outs/spatial/tissue_lowres_image.png",
                 "S5/outs/spatial/tissue_lowres_image.png",
                 "S6/outs/spatial/tissue_lowres_image.png",
                 "S7/outs/spatial/tissue_lowres_image.png",
                 "S8/outs/spatial/tissue_lowres_image.png")

scalefactor_paths <- c("S1/outs/spatial/scalefactors_json.json",
                       "S2/outs/spatial/scalefactors_json.json",
                       "S3/outs/spatial/scalefactors_json.json",
                       "S4/outs/spatial/scalefactors_json.json",
                       "S5/outs/spatial/scalefactors_json.json",
                       "S6/outs/spatial/scalefactors_json.json",
                       "S7/outs/spatial/scalefactors_json.json",
                       "S8/outs/spatial/scalefactors_json.json")

tissue_paths <- c("S1/outs/spatial/tissue_positions_list.csv",
                       "S2/outs/spatial/tissue_positions_list.csv",
                       "S3/outs/spatial/tissue_positions_list.csv",
                       "S4/outs/spatial/tissue_positions_list.csv",
                       "S5/outs/spatial/tissue_positions_list.csv",
                       "S6/outs/spatial/tissue_positions_list.csv",
                       "S7/outs/spatial/tissue_positions_list.csv",
                       "S8/outs/spatial/tissue_positions_list.csv") 

matrix_paths <- c("S1/outs/filtered_feature_bc_matrix.h5",
                  "S2/outs/filtered_feature_bc_matrix.h5",
                  "S3/outs/filtered_feature_bc_matrix.h5",
                  "S4/outs/filtered_feature_bc_matrix.h5",
                  "S5/outs/filtered_feature_bc_matrix.h5",
                  "S6/outs/filtered_feature_bc_matrix.h5",
                  "S7/outs/filtered_feature_bc_matrix.h5",
                  "S8/outs/filtered_feature_bc_matrix.h5") 


cluster_paths <- c("S1/outs/clusters.csv",
                   "S2/outs/clusters.csv",
                   "S3/outs/clusters.csv",
                   "S4/outs/clusters.csv",
                   "S5/outs/clusters.csv",
                   "S6/outs/clusters.csv",
                   "S7/outs/clusters.csv",
                   "S8/outs/clusters.csv")

images_cl <- list()

for (i in 1:length(sample_names)) {
  images_cl[[i]] <- read.bitmap(image_paths[i])
}

height <- list()

for (i in 1:length(sample_names)) {
  height[[i]] <-  data.frame(height = nrow(images_cl[[i]]))
}

height <- bind_rows(height)

width <- list()

for (i in 1:length(sample_names)) {
  width[[i]] <- data.frame(width = ncol(images_cl[[i]]))
}

width <- bind_rows(width)

grobs <- list()
for (i in 1:length(sample_names)) {
  grobs[[i]] <- rasterGrob(images_cl[[i]], width=unit(1,"npc"), height=unit(1,"npc"))
}

images_tibble <- tibble(sample=factor(sample_names), grob=grobs)
images_tibble$height <- height$height
images_tibble$width <- width$width
scales <- list()

for (i in 1:length(sample_names)) {
  scales[[i]] <- rjson::fromJSON(file = scalefactor_paths[i])
}

clusters <- list()
for (i in 1:length(sample_names)) {
  clusters[[i]] <- read.csv(cluster_paths[i])
}

bcs <- list()

for (i in 1:length(sample_names)) {
  bcs[[i]] <- read.csv(tissue_paths[i],col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = TRUE)
  bcs[[i]]$imagerow <- bcs[[i]]$imagerow * scales[[i]]$tissue_lowres_scalef    # scale tissue coordinates for lowres image
  bcs[[i]]$imagecol <- bcs[[i]]$imagecol * scales[[i]]$tissue_lowres_scalef
  bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue)
  bcs[[i]] <- merge(bcs[[i]], clusters[[i]], by.x = "barcode", by.y = "Barcode", all = TRUE)
  bcs[[i]]$height <- height$height[i]
  bcs[[i]]$width <- width$width[i]
}

names(bcs) <- sample_names

matrix <- list()

for (i in 1:length(sample_names)) {
  matrix[[i]] <- as.data.frame(t(Read10X_h5(matrix_paths[i])))
}

bcs_merge <- bind_rows(bcs, .id = "sample")

umi_sum <- list()

for (i in 1:length(sample_names)) {
  umi_sum[[i]] <- data.frame(barcode =  row.names(matrix[[i]]),
                             sum_umi = Matrix::rowSums(matrix[[i]]))
}

names(umi_sum) <- sample_names

umi_sum <- bind_rows(umi_sum, .id = "sample")

bcs_merge <- merge(bcs_merge,umi_sum, by = c("barcode", "sample"))

gene_sum <- list()

for (i in 1:length(sample_names)) {
  gene_sum[[i]] <- data.frame(barcode =  row.names(matrix[[i]]),
                              sum_gene = Matrix::rowSums(matrix[[i]] != 0))
  
}

names(gene_sum) <- sample_names

gene_sum <- bind_rows(gene_sum, .id = "sample")

bcs_merge <- merge(bcs_merge,gene_sum, by = c("barcode", "sample"))

xlim(0,max(bcs_merge %>%
             filter(sample ==sample_names[i]) %>%
             select(width)))

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

#Figure 2
#Substitute Cartpt, Th, Snca, and Sncg in for Adcyap1 in BOTH places below
plots <- list()

for (i in 1:length(sample_names)) {
  plots[[i]] <- bcs_merge %>% 
    filter(sample ==sample_names[i]) %>% 
    bind_cols(as.data.table(matrix[i])[, "Adcyap1", with=FALSE]) %>% 
    ggplot(aes(x=imagecol,y=imagerow,fill=Adcyap1)) +
    geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.50, y=0.50)+
    geom_point(shape = 21, colour = "black", size = 5, stroke = 0.25)+
    coord_cartesian(expand=FALSE)+
    #scale_fill_gradientn(colours = c("lightgray", "red"))+
    scale_fill_gradientn(colours = myPalette(100))+
    #scale_fill_gradientn(colours = colors)+
    xlim(0,max(bcs_merge %>% 
                 filter(sample ==sample_names[i]) %>% 
                 select(width)))+
    ylim(max(bcs_merge %>% 
               filter(sample ==sample_names[i]) %>% 
               select(height)),0)+
    xlab("") +
    ylab("") +
    ggtitle(sample_names[i])+
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

plot_grid(plotlist = plots)

#Figure 3A

data<-readRDS("10xsc_SeuratObject.rds")
data$seurat_clusters<-mapvalues(data$seurat_clusters, c(0:9), c(1:10))
Idents(data)<-data$seurat_clusters
DimPlot(data, reduction="umap")

#Figure 3B
sample_annot<-read.table("10xsc-Sample_Annotations.txt", sep ="\t", header= T, row.names = 1)
sample_annot<-sample_annot[which(rownames(sample_annot) %in% gsub("\\.", "-",colnames(data))),]
Idents(data)<-sample_annot$Cell_Type
DimPlot(data, reduction="umap")

#Figure 3C
Neuron<-readRDS("10xsc_NeuronalSubset_SeuratObject.rds")
DimPlot(Neuron, reduction="umap")

#Figure 3D
Neuron <- SCTransform(Neuron, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

anchors <- FindTransferAnchors(reference = Neuron, query = brain.merge, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Neuron$seurat_clusters, prediction.assay = TRUE,
                                  weight.reduction = brain.merge[["pca"]], dims = 1:30)

brain.merge[["predictions"]] <- predictions.assay
DefaultAssay(brain.merge) <- "predictions"
SpatialFeaturePlot(brain.merge, features = c("0"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

#Figure 3F
DE<-FindMarkers(Neuron, ident.1 = "0", ident.2 = "1")
VlnPlot(Neuron, features = c("Btg2", "Ier2", "Egr1", "Nfkbiz", "Tnfaip3", "Nfkbia"))

#Figure 4
library(pheatmap)
library(pcaMethods)
library(ggplot2)
library(cluster)
library(gtools)
library(ggbeeswarm)
library(rgl)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)

#Figure 4A
data<-read.table("LCM-qPCR_Negddct.txt", sep = "\t", header = T, row.names = 1)
colnames(data)<-gsub("\\.", "-", colnames(data))

annotations<-read.table("LCM-qPCR_Sample_Annotations.txt", sep = "\t", header = T, row.names = 1)
annotations<-annotations[match(colnames(data), rownames(annotations)),]

data<-data[,which(annotations$Region %in% "DMV")]
annotations<-annotations[which(annotations$Region %in% "DMV"),]

Connectivity<-rep(NA,length=length(colnames(data)))
Connectivity[grep("FB", annotations$Connectivity)]<-"FB"
Connectivity[grep("Non-FB", annotations$Connectivity)]<-"Non-FB"

library(dataVisEasy)
initiate_params()
set_annotations(annotations)
set_annot_samps(c("Connectivity"))

aov.results <- AOV1way(data,  "Connectivity")
head(aov.results$AOV.Results)

myHeatmapByAnnotation(data, list = aov.results$Sig.Genes, groupings = "Connectivity")

#Figure 4B
data<-read.table("LCM-qPCR_Negddct.txt", sep = "\t", header = T, row.names = 1)
colnames(data)<-gsub("\\.", "-", colnames(data))

setwd("/Annotations")
annotations<-read.table("LCM-qPCR_Sample_Annotations.txt", sep = "\t", header = T, row.names = 1)
annotations<-annotations[match(colnames(data), rownames(annotations)),]

data<-data[,which(annotations$Connectivity %in% "FB")]
annotations<-annotations[which(annotations$Connectivity %in% "FB"),]

Region<-rep(NA,length=length(colnames(data)))
Region[grep("DMV", annotations$Region)]<-"DMV"
Region[grep("NAmb", annotations$Region)]<-"NAmb"

library(dataVisEasy)
initiate_params()
set_annotations(annotations)
set_annot_samps("Region")

aov.results <- AOV1way(data, "Region")
head(aov.results$AOV.Results)

aov.results
myHeatmapByAnnotation(data, list = aov.results$Sig.Genes, groupings = "Region")

#Figure 4G
Pos<-read.csv("Male_BS-127_111_Left_Medial_DMV_PACAPColocalizedObjects.csv", sep = ",")
length(Pos$Intensity_MeanIntensity_PACAP)

Neg<-read.csv("Male_BS-127_111_Left_Medial_DMV_PACAPRemovedObjects.csv", sep = ",")
length(Neg$Intensity_MeanIntensity_PACAP)

df<-data.frame(Intensity = c(Pos$Intensity_MeanIntensity_PACAP, Neg$Intensity_MeanIntensity_PACAP), 
               FB = c(rep("Cardiac", nrow(Pos)), rep("Other", nrow(Neg))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = FB)) + geom_density(size =2) +
  scale_x_continuous(breaks=seq(0.01,0.04,by = 0.005),limits = c(0.01,0.04))+
  scale_color_manual(values=c("turquoise3", "grey69"), name = "FB") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure 4H
DMV<-read.csv("Male_BS-127_111_Left_Medial_DMV_PACAPColocalizedObjects.csv", sep = ",")
length(DMV$Intensity_MeanIntensity_PACAP)

NAm<-read.csv("Male_BS-127_111_NA_PACAPColocalizedObjects.csv", sep = ",")
length(NAm$Intensity_MeanIntensity_PACAP)

df<-data.frame(Intensity = c(DMV$Intensity_MeanIntensity_PACAP, NAm$Intensity_MeanIntensity_PACAP), 
               Region = c(rep("DMV", nrow(DMV)), rep("NA", nrow(NAm))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Region)) + geom_density(size =2) +
  scale_x_continuous(breaks=seq(0.01,0.04,by = 0.005),limits = c(0.01,0.04))+
  scale_color_manual(values=c("green","orange"), name = "Region") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure 4I
Pos<-read.csv("Female_BS-124_132_Left__DMV_PACAPColocalizedObjects.csv", sep = ",")
length(Pos$Intensity_MeanIntensity_PACAP)

Neg<-read.csv("Female_BS-124_132_Left__DMV_PACAPRemovedObjects.csv", sep = ",")
length(Neg$Intensity_MeanIntensity_PACAP)

df<-data.frame(Intensity = c(Pos$Intensity_MeanIntensity_PACAP, Neg$Intensity_MeanIntensity_PACAP), 
               FB = c(rep("Cardiac", nrow(Pos)), rep("Other", nrow(Neg))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = FB)) + geom_density(size =2) +
  scale_x_continuous(breaks=seq(0.02,0.07,by = 0.005),limits = c(0.02,0.07))+
  scale_color_manual(values=c("turquoise3", "grey69"), name = "FB") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure 4J
DMV<-read.csv("Female_BS-124_132_Left__DMV_PACAPColocalizedObjects.csv", sep = ",")
length(DMV$Intensity_MeanIntensity_PACAP)

NAm<-read.csv("Female_BS-124_132_Left__NA_PACAPColocalizedObjects.csv", sep = ",")
length(NAm$Intensity_MeanIntensity_PACAP)


df<-data.frame(Intensity = c(DMV$Intensity_MeanIntensity_PACAP, NAm$Intensity_MeanIntensity_PACAP), 
               Region = c(rep("DMV", nrow(DMV)), rep("NA", nrow(NAm))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Region)) + geom_density(size =2) +
  scale_x_continuous(breaks=seq(0.02,0.06,by = 0.005),limits = c(0.02,0.06))+
  scale_color_manual(values=c("green","orange"), name = "Region") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure 5
data<-read.table("LCM-qPCR_Negddct.txt", sep = "\t", header = T, row.names = 1)
colnames(data)<-gsub("\\.", "-", colnames(data))

setwd("/Annotations")
annotations<-read.table("LCM-qPCR_Sample_Annotations.txt", sep = "\t", header = T, row.names = 1)
gene_annot<-read.table("LCM-qPCR_Gene_Annotations.txt", sep = "\t", header = T, row.names = 1)
annotations<-annotations[match(colnames(data), rownames(annotations)),]

library(dataVisEasy)
initiate_params()
set_annotations(annotations)
set_annot_samps(c("Connectivity", "Cell_State", "Region"))
gene_annot$Order<-gene_annot$Gene
set_annotations.genes(gene_annot)
set_annot_genes(c("Rationale","Order"))
annot_cols <- list('Connectivity'=c('FB'='turquoise3', 'Non-FB'='grey69'))
set_annot_cols(annot_cols)

GOI<-gene_annot[which(gene_annot$Rationale %in% c("Neuropeptide", "Neurotransmitter synthesis enzymes")),]
GOI<-GOI$Gene
GOI<-GOI[order(GOI)]
#Figure 5A
myHeatmapByAnnotation(data, list = GOI, clust.rows = F, groupings.gaps = c(0,1,6), groupings.genes.gaps = c(0,3), groupings = c("Connectivity", "Cell_State", "Region"), groupings.genes = c("Order","Rationale"))

#Figure 5B
DMVannot<-annotations[which(annotations$Region %in% "DMV"),]

dat<- c(length(intersect(which(DMVannot$Cell_State %in% "A"),which(DMVannot$Connectivity %in% "FB")))/ length(which(DMVannot$Connectivity %in% "FB")),
        length(intersect(which(DMVannot$Cell_State %in% "A"),which(DMVannot$Connectivity %in% "Non-FB")))/ length(which(DMVannot$Connectivity %in% "Non-FB")),
        length(intersect(which(DMVannot$Cell_State %in% "B"),which(DMVannot$Connectivity %in% "FB")))/ length(which(DMVannot$Connectivity %in% "FB")),
        length(intersect(which(DMVannot$Cell_State %in% "B"),which(DMVannot$Connectivity %in% "Non-FB")))/ length(which(DMVannot$Connectivity %in% "Non-FB")),
        length(intersect(which(DMVannot$Cell_State %in% "C"),which(DMVannot$Connectivity %in% "FB")))/ length(which(DMVannot$Connectivity %in% "FB")),
        length(intersect(which(DMVannot$Cell_State %in% "C"),which(DMVannot$Connectivity %in% "Non-FB")))/ length(which(DMVannot$Connectivity %in% "Non-FB")),
        length(intersect(which(DMVannot$Cell_State %in% "D"),which(DMVannot$Connectivity %in% "FB")))/ length(which(DMVannot$Connectivity %in% "FB")),
        length(intersect(which(DMVannot$Cell_State %in% "D"),which(DMVannot$Connectivity %in% "Non-FB")))/ length(which(DMVannot$Connectivity %in% "Non-FB")),
        length(intersect(which(DMVannot$Cell_State %in% "E"),which(DMVannot$Connectivity %in% "FB")))/ length(which(DMVannot$Connectivity %in% "FB")),
        length(intersect(which(DMVannot$Cell_State %in% "E"),which(DMVannot$Connectivity %in% "Non-FB")))/ length(which(DMVannot$Connectivity %in% "Non-FB")),
        length(intersect(which(DMVannot$Cell_State %in% "F"),which(DMVannot$Connectivity %in% "FB")))/ length(which(DMVannot$Connectivity %in% "FB")),
        length(intersect(which(DMVannot$Cell_State %in% "F"),which(DMVannot$Connectivity %in% "Non-FB")))/ length(which(DMVannot$Connectivity %in% "Non-FB"))
)

colors <- rep(c("FB", "Non-FB"), 6)
Names <- rep(c("A","B","C","D","E","F"), each = 2)
df <- data.frame("Proportion" = dat, "Cell_State" = Names, "Connectivity" = colors)

ggplot(df, aes(x=Cell_State, y=Proportion, fill=Connectivity)) + 
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  scale_x_discrete(limits = c("A", "B","C","D","E","F")) +
  scale_y_continuous(breaks=seq(0,0.7,by =0.1),limits = c(0,0.7), expand = c(0,0))+
  scale_fill_manual(values=c("turquoise3","grey69"), name = "Connectivity") + 
  labs(x="Neuronal State", y = "Proportion")+
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"))

NAmbannot<-annotations[which(annotations$Region %in% "NAmb"),]

dat<- c(length(intersect(which(NAmbannot$Cell_State %in% "A"),which(NAmbannot$Connectivity %in% "FB")))/ length(which(NAmbannot$Connectivity %in% "FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "A"),which(NAmbannot$Connectivity %in% "Non-FB")))/ length(which(NAmbannot$Connectivity %in% "Non-FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "B"),which(NAmbannot$Connectivity %in% "FB")))/ length(which(NAmbannot$Connectivity %in% "FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "B"),which(NAmbannot$Connectivity %in% "Non-FB")))/ length(which(NAmbannot$Connectivity %in% "Non-FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "C"),which(NAmbannot$Connectivity %in% "FB")))/ length(which(NAmbannot$Connectivity %in% "FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "C"),which(NAmbannot$Connectivity %in% "Non-FB")))/ length(which(NAmbannot$Connectivity %in% "Non-FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "D"),which(NAmbannot$Connectivity %in% "FB")))/ length(which(NAmbannot$Connectivity %in% "FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "D"),which(NAmbannot$Connectivity %in% "Non-FB")))/ length(which(NAmbannot$Connectivity %in% "Non-FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "E"),which(NAmbannot$Connectivity %in% "FB")))/ length(which(NAmbannot$Connectivity %in% "FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "E"),which(NAmbannot$Connectivity %in% "Non-FB")))/ length(which(NAmbannot$Connectivity %in% "Non-FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "F"),which(NAmbannot$Connectivity %in% "FB")))/ length(which(NAmbannot$Connectivity %in% "FB")),
        length(intersect(which(NAmbannot$Cell_State %in% "F"),which(NAmbannot$Connectivity %in% "Non-FB")))/ length(which(NAmbannot$Connectivity %in% "Non-FB"))
)

colors <- rep(c("FB", "Non-FB"), 6)
Names <- rep(c("A","B","C","D","E","F"), each = 2)
df <- data.frame("Proportion" = dat, "Cell_State" = Names, "Connectivity" = colors)

ggplot(df, aes(x=Cell_State, y=Proportion, fill=Connectivity)) + 
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  scale_x_discrete(limits = c("A", "B","C","D","E","F")) +
  scale_y_continuous(breaks=seq(0,0.7,by =0.1),limits = c(0,0.7), expand = c(0,0))+
  scale_fill_manual(values=c("turquoise3","grey69"), name = "Connectivity") + 
  labs(x="Neuronal State", y = "Proportion")+
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"))

#Figure 5C
library(dataVisEasy)
initiate_params()

DMV<-data[,which(annotations$Region %in% "DMV")]
set_annotations(DMVannot)
set_annot_samps(c("Connectivity"))

annot_cols <- list('Connectivity'=c('FB'='turquoise3', 'Non-FB'='grey69'))
set_annot_cols(annot_cols)

beeswarmGenes(DMV, c("Adcyap1", "Cartpt", "Chat", "Th"), groupby.x = "Connectivity", color.by = "Connectivity", facet.wrap = T, scales = "fixed", na.fix = T) + geom_boxplot(alpha = 0.5)

NAmb<-data[,which(annotations$Region %in% "NAmb")]
set_annotations(NAmbannot)
set_annot_samps(c("Connectivity"))
beeswarmGenes(NAmb, c("Adcyap1", "Cartpt", "Chat", "Th"), groupby.x = "Connectivity", color.by = "Connectivity", facet.wrap = T, scales = "fixed", na.fix = T) + geom_boxplot(alpha = 0.5)

#Figure 5D
set_annotations(DMVannot)
set_annot_samps(c("Connectivity"))
scatterGenes(DMV,"Chat","Adcyap1", na.fix =  T, color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)
scatterGenes(DMV,"Chat","Th", na.fix =  T, color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)

set_annotations(NAmbannot)
set_annot_samps(c("Connectivity"))
scatterGenes(NAmb,"Chat","Adcyap1", na.fix =  T, color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)
scatterGenes(NAmb,"Chat","Th", na.fix =  T, color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)

#Figure 6A
setwd("Data/LCM-RNA-seq")
data<-readRDS("LCMseq_SeuratObject.rds")
DimPlot(data, reduction="umap")

#Figure 6B
all.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers<-all.markers[order(all.markers$avg_log2FC, decreasing = T),]
zero<-all.markers$gene[which(all.markers$cluster %in% "0")][1:10]
one<-all.markers$gene[which(all.markers$cluster %in% "1")][1:10]
two<-all.markers$gene[which(all.markers$cluster %in% "2")][1:10]
three<-all.markers$gene[which(all.markers$cluster %in% "3")][1:10]
four<-all.markers$gene[which(all.markers$cluster %in% "4")][1:10]
five<-all.markers$gene[which(all.markers$cluster %in% "5")][1:10]
genes<-c(zero,one,two,three,four,five)
gene_annot<-data.frame(genes = genes, cluster = rep(0:5, each = 10))
rownames(gene_annot)<-genes
#top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10)
data<-GetAssayData(data, slot="scale.data")

setwd("/Annotations")
annotations<-read.table("LCMseq-Sample_Annotations.txt", sep = "\t", header =T, row.names = 1)
annotations<-annotations[match(colnames(data), rownames(annotations)),]

initiate_params()
set_annotations.genes(gene_annot)
set_annot_genes("cluster")
set_annotations(annotations)
set_annot_samps("Cell_State")

myHeatmapByAnnotation(data, list = genes, groupings.genes = "cluster",groupings = "Cell_State")

#Figure 6C
dat<- c(length(intersect(which(annotations$Cell_State %in% "0"),which(annotations$Connectivity %in% "FB")))/ length(which(annotations$Connectivity %in% "FB")),
        length(intersect(which(annotations$Cell_State %in% "0"),which(annotations$Connectivity %in% "Non-FB")))/ length(which(annotations$Connectivity %in% "Non-FB")),
        length(intersect(which(annotations$Cell_State %in% "1"),which(annotations$Connectivity %in% "FB")))/ length(which(annotations$Connectivity %in% "FB")),
        length(intersect(which(annotations$Cell_State %in% "1"),which(annotations$Connectivity %in% "Non-FB")))/ length(which(annotations$Connectivity %in% "Non-FB")),
        length(intersect(which(annotations$Cell_State %in% "2"),which(annotations$Connectivity %in% "FB")))/ length(which(annotations$Connectivity %in% "FB")),
        length(intersect(which(annotations$Cell_State %in% "2"),which(annotations$Connectivity %in% "Non-FB")))/ length(which(annotations$Connectivity %in% "Non-FB")),
        length(intersect(which(annotations$Cell_State %in% "3"),which(annotations$Connectivity %in% "FB")))/ length(which(annotations$Connectivity %in% "FB")),
        length(intersect(which(annotations$Cell_State %in% "3"),which(annotations$Connectivity %in% "Non-FB")))/ length(which(annotations$Connectivity %in% "Non-FB")),
        length(intersect(which(annotations$Cell_State %in% "4"),which(annotations$Connectivity %in% "FB")))/ length(which(annotations$Connectivity %in% "FB")),
        length(intersect(which(annotations$Cell_State %in% "4"),which(annotations$Connectivity %in% "Non-FB")))/ length(which(annotations$Connectivity %in% "Non-FB")),
        length(intersect(which(annotations$Cell_State %in% "5"),which(annotations$Connectivity %in% "FB")))/ length(which(annotations$Connectivity %in% "FB")),
        length(intersect(which(annotations$Cell_State %in% "5"),which(annotations$Connectivity %in% "Non-FB")))/ length(which(annotations$Connectivity %in% "Non-FB"))
)

colors <- rep(c("FB", "Non-FB"), 6)
Names <- rep(c("1","2","3","4","5","6"), each = 2)
df <- data.frame("Proportion" = dat, "Cell_State" = Names, "Connectivity" = colors)

ggplot(df, aes(x=Cell_State, y=Proportion, fill=Connectivity)) + 
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  scale_x_discrete(limits = c("1", "2","3","4","5","6")) +
  scale_y_continuous(breaks=seq(0,0.5,by =0.1),limits = c(0,0.5), expand = c(0,0))+
  scale_fill_manual(values=c("turquoise3","grey69"), name = "Connectivity") + 
  labs(x="Neuronal State", y = "Proportion")+
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"))

dat<- c(length(intersect(which(annotations$Cell_State %in% "0"),which(annotations$Side %in% "Left")))/ length(which(annotations$Side %in% "Left")),
        length(intersect(which(annotations$Cell_State %in% "0"),which(annotations$Side %in% "Right")))/ length(which(annotations$Side %in% "Right")),
        length(intersect(which(annotations$Cell_State %in% "1"),which(annotations$Side %in% "Left")))/ length(which(annotations$Side %in% "Left")),
        length(intersect(which(annotations$Cell_State %in% "1"),which(annotations$Side %in% "Right")))/ length(which(annotations$Side %in% "Right")),
        length(intersect(which(annotations$Cell_State %in% "2"),which(annotations$Side %in% "Left")))/ length(which(annotations$Side %in% "Left")),
        length(intersect(which(annotations$Cell_State %in% "2"),which(annotations$Side %in% "Right")))/ length(which(annotations$Side %in% "Right")),
        length(intersect(which(annotations$Cell_State %in% "3"),which(annotations$Side %in% "Left")))/ length(which(annotations$Side %in% "Left")),
        length(intersect(which(annotations$Cell_State %in% "3"),which(annotations$Side %in% "Right")))/ length(which(annotations$Side %in% "Right")),
        length(intersect(which(annotations$Cell_State %in% "4"),which(annotations$Side %in% "Left")))/ length(which(annotations$Side %in% "Left")),
        length(intersect(which(annotations$Cell_State %in% "4"),which(annotations$Side %in% "Right")))/ length(which(annotations$Side %in% "Right")),
        length(intersect(which(annotations$Cell_State %in% "5"),which(annotations$Side %in% "Left")))/ length(which(annotations$Side %in% "Left")),
        length(intersect(which(annotations$Cell_State %in% "5"),which(annotations$Side %in% "Right")))/ length(which(annotations$Side %in% "Right"))
)

colors <- rep(c("Left", "Right"), 6)
Names <- rep(c("1","2","3","4","5","6"), each = 2)
df <- data.frame("Proportion" = dat, "Cell_State" = Names, "Side" = colors)

ggplot(df, aes(x=Cell_State, y=Proportion, fill=Side)) + 
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  scale_x_discrete(limits = c("1", "2","3","4","5","6")) +
  scale_y_continuous(breaks=seq(0,0.5,by =0.1),limits = c(0,0.5), expand = c(0,0))+
  scale_fill_manual(values=c("tan4","violet"), name = "Side") + 
  labs(x="Neuronal State", y = "Proportion")+
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"))

dat<- c(length(intersect(which(annotations$Cell_State %in% "0"),which(annotations$Position %in% "Rostral")))/ length(which(annotations$Position %in% "Rostral")),
        length(intersect(which(annotations$Cell_State %in% "0"),which(annotations$Position %in% "Intermediate")))/ length(which(annotations$Position %in% "Intermediate")),
        length(intersect(which(annotations$Cell_State %in% "0"),which(annotations$Position %in% "Caudal")))/ length(which(annotations$Position %in% "Caudal")),
        length(intersect(which(annotations$Cell_State %in% "1"),which(annotations$Position %in% "Rostral")))/ length(which(annotations$Position %in% "Rostral")),
        length(intersect(which(annotations$Cell_State %in% "1"),which(annotations$Position %in% "Intermediate")))/ length(which(annotations$Position %in% "Intermediate")),
        length(intersect(which(annotations$Cell_State %in% "1"),which(annotations$Position %in% "Caudal")))/ length(which(annotations$Position %in% "Caudal")),
        length(intersect(which(annotations$Cell_State %in% "2"),which(annotations$Position %in% "Rostral")))/ length(which(annotations$Position %in% "Rostral")),
        length(intersect(which(annotations$Cell_State %in% "2"),which(annotations$Position %in% "Intermediate")))/ length(which(annotations$Position %in% "Intermediate")),
        length(intersect(which(annotations$Cell_State %in% "2"),which(annotations$Position %in% "Caudal")))/ length(which(annotations$Position %in% "Caudal")),
        length(intersect(which(annotations$Cell_State %in% "3"),which(annotations$Position %in% "Rostral")))/ length(which(annotations$Position %in% "Rostral")),
        length(intersect(which(annotations$Cell_State %in% "3"),which(annotations$Position %in% "Intermediate")))/ length(which(annotations$Position %in% "Intermediate")),
        length(intersect(which(annotations$Cell_State %in% "3"),which(annotations$Position %in% "Caudal")))/ length(which(annotations$Position %in% "Caudal")),
        length(intersect(which(annotations$Cell_State %in% "4"),which(annotations$Position %in% "Rostral")))/ length(which(annotations$Position %in% "Rostral")),
        length(intersect(which(annotations$Cell_State %in% "4"),which(annotations$Position %in% "Intermediate")))/ length(which(annotations$Position %in% "Intermediate")),
        length(intersect(which(annotations$Cell_State %in% "4"),which(annotations$Position %in% "Caudal")))/ length(which(annotations$Position %in% "Caudal")),
        length(intersect(which(annotations$Cell_State %in% "5"),which(annotations$Position %in% "Rostral")))/ length(which(annotations$Position %in% "Rostral")),
        length(intersect(which(annotations$Cell_State %in% "5"),which(annotations$Position %in% "Intermediate")))/ length(which(annotations$Position %in% "Intermediate")),
        length(intersect(which(annotations$Cell_State %in% "5"),which(annotations$Position %in% "Caudal")))/ length(which(annotations$Position %in% "Caudal"))
        
)

colors <- rep(c("Rostral", "Intermediate", "Caudal"), 6)
Names <- rep(c("1","2","3","4","5","6"), each = 3)
df <- data.frame("Proportion" = dat, "Cell_State" = Names, "Position" = colors)
df$Position<-factor(x = df$Position, levels = c("Rostral", "Intermediate", "Caudal"))
ggplot(df, aes(x=Cell_State, y=Proportion, fill=Position)) + 
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  scale_x_discrete(limits = c("1", "2","3","4","5","6")) +
  scale_y_continuous(breaks=seq(0,0.5,by =0.1),limits = c(0,0.5), expand = c(0,0))+
  scale_fill_manual(values=c("Rostral" = "green","Intermediate"="orange", "Caudal"="purple"), name = "Position") + 
  labs(x="Neuronal State", y = "Proportion")+
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"))

#Figure 6D
data<-readRDS("LCMseq_SeuratObject.rds")
DE<-FindMarkers(data, ident.1 = "1")
FeaturePlot(data, features =c("Adcyap1", "Snca", "Sncg"))

#Figure 6E
initiate_params()
set_annotations(annotations)
set_annot_samps(c("Connectivity"))

annot_cols <- list('Connectivity'=c('FB'='turquoise3', 'Non-FB'='grey69'))
set_annot_cols(annot_cols)

data<-GetAssayData(data, slot = "scale.data")
scatterGenes(data, "Chat","Adcyap1", na.fix =  T, color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)
scatterGenes(data, "Chat","Snca", na.fix =  T,  color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)
scatterGenes(data, "Chat","Sncg", na.fix =  T,  color.by = "Connectivity", point.size = 5) + geom_hline(yintercept=0) +geom_vline(xintercept = 0)

#Figure 6F
data<-readRDS("LCMseq_SeuratObject.rds")
data <- SCTransform(data, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

anchors <- FindTransferAnchors(reference = data, query = brain.merge, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = data$seurat_clusters, prediction.assay = TRUE,
                                  weight.reduction = brain.merge[["pca"]], dims = 1:30)
brain.merge[["predictions"]] <- predictions.assay
DefaultAssay(brain.merge) <- "predictions"
SpatialFeaturePlot(brain.merge, features = c("1"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

#Figure 7
#Left
Leftbp<-read.table("Male_BS-127_111_Left_LateralDMVImage.csv", sep = ",", header = T)

#Total FB
L_TFB<-Leftbp$Count_IdentifyPrimaryObjectsFB
#Total PACAP
L_TP<-Leftbp$Count_IdentifyPrimaryObjectsPACAP
#Within FB+, how many are PACAP+?
L_FB_P<-Leftbp$Count_ColocalizedRegion / L_TFB
#Within FB+, how many are PACAP-?
L_Pneg<-Leftbp$Count_IdentifyPrimaryObjectsFB - Leftbp$Count_ColocalizedRegion
L_FB_Pneg<-L_Pneg/L_TFB

#Right
Rightbp<-read.table("Male_BS-127_111_Right_Lateral_DMV_Image", sep = ",", header = T)

#Total FB
R_TFB<-Rightbp$Count_IdentifyPrimaryObjectsFB
#Total PACAP
R_TP<-Rightbp$Count_IdentifyPrimaryObjectsPACAP
#Within FB+, how many are PACAP+?
R_FB_P<-Rightbp$Count_ColocalizedObjects / R_TFB
#Within FB+, how many are PACAP-?
R_Pneg<-Rightbp$Count_IdentifyPrimaryObjectsFB - Rightbp$Count_ColocalizedRegion
R_FB_Pneg<-R_Pneg/R_TFB

dfbp<-data.frame(Percentage = c(L_FB_P, L_FB_Pneg, R_FB_P, R_FB_Pneg), Side = c("Left", "Left", "Right", "Right"), Detection = c("PACAP+", "PACAP-", "PACAP+", "PACAP-") )
dfbp$Detection<-factor(dfbp$Detection, levels = c("PACAP+", "PACAP-"))
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

Left<-read.csv("Male_BS-127_111_Left_LateralDMVColocalizedObjects.csv", sep = ",")
length(Left$Intensity_MeanIntensity_PACAP)
Right<-read.csv("Male_BS-127_111_Right_Lateral_DMV_PACAPColocalizedObjects.csv", sep = ",")
length(Right$Intensity_MeanIntensity_PACAP)

df<-data.frame(Intensity = c(Left$Intensity_MeanIntensity_PACAP, Right$Intensity_MeanIntensity_PACAP), 
               Side = c(rep("Left", nrow(Left)), rep("Right", nrow(Right))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Side)) + geom_density(size =2) +
  #+ xlim(2,200) 
  scale_x_continuous(breaks=seq(0.02,0.07,by = 0.01),limits = c(0.02,0.07))+
  scale_color_manual(values=c("tan4","violet")) + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure 7D
Caudalbp<-read.table("Male_BS-127_85_Caudal_Left_DMV_PACAPImage.csv", sep = ",", header = T)
#Total FB
mC_TFB<-Caudalbp$Count_IdentifyPrimaryObjectsFB
#Total PACAP
mC_TP<-Caudalbp$Count_IdentifyPrimaryObjectsPACAP
#Within FB+, how many are PACAP+?
mC_FB_P<-Caudalbp$Count_ColocalizedRegion / mC_TFB

#Within FB+, how many are PACAP-?
mC_Pneg<-Caudalbp$Count_IdentifyPrimaryObjectsFB - Caudalbp$Count_ColocalizedRegion
mC_FB_Pneg<-mC_Pneg/mC_TFB


Intermediatebp<-read.table("Male_BS-127_111_Left_Lateral_DMV_PACAPImage.csv", sep = ",", header = T)
#Total FB
mI_TFB<-Intermediatebp$Count_IdentifyPrimaryObjectsFB
#Total PACAP
mI_TP<-Intermediatebp$Count_IdentifyPrimaryObjectsPACAP
#Within FB+, how many are PACAP+?
mI_FB_P<-Intermediatebp$Count_ColocalizedRegion / mI_TFB
#Within FB+, how many are PACAP-?
mI_Pneg<-Intermediatebp$Count_IdentifyPrimaryObjectsFB - Intermediatebp$Count_ColocalizedRegion
mI_FB_Pneg<-mI_Pneg/mI_TFB


Rostralbp<-read.table("Male_BS-127_136_Rostral_Left_Lateral_DMV_PACAPImage", sep = ",", header = T)
#Right
#Total FB
mR_TFB<-Rostralbp$Count_IdentifyPrimaryObjectsFB
#Total PACAP
mR_TP<-Rostralbp$Count_IdentifyPrimaryObjectsPACAP
#Within FB+, how many are PACAP+?
mR_FB_P<-Rostralbp$Count_ColocalizedRegion / mR_TFB
#Within FB+, how many are PACAP-?
mR_Pneg<-Rostralbp$Count_IdentifyPrimaryObjectsFB - Rostralbp$Count_ColocalizedRegion
mR_FB_Pneg<-mR_Pneg/mR_TFB

dfbp<-data.frame(Percentage = c(mC_FB_P, mC_FB_Pneg, mI_FB_P, mI_FB_Pneg, mR_FB_P, mR_FB_Pneg), Region = c("Caudal", "Caudal", "Intermediate", "Intermediate", "Rostral", "Rostral"), Detection = c("PACAP+", "PACAP-","PACAP+", "PACAP-", "PACAP+", "PACAP-") )
dfbp$Detection<-factor(dfbp$Detection, levels = c("PACAP+", "PACAP-"))
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

Caudal<-read.csv("Male_BS-127_85_Caudal_Left_DMV_PACAPColocalizedRegion.csv", sep = ",")
#Caudal <- Caudal[-3,]
length(Caudal$Intensity_MeanIntensity_PACAP)

Intermediate<-read.csv("Male_BS-127_111_Left_Lateral_DMV_PACAPColocalizedRegion.csv", sep = ",")
#Intermediate<-Intermediate[c(6:9),]
length(Intermediate$Intensity_MeanIntensity_PACAP)

Rostral<-read.csv("Male_BS-127_136_Rostral_Left_Lateral_DMV_PACAPColocalizedRegion.csv", sep = ",")
length(Rostral$Intensity_MeanIntensity_PACAP)

df<-data.frame(Intensity = c(Caudal$Intensity_MeanIntensity_PACAP, Intermediate$Intensity_MeanIntensity_PACAP, Rostral$Intensity_MeanIntensity_PACAP), 
               Region = c(rep("Caudal", nrow(Caudal)), rep("Intermediate", nrow(Intermediate)), rep("Rostral", nrow(Rostral))))

range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Region)) + geom_density(size = 2) +
  #+ xlim(2,200) 
  scale_x_continuous(breaks=seq(0.00,0.05,by = 0.01),limits = c(0.00,0.05))+
  scale_color_manual(values=c("purple","orange","green"), name = "Region") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure7F
#Left
Leftbp<-read.table("Female_BS124-90_Left_DMVImage.csv", sep = ",", header = T)

#Total FB
L_TFB<-Leftbp$Count_IdentifyPrimaryObjectsFB
#Total PACAP
L_TP<-Leftbp$Count_IdentifyPrimaryObjectsPACAP
#Within FB+, how many are PACAP+?
L_FB_P<-Leftbp$Count_ColocalizedObjects / L_TFB
#Within FB+, how many are PACAP-?
L_Pneg<-Leftbp$Count_IdentifyPrimaryObjectsFB - Leftbp$Count_ColocalizedObjects
L_FB_Pneg<-L_Pneg/L_TFB

#Right
Rightbp<-read.table("Female_BS-124_90_RightDMVImage.csv", sep = ",", header = T)

#Total FB
R_TFB<-Rightbp$Count_IdentifyPrimaryObjectsFB
#Total PACAP
R_TP<-Rightbp$Count_IdentifyPrimaryObjectsPACAP
#Within FB+, how many are PACAP+?
R_FB_P<-Rightbp$Count_ColocalizedObjects / R_TFB
#Within FB+, how many are PACAP-?
R_Pneg<-Rightbp$Count_IdentifyPrimaryObjectsFB - Rightbp$Count_ColocalizedObjects
R_FB_Pneg<-R_Pneg/R_TFB

dfbp<-data.frame(Percentage = c(L_FB_P, L_FB_Pneg, R_FB_P, R_FB_Pneg), Side = c("Left", "Left", "Right", "Right"), Detection = c("PACAP+", "PACAP-", "PACAP+", "PACAP-") )
dfbp$Detection<-factor(dfbp$Detection, levels = c("PACAP+", "PACAP-"))
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

Left<-read.csv("Female_BS-124_90_Left_DMVColocalizedObjects.csv", sep = ",")
length(Left$Intensity_MeanIntensity_PACAP)
Right<-read.csv("Female_BS-124_90_RightDMVColocalizedObjects.csv", sep = ",")
length(Right$Intensity_MeanIntensity_PACAP)

df<-data.frame(Intensity = c(Left$Intensity_MeanIntensity_PACAP, Right$Intensity_MeanIntensity_PACAP), 
               Side = c(rep("Left", nrow(Left)), rep("Right", nrow(Right))))
range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Side)) + geom_density(size =2) +
  #+ xlim(2,200) 
  scale_x_continuous(breaks=seq(0.02,0.06,by = 0.01),limits = c(0.02,0.06))+
  scale_color_manual(values=c("tan4","violet")) + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))

#Figure 7H
setwd("Data/IF_quantification/Figure 7/Rostral_Int_Caudal")

Caudalbp<-read.table("Female_BS-124_24_Caudal_Left_DMV_PACAPImage.csv", sep = ",", header = T)
#Total FB
mC_TFB<-Caudalbp$Count_IdentifyPrimaryObjectsFB
#Total PACAP
mC_TP<-Caudalbp$Count_IdentifyPrimaryObjectsPACAP
#Within FB+, how many are PACAP+?
mC_FB_P<-Caudalbp$Count_ColocalizedRegion / mC_TFB

#Within FB+, how many are PACAP-?
mC_Pneg<-Caudalbp$Count_IdentifyPrimaryObjectsFB - Caudalbp$Count_ColocalizedRegion
mC_FB_Pneg<-mC_Pneg/mC_TFB

Intermediatebp<-read.table("Female_BS-124_80_Left_Medial_DMV_PACAPImage.csv", sep = ",", header = T)
#Total FB
mI_TFB<-Intermediatebp$Count_IdentifyPrimaryObjectsFB
#Total PACAP
mI_TP<-Intermediatebp$Count_IdentifyPrimaryObjectsPACAP
#Within FB+, how many are PACAP+?
mI_FB_P<-Intermediatebp$Count_ColocalizedRegion / mI_TFB
#Within FB+, how many are PACAP-?
mI_Pneg<-Intermediatebp$Count_IdentifyPrimaryObjectsFB - Intermediatebp$Count_ColocalizedRegion
mI_FB_Pneg<-mI_Pneg/mI_TFB

Rostralbp<-read.table("Female_BS-124_132_Rostral_Left_DMV__PACAPImage.csv", sep = ",", header = T)
#Right
#Total FB
mR_TFB<-Rostralbp$Count_IdentifyPrimaryObjectsFB
#Total PACAP
mR_TP<-Rostralbp$Count_IdentifyPrimaryObjectsPACAP
#Within FB+, how many are PACAP+?
mR_FB_P<-Rostralbp$Count_ColocalizedRegion / mR_TFB
#Within FB+, how many are PACAP-?
mR_Pneg<-Rostralbp$Count_IdentifyPrimaryObjectsFB - Rostralbp$Count_ColocalizedRegion
mR_FB_Pneg<-mR_Pneg/mR_TFB

dfbp<-data.frame(Percentage = c(mC_FB_P, mC_FB_Pneg, mI_FB_P, mI_FB_Pneg, mR_FB_P, mR_FB_Pneg), Region = c("Caudal", "Caudal", "Intermediate", "Intermediate", "Rostral", "Rostral"), Detection = c("PACAP+", "PACAP-","PACAP+", "PACAP-", "PACAP+", "PACAP-") )
dfbp$Detection<-factor(dfbp$Detection, levels = c("PACAP+", "PACAP-"))
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

Caudal<-read.csv("Female_BS-124_24_Caudal_Left_DMV_PACAPColocalizedRegion.csv", sep = ",")
length(Caudal$Intensity_MeanIntensity_PACAP)

Intermediate<-read.csv("Female_BS-124_80_Left_Medial_DMV_PACAPColocalizedRegion.csv", sep = ",")
length(Intermediate$Intensity_MeanIntensity_PACAP)

Rostral<-read.csv("Female_BS-124_132_Rostral_Left_DMV_PACAPColocalizedRegion.csv", sep = ",")
length(Rostral$Intensity_MeanIntensity_PACAP)

df<-data.frame(Intensity = c(Caudal$Intensity_MeanIntensity_PACAP, Intermediate$Intensity_MeanIntensity_PACAP, Rostral$Intensity_MeanIntensity_PACAP), 
               Region = c(rep("Caudal", nrow(Caudal)), rep("Intermediate", nrow(Intermediate)), rep("Rostral", nrow(Rostral))))

range(df$Intensity)

ggplot(df, aes(x=Intensity, col = Region)) + geom_density(size = 2) +
  #+ xlim(2,200) 
  scale_x_continuous(breaks=seq(0.02,0.09,by = 0.01),limits = c(0.02,0.09))+
  scale_color_manual(values=c("purple","orange","green"), name = "Region") + 
  theme_minimal() +
  theme(axis.line = element_line(color="black"),
        text = element_text(size=26),
        axis.ticks = element_line(color="black"))