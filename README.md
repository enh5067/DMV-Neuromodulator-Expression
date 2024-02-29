# DMV-Neuromodulator-Expression
## This repository contains all the code and data needed to reproduce the results of the paper “Neuromodulatory co-expression in cardiac vagal motor neurons of the Dorsal Motor Nucleus of the Vagus".

### This repository contains 6 files.

•	ReplicateFigures.R—R script for replicating the main figures present in the manuscript, including visualizing Visium sections colored according to spatial transcriptomic phenotype, uniform manifold approximation and projection (UMAP) plots, raw spatial transcriptomic 
  gene expression counts shown in the context of the tissue, kernel density estimates of mean protein intensity level, bar plots of proportions of neurons, as well as violin plots, heat maps, boxplots, and scatterplots of relative expression.

•	ReplicateFigures_Supplemental.R—R script for replicating the supplemental figures present in the manuscript.

•	Annotations—folder containing all annotation files necessary to replicate the figures.
  
    o	10xsc-Sample_Annotations.txt—sample annotations for the 10x genomics single cell RNAseq data stored in 10xsc_SeuratObject.rds.
  
    o	LCM-qPCR_Gene_Annotations.txt—gene annotations for the HT-qPCR data, LCM-qPCR_Negddct.txt.
  
    o	LCM-qPCR_Sample_Annotations.txt—sample annotations for the HT-qPCR data, LCM-qPCR_Negddct.txt.
  
    o	LCMseq-Sample_Annotations.txt—sample annotations for the LCM-RNAseq data stored in LCMseq_SeuratObject.rds.


