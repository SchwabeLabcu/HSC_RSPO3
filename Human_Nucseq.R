library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

setwd("D:/projects/Human_Nucseq")

## preprocessing from scrublet doublet removed file

# reading each sample, creating a seurat object and then a list of seurat objects 

loc <- "CellBender/"
dir <- list.dirs(path = loc,full.names = F, recursive = F)

# read the dataset from a csv file
df <- data.frame(read.csv(file = "humNucSeqSampleInfo.csv"))

seurat.list <- vector(mode = "list", length = nrow(df)) 
names(seurat.list) <- df$Name

#create seurat list with all samples
i=1
for (i in 1:nrow(df)) { 
  
    sampleName <- df$Name[i]
    datLoc <- file.path("CellBender",sampleName,"CB_raw_feature_bc_matrix_filtered_scrublet.h5")
    dat <- Read10X_h5(datLoc)
    seurat <- CreateSeuratObject(counts = dat, project = sampleName, min.cells = 0, min.features = 0)
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
    seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RPL|^RPS|MRPS|^MRPL")
    seurat <- subset(seurat, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 6500 & nCount_RNA < 40000 )
    data <- rbind(data,data_seurat)
    seurat.list[i] <- seurat
}

#select reference samples for RPCA
tmp <- names(seurat.list)
req_samples <-df$Name
req_sample_pos <- which(tmp %in% req_samples)

# reference samples
ref_samples <- c("S4","S5","S12","S15","S21","S23","S28","S33","S35","S36")
ref_sample_pos <- which(req_samples %in% ref_samples)

# normalize and identify variable features for each dataset independently
log.data <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

features <- SelectIntegrationFeatures(object.list = log.data)

log.data <- lapply(X = log.data, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T)
})

#perform anchor
anchors <- FindIntegrationAnchors(object.list = log.data, reference = ref_sample_pos, anchor.features = features, reduction = "rpca", k.anchor = 10)

log.data <- IntegrateData(anchorset = anchors)

DefaultAssay(log.data) <- "integrated"

# Run the standard workflow for visualization and clustering
log.data <- ScaleData(log.data, verbose = FALSE)
log.data <- RunPCA(log.data, npcs = 50, verbose = FALSE)
ElbowPlot(log.data, ndims = 50)
log.data <- RunUMAP(log.data, reduction = "pca", dims = 1:30)
log.data <- FindNeighbors(log.data, reduction = "pca", dims = 1:30)
log.data <- FindClusters(log.data, resolution = 0.5)

DimPlot(log.data, reduction = "umap", label = F, group.by = "orig.ident", raster = F)

# Visualization
DimPlot(log.data, reduction = "umap", label = TRUE, repel = TRUE, raster = F)

DefaultAssay(log.data) <- "RNA"

## Cell clustering 
log.data <- SetIdent(log.data, value = "orig.ident")
log.data[["Patient.ID"]] <- Idents(object = log.data)

log.data <- RenameIdents(object = log.data, 
                         "S1" = "Alcoholic Cirrhosis",
                         "S2" = "Alcoholic Cirrhosis",
                         "S3" = "Alcoholic Cirrhosis",
                         "S4" = "Alcoholic Cirrhosis",
                         
                         "S5" = "Alcoholic Hepatitis",
                         "S6" = "Alcoholic Hepatitis",
                         "S7" = "Alcoholic Hepatitis",
                         "S8" = "Alcoholic Hepatitis",
                         "S9" = "Alcoholic Hepatitis",
                         
                         "S12" = "Normal",
                         "S13" = "Normal",
                         "S14" = "Normal",
                         "S23" = "Normal",
                         "S27" = "Normal",
                         "S28" = "Normal",
                         
                         "S15" = "Cirhosis (MASH)",
                         "S16" = "Cirhosis (MASH)",
                         "S17" = "Cirhosis (MASH)",
                         "S18" = "Cirhosis (MASH)",
                         
                         "S19" = "Normal (with MASH)",
                         "S20" = "Normal (with MASH)",
                         "S21" = "Normal (with MASH)",
                         
                         "S33" = "Intermediate (MASH)",
                         "S35" = "Intermediate (MASH)",
                         "S36" = "Intermediate (MASH)",
                         "S38" = "Intermediate (MASH)"

)

log.data[["Fibrosis.Group"]] <- Idents(object = log.data)

log.data <- JoinLayers(log.data)

DefaultAssay(log.data) <- "RNA"
# Run the standard workflow for visualization and clustering on not batch corrected data
log.data <- ScaleData(log.data, verbose = FALSE)
log.data <- RunPCA(log.data, npcs = 50, verbose = FALSE)
ElbowPlot(log.data, ndims = 50)
log.data <- RunUMAP(log.data, reduction = "pca", dims = 1:40)
log.data <- FindNeighbors(log.data, reduction = "pca", dims = 1:30)
log.data <- FindClusters(log.data, resolution = 0.5)

DimPlot(log.data, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)
DimPlot(log.data, reduction = "umap", label = T, group.by = "orig.ident", raster = F)


Cell.type.new <- c("VSMC","HSC","PF","Doublets", "Resident NK", "Circulating NK/NKT",
                   "T Cell", "plasma cell", "B cell", "cDC1", "cDC2", "pDC",
                   "Kupffer Cell", "Monocyte", "Neutrophil","Basophil",
                   "Endothelial Cell", "Cholangiocytes","Cholangiocytes with Hepatocytes Markers", 
                   "Hepatocytes", "Hepatocytes with Cholangiocytes Markers", "Proliferating Cells")

markerGenes <- c("ALB","APOE", #hepatocyte
                 "ANXA4","EPCAM", "CFTR",#cholangiocyte
                 "COL4A1","LAMA2","CCBE1","HGF","RSPO3", #HSC
                 "PTPRB","PECAM1","FLT1", #endothelial
                 'UPK3B','MSLN', #methothelial
                 "IL7R","THEMIS","FYN", #T cell
                 "TOP2A","MKI67", #proliferating cell
                 "DERL3","TNFRSF17","MZB1", #plasma cell
                 "CD22","FCRL1", #B cell
                 "CLEC9A","XCR1", #cDC1
                 "FCER1A", #cDC2
                 "LILRA4", #pDC
                 "XCL1","NCR1", # Resident NK
                 "FGFBP2","GNLY", # Circulating NK/NKT
                 "VSIG4","CD86", # Kupffer cell
                 "FCN1","APOBEC3A" # Monocyte
) 


markerGenesFib <- c( "EBF1","SLIT3","RGS5",#VSMC
                     "COL25A1","SYT1","RELN","ADAMTSL1", #qHSC
                     "NR1H4","HGF","CCBE1","COL4A1","LAMA2", #InterHSC
                     "VCAN","KCND2","FN1","TIMP1","COL1A1", #aHSC
                     "PER2","NFATC1","ADAMTS3","NFATC2","ICAM1", #PF
                     "C3","PRG4","CLDN1")# Meso

markerGenesmyeloid <- c(
  "TNFRSF17","MZB1", #plasma cell
  "CD22","FCRL1", #B cell
  "CLEC9A","XCR1", #cDC1
  "FCER1A", "CD1C",#cDC2
  "LILRA4","SCT", #pDC
  "VSIG4","NDST3",  # Kupffer cell
  "VCAN","CD300E",# Monocyte
  "APOBEC3A", "FCGR3B", # Neutrophil
  "TPSAB1","CPA3" #Basophil
) 

markerGenesTCell <- c(
  "IL7R","THEMIS", #T cell
  "XCL1","NCR1", # Resident NK
  "FGFBP2","GNLY"# Circulating NK/NKT
) 

# CELL TYPES DOT PLOT
DotPlot(log.data, features = markerGenes) + RotatedAxis() 

## Cell clustering 
log.data <- SetIdent(log.data, value = "seurat_clusters")
log.data <- RenameIdents(object = log.data, 
                         "3" = "T Cell/NK",
                         "23" = "T Cell/NK",
                         
                         "2" = "Myeloid Cell",
                         "21" = "Myeloid Cell",
                         "22" = "Myeloid Cell",
                         
                         "4" = "HSC/Fib",     
                         "19" = "HSC/Fib",  
                         "24" = "HSC/Fib",   
                         
                         "11" = "Endothelial Cell",  
                         "13" = "Endothelial Cell",  
                         
                         "0" = "Cholangiocytes",  
                         "18" = "Cholangiocytes" ,
                         
                         "1" = "Hepatocytes" ,
                         "5" = "Hepatocytes" ,
                         "6" = "Hepatocytes" ,
                         "8" = "Hepatocytes" ,
                         "9" = "Hepatocytes" ,
                         "10" = "Hepatocytes" ,
                         "12" = "Hepatocytes" ,
                         "15" = "Hepatocytes" ,
                         "16" = "Hepatocytes" ,
                         "17" = "Hepatocytes" ,
                         "20" = "Hepatocytes" ,
                         
                         "7" = "Hepatocytes" ,
                         "14" = "Hepatocytes" ,
                         
                         "25" = "Doublets" ,
                         "26" = "Proliferating Cells" 
                         
)

log.data[["Cell.type"]] <- Idents(object = log.data)

DimPlot(log.data, reduction = "umap", label = F, group.by = "Cell.type", raster = F)
DotPlot(log.data, features = markerGenes) + RotatedAxis() 


log.data <- SetIdent(log.data, value = "Cell.type")

#subset, assign cell type and transfer label back to the old population

#myeloid cells
subset <- subset(log.data,idents = c("Myeloid Cell"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:30)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:30)
subset <- FindClusters(subset, resolution = 0.5)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)

## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       "7" = "plasma cell",
                       
                       "10" = "B cell",
                       
                       "12" = "cDC1",
                       
                       "5" = "cDC2",
                       
                       "13" = "pDC",
                       
                       "1" = "Kupffer Cell",
                       "3" = "Kupffer Cell",
                       "8" = "Kupffer Cell",
                       "9" = "Kupffer Cell",
                       
                       "2" = "Monocyte",
                       "6" = "Monocyte",
                       "0" = "Monocyte",
                       
                       "4" = "Neutrophil",
                       "17" = "Neutrophil",
                       
                       "15" = "Basophil",
                       
                       "11" = "Doublets",                        
                       "14" = "Doublets" ,
                       
                       "16" = "Doublets",
                       "18" = "Doublets"
)

subset[["Cell.type.new"]] <- Idents(object = subset)

Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)


#T/NK
subset <- subset(log.data,idents = c("T Cell/NK"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:30)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:30)
subset <- FindClusters(subset, resolution = 0.5)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)

## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       "5" = "Resident NK",
                       "11" = "Circulating NK/NKT",
                       
                       "0" = "T Cell",
                       "1" = "T Cell",
                       "2" = "T Cell",
                       "3" = "T Cell",
                       "4" = "T Cell",
                       "6" = "T Cell",
                       "7" = "T Cell",
                       "8" = "T Cell",
                       "9" = "T Cell",
                       "10" = "T Cell",
                       
                       "14" = "Doublets",
                       "12" = "Doublets",
                       "13" = "Doublets",
                       "15" = "Doublets"
                       
)

subset[["Cell.type.new"]] <- Idents(object = subset)

Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)


#HSC
subset <- subset(log.data,idents = c("HSC/Fib"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:13)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:13)
subset <- FindClusters(subset, resolution = 0.5)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)

## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       
                       "8" = "VSMC",
                       
                       "0" = "HSC",
                       "1" = "HSC",
                       "2" = "HSC",
                       "3" = "HSC",
                       "4" = "HSC",
                       "5" = "HSC",
                       "6" = "HSC",
                       "10" = "HSC",
                       "12" = "HSC",
                       
                       "7" = "HSC",
                       
                       "9" = "Doublets",
                       "11" = "Doublets",
                       "13" = "Doublets",
                       "14" = "Doublets"
                       
)

subset[["Cell.type.new"]] <- Idents(object = subset)

Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)

#Endothelial
subset <- subset(log.data,idents = c("Endothelial Cell"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:13)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:13)
subset <- FindClusters(subset, resolution = 0.3)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)

## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       
                       "0" = "Endothelial Cell",
                       "1" = "Endothelial Cell",
                       "2" = "Endothelial Cell",
                       "3" = "Endothelial Cell",
                       "4" = "Endothelial Cell",
                       "5" = "Endothelial Cell",
                       "6" = "Endothelial Cell",
                       "7" = "Endothelial Cell",
                       "10" = "Endothelial Cell",
                       
                       "8" = "Doublets",
                       "9" = "Doublets",
                       "11" = "Doublets",
                       "12" = "Doublets"
                       
)

subset[["Cell.type.new"]] <- Idents(object = subset)

Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)

#Cholangiocytes
subset <- subset(log.data,idents = c("Cholangiocytes"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:13)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:13)
subset <- FindClusters(subset, resolution = 0.3)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)

## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       
                       "8" = "Cholangiocytes with Hepatocytes Markers",

                       "0" = "Cholangiocytes",
                       "1" = "Cholangiocytes",
                       "2" = "Cholangiocytes",
                       "3" = "Cholangiocytes",
                       "4" = "Cholangiocytes",
                       "5" = "Cholangiocytes",
                       "6" = "Cholangiocytes",
                       "7" = "Cholangiocytes",
                       "9" = "Cholangiocytes",
                       "10" = "Cholangiocytes"
                       
)

subset[["Cell.type.new"]] <- Idents(object = subset)

Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)


#hep
subset <- subset(log.data,idents = c("Hepatocytes"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:15)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:15)
subset <- FindClusters(subset, resolution = 0.3)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)

## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       
                       "11" = "Doublets",
                       
                       "0" = "Hepatocytes",
                       "1" = "Hepatocytes",
                       "2" = "Hepatocytes",
                       "3" = "Hepatocytes",
                       "4" = "Hepatocytes",
                       "5" = "Hepatocytes",
                       "6" = "Hepatocytes",
                       "7" = "Hepatocytes",
                       "8" = "Hepatocytes",
                       "9" = "Hepatocytes",
                       "10" = "Hepatocytes",
                       "12" = "Hepatocytes",
                       "13" = "Hepatocytes"
                       
)

subset[["Cell.type.new"]] <- Idents(object = subset)


subset2 <- subset(subset, idents = "Doublets")

subset2 <- ScaleData(subset2, verbose = FALSE)
subset2 <- RunPCA(subset2, npcs = 50, verbose = FALSE)
subset2 <- RunUMAP(subset2, reduction = "pca", dims = 1:13)
subset2 <- FindNeighbors(subset2, reduction = "pca", dims = 1:13)
subset2 <- FindClusters(subset2, resolution = 0.3)

DotPlot(subset2, features = markerGenes) + RotatedAxis() 

DimPlot(subset2, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)
FeaturePlot(subset2, features =c("ALB","RSPO3","CFTR","THEMIS"),order = T, raster = F)


## Cell clustering 
subset2 <- SetIdent(subset2, value = "seurat_clusters")
subset2 <- RenameIdents(object = subset2, 
                        
                        "2" = "Hepatocytes with Cholangiocytes Markers",
                        
                        "0" = "Doublets",
                        "1" = "Doublets",
                        "3" = "Doublets",
                        "4" = "Doublets",
                        "5" = "Doublets",
                        "6" = "Doublets",
                        "7" = "Doublets"
                        
)

subset2[["Cell.type.new"]] <- Idents(object = subset2)

#remove doublet
Idents(subset) <- Idents(subset2)
subset[["Cell.type.new"]] <- Idents(object = subset)


Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)

DimPlot(log.data,label = F, group.by = "Cell.type.new",raster = F)

log.data <- subset(log.data, ident = c("VSMC","HSC", "Resident NK", "Circulating NK/NKT",
                                       "T Cell", "plasma cell", "B cell", "cDC1", "cDC2", "pDC",
                                       "Kupffer Cell", "Monocyte", "Neutrophil","Basophil",
                                       "Endothelial Cell", "Cholangiocytes","Cholangiocytes with Hepatocytes Markers", 
                                       "Hepatocytes", "Hepatocytes with Cholangiocytes Markers"))


# plot with cell numbers
cell_counts <- table(log.data$Cell.type)
clusterlabels = paste(names(cell_counts),"(",cell_counts,sep = "")
DimPlot(log.data, label = T, group.by = "Cell.type",raster = F) + scale_color_discrete(label = clusterlabels)

saveRDS(log.data, file = "humnucseq_dim30.rds")

##
## generating mean expression for analysis
# read the dataset from a csv file
df <- data.frame(read.csv(file = "humNucSeqSampleInfo.csv"))
data_all <- data.frame(row.names = c("sampleName", "RSPO3_HSC","CYP2E1_Hep","TBX3_Hep","CYP1A2_Hep"))

log.data <- JoinLayers(log.data)
#create seurat list with all samples
i=1
for (i in 1:nrow(df)) { 
  sampleName <- df$Name[i]
  datLoc <- file.path(sampleName)
  patient <- subset(log.data,idents = datLoc)
  patient <- SetIdent(patient, value = "Cell.type.new")
  hepatocyte <- subset(patient, idents = "Hepatocyte")
  HSC <- subset(patient, idents = "HSC")
  RSPO3_HSC <- mean(HSC[["RNA"]]$data['RSPO3',])
  CYP2E1_Hep <- mean(hepatocyte[["RNA"]]$data['CYP2E1',])
  TBX3_Hep <- mean(hepatocyte[["RNA"]]$data['TBX3',])
  CYP1A2_Hep <- mean(hepatocyte[["RNA"]]$data['CYP1A2',])  
  data <- c(sampleName, RSPO3_HSC,CYP2E1_Hep,TBX3_Hep,CYP1A2_Hep)
  data_all <- rbind(data_all, data)
  print(data)
}

write.csv(data_all, "data.csv")

## identify cymy hsc

library(readxl)
dataset <- data.frame(read_excel("cymysigs.xlsx"))
head(dataset)

cySig <- na.omit(dataset$cysigs_human)
mySig <- na.omit(dataset$mysigs_human)


sig <- cySig
sigName <- "cySig"
callName <- paste0(sigName,"1")
req_subset <- AddModuleScore(log.data,features = list(sig), name = sigName)
p1 <- FeaturePlot(req_subset, callName, max.cutoff = "q99", pt.size = 1) + scale_color_viridis_c(option = "viridis")
VlnPlot(req_subset, callName, sort = T)
# VlnPlot(req_subset, callName, group.by = "Cell.clusters", sort = T)

sig <- mySig
sigName <- "mySig"
callName <- paste0(sigName,"1")
req_subset <- AddModuleScore(req_subset,features = list(sig), name = sigName)
p2 <- FeaturePlot(req_subset, callName, max.cutoff = "q99", pt.size = 1) + scale_color_viridis_c(option = "viridis")
VlnPlot(req_subset, callName, sort = T)
# VlnPlot(req_subset, callName, group.by = "Cell.clusters", sort = T)


VlnPlot(req_subset, features = c("cySig1","mySig1"), sort = F)
DotPlot(req_subset, features = c("cySig1","mySig1")) + RotatedAxis()
p3 <- FeaturePlot(req_subset, "RSPO3", max.cutoff = "q99", pt.size = 1, order = F) 
p4 <- FeaturePlot(req_subset, "HHIP", max.cutoff = "q99", pt.size = 1, order = F) 
p5 <- FeaturePlot(req_subset, "COL15A1", max.cutoff = "q99", pt.size = 1, order = F) 
p6 <- FeaturePlot(req_subset, c("TIMP1"),  max.cutoff = "q99", pt.size = 1,order = F) 

p <- p1 + p2 + p3 + p4 + p5 + p6

ggsave("RSPO3_cymy.tiff", width = 5, height = 25, units = "in", dpi = 300)

req_subset <- SetIdent(req_subset, value = "seurat_clusters")
df <- data.frame(clust=req_subset@active.ident, cy=req_subset@meta.data$cySig1, my=req_subset@meta.data$mySig1)

library(tidyverse)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE) )
  } 
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
  return(data_sum)
}

tmp <- data_summary(df, varname = "cy", groupnames = "clust")
tmp
cymyscore <- tmp[,c(1,2)]
tmp <- data_summary(df, varname = "my", groupnames = "clust")
tmp
cymyscore <- cbind(cymyscore, my=tmp[,c(2)], cy_my=(cymyscore[,2]-tmp[,c(2)]))
cymyscore

# req_subset[["Cell.clusters.HSC"]] <- Idents(object = req_subset)
DimPlot(req_subset, reduction = "umap",group.by = "seurat_clusters", label = T)
req_subset <- SetIdent(req_subset, value = "seurat_clusters")
req_subset <- RenameIdents(object = req_subset,
                           '0' = "cyHSC",
                           '2' = "cyHSC",
                           '3' = "cyHSC",
                           '6' = "cyHSC",
                           '8' = "cyHSC",
                           
                           
                           '1' = "myHSC",
                           '4' = "myHSC",
                           '5' = "myHSC",
                           '7' = "myHSC",
                           '12' = "myHSC"
                           )

req_subset[["Cell.types.cymy"]] <- Idents(object = req_subset)
p6 <- DimPlot(req_subset,label = T, group.by = "Cell.types.cymy")
table(req_subset@active.ident)

cols_cymy <- c("#FF9900","#000099")
p6 <- DimPlot(req_subset, reduction = "umap",group.by = "Cell.types.cymy", pt.size = 1,cols = cols_cymy, label = F)

p <- p6 + p3 + p4 + p5 
ggsave("RSPO3_cymy.tiff", width = 10, height = 9, units = "in", dpi = 300)


saveRDS(req_subset, "humnucseq_HSC_cymy.rds")

#'
#'
#'pseudobulk analysis


log.data <- readRDS("RDS_File/02_humnucseq_dim30_doublet_removed.rds")
log.data <- SetIdent(log.data, value = "orig.ident")

# read the dataset from a csv file
df <- data.frame(read.csv(file = "humNucSeqSampleInfo.csv"))
HSC_bulk_all <- data.frame(row.names = row.names(log.data))
Hep_bulk_all <- data.frame(row.names = row.names(log.data))

log.data <- JoinLayers(log.data)
#create seurat list with all samples
i=1
for (i in 1:nrow(df)) { 
  print(sampleName)
  sampleName <- df$Name[i]
  datLoc <- file.path(sampleName)
  patient <- subset(log.data,idents = datLoc)
  patient <- SetIdent(patient, value = "Cell.type.new")
  hepatocyte <- subset(patient, idents = "Hepatocyte")
  HSC <- subset(patient, idents = "HSC")
  
  hep_mat <- hepatocyte[["RNA"]]$counts
  hep_mat <- as.matrix(hep_mat)
  hep_mat_norm <- t(t(hep_mat)/(colSums(hep_mat)/1e6))
  hep_bulk <- rowMeans(hep_mat_norm)
  
  HSC_mat <- HSC[["RNA"]]$counts
  HSC_mat <- as.matrix(HSC_mat)
  HSC_mat_norm <- t(t(HSC_mat)/(colSums(HSC_mat)/1e6))
  HSC_bulk <- rowMeans(HSC_mat_norm)
  
  Hep_bulk_all[[sampleName]] <- hep_bulk
  HSC_bulk_all[[sampleName]] <- HSC_bulk
}

write.csv(Hep_bulk_all, "Hep_bulk_norm.csv")
write.csv(HSC_bulk_all, "HSC_bulk_norm.csv")