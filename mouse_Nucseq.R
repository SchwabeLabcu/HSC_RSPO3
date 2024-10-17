library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scCustomize)

setwd("D:/projects/mousenucseq")

# read the dataset from a csv file
df <- data.frame(read.csv("sampleinfo.csv"))

seurat.list <- vector(mode = "list", length = nrow(df)) 
names(seurat.list) <- df$Name

#create seurat list with all samples
i=1
for (i in 1:nrow(df)) { 
  
    sampleName <- df$Name[i]
    
    datLoc <- file.path("CellBender",sampleName,"CB_raw_feature_bc_matrix_filtered.h5")
    dat <- Read_CellBender_h5_Mat(datLoc)
    
    seurat <- CreateSeuratObject(counts = dat, project = sampleName, min.cells = 0, min.features = 0)
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
    seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^Rpl|^Rps|Mrps|^Mrpl")

    seurat <- subset(seurat, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA < 60000 )
    seurat.list[i] <- seurat
}

# normalize and identify variable features for each dataset independently
log.data <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = log.data)

log.data <- lapply(X = log.data, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T)
})

#perform anchor
anchors <- FindIntegrationAnchors(object.list = log.data, reference = NULL, anchor.features = features, reduction = "rpca", k.anchor = 10)

# this command creates an 'integrated' data assay
log.data <- IntegrateData(anchorset = anchors)

DefaultAssay(log.data) <- "integrated"

# Run the standard workflow for visualization and clustering
log.data <- ScaleData(log.data, verbose = FALSE)
log.data <- RunPCA(log.data, npcs = 50, verbose = FALSE)
ElbowPlot(log.data, ndims = 50)
log.data <- RunUMAP(log.data, reduction = "pca", dims = 1:30)
log.data <- FindNeighbors(log.data, reduction = "pca", dims = 1:30)
log.data <- FindClusters(log.data, resolution = 0.5)

DimPlot(log.data, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)

cell_counts <- table(log.data$seurat_clusters)
clusterlabels = paste(names(cell_counts),"(",cell_counts,sep = "")
DimPlot(log.data, label = F, group.by = "seurat_clusters",raster = F) + scale_color_discrete(label = clusterlabels)

DefaultAssay(log.data) <- "RNA"

saveRDS(log.data, "01_musnucseq_merged_CB.rds")

VlnPlot(log.data, features = "nCount_RNA", pt.size = 0)

markerGenes <- c("Apoe","Alb", #hepatocyte
                  "Cftr",#cholangiocyte
                 "Lama2",#HSC
                 "Ebf1", #VSMC
                 "Ptprb", #endothelial
                 'Upk3b', #methothelial
                 "Il7r", #T cell
                 "Mzb1", #plasma cell
                 "Cd22", #B cell
                 "Clec9a", #cDCs
                 "Xcl1", # NK
                 "Vsig4" # Kupffer cell
) 

markerGenes1 <- c("Rspo3","Cyp2e1","Cyp2f2","Lgr4",
                  "Lgr5","Lgr6"
) 


# CELL TYPES DOT PLOT
DotPlot(log.data, features = markerGenes) + RotatedAxis() 

## Cell clustering 
log.data <- SetIdent(log.data, value = "seurat_clusters")
log.data <- RenameIdents(object = log.data, 
                         "15" = "T Cell/NK",
                         "21" = "T Cell/NK",
                         
                         "18" = "Myeloid Cell",
                         "19" = "Myeloid Cell",
                         "9" = "Myeloid Cell",
                         "13" = "Myeloid Cell",
                         
                         "11" = "HSC/Fib",     
                         "17" = "HSC/Fib",  
                         
                         "2" = "Endothelial Cell",  
                         "20" = "Endothelial Cell",  
                         
                         "14" = "Cholangiocytes",  
                         
                         "0" = "Hepatocytes" ,
                         "1" = "Hepatocytes" ,
                         "3" = "Hepatocytes" ,
                         "4" = "Hepatocytes" ,
                         "5" = "Hepatocytes" ,
                         
                         "7" = "Hepatocytes" ,
                         "8" = "Hepatocytes" ,
                         
                         "6" = "low quality Hepatocytes" ,
                         "12" = "low quality Hepatocytes" ,
                         
                         "10" = "Doublets" ,
                         "16" = "Unknown" ,
                         "22" = "Proliferating Cells" 
                         
)

log.data[["Cell.type"]] <- Idents(object = log.data)

DimPlot(log.data, reduction = "umap", label = F, group.by = "Cell.type", raster = F)
ggsave("Figures/cell.type.jpeg",dpi = 300, units = "in", width = 7, height = 4.8)


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

cell_counts <- table(subset$seurat_clusters)
clusterlabels = paste(names(cell_counts),"(",cell_counts,sep = "")
DimPlot(subset, label = F, group.by = "seurat_clusters",raster = F) + scale_color_discrete(label = clusterlabels)

markerGenesmyeloid <- c(
  "Tnfrsf17","Mzb1", #plasma cell
  "Cd22","Fcrl1",#B cell
  "Clec9a","Xcr1", #cDCs
  "Vsig4","Cd5l",  # Kupffer cell
  "Vcan","Chil3",# Monocyte
  "Arhgap22","H2-Ab1" #other myeloid
) 

subset <- SetIdent(subset, value = "Cell.type.new")


DotPlot(subset, features = markerGenes) + RotatedAxis() 
DotPlot(subset, features = markerGenesmyeloid) + RotatedAxis() 
VlnPlot(subset, features = "nCount_RNA", pt.size = 0)
FeaturePlot(subset, features = c("Arhgap22","Chrm3"))

## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       "14" = "plasma cell",
                       
                       "6" = "B cell",
                       
                       "12" = "cDCs",
                       
                       
                       "0" = "Kupffer Cell",
                       "1" = "Kupffer Cell",
                       "2" = "Kupffer Cell",
                       "3" = "Kupffer Cell",
                       
                       "13" = "Monocyte",
                       
                       "9" = "Other Myeloid cell",

                       "4" = "Doublets",                        
                       "5" = "Doublets" ,
                       
                       "10" = "Doublets",                        
                       "11" = "Doublets" ,
                       
                       "7" = "Doublets",
                       "8" = "Doublets"
)

subset[["Cell.type.new"]] <- Idents(object = subset)

DimPlot(subset, reduction = "umap", label = F, group.by = "Cell.type.new", raster = F)
ggsave("Figures/cell.type.new.myeloid.jpeg",dpi = 300, units = "in", width = 6, height = 4.8)

DotPlot(subset, features = markerGenesmyeloid) + RotatedAxis() 
ggsave("Figures/cell.type.dotplot.myeloid.jpeg",dpi = 300, units = "in", width = 12, height = 4.8)

FeaturePlot(subset, features =c(markerGenesmyeloid),order = T, raster = F)
ggsave("Figures/myeloidexpression.jpeg",dpi = 300, units = "in", width = 16, height = 12)

plot <- DimPlot(subset, reduction = "umap", cells = unlist(CellsByIdentities(subset)))
select.cells <- CellSelector(plot = plot)
head(select.cells)
Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)

DimPlot(log.data,label = F, group.by = "Cell.type.new",raster = F)
ggsave("Figures/cell.type.new.jpeg",dpi = 300, units = "in",  width = 6, height = 4.8)

saveRDS(log.data, file = "02_musnucseq_Celltype.rds")
saveRDS(subset, file = "11_musnucseq_myeloid.rds")


#T/NK
subset <- subset(log.data,idents = c("T Cell/NK"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:30)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:30)
subset <- FindClusters(subset, resolution = 0.5)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)

markerGenesTCell <- c(
  "Klra4","Cma1", "Ncr1",#NK
  "Cd3e","Il7r","Themis" # T
) 

DotPlot(subset, features = markerGenes) + RotatedAxis() 
DotPlot(subset, features = markerGenesTCell) + RotatedAxis() 
VlnPlot(subset, features = "nCount_RNA", pt.size = 0)
FeaturePlot(subset, features = c("Arhgap22","Chrm3"))


## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       "1" = "NK",
                       
                       "0" = "T Cell",
                       "3" = "T Cell",
                       
                       "4" = "cDCs",
 
                       "6" = "Doublets",
                       "2" = "Doublets",
                       "5" = "Doublets"
                       
)

subset[["Cell.type.new"]] <- Idents(object = subset)

DimPlot(subset, reduction = "umap", label = F, group.by = "Cell.type.new", raster = F)
ggsave("Figures/cell.type.new.TNK.jpeg",dpi = 300, units = "in", width = 6, height = 4.8)

DotPlot(subset, features = markerGenesTCell) + RotatedAxis() 
ggsave("Figures/cell.type.dotplot.TNK.jpeg",dpi = 300, units = "in", width = 12, height = 4.8)

FeaturePlot(subset, features =c(markerGenesTCell),order = T, raster = F)
ggsave("Figures/NKandTexpression.jpeg",dpi = 300, units = "in", width = 6, height = 8)

Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)

DimPlot(log.data,label = F, group.by = "Cell.type.new",raster = F)
ggsave("Figures/cell.type.new.jpeg",dpi = 300, units = "in",  width = 6, height = 4.8)

saveRDS(log.data, file = "02_musnucseq_Celltype.rds")
saveRDS(subset, file = "11_musnucseq_TNK.rds")

log.data <- readRDS(file = "RDS_file/02_musnucseq_Celltype.rds")

#HSC
subset <- subset(log.data,idents = c("HSC"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:13)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:13)
subset <- FindClusters(subset, resolution = 0.5)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)

markerGenesFib <- c( "Ebf1","Slit3","Rgs5",#VSMC
                     "Col25a1","Reln", #qHSC
                     "Hgf","Ccbe1", #InterHSC
                     "Lama2","Col1a1", #aHSC
                     "Msln","Upk3b")# Meso

DotPlot(subset, features = markerGenes) + RotatedAxis() 
DotPlot(subset, features = markerGenesFib) + RotatedAxis() 
VlnPlot(subset, features = "nCount_RNA", pt.size = 0)
FeaturePlot(subset, features = c("Rspo3","Apoe"))


## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       
                       
                       "0" = "HSC",
                       "1" = "HSC",
                       "2" = "HSC",
                       "3" = "HSC",
                       
                       "4" = "HSC"
                       
)

subset[["Cell.type.new"]] <- Idents(object = subset)

DimPlot(subset, reduction = "umap", label = F, group.by = "Cell.type.new", raster = F)
ggsave("Figures/cell.type.new.FibHSC.jpeg",dpi = 300, units = "in", width = 6, height = 4.8)

DotPlot(subset, features = markerGenesFib) + RotatedAxis() 
ggsave("Figures/cell.type.dotplot.FibHSC.jpeg",dpi = 300, units = "in", width = 12, height = 4.8)

FeaturePlot(subset, features =c(markerGenesFib),order = T, raster = F)
ggsave( "Figures/FibHSCexpression.jpeg",dpi = 300, units = "in", width = 12, height =8)

Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)

DimPlot(log.data,label = F, group.by = "Cell.type.new",raster = F)
ggsave("Figures/cell.type.new.jpeg",dpi = 300, units = "in",  width = 6, height = 4.8)

saveRDS(log.data, file = "02_musnucseq_Celltype.rds")
saveRDS(subset, file = "11_musnucseq_Hscfib.rds")


#Endothelial
subset <- subset(log.data,idents = c("Endothelial Cell"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:13)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:13)
subset <- FindClusters(subset, resolution = 0.3)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)


DotPlot(subset, features = markerGenes) + RotatedAxis() 
DotPlot(subset, features = markerGenesFib) + RotatedAxis() 
VlnPlot(subset, features = "nCount_RNA", pt.size = 0)
FeaturePlot(subset, features = c("Rspo3","Ptprb"))

## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       
                       "0" = "Endothelial Cell",
                       "1" = "Endothelial Cell",
                       "3" = "Endothelial Cell",
                       "4" = "Endothelial Cell",
                       "5" = "Endothelial Cell",
                       "6" = "Endothelial Cell",
                       
                       "2" = "Doublets",
                       "7" = "Doublets"
                       
)

subset[["Cell.type.new"]] <- Idents(object = subset)

DimPlot(subset, reduction = "umap", label = F, group.by = "Cell.type.new", raster = F)

DotPlot(subset, features = markerGenes) + RotatedAxis() 

Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)

DimPlot(log.data,label = F, group.by = "Cell.type.new",raster = F)

ggsave("Figures/cell.type.new.jpeg",dpi = 300, units = "in",  width = 6, height = 4.8)

saveRDS(log.data, file = "02_musnucseq_Celltype.rds")
saveRDS(subset, file = "11_musnucseq_endo.rds")


#Cholangiocytes
subset <- subset(log.data,idents = c("Cholangiocytes"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:13)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:13)
subset <- FindClusters(subset, resolution = 0.3)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)

DotPlot(subset, features = markerGenes) + RotatedAxis() 
DotPlot(subset, features = markerGenesFib) + RotatedAxis() 
VlnPlot(subset, features = "nCount_RNA", pt.size = 0)
FeaturePlot(subset, features = c("Krt19","Anxa4"))


## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       
                       
                       "0" = "Cholangiocytes",
                       
                       "1" = "Doublets",
                       "2" = "Doublets"
  
                       
)

subset[["Cell.type.new"]] <- Idents(object = subset)

DimPlot(subset, reduction = "umap", label = F, group.by = "Cell.type.new", raster = F)
ggsave("Figures/cell.type.new.Chol.jpeg",dpi = 300, units = "in", width = 6, height = 4.8)

DotPlot(subset, features = markerGenes) + RotatedAxis() 

FeaturePlot(subset, features =c("ALB","APOE","CFTR"),order = T, raster = F)
ggsave("Figures/expression.Chol.jpeg",dpi = 300, units = "in", width = 6, height = 4.8)

Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)

DimPlot(log.data,label = F, group.by = "Cell.type.new",raster = F)

ggsave("Figures/cell.type.new.jpeg",dpi = 300, units = "in",  width = 10, height = 4.8)

saveRDS(log.data, file = "02_musnucseq_Celltype.rds")
saveRDS(subset, file = "11_musnucseq_chol.rds")


#hep
subset <- subset(log.data,idents = c("Hepatocytes"))

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:15)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:15)
subset <- FindClusters(subset, resolution = 0.3)

DimPlot(subset, reduction = "umap", label = T, group.by = "seurat_clusters", raster = F)
DimPlot(subset, reduction = "umap", label = T, group.by = "Cell.type", raster = F)


DotPlot(subset, features = markerGenes) + RotatedAxis() 
DotPlot(subset, features = markerGenesmyeloid) + RotatedAxis() 
VlnPlot(subset, features = "nCount_RNA", pt.size = 0)
FeaturePlot(subset, features = c("Alb","Apoe"))


## Cell clustering 
subset <- SetIdent(subset, value = "seurat_clusters")
subset <- RenameIdents(object = subset, 
                       
                       "7" = "Doublets",
                       
                       "0" = "Hepatocytes",
                       "1" = "Hepatocytes",
                       "2" = "Hepatocytes",
                       "3" = "Hepatocytes",
                       "4" = "Hepatocytes",
                       "5" = "Hepatocytes",
                       "6" = "Hepatocytes",
                       "8" = "Hepatocytes",
                       "9" = "Hepatocytes"

)

subset[["Cell.type.new"]] <- Idents(object = subset)

DimPlot(subset, reduction = "umap", label = F, group.by = "Cell.type", raster = F)
ggsave("Figures/Cell.type.Hep.jpeg",dpi = 300, units = "in", width = 8, height = 4.8)

DotPlot(subset, features = markerGenes) + RotatedAxis() 
ggsave("Figures/Dotplot.Hep.jpeg",dpi = 300, units = "in", width = 6, height = 4.8)

FeaturePlot(subset, features =c(markerGenes),order = T, raster = F)
ggsave("Figures/expression.Hep.jpeg",dpi = 300, units = "in", width = 8, height = 4.8)


Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)

DimPlot(subset,label = F, group.by = "orig.ident",raster = F)

ggsave("Figures/cell.type.new.jpeg",dpi = 300, units = "in",  width = 7, height = 4.8)

saveRDS(log.data, file = "02_musnucseq_Celltype.rds")
saveRDS(subset, file = "11_musnucseq_hep.rds")

log.data <- SetIdent(log.data, value = "Cell.type")
subset <- subset(log.data,idents = c("low quality Hepatocytes"))
log.data <- SetIdent(log.data, value = "Cell.type.new")
subset[["Cell.type.new"]] <- Idents(object = subset)
Idents(log.data) <- Idents(subset)
log.data[["Cell.type.new"]] <- Idents(object = log.data)

#
cell_counts <- table(log.data$Cell.type)
clusterlabels = paste(names(cell_counts),"(",cell_counts,")",sep = "")
DimPlot(log.data, label = F, group.by = "Cell.type",raster = F) + scale_color_discrete(label = clusterlabels)
ggsave("Figures/cell.type.count.jpeg",dpi = 300, units = "in",  width = 8, height = 4.8)

FeaturePlot(log.data, features =c(markerGenes),order = F, raster = F)
ggsave("Figures/expression_features.jpeg",dpi = 300, units = "in", width = 15, height = 12)

FeaturePlot(log.data, features =c(markerGenes1),order = T, raster = F)
ggsave("Figures/expression_rspo3_cyps.jpeg",dpi = 300, units = "in", width = 8, height = 10)

subset <- subset(log.data, ident = c("VSMC","HSC", "NK", "Mesothelial", 
                                       "Other Myeloid cell",
                                       "T Cell", "plasma cell", "B cell", "cDCs",
                                       "Kupffer Cell", "Monocyte", 
                                       "Endothelial Cell", "Cholangiocytes",
                                       "Hepatocytes"))


DefaultAssay(subset) <- "integrated"

subset <- ScaleData(subset, verbose = FALSE)
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:29)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:30)
subset <- FindClusters(subset, resolution = 0.5)

subset <- SetIdent(subset, value = "Cell.type.new")
subset[["Cell.type.new"]] <- Idents(object = subset)

subset <- SetIdent(subset, value = "Cell.type")


#
cell_counts <- table(subset$Cell.type.new)
clusterlabels = paste(names(cell_counts),"(",cell_counts,")",sep = "")
DimPlot(subset, label = F, group.by = "Cell.type.new",raster = F) + scale_color_discrete(label = clusterlabels)
ggsave("Figures/control_umap.jpeg",dpi = 300, units = "in",  width = 8, height = 4.8)

DefaultAssay(subset) <- "RNA"

saveRDS(log.data, file = "RDS_file/03_musnucseq_Celltype_doublet_Removed_2.0.rds")

log.data <- readRDS("RDS_file/03_musnucseq_Celltype_doublet_Removed_2.0.rds")
VlnPlot(log.data, features = "Cyp2e1", pt.size = 0)

log.data <- JoinLayers(log.data)

log.data_V3 <- Convert_Assay(seurat_object = log.data, convert_to = "V3", assay = "RNA")

saveRDS(log.data_V3, file = "RDS_file/03_musnucseq_Celltype_doublet_Removed_V3.rds")

celltype <- c("HSC", "NK", "Mesothelial", 
          "Other Myeloid cell",
          "T Cell", "B cell", "cDCs",
          "Kupffer Cell", "Monocyte", 
          "Endothelial Cell", "Cholangiocytes",
          "Hepatocytes")
