library(Seurat)
library(cowplot)
library(limma)
library(Matrix)
library(umap)
library(dplyr)
library('scran')
library(RColorBrewer)
library(Rtsne)
rm(list=ls())
graphics.off()
#memory.limit(size=56000)

setwd(workpath)
# load the data after quality control, removing the batch effect, Dimensionality reduction and clustering
seuratobj.data <- readRDS("processed.rds")  

# annotation and visualization
DefaultAssay(seuratobj.data) <- "RNA"
tiff('NAME.tiff', units="in", width=9, height=7, res=300, compression = 'lzw')
FeaturePlot(seuratobj.data, features = c(marker_gene), min.cutoff = "q9", label.size = 10)
dev.off()

DefaultAssay(seuratobj.data) <- "integrated"
cluster <- seuratobj.data$seurat_clusters
ident <- matrix("0", dim(seuratobj.data)[2],1)
cellTypeLabel <- ident

idx_T <- c(which(cluster==2), which(cluster==3), which(cluster==6), which(cluster==7),
           which(cluster==8), which(cluster==9), which(cluster==13), which(cluster==14), which(cluster==16),
           which(cluster==17), which(cluster==18), which(cluster==20), which(cluster==21), which(cluster==25), 
           which(cluster==27), which(cluster==28), which(cluster==31))
ident[idx_T] <- "1"
cellTypeLabel[idx_T] <- "T"
# get Subgroup data
seuratobj.data_T <- seuratobj.data[,idx_T]
saveRDS(seuratobj.data_T, file = "seuratobj.data_T.rds")

# after re-Dimensionality reduction and re-clustering
tiff('T_recluster.tiff', units="in", width=15, height=5, res=300, compression = 'lzw')
DimPlot(seuratobj.data_T, reduction = "umap", label = FALSE)
dev.off()

# visualization_heatmap
DefaultAssay(seuratobj.data_T) <- "RNA"
cluster <- seuratobj.data_T$seurat_clusters
cl <- c(cluster_order)
IDX <- list()
for (i in 1:length(cl)){
  IDX[[i]] <- which(cluster==cl[i])
}

geneList_selected <- c(function_signature_gene)
idx <- matrix(1, 1, length(geneList_selected))
for (i in 1:length(geneList_selected)){
  idx[i] <- which(rownames(seuratobj.data_T)  %in%  geneList_selected[i])
}

expression <- matrix(0, length(geneList_selected), length(cl))
for (i in 1:length(cl)){
  buffer <- as.matrix(seuratobj.data_T@assays$RNA@data[,IDX[[i]]])
  buffer <- buffer[idx,]
  expression[,i] <- rowMeans(buffer, na.rm = TRUE)
}
rownames(expression) <- geneList_selected
colnames(expression) <- c(cluster_name_order)

tiff('NAME.tiff', units="in", width=6.5, height=5.5, res=300, compression = 'lzw')
pheatmap(expression, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, 
         show_colnames = TRUE, show_rownames = TRUE, gaps_row = c(4,16,21, 26, 34, 37), gaps_col = c(4, 6),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
dev.off()