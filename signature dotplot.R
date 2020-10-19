library(Seurat)
library(cowplot)
library(limma)
library(Matrix)
library(umap)
library(dplyr)
library('scran')
library(readxl)
library(pheatmap)
library(RColorBrewer) 
library(monocle)
library(ggpubr)

rm(list=ls())
graphics.off()
#memory.limit(size=56000)

setwd(workpath)
# load the data after quality control, removing the batch effect, Dimensionality reduction and clustering
seuratobj.data <- readRDS("processed.rds") 

# Select the clusters want to map
idx_cluster <- c(which(cluster==1), which(cluster==2), which(cluster==3))
seuratobj.sel <- seuratobj.data[,idx_cluster]

# cell cycle, dotplot
DefaultAssay(seuratobj.data) <- "RNA"
seuratobj.sel <- seuratobj.data[,idx_cluster]

# G1 list
T <- read_excel("cell cycle genes.xlsx", sheet = 1, col_names = FALSE)
geneList <- T$...1
geneList <- geneList[-1]
idx_G1 <- match(geneList, rownames(seuratobj.data))
idx_G1 <- idx_G1[ !is.na(idx_G1) ]
# G2 list
T <- read_excel("cell cycle genes.xlsx", sheet = 2, col_names = FALSE)
geneList <- T$...1
geneList <- geneList[-1]
idx_G2 <- match(geneList, rownames(seuratobj.data))
idx_G2 <- idx_G2[ !is.na(idx_G2) ]


expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_G1,])
G1_scores <- colMeans(expression, na.rm = FALSE)
expression <- as.matrix(seuratobj.sel@assays$RNA@data[idx_G2,])
G2_scores <- colMeans(expression, na.rm = FALSE)
expr <- rbind(G1_scores, G2_scores)
rownames(expr) <- c("G1", "G2")

sd <- colSds(t(expr))*sqrt((dim(expr)[2]-1)/(dim(expr)[2]))
mean <- colMeans(t(expr), na.rm = FALSE, dims = 1)
z_score <- matrix(0, dim(expr)[1], dim(expr)[2])
for (i in 1:2){
  z_score[i,] <- (expr[i,]-mean[i])/sd[i]
}
rownames(z_score) <- rownames(expr)
colnames(z_score) <- colnames(expr)
temp <- rbind(seuratobj.sel@assays$RNA@data, z_score)
seuratobj.sel@assays$RNA@data <- temp



markers.to.plot <- c("G1", "G2")
tiff('Fig3m.tiff', units="in", width=6, height=4.5, res=300, compression = 'lzw')
DotPlot(seuratobj.sel, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 16) + RotatedAxis() +  ggtitle("NK+ cell cycle")
dev.off()


