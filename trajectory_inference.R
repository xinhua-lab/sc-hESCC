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
#load data
cds_for_pseudo <- readRDS(file = "seuratobjdata_T_traj.rds") # 10000 random picked cells


tiff('NAME.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
plot_cell_trajectory(cds_for_pseudo , markers_linear = TRUE, color_by = "cluster", cell_size = 0.5) + scale_color_brewer(palette="Dark2")
dev.off()

DefaultAssay(seurat.objdata) <- "RNA"
seuratobjdata_T <- seurat.objdata[,idx_T]


#T <- read_excel("T signatures.xlsx", col_names = FALSE, sheet = 3) # self-defined gene list
T <- read_excel("function gene set.xlsx", col_names = FALSE) # published gene list
geneList <- T$...1
idx_function <- match(geneList, rownames(seurat.objdata))
idx_function <- idx_Cytotoxic[ !is.na(idx_function) ]

expression <- as.matrix(expression <- as.matrix(seurat.objdata_T@assays$RNA@data[idx_function,])
                        function_scores <- colMeans(expression, na.rm = FALSE)@assays$RNA@data[idx_function,])
function_scores <- colMeans(expression, na.rm = FALSE)
expr <- rbind(several function_scores)
rownames(expr) <- c(function_name)

sd <- colSds(t(expr))*sqrt((dim(expr)[2]-1)/(dim(expr)[2]))
mean <- colMeans(t(expr), na.rm = FALSE, dims = 1)
z_score <- matrix(0, dim(expr)[1], dim(expr)[2])
for (i in 1:2){
  z_score[i,] <- (expr[i,]-mean[i])/sd[i]
}
rownames(z_score) <- rownames(expr)
colnames(z_score) <- colnames(expr)

temp <- rbind(seuratobjdata_T@assays$RNA@data, z_score)
seuratobjdata_T@assays$RNA@data <- temp
markers.to.plot <- c(function_name)

tiff('NAME.tiff', units="in", width=6, height=5, res=300, compression = 'lzw')
DotPlot(function_name, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 16) + RotatedAxis() +  ggtitle("CD4+")
dev.off()