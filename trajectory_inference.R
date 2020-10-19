library(Seurat)
library(cowplot)
library(limma)
library(Matrix)
library(umap)
library(dplyr)
library('scran')
library(readxl)
library(monocle)
library(xlsx)
rm(list=ls())
graphics.off()
#memory.limit(size=56000)

setwd(workpath)
#fig 5c
#load data
seuratobj.data <- readRDS("processed.rds")
cluster <- seuratobj.data$seurat_clusters
# Select the clusters want to map
idx_sel <- c(which(cluster==1), which(cluster==2), which(cluster==3)
             )

cluster <- cluster[idx_sel]

expression <- as.matrix(seuratobj.data@assays$integrated@data[,idx_sel])
gene_feature <- rownames(expression)
gene_feature <- as.matrix(gene_feature)
rownames(gene_feature) <- gene_feature
colnames(gene_feature) <- "gene_short_name"
gene_feature <- data.frame(gene_feature)
cluster <- data.frame(cluster)
pd <- new("AnnotatedDataFrame", data = cluster);
fd <- new("AnnotatedDataFrame", data = gene_feature);
cds <- newCellDataSet(expression, phenoData = pd, featureData = fd);

cds<- estimateSizeFactors(cds);
cds <- estimateDispersions(cds);

disp_table <- dispersionTable(cds);
ordering_genes <- subset(disp_table, mean_expression >= 0.1);
cds <- setOrderingFilter(cds, ordering_genes$gene_id)

#cds_reduced <- reduceDimension(cds);


cds_reduced <- reduceDimension(cds)#, max_components = 2, num_dim = 6, reduction_method = 'DDRTree', verbose = T)

cds_for_pseudo <- orderCells(cds_reduced)


tiff('Fig5c.tiff', units="in", width=6, height=4, res=300, compression = 'lzw')
plot_cell_trajectory(cds_for_pseudo , markers_linear = TRUE, cell_size = 0.5, color_by = "cluster") + scale_color_brewer(palette="Dark2")
dev.off()
