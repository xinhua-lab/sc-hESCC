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
# load the data after quality control, removing the batch effect, Dimensionality reduction and clustering
seuratobj.data <- readRDS("processed.rds") 

# Select the clusters want to map
idx_cluster <- c(which(cluster==1), which(cluster==2), which(cluster==3))
seuratobj.sel <- seuratobj.data[,idx_cluster]

# WGCNA analysis
cluster <- seuratobj.data$seurat_clusters
idx <- c(which(cluster==1), which(cluster==2), which(cluster==3), which(cluster==4) 
         )
seuratobj.data <- seuratobj.data[,idx]
# select 300 cells
x <- 1:length(colnames(seuratobj.data@assays$RNA@counts))
id <- sample(x, size = 500, replace = F)
saveRDS(id, "id.rds")
id <- readRDS(file = "id.rds")

datExpr <- as.matrix(seuratobj.data@assays$RNA@counts)[,id]
# OGFSC
log2Data <- log2(datExpr +1)
## gene filtering by OGFSC
OGF <- OGFSC(log2Data, plot_option = 1, nBins = 30, minBinSize=100, LR_p=0.01,
             alpha=c(0.5), TW_threshold=0.0001) 
OGFSC_idx <- OGF$OGFSC_idx 
datExpr <- as.matrix(datExpr)[OGFSC_idx,]
datExpr <- t(datExpr)
datTraits <- as.matrix(seuratobj.data$seurat_clusters)
datTraits <- cbind(datTraits, seuratobj.data$orig.ident)
datTraits <- as.data.frame(datTraits[id,])
colnames(datTraits) <- c("cluster","ident")
for (i in 1:length(id)){
  datTraits$cluster[i] <- paste0("cluster", datTraits$cluster[i])
}
datTraits$cluster <- as.factor(datTraits$cluster)
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

# confirm beta
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft$powerEstimate
par(mfrow = c(1,2))
cex1 = 0.9
# 
net = blockwiseModules(datExpr, power = 1,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, maxBlockSize = 5000,
                       saveTOMFileBase = "300genes",
                       verbose = 3)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
pdf('Fig5d.pdf', width=6, height=6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# save
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

# fig 5e
nSamples <- nrow(datExpr)
nGenes = ncol(datExpr)
design=model.matrix(~0+ datTraits$cluster)
colnames(design)=levels(datTraits$cluster)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
sta <- array(data = NA)
for (i in 1:nrow(moduleTraitCor)){
  sta[i] <- max(moduleTraitCor[i,])-min(moduleTraitCor[i,])
}
moduleTraitCor <- cbind(moduleTraitCor, sta)
moduleTraitCor <- moduleTraitCor[order(moduleTraitCor[,10], decreasing = T),]
id <- which(moduleTraitCor[,10] < -0.13 | moduleTraitCor[,10] > 0.13)
id <- id[-2]
moduleTraitCor <- moduleTraitCor[id,-10]
textMatrix <- cbind(textMatrix, sta)
textMatrix <- textMatrix[order(textMatrix[,10], decreasing = T),]
textMatrix <- textMatrix[,-10]
textMatrix <- textMatrix[id,]
# Display the correlation values within a heatmap plot
pdf('Fig5e.pdf', width=6, height=6)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = T,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               cex.lab.y = 0.3,
               cex.lab.x = 0.5,
               yLabelsAngle = 45,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# fig 5f
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
cluster5 = as.data.frame(design[,5])
names(cluster5) = "cluster5"
geneTraitSignificance = as.data.frame(cor(datExpr, cluster5, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(cluster5), sep="");
names(GSPvalue) = paste("p.GS.", names(cluster5), sep="")
module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1))
module_gene <- abs(geneModuleMembership[moduleGenes, column])
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", colnames(cluster5)),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()


x <- abs(geneModuleMembership[moduleGenes, column])
y <- abs(geneTraitSignificance[moduleGenes, 1])
data <- data.frame(x, y, row.names = rownames(geneModuleMembership)[moduleGenes])
data$name <- rownames(data)
data$name[which(data$x < 0.6 | data$y < 0.4)] <- NA
write.csv(data, "dot_data.csv")
data <- read.csv("dot_data.csv", row.names = 1)

tiff('Fig5f.tiff', units="in", width=9, height=8, res=300, compression = 'lzw')
ggplot(data)+
  geom_point(aes(x=x,y=y),size = 2, alpha=1, color="DarkCyan")+
  geom_text_repel(aes(x=x,y=y,label=name))+ xlab(paste("Module Membership in", "turquoise", "module"))+
  ylab(paste("Gene significance for","cluster5"))+ theme_bw()+
  labs(title = paste("Module membership vs. gene significance\n", "cor=0.55,p=8.9e-24"))+ # need change
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("DarkCyan"))
dev.off()