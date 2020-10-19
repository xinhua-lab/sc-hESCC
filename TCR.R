library(Seurat)
library(readxl)
library(pheatmap)
library(RColorBrewer) 
library(ggplot2)

rm(list=ls())
graphics.off()

### load clustering results
ZYX.combined <- readRDS(file = "ZYX.NKTcell.rds") # NKT has 0-16 clusters
cluster <- ZYX.combined$seurat_clusters

idx_C0 <- which(cluster==0)
ZYX.combined.C0 <- ZYX.combined[,idx_C0]
idx_C4 <- which(cluster==4)
ZYX.combined.C4 <- ZYX.combined[,idx_C4]
idx_C9 <- which(cluster==9)
ZYX.combined.C9 <- ZYX.combined[,idx_C9]
idx_C10 <- which(cluster==10)
ZYX.combined.C10 <- ZYX.combined[,idx_C10]
idx_C13 <- which(cluster==13)
ZYX.combined.C13 <- ZYX.combined[,idx_C13]
idx_C14 <- which(cluster==14)
ZYX.combined.C14 <- ZYX.combined[,idx_C14]
idx_C1 <- which(cluster==1)
ZYX.combined.C1 <- ZYX.combined[,idx_C1]
idx_C3 <- which(cluster==3)
ZYX.combined.C3 <- ZYX.combined[,idx_C3]
idx_C5 <- which(cluster==5)
ZYX.combined.C5 <- ZYX.combined[,idx_C5]
idx_C8 <- which(cluster==8)
ZYX.combined.C8 <- ZYX.combined[,idx_C8]
idx_C16 <- which(cluster==16)
ZYX.combined.C16 <- ZYX.combined[,idx_C16]
idx_C6 <- which(cluster==6)
ZYX.combined.C6 <- ZYX.combined[,idx_C6]
idx_C7 <- which(cluster==7)
ZYX.combined.C7 <- ZYX.combined[,idx_C7]
idx_C12 <- which(cluster==12)
ZYX.combined.C12 <- ZYX.combined[,idx_C12]

C <- c("C0", "C4","C9", "C10","C13", "C14","C1", "C3","C5", "C8","C16", "C6","C7", "C12")
Cn <- c(0, 4, 9, 10, 13, 14, 1, 3, 5, 8, 16, 6, 7, 12)

ZYX.combined.sub <- list()
for (i in 1:length(Cn)){
  idx <- which(cluster==Cn[i])
  ZYX.combined.sub[[i]] <- ZYX.combined[,idx]
  #  assign(paste0("TCR_barcodes_C", Cn[i]), list())
  #  assign(paste0("TCR_cloneTypes_C", Cn[i]), list())
}

TCR_barcode_buffer = NULL
TCR_cloneType_buffer = NULL

### load  sampleA-tumor TCR results
data <- read.table("./filtered_contig_annotations.csv", sep=",", header=TRUE)
buffer_barcode <- data[["barcode"]]
buffer_cloneType <- data[["raw_clonotype_id"]]
buffer_barcode_unique <- unique(buffer_barcode)
idx_selected <- match(buffer_barcode_unique, buffer_barcode)
buffer_barcode <- buffer_barcode[idx_selected]
buffer_cloneType <- buffer_cloneType[idx_selected]
idx_none <- which(buffer_cloneType=="None")
buffer_barcode <- buffer_barcode[-idx_none]
buffer_cloneType <- buffer_cloneType[-idx_none]
buffer_barcode <- sub("-.*", "", buffer_barcode) # extract the barcode sequence
TCR_barcode_buffer = c(TCR_barcode_buffer, buffer_barcode)
TCR_cloneType_buffer = c(TCR_cloneType_buffer, buffer_cloneType)

# for each cluster

for (i in 1:length(Cn)){
  idx1 <- which(ZYX.combined.sub[[i]]$sample == "tumor")
  idx2 <- which(ZYX.combined.sub[[i]]$group == "sampleA")
  idx_sample <- intersect(idx1, idx2)
  Cell_barcodes_buffer <- colnames(ZYX.combined.sub[[i]])[idx_sample]
  Cell_barcodes_buffer <- sub("_.*", "", Cell_barcodes_buffer)
  idx_TCRinRNAseq <- match(Cell_barcodes_buffer, buffer_barcode)
  idx_NA <- which(is.na(idx_TCRinRNAseq))
  idx_TCRinRNAseq <- idx_TCRinRNAseq[-idx_NA]
  assign(paste0("TCR_barcodes_C", Cn[i], "_sampleA-tumor"), buffer_barcode[idx_TCRinRNAseq])
  assign(paste0("TCR_cloneTypes_C", Cn[i], "_sampleA-tumor"), buffer_cloneType[idx_TCRinRNAseq])
}

# fig 4a
sample <- c("sampleA-tumor","sampleA-normal")
buffer <- factor()
for (i in 1:length(sample)){
  for (j in 1:length(Cn)){
    temp <- eval(parse(text = paste0("TCR_cloneTypes_C", Cn[j],"_", sample[i])))
    buffer <- c(buffer, as.character(temp))
  }
}

a <- table(buffer)
b <- matrix(0, length(a),1)
for (i in 1:length(a)){
  idx <- which(a==a[i])
  b[i] <- length(idx)
}
a <- log2(a)
b <- log2(b)
X <- data.frame("Ncells" = a, "Nclones" = b)
tiff('Fig4a.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
ggplot(X, aes(x=Ncells.Freq, y=Nclones)) + geom_point(size=2)+geom_smooth()+
  labs(x="log2(number of cells per clonetypes)", y = "log2(number of clonetypes)") + 
  theme(text = element_text(size=14))+geom_vline(xintercept=1, linetype="dashed", color = "red", size=0.5)
dev.off()

# fig 4b
sample <- c("sampleA-tumor","sampleA-normal")
buffer <- list()
buffer2 <- matrix(0, 3, length(sample))
for (i in 1:length(sample)){
  buffer[[i]] <- factor()
  for (j in 1:length(Cn)){
    temp <- eval(parse(text = paste0("TCR_cloneTypes_C", Cn[j],"_", sample[i])))
    buffer[[i]] <- c(buffer[[i]], as.character(temp))
  }
  a <- table(buffer[[i]])
  idx <- which(a==1)
  buffer2[1,i] <- sum(a[which(a==1)])/sum(a)
  buffer2[2,i] <- sum(a[which(a==2)])/sum(a)
  buffer2[3,i] <- sum(a[which(a>2)])/sum(a)
}
colnames(buffer2) <- sample
temp <- c(buffer2[1,], buffer2[2,], buffer2[3,])
Clones <- c(rep("n=1", 14), rep("n=2", 14), rep("n>2", 14))
sampleTypes <- c(rep(sample, 3))
buffer3 <- data.frame("PC" = temp, "noClones" = Clones, "Samples" = sampleTypes)

tiff('Fig4b.tiff', units="in", width=6, height=4.5, res=300, compression = 'lzw')
ggplot(data=buffer3, aes(x=Samples, y=PC, fill=Clones)) + geom_col(aes(y = PC, fill = noClones), 
                                                                   position = position_stack(reverse = TRUE),
                                                                   width = 0.8)+theme_minimal()+
  labs(y="Percentage of TCs (%)")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
  theme(text = element_text(size=14))+geom_text(aes(label = paste0(floor(100*temp),"%")),
                                                position = position_stack(vjust = 0.5, reverse = TRUE), size = 3)
dev.off()

# fig 4d
library(gridExtra)
library(dplyr)
Cn <- c(0, 4, 9, 10, 13, 14, 1, 3, 5, 8, 16, 6, 7, 12)
C <- c("C0", "C4","C9", "C10","C13", "C14","C1", "C3","C5", "C8","C16", "C6","C7", "C12")
sample <- c("sampleA-tumor","sampleA-normal")
buffer <- matrix(0, length(sample), length(Cn))
for (i in 1:length(sample)){
  for (j in 1:length(Cn)){
    temp <- eval(parse(text = paste0("TCR_cloneTypes_C", Cn[j],"_", sample[i])))
    a <- table(temp)
    buffer[i,j] <- sum(a[which(a>1)])/sum(a)
  }
}
colnames(buffer) <- C
rownames(buffer) <- sample
#buffer[is.na(buffer)] <- 0
buffer <- as.vector(t(buffer))
sampleTypes <- c(rep("sampleA-tumor", length(C)), rep("sampleA-normal", length(C)),  
                 )
Clusters <- c(rep(C, length(sample)))
Group <- rep(c(rep("CD4+", 6), rep("CD8+", 7)), length(sample))
X <- data.frame(clonalCells = buffer, Clusters = Clusters, Samples = sampleTypes, Group = Group)
colors <- c("#999999", "#E69F00","#999999", "#56B4E9", "#999999", "#999999", "#E69F00", 
            "#E69F00", "#999999", "#E69F00", "#E69F00", "#E69F00", "#E69F00")
my_dat <- summarise(group_by(X, Clusters), my_mean=mean(clonalCells, na.rm = TRUE), 
                    my_se = sd(clonalCells,na.rm = TRUE)/sqrt(n()))

wilcox.test(my_dat$my_mean[c(1,9,14,3,5,6)], my_dat$my_mean[c(2,8,10,13,7,11,12)])

tiff('Fig4d.tiff', units="in", width=6, height=4.5, res=300, compression = 'lzw')
ggplot() +  geom_bar(data = my_dat, aes(y = my_mean, x = Clusters, fill=Clusters), stat="identity", width=0.5)+
  geom_errorbar(data = my_dat, aes(x = Clusters, ymin = my_mean - my_se, ymax = my_mean + my_se), width=0.2)+
  scale_x_discrete(name ="Clusters", limits=C)+
  scale_fill_manual(values=colors, guide=FALSE)+
  geom_point(data = X, aes(x = Clusters, y = clonalCells, color = Samples, group = Samples))+ theme_bw() +
  labs(y="Clonal cells in each cluster (%)")
dev.off()

# fig 4f

Cn <- c(0)
C <- c("C0")
sample <- c("sampleA-tumor","sampleA-normal","sampleB-tumor","sampleB-normal")
I <- seq(from = 1 , to = 4 ,by = 2)
buffer1 <- matrix(0, length(I), length(C))
buffer2 <- matrix(0, length(I), length(C))
for (i in 1:length(I)){
  for (j in 1:length(C)){
    temp1 <- eval(parse(text = paste0("TCR_cloneTypes_C", Cn[j],"_", sample[I[i]])))
    a1 <- table(temp1)
    temp2 <- eval(parse(text = paste0("TCR_cloneTypes_C", Cn[j],"_", sample[I[i]+1])))
    a2 <- table(temp2)
    buffer1[i,j] <- sum(a1[which(a1>1)])/sum(a1)
    buffer2[i,j] <- sum(a2[which(a2>1)])/sum(a2)
  }
}
colnames(buffer1) <- C
colnames(buffer2) <- C
buffer1 <- rowMeans(buffer1, na.rm = TRUE)
buffer2 <- rowMeans(buffer2, na.rm = TRUE)
buffer <- c(buffer1, buffer2)
Samples <- c(rep("tumor", 2), rep("normal", 2))
Group <- rep(c("sampleA", "sampleB"), 2)
X <- data.frame(clonal = buffer, Samples = Samples, Group = Group)

tiff('Fig4F.tiff', units="in", width=6, height=4.5, res=300, compression = 'lzw')
ggplot(X, aes(x = Samples, y = clonal, color = Group, group = Group)) + 
  geom_point(size=5)+ geom_line(size=2)+ ggtitle("T Cells C0") +
  labs(y = "Percentage of clonal cells (%)")+ 
  scale_color_manual(values=c("#D9D9D9", "#8DD3C7"))
dev.off()
