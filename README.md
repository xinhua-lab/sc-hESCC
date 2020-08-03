# sc-hESCC
Example scripts to process and analyze data for the ESCC scRNA-seq article. 
Detailed information will be available from authors upon reasonable request.

Key software and algorithms：
10x Genomics Cell Ranger 3.0.1 version
R 3.6

R package：
Seurat (v3)
cowplot
limma
Matrix
umap
dplyr
'scran'
Rtsne
RColorBrewer
monocle
scTHI
SCENIC

Basic information of preprocessing

quality control, removing the batch effect

The 10x Genomics Cell Ranger pipeline was used to demultiplex raw files into FASTQ files, extract barcodes and UMI, filter, and map reads to the GRCh38 reference genome, and generate a matrix containing normalized gene counts versus cells per sample. 
Counts data was then imported into the Seurat object for quality control. Low-quality cells (<400 genes/cell and >10% mitochondrial genes) were excluded. To remove the batch effect, the datasets collected from different samples were integrated with default parameters.

Dimensionality reduction and clustering

The Seurat function ‘FindVariableFeatures’ was applied to identify the highly variable genes (HVG). The top 2000 HVGs were used for data integration. The data were scaled using ‘ScaleData’ and the first 20 principle components were adopted for auto-clustering analyses using ‘FindNeighbors’ and ‘FindClusters’ functions. For all cells, we identified clusters setting the resolution parameter as 1.5, and the clustering results were visualized with the UMAP scatter plot. The marker genes of each cell cluster were identified using the ROC analysis function provided by the Seurat ‘FindAllMarkers’ function for the top genes with the largest AUC (area under curve).
