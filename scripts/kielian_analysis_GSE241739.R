library(Seurat)
library(dplyr)
library(Matrix)
library(BiocManager)
library(glmGamPoi)
library(celldex)
library(SingleR)
library(CellChat)
library(ggplot2)
library(dplyr)    # <- provides %>% and group_by/summarise
library(rlang)    # for !!sym()
library(ggplot2)
library(ggpubr)
library(edgeR)
library(fgsea)
library(msigdbr)
library(tidyr)
library(RColorBrewer)

# loading in the data, each sample individually, select that I'm looking at the gene expression, then add metadata
Tissue_1 <- Read10X_h5("/Users/kmcdon27/Downloads/GSM7734668_PJI_1_Tissue.h5")
gene_expression_4 <- Tissue_1[["Gene Expression"]]
Tissue_1_Seurat <- CreateSeuratObject(counts = gene_expression_4, project = "GSE241739")
Tissue_1_Seurat[['Sample']] <- "Sample 1"

Tissue_2 <- Read10X_h5("/Users/kmcdon27/Downloads/GSM7734670_PJI_2_Tissue.h5")
gene_expression_5 <- Tissue_2[["Gene Expression"]]
Tissue_2_Seurat <- CreateSeuratObject(counts = gene_expression_5, project = "GSE241739")
Tissue_2_Seurat[['Sample']] <- "Sample 2"

Tissue_3 <- Read10X_h5("/Users/kmcdon27/Downloads/GSE241739_PJI_5_Tissue.h5")
gene_expression_6 <- Tissue_3[["Gene Expression"]]
Tissue_3_Seurat <- CreateSeuratObject(counts = gene_expression_6, project = "GSE241739")
Tissue_3_Seurat[['Sample']] <- "Sample 3"

# merging into one object
mydata <- merge(Tissue_1_Seurat, y = c( Tissue_2_Seurat, Tissue_3_Seurat), all = T, add.cell.ids = c("Tissue_1_Seurat", "Tissue_2_Seurat", "Tissue_3_Seurat"), project = "GSE241739")
# looking at mitochondrial reads
mydata <- PercentageFeatureSet(mydata, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#filtering out bad reads
mydata <- subset(mydata, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA < 50000 & percent.mt < 5)
#plotting it again
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# saving this so I can just have it in my enivornment
saveRDS(mydata, file = "merged_seurat_kielian_tissue.rds")

# loading it in after clearing the environment 
mydata <- readRDS("merged_seurat_kielian_tissue.rds")
# joining the layers together
mydata <- JoinLayers(mydata)

#Splitting the object into a list of objects, by sample
mydata.list <- SplitObject(mydata, split.by = "Sample")

# Running SCTransform on each individual object from the list to normalize and scale the data, and identify variable features.
mydata.list <- lapply(X = mydata.list, FUN = function(x) {x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = FALSE)})

# selecting features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = mydata.list, nfeatures = 3000)

#normalizing and filtering for integration purposes
mydata.list <- PrepSCTIntegration(object.list = mydata.list, anchor.features = features)

#Finding anchors 
immune.anchors <- FindIntegrationAnchors(object.list = mydata.list, anchor.features = features, normalization.method = "SCT")#, reduction = "rpca")

# integration across Donors and creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

#See a summary of the object
immune.combined
DefaultAssay(immune.combined)

saveRDS(immune.combined, file = "merged_seurat_kielian_int_tissue.rds")
immune.combined <- readRDS("merged_seurat_kielian_int_tissue.rds")

# Running principal component analysis + other analysis
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.4)

# visualizing
DimPlot(immune.combined, reduction = "umap")
DimPlot(immune.combined, group.by = "Sample", reduction = "umap")
FeaturePlot(immune.combined, features = "CD3E", reduction = "umap")

# next, doing automatic annotations based on human reference dataset

#Formatting the normalized expression values in preparation for finding cluster marker genes
immune.combined <- PrepSCTFindMarkers(immune.combined)

#Finding all positive cluster marker genes using the SCTransform normalized values
all.markers <- FindAllMarkers(object = immune.combined, assay = "SCT", only.pos = TRUE)

#Switching from using the "integrated" assay to using the "RNA" assay.  
DefaultAssay(immune.combined) <- "RNA"

#Performning basic log-normalization 
immune.combined <- NormalizeData(immune.combined)

#adding the human cell atlas for a reference genome
hpca.se <- HumanPrimaryCellAtlasData()

#Using the reference, automatically determining each cluster's course-grained cell type
preds <- SingleR(test = GetAssayData(immune.combined, assay = "RNA", slot = "data"), clusters = immune.combined$seurat_clusters, ref = hpca.se, labels = hpca.se$label.main)
#Adding the cell type labels to a variable in the seurat object called SingleR.cluster.labels.main
immune.combined[["SingleR.cluster.labels.main"]] <- preds$labels[match(immune.combined[[]][["seurat_clusters"]], rownames(preds))]

#Repeating the above except for fine-grained cell types
preds <- SingleR(test = GetAssayData(immune.combined, assay = "RNA", slot = "data"), clusters = immune.combined$seurat_clusters, ref = hpca.se, labels = hpca.se$label.fine)
immune.combined[["SingleR.cluster.labels.fine"]] <- preds$labels[match(immune.combined[[]][["seurat_clusters"]], rownames(preds))]
# plotting them
DimPlot(immune.combined, reduction = "umap", repel = TRUE, group.by = "SingleR.cluster.labels.main")
ggsave("All_cellsUMAP.pdf", width=7, height=6)
DimPlot(immune.combined, reduction = "umap", repel = TRUE, group.by = "SingleR.cluster.labels.fine")
ggsave("UMAP_FineLabels.pdf", width=16, height=12)

#saving
saveRDS(immune.combined, file = "singler_processed_integrated_seurat_kielian_tissue.rds")

# plotting it
FeaturePlot(immune.combined, c("CD3E", "CD4", "CD8A"))

#Pulling out just the t-cells into a separate object and repeat some of the visualizations
tcell.subset <- subset(immune.combined, subset = SingleR.cluster.labels.main == "T_cells")

#Redoing the diagnostic plots and filtering, as before, except this time just on the T-cells
tcell.subset <- PercentageFeatureSet(tcell.subset, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(tcell.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "Sample", ncol = 3)
tcell.subset <- subset(tcell.subset, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA < 50000 & percent.mt < 5)
VlnPlot(tcell.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "Sample", ncol = 3)
saveRDS(tcell.subset, file = "Tcells_filt_seurat_tissue.rds")

#Redoing the integration except just on the T-cells
mydata <- readRDS("Tcells_filt_seurat_tissue.rds")
mydata.list <- SplitObject(mydata, split.by = "Sample")
mydata.list <- lapply(X = mydata.list, FUN = function(x) {
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = FALSE)})

# selecting features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = mydata.list, nfeatures = 3000)
mydata.list <- PrepSCTIntegration(object.list = mydata.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = mydata.list, anchor.features = features, normalization.method = "SCT")#, reduction = "rpca")
immune.combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", k.weight = 20)

#Redoing the clustering process, except just on the T-cells
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = c(0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0))

#Picking some appropriate resolution
immune.combined[["seurat_clusters"]] <- immune.combined[["integrated_snn_res.1"]]
Idents(object = immune.combined) <- "seurat_clusters"
DimPlot(immune.combined, reduction = "umap")

saveRDS(immune.combined, file = "Tcells_processed_interegrated_seurat_tissue.rds")
immune.combined <- readRDS("Tcells_processed_interegrated_seurat_tissue.rds")

#finding the cd4 T cells
FeaturePlot(immune.combined, c("CD4"))
ggsave("cd4.pdf", width=7, height=6)
FeaturePlot(immune.combined, c("CD8A"))
ggsave("cd8.pdf", width=7, height=6)

# taking out the cd8s 
DefaultAssay(immune.combined) <- "RNA"
cd4.tcells <- subset(immune.combined, subset = CD4 > 0 & CD8A == 0 & CD8B == 0)

# saving the cd4 files
saveRDS(cd4.tcells, file = "CD4_processed_integrated_seurat_kielian_tissue.rds")

# looking at the cd4s
DimPlot(cd4.tcells, reduction = "umap")

# making gene expression dot plots of checkpoint proteins
genes <- c("CTLA4", "HAVCR2", "LAG3", "TIGIT", "PDCD1")
Idents(cd4.tcells) <- "Sample"
DotPlot(cd4.tcells, features = genes, cols = "BuPu") + RotatedAxis()
ggsave("checkpointtissuenew.pdf", width=7, height=6)

# doing the same now for exhaustion-related transcription factors
tfs <- c("EOMES", "TCF7", "TOX", "TOX2")
DotPlot(cd4.tcells, features = tfs, cols = "BuPu") + RotatedAxis()
ggsave("tftissuenew.pdf", width=7, height=6)

# doing the same now for interesting cytokines
cytokines <- c("CCR7", "IFNG", "IL17A", "MKI67", "TNF")
DotPlot(cd4.tcells, features = cytokines, cols = "BuPu") + RotatedAxis()
ggsave("cytokinetissuenew.pdf", width=7, height=6)
