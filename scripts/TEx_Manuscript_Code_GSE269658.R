library(Seurat)
library(clustree)
library(SingleR)
library(celldex)
library(sctransform)
library(glmGamPoi)
library(ggplot2)

#Load the data for each sample, where STE_057 is a directory in the current working directory and it contains the files:
#barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
expression_matrix <- Read10X(data.dir = 'STE_057', gene.column = 2)
#Convert the data to a seurat object, removing genes present in <3 cells and cells with <200 genes detected
STE_057 <- CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 200, project = "STE_057")

#Apply metadata to the seurat object for this sample
STE_057[['Sample']] <- "STE_057"
STE_057[['Condition']] <- "Sterile"
STE_057[['Donor']] <- "D057"

#Repeat the above for subsequent samples
expression_matrix <- Read10X(data.dir = 'STE_058', gene.column = 2)
STE_058 <- CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 200, project = "STE_058")

STE_058[['Sample']] <- "STE_058"
STE_058[['Condition']] <- "Sterile"
STE_058[['Donor']] <- "D058"

expression_matrix <- Read10X(data.dir = 'STE_059', gene.column = 2)
STE_059 <- CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 200, project = "STE_059")

STE_059[['Sample']] <- "STE_059"
STE_059[['Condition']] <- "Sterile"
STE_059[['Donor']] <- "D059"

expression_matrix <- Read10X(data.dir = 'INF_057', gene.column = 2)
INF_057 <- CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 200, project = "INF_057")

INF_057[['Sample']] <- "INF_057"
INF_057[['Condition']] <- "Infected"
INF_057[['Donor']] <- "D057"

expression_matrix <- Read10X(data.dir = 'INF_058', gene.column = 2)
INF_058 <- CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 200, project = "INF_058")

INF_058[['Sample']] <- "INF_058"
INF_058[['Condition']] <- "Infected"
INF_058[['Donor']] <- "D058"

expression_matrix <- Read10X(data.dir = 'INF_059', gene.column = 2)
INF_059 <- CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 200, project = "INF_059")

INF_059[['Sample']] <- "INF_059"
INF_059[['Condition']] <- "Infected"
INF_059[['Donor']] <- "D059"

#Merge the sample-specific seurat objects into one object
mydata <- merge(STE_057, y = c(STE_058, STE_059, INF_057, INF_058, INF_059), all = T, add.cell.ids = c("STE_057", "STE_058", "STE_059", "INF_057", "INF_058", "INF_059"), project = "TEx")

#Calculate the percentage of all reads coming from mitochonrial genes in each cell and store it as a variable called percent.mt
mydata <- PercentageFeatureSet(mydata, pattern = "^MT-", col.name = "percent.mt")

#Plot the number of genes detected (nFeatures_RNA), the number of reads (nCount_RNA), and the percent mitochondrial reads for each cell
#These are useful diagnostics
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Save the plot
ggsave("Pre-filtering_Diagnostic_Violin.pdf", width=12, height=8)

###

VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha = 0.1)
#Save the plot
ggsave("LowAlpha-Pre-filtering_Diagnostic_Violin.pdf", width=12, height=8)

VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0)
#Save the plot
ggsave("NoDots-Pre-filtering_Diagnostic_Violin.pdf", width=12, height=8)

#

VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha = 0.1, log = TRUE)
#Save the plot
ggsave("Log-LowAlpha-Pre-filtering_Diagnostic_Violin.pdf", width=12, height=8)

VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0, log = TRUE)
#Save the plot
ggsave("Log-NoDots-Pre-filtering_Diagnostic_Violin.pdf", width=12, height=8)

###

#Examine the above plot and pick reasonable filtering values.  Extreme outliers and high percent.mt should be excluded.
#This command filters out cells based on the "subset =" criteria
mydata <- subset(mydata, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA < 50000 & percent.mt < 5)

#Remake the plots after having removed the offending cells
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("Post-filtering_Diagnostic_Violin.pdf", width=12, height=8)

###

#Remake the plots after having removed the offending cells
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha = 0.1)
ggsave("LowAlpha-Post-filtering_Diagnostic_Violin.pdf", width=12, height=8)

VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0)
ggsave("NoDots-Post-filtering_Diagnostic_Violin.pdf", width=12, height=8)

#

VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha = 0.1, log = TRUE)
ggsave("Log-LowAlpha-Post-filtering_Diagnostic_Violin.pdf", width=12, height=8)

VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0, log = TRUE)
ggsave("Log-NoDots-Post-filtering_Diagnostic_Violin.pdf", width=12, height=8)

###

#Split the object into a list of objects, by Donor
mydata.list <- SplitObject(mydata, split.by = "Donor")

#Run SCTransform on each individual object from the list.
#SCTransform normalizes and scales the data, and identifies variable features (important later)
#SCTransform can also be used to adjust for confounding factors.  Here we're adjusting for percent.mt
#More about SCTransform here: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html and here: https://satijalab.org/seurat/articles/sctransform_vignette.html
mydata.list <- lapply(X = mydata.list, FUN = function(x) {
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = FALSE)
})

# select features that are repeatedly variable across datasets for integration
#More about integration here: https://satijalab.org/seurat/articles/integration_introduction.html
#Also see the first SCTransform link above for integration with SCTransform
features <- SelectIntegrationFeatures(object.list = mydata.list, nfeatures = 3000)

#Some normalization and filtering for integration purposes
mydata.list <- PrepSCTIntegration(object.list = mydata.list, anchor.features = features)

#Find anchors - this is an important part of integration
immune.anchors <- FindIntegrationAnchors(object.list = mydata.list, anchor.features = features, normalization.method = "SCT")#, reduction = "rpca")

# this command does the integration across Donors and creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

#See a summary of the object
immune.combined
DefaultAssay(immune.combined)

#Run principal component analysis.  Read more here: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)

#Plot the % variance explained by each PC
ElbowPlot(immune.combined, ndims = 50)
ggsave("PCA_Elbow_Plot.pdf", width=12, height=8)

#Based on the above plot, pick some number of PCs (ideally the ones before it flattens out, see more in the link above)
#Run UMAP (dimension reduction) using the select # of PCs (30 here)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)

#FindNeighbors is a step in the clustering processing
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)

#Run the final step in the clustering process at different resolutions (higher resolutions will give more clusters)
immune.combined <- FindClusters(immune.combined, resolution = c(1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4))

#This will create a diagnostic plot of the clusters at each resolution.  More info here: https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html
clustree(immune.combined, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
ggsave("Fine_Clustree_Plot.pdf", width=12, height=12)

#Pick the resolution you want to use and set those as the clusters
immune.combined[["seurat_clusters"]] <- immune.combined[["integrated_snn_res.2.4"]]
Idents(object = immune.combined) <- "seurat_clusters"

#

#Plot a UMAP colored by cluster, with no labels
DimPlot(immune.combined, reduction = "umap", label = FALSE, repel = TRUE)
ggsave("2B_UMAP_Clusters_Unlabeled.pdf", width=10, height=8)

###

#Load a reference dataset, in this case the Human Primary Cell Atlas: https://bioconductor.org/packages/3.18/data/experiment/vignettes/celldex/inst/doc/userguide.html#21_Human_primary_cell_atlas_(HPCA)
hpca.se <- HumanPrimaryCellAtlasData()

#Using the reference, automatically determine each cluster's course-grained cell type
preds <- SingleR(test = GetAssayData(immune.combined, assay = "RNA", slot = "data"), clusters = immune.combined$seurat_clusters, ref = hpca.se,
                 labels = hpca.se$label.main)

#Add the cell type labels to a variable in the seurat object called SingleR.cluster.labels.main
immune.combined[["SingleR.cluster.labels.main"]] <- 
  preds$labels[match(immune.combined[[]][["seurat_clusters"]], rownames(preds))]

#Repeat the above except for fine-grained cell types
preds <- SingleR(test = GetAssayData(immune.combined, assay = "RNA", slot = "data"), clusters = immune.combined$seurat_clusters, ref = hpca.se,
                 labels = hpca.se$label.fine)

immune.combined[["SingleR.cluster.labels.fine"]] <- 
  preds$labels[match(immune.combined[[]][["seurat_clusters"]], rownames(preds))]

###
#############################################################################################################################################################
###

FeaturePlot(immune.combined,
            reduction = "umap",
            features = "CD19", label = FALSE, repel = TRUE)

ggsave("2C2_UMAP_CD19.pdf", width=10, height=8)

###
#############################################################################################################################################################
###

FeaturePlot(immune.combined,
            reduction = "umap",
            features = "CD3E", label = FALSE, repel = TRUE)

ggsave("2C1_UMAP_CD3E.pdf", width=10, height=8)

###
#############################################################################################################################################################
###

#Pull out (into a new seurat object) just the cells that were automatically annotated as T-cells AND that are not expressing CD19
mydata <- subset(immune.combined, subset = (SingleR.cluster.labels.main == "T_cells") & (CD19 == 0))

#Redo the diagnostic plots and filtering, as before, except this time just on the T-cells
mydata <- PercentageFeatureSet(mydata, pattern = "^MT-", col.name = "percent.mt")

VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "Sample", ncol = 3)
ggsave("Tcells_Pre-filtering_Diagnostic_Violin.pdf", width=12, height=8)

mydata <- subset(mydata, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA < 50000 & percent.mt < 5)

VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "Sample", ncol = 3)

ggsave("Tcells_Post-filtering_Diagnostic_Violin.pdf", width=12, height=8)

mydata.list <- SplitObject(mydata, split.by = "Donor")

mydata.list <- lapply(X = mydata.list, FUN = function(x) {
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = FALSE)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = mydata.list, nfeatures = 3000)

mydata.list <- PrepSCTIntegration(object.list = mydata.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = mydata.list, anchor.features = features, normalization.method = "SCT")#, reduction = "rpca")

immune.combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

#Redo the clustering process, as in Go1, except just on the T-cells
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
ElbowPlot(immune.combined, ndims = 50)
ggsave("Tcells_PCA_Elbow_Plot.pdf", width=12, height=8)


immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = c(0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0))
#immune.combined <- FindClusters(immune.combined, resolution = c(1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4))

clustree(immune.combined, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
ggsave("Tcells_Fine_Clustree_Plot.pdf", width=12, height=12)

clustree(immune.combined, prefix = "integrated_snn_res.", node_colour = "sc3_stability", node_text_size = 0)
ggsave("NoLabels_Tcells_Fine_Clustree_Plot.pdf", width=12, height=12)

#Pick some appropriate resolution
immune.combined[["seurat_clusters"]] <- immune.combined[["integrated_snn_res.1"]]
Idents(object = immune.combined) <- "seurat_clusters"

#Set the levels of the Condition variable.  This determines the order (Sterile vs. Infected or Infected vs. Sterile)
#of the differential expression comparison and the order in which the groups show up in figure legends.
immune.combined[["Condition"]] <- factor(immune.combined$Condition, levels = c("Sterile", "Infected"))

#UMAP split by Condition
DimPlot(immune.combined, reduction = "umap", split.by = "Condition", label = TRUE, repel = TRUE)
ggsave("2D_Tcells_UMAP_Condition.pdf", width=16, height=10)

#Plot UMAPs, colored by the various metadata variables
DimPlot(immune.combined, reduction = "umap", group.by = "Condition")
ggsave("Tcells_UMAP_Condition.pdf", width=10, height=8)
DimPlot(immune.combined, reduction = "umap", group.by = "Donor")
ggsave("Tcells_UMAP_Donor.pdf", width=10, height=8)
DimPlot(immune.combined, reduction = "umap", group.by = "Sample")
ggsave("Tcells_UMAP_Sample.pdf", width=10, height=8)

#As before, we're normalizing/reformatting expression values and finding all positive cluster markers.
immune.combined <- PrepSCTFindMarkers(immune.combined)
all.markers <- FindAllMarkers(object = immune.combined, assay = "SCT", only.pos = TRUE)
write.csv(all.markers, "Tcells_Cluster_Markers.csv")

###

DefaultAssay(immune.combined) <- "RNA"

#Do basic log normalization of the expression values.
immune.combined <- NormalizeData(immune.combined)

# Set a list of features to plot.
features.integrated.t <- c("CD4", "CD8A", "CCR7", "IL7R", 
                           "IFNG", "RORC", "GZMK", "MKI67", 
                           "FOXP3", "IL2RA","TNFRSF18", "PDCD1", 
                           "CTLA4", "HAVCR2", "TIGIT","LAG3")

# Using StackedVlnPlot function in "CellChat" library
#Plot the expression of each of these genes in each cluster.
StackedVlnPlot(immune.combined, features = features.integrated.t)
ggsave("2E_Tcells_Stacked_Violin.pdf", width=16, height=12)

#Make a similar plot, except split each violin by Condition
StackedVlnPlot(immune.combined, features = features.integrated.t, split.by = "Condition", color.use = c("royalblue", "red"))
ggsave("Tcells_Condition_Stacked_Violin.pdf", width=24, height=12)

###
#############################################################################################################################################################
###

immune.combined <- subset(immune.combined, idents = c(8, 20))

DefaultAssay(immune.combined) <- "RNA"

immune.combined <- NormalizeData(immune.combined)
immune.combined <- ScaleData(immune.combined, features = rownames(GetAssayData(object = immune.combined, slot = "counts")))
immune.combined <- FindVariableFeatures(immune.combined)

###

goisA <- c("CTLA4", "PDCD1", "LAG3", "HAVCR2", "CD244", "TIGIT", "BTLA", "CD160", #Checkpoint proteins
           "EOMES", # Point 2 - terminally exhausted
           "IL10", "TGFB1", "NR4A1", "NR4A2", "NR4A3", "NFATC1", "NFATC2", "NFATC3", "NFATC4", "NFAT5", "HNF1A") #Increased production, TCF1 = HNF1A

goisB <- c("CD4", "CXCR5", "ICOS", #Surface markers
           "IL6", "IL21", #Secreted cytokines
           "BCL6", "IRF4", "STAT4", #Transcription factors
           "TRBV20OR9-2", "TOX", #Epigenetic commitment, TRBV20OR9-2 = TCR
           "TRAF1", "IRF9", # Downregulated
           "IFNG", "IL2", "TNF", "IL4", "IL17A", "IL22", "IL27") # Lower production

goisC <- c("CCL1", "CCL2", "CCL13", "CCL22", "CCL3", "CCL23", "CCL4", # Beta chemokines - continued on next line
           "CCL15", "CCL24", "CCL5", "CCL25")

#CCL12, CCL6 in humans is no longer a thing
#CCL21, CCL14, CCL16 not present in dataset

goisD <- c("CCL17", "CCL18",
           "CCL27", "CCL19", "CCL28", "CCL20",
           "CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CX3CR1", "TAFA2")

#CCL26, CCL7 not present in dataset
#CCRL1 = CX3CR1

goisE <- c("IL2", "IL15", "IL21", "IL4", "IL18", "TNF", "IL1B", "IFNG", "CSF2", "CD40LG", #Cytokines inhibited - Continued next line
           "IL27", "TNFSF15", "IL12B", "IL7", "MIF", "CNTF", "TNFSF13B", "IL1A", "IL3", "CXCL12", "CCR2")

#IL33, IFNA2 not present in dataset

goisA <- unique(goisA)
goisB <- unique(goisB)
goisC <- unique(goisC)
goisD <- unique(goisD)
goisE <- unique(goisE)

gois <- c(goisA, goisB, goisC, goisD, goisE)

gois <- unique(gois)

###

DefaultAssay(immune.combined) <- "integrated"

useless <- c("CD160", "IL10", "NFATC4", "HNF1A", "CXCR5", "IL6", "TRBV20OR9-2", "IL2", "IL4", "IL22", "IL27", "CCL1", "CCL2", "CCL13", "CCL22",
             "CCL23", "CCL15", "CCL24", "CCL25", "CCL17", "CCL18", "CCL27", "CCL19", "CCR3", "IL18", "IL1B", "TNFSF15", "IL12B", "CNTF", "IL1A", "IL3", "CXCL12")

tex_gois <- gois[!(gois %in% useless)]

immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE, features = tex_gois)

ElbowPlot(immune.combined, ndims = 30)
ggsave("Reclustered_Th1-Th17_PCA_Elbow_Plot.pdf", width=12, height=8)

#SELECT NUMBER OF PCs and Genes and Assay! Integrated, Non-GOIs, 30 PCs
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:6)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:6)
immune.combined <- FindClusters(immune.combined, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

clustree(immune.combined, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
ggsave("Reclustered_Th1-Th17_Clustree_Plot.pdf", width=12, height=12)

#Pick some appropriate resolution
immune.combined[["reclustered_clusters"]] <- immune.combined[["integrated_snn_res.0.6"]]
Idents(object = immune.combined) <- "reclustered_clusters"

DefaultAssay(immune.combined) <- "RNA"

#Note that default assay is integrated at this point
saveRDS(immune.combined, file = "Reclustered_Th1-Th17_seurat.rds")

DefaultAssay(immune.combined) <- "integrated"

###

DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
ggsave("Reclustered_Th1-Th17_UMAP_Clusters.pdf", width=8, height=8)

#UMAP split by Condition
DimPlot(immune.combined, reduction = "umap", split.by = "Condition", label = TRUE, repel = TRUE)
ggsave("Reclustered_Th1-Th17_UMAP_Split_Condition.pdf", width=12, height=8)

DimPlot(immune.combined, reduction = "umap", group.by = "Condition", label = TRUE, repel = TRUE)
ggsave("Reclustered_Th1-Th17_UMAP_Colored_Condition.pdf", width=8, height=8)

DimPlot(immune.combined, reduction = "umap", group.by = "Sample", label = TRUE, repel = TRUE)
ggsave("Reclustered_Th1-Th17_UMAP_Colored_Sample.pdf", width=8, height=8)

#As before, we're normalizing/reformatting expression values and finding all positive cluster markers.
DefaultAssay(immune.combined) <- "SCT"
immune.combined <- PrepSCTFindMarkers(immune.combined)
all.markers <- FindAllMarkers(object = immune.combined, assay = "SCT", only.pos = TRUE)
write.csv(all.markers, "Reclustered_Th1-Th17_Cluster_Markers.csv")

#

#Get the proportions of each cell type that are Sterile and Infected
tcell.integrated.cellcounts <- table(Idents(immune.combined), immune.combined$Condition)
tcell.integrated.cellcounts <- as.data.frame(tcell.integrated.cellcounts)
tcell.integrated.cellcounts$Var1 <- as.character(tcell.integrated.cellcounts$Var1)
tcell.integrated.cellcounts <- cbind(" "=rownames(tcell.integrated.cellcounts), tcell.integrated.cellcounts)

#Plot these proporations for manual_clusters
ggplot(tcell.integrated.cellcounts, aes(x = Var1, y = Freq, fill = Var2)) + 
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Cluster") +
  ylab("Proportion") + 
  #scale_fill_manual(values = brewer.pal(12, "Paired")) + 
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c("royalblue", "red"))

ggsave("Reclustered_Th1-Th17_Cluster_Frequencies.pdf", width=6, height=8)


###
#############################################################################################################################################################
###

DefaultAssay(immune.combined) <- "RNA"

# Using StackedVlnPlot function in "CellChat" library
#Plot the expression of each of these genes in each cluster.
StackedVlnPlot(immune.combined, features = goisA)
ggsave("Reclustered_Th1-Th17_Stacked_Violin-IncreasedProduction.pdf", width=8, height=16)

#Make a similar plot, except split each violin by Condition
StackedVlnPlot(immune.combined, features = goisA, split.by = "Condition", color.use = c("royalblue", "red"))
ggsave("Reclustered_Th1-Th17_Condition_Stacked_Violin-IncreasedProduction.pdf", width=8, height=16)

#

StackedVlnPlot(immune.combined, features = goisB)
ggsave("Reclustered_Th1-Th17_Stacked_Violin-Mixed.pdf", width=8, height=16)

#Make a similar plot, except split each violin by Condition
StackedVlnPlot(immune.combined, features = goisB, split.by = "Condition", color.use = c("royalblue", "red"))
ggsave("Reclustered_Th1-Th17_Condition_Stacked_Violin-Mixed.pdf", width=8, height=16)

#

StackedVlnPlot(immune.combined, features = goisC)
ggsave("Reclustered_Th1-Th17_Stacked_Violin-ChemokinesA.pdf", width=8, height=16)

#Make a similar plot, except split each violin by Condition
StackedVlnPlot(immune.combined, features = goisC, split.by = "Condition", color.use = c("royalblue", "red"))
ggsave("Reclustered_Th1-Th17_Condition_Stacked_Violin-ChemokinesA.pdf", width=8, height=16)

#

StackedVlnPlot(immune.combined, features = goisD)
ggsave("Reclustered_Th1-Th17_Stacked_Violin-ChemokinesB.pdf", width=8, height=16)

#Make a similar plot, except split each violin by Condition
StackedVlnPlot(immune.combined, features = goisD, split.by = "Condition", color.use = c("royalblue", "red"))
ggsave("Reclustered_Th1-Th17_Condition_Stacked_Violin-ChemokinesB.pdf", width=8, height=16)

#

StackedVlnPlot(immune.combined, features = goisE)
ggsave("Reclustered_Th1-Th17_Stacked_Violin-Cytokines.pdf", width=8, height=16)

#Make a similar plot, except split each violin by Condition
StackedVlnPlot(immune.combined, features = goisE, split.by = "Condition", color.use = c("royalblue", "red"))
ggsave("Reclustered_Th1-Th17_Condition_Stacked_Violin-Cytokines.pdf", width=8, height=16)

##

DoHeatmap(immune.combined, features = gois, raster = FALSE)
ggsave("Reclustered_Th1-Th17_Heatmap.pdf", width=8, height=24)

DoHeatmap(immune.combined, features = gois, group.by = "Condition", raster = FALSE)
ggsave("Reclustered_Th1-Th17_Condition_Heatmap.pdf", width=8, height=24)

DotPlot(immune.combined, features = gois) + RotatedAxis() + coord_flip()
ggsave("Reclustered_Th1-Th17_Dotplot.pdf", width=6, height=24)

DotPlot(immune.combined, features = gois) + RotatedAxis() + coord_flip()
ggsave("Reclustered_Th1-Th17_Dotplot.pdf", width=6, height=24)

DotPlot(immune.combined, features = gois, split.by = "Condition") + RotatedAxis()
ggsave("Reclustered_Th1-Th17_Condition_Dotplot.pdf", width=24, height=6)

#

markers <- FindAllMarkers(immune.combined, only.pos = TRUE)

markers %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(immune.combined, features = top10$gene, raster = FALSE, group.by = 'reclustered_clusters') + NoLegend()

ggsave("Reclustered_Th1-Th17_ClusterMarker_Heatmap.pdf", width=8, height=24)

###

#Try to define exhausted population

immune.combined[['Exhaustion']] <- immune.combined[['reclustered_clusters']]

#Set this new variable as the Idents of the object.  Idents are always a metadata variable, and a variable being set as Idents
#gives is a special role in certain functions.
Idents(immune.combined) <- 'Exhaustion'

#Rename the Idents (the variable called manual_clusters) from cluster number to cell type annotation.
#Note that we set multiple clusters to the same annotation sometimes.
#I figured out the mapping of cluster number to annotation by looking at the stacked violin plot.
immune.combined <- RenameIdents(immune.combined, `0` = 'Exhausted', 
                                `1` = 'Unexhausted',
                                `2` = 'Unexhausted',
                                `3` = 'Exhausted',
                                `4` = 'Unexhausted',
                                `5` = 'Unexhausted',
                                `6` = 'Exhausted')

#For some reason you need to do this after renaming them in order for it to actually show up in the manual_clusters column/variable.
immune.combined[['Exhaustion']] <- Idents(immune.combined)

#

Idents(immune.combined) <- 'Exhaustion'

DoHeatmap(immune.combined, features = gois, raster = FALSE)
ggsave("Exhaustion_Reclustered_Th1-Th17_Heatmap.pdf", width=8, height=24)

DoHeatmap(immune.combined, features = gois, group.by = "Condition", raster = FALSE)
ggsave("Exhaustion_Reclustered_Th1-Th17_Condition_Heatmap.pdf", width=8, height=24)

DotPlot(immune.combined, features = gois) + RotatedAxis() + coord_flip()
ggsave("Exhaustion_Reclustered_Th1-Th17_Dotplot.pdf", width=6, height=24)

DotPlot(immune.combined, features = gois) + RotatedAxis() + coord_flip()
ggsave("Exhaustion_Reclustered_Th1-Th17_Dotplot.pdf", width=6, height=24)

DotPlot(immune.combined, features = gois, split.by = "Condition") + RotatedAxis()
ggsave("Exhaustion_Reclustered_Th1-Th17_Condition_Dotplot.pdf", width=24, height=6)

#

DE <- FindMarkers(immune.combined, ident.1 = 'Unexhausted', ident.2 = 'Exhausted')

de <- DE[DE$p_val_adj <= 0.05 & DE$avg_log2FC > 0,]
de <- de[order(de$p_val), ]

if (nrow(de) > 50) {
  diffs <- rownames(de)[1:50]
} else {
  diffs <- rownames(de)
}

de <- DE[DE$p_val_adj <= 0.05 & DE$avg_log2FC < 0,]
de <- de[order(de$p_val), ]

if (nrow(de) > 50) {
  diffs <- append(diffs, de$gene[1:50])
} else {
  diffs <- append(diffs, de$gene)
}

write.csv(DE, "Reclustered_Th1-Th17_ExhaustionMarkers.csv")

DotPlot(immune.combined, features = diffs) + RotatedAxis()
ggsave("DE_ExhaustedClusters_Dotplot.pdf", width=24, height=6)

DoHeatmap(immune.combined, features = diffs, raster = FALSE)
ggsave("DE_ExhaustedClusters_Heatmap.pdf", width=8, height=24)

DoHeatmap(immune.combined, features = diffs, raster = FALSE, group.by = 'reclustered_clusters')
ggsave("Clusters_DE_ExhaustedClusters_Heatmap.pdf", width=8, height=24)

###

Idents(immune.combined) <- immune.combined$reclustered_clusters

# Set a list of features to plot.
features.integrated.t <- c("CD4", "CD8A", "CCR7", "IL7R", 
                           "IFNG", "RORC", "GZMK", "MKI67", 
                           "FOXP3", "IL2RA","TNFRSF18", "PDCD1", 
                           "CTLA4", "HAVCR2", "TIGIT","LAG3")

# Using StackedVlnPlot function in "CellChat" library
#Plot the expression of each of these genes in each cluster.
StackedVlnPlot(immune.combined, features = features.integrated.t)
ggsave("Original-GOIs_Stacked_Violin.pdf", width=16, height=12)

#Make a similar plot, except split each violin by Condition
StackedVlnPlot(immune.combined, features = features.integrated.t, split.by = "Condition", color.use = c("royalblue", "red"))
ggsave("Original-GOIs_Condition_Stacked_Violin.pdf", width=24, height=12)

###
#############################################################################################################################################################
###

Idents(immune.combined) <- immune.combined$reclustered_clusters

DefaultAssay(immune.combined) <- "RNA"

# Set a list of features to plot.

features.integrated.t <- c("CD4", "CD8A", "CCR7", "IL7R", 
                           "IFNG", "RORC", "GZMK", "MKI67", 
                           "FOXP3", "IL2RA","TNFRSF18", "PDCD1", 
                           "CTLA4", "HAVCR2", "TIGIT","LAG3",
                           "TCF7", "TOX2", "EOMES", "CD40LG",
                           "IL-17A", "CXCL13", "TOX", "TCF7",
                           "TOX2", "EOMES", "CD40LG", "TGFB1",
                           "IL-17A", "CXCL13", "IRF9", "NR4A1")

features.integrated.t <- c("CD4", "CD8A", "CCR7", "IL7R", 
                           "IFNG", "RORC", "GZMK", 
                           "FOXP3", "IL2RA","TNFRSF18", "PDCD1", 
                           "CTLA4", "HAVCR2", "TIGIT","LAG3")

# Using StackedVlnPlot function in "CellChat" library
#Plot the expression of each of these genes in each cluster.
StackedVlnPlot(immune.combined, features = features.integrated.t)
ggsave("WithTFs-GOIs_Stacked_Violin-A.pdf", width=16, height=12)

#Make a similar plot, except split each violin by Condition
StackedVlnPlot(immune.combined, features = features.integrated.t, split.by = "Condition", color.use = c("royalblue", "red"))
ggsave("WithTFs-GOIs_Condition_Stacked_Violin-A.pdf", width=24, height=12)

features.integrated.t <- c("TCF7", "TOX", "TOX2","EOMES", "CD40LG",
                           "TGFB1", "NR4A1", "MKI67", "IL1A", "IL1B",
                           "IL17A", "CXCL13", "CXCR5")

# Using StackedVlnPlot function in "CellChat" library
#Plot the expression of each of these genes in each cluster.
StackedVlnPlot(immune.combined, features = features.integrated.t)
ggsave("WithTFs-GOIs_Stacked_Violin-B.pdf", width=16, height=12)

#Make a similar plot, except split each violin by Condition
StackedVlnPlot(immune.combined, features = features.integrated.t, split.by = "Condition", color.use = c("royalblue", "red"))
ggsave("WithTFs-GOIs_Condition_Stacked_Violin-B.pdf", width=24, height=12)

features.integrated.t <- c("TCF7", "TOX", "TOX2","EOMES", "CD40LG",
                           "TGFB1", "NR4A1", "MKI67", "IL1A", "IL1B",
                           "IL17A")

# Using StackedVlnPlot function in "CellChat" library
#Plot the expression of each of these genes in each cluster.
StackedVlnPlot(immune.combined, features = features.integrated.t)
ggsave("WithTFs-GOIs_Stacked_Violin-C.pdf", width=16, height=12)

#Make a similar plot, except split each violin by Condition
StackedVlnPlot(immune.combined, features = features.integrated.t, split.by = "Condition", color.use = c("royalblue", "red"))
ggsave("WithTFs-GOIs_Condition_Stacked_Violin-C.pdf", width=24, height=12)

###
#############################################################################################################################################################
###

#3D - The input to this is a subset of the output of IPA.  There is no code to provide for the IPA analysis as it is a GUI.

df <- read.csv("IPA_Th1Th17_SelectRegulators.tsv", sep = "\t", header = TRUE, row.names = 1, comment.char = "")

df$Gene <- rownames(df)

types <- c("cytokine", "transcription regulator", "transmembrane receptor")

for (type in types) {
  
  df.sub <- df[df$Molecule.Type == type, ]
  df.sub <- df.sub[order(df.sub$Activation.z.score), ]
  
  df.sub$Gene <- factor(df.sub$Gene, levels = rownames(df.sub))
  
  if (type == "cytokine") {
    
    ggplot(data = df.sub) +
      geom_col(aes(Activation.z.score, Gene, fill = Predicted.Activation.State), width = 0.6) +
      xlab("Activation z-score") +
      ylab("Upstream Regulator") +
      scale_fill_manual(values = c("blue")) +
      guides(fill=guide_legend(title="Predicted Activation State")) +
      theme(text = element_text(size = 14), axis.text.y = element_text(face = "bold")) + theme_bw()
    #theme(axis.text.y = element_text(size = 14, face = "bold")) 
    
  } else {
    
    ggplot(data = df.sub) +
      geom_col(aes(Activation.z.score, Gene, fill = Predicted.Activation.State), width = 0.6) +
      xlab("Activation z-score") +
      ylab("Upstream Regulator") +
      scale_fill_manual(values = c("red", "blue")) +
      guides(fill=guide_legend(title="Predicted Activation State")) +
      theme(text = element_text(size = 14), axis.text.y = element_text(face = "bold")) + theme_bw()
    #theme(axis.text.y = element_text(size = 14, face = "bold")) 
    
  }
  
  ggsave(paste("UpstreamRegulator_Barplot_", gsub(" ", "-", type), ".pdf", sep = ""), width=6, height=4+nrow(df.sub)*0.25)
  
}

###

df <- read.csv("IPA_Th1Th17_SelectPathways.tsv", sep = "\t", header = TRUE, row.names = 1, comment.char = "")

df <- df[order(df$NegLogP), ]

df$Pathway <- factor(df$Pathway, levels = df$Pathway)

ggplot(data = df) +
  geom_col(aes(NegLogP, Pathway), fill = "forestgreen", width = 0.6) +
  xlab("-log(p-value)") +
  ylab("Pathway") +
  #scale_fill_manual(values = c("green")) +
  #guides(fill=guide_legend(title="Predicted Activation State")) +
  theme_bw() + theme(text = element_text(size = 16, face = "bold"))
#theme(text = element_text(size = 14), axis.text.y = element_text(size = 16, face = "bold"), axis.text.x = element_text(size = 16, face = "bold"))
#theme(text = element_text(size = 14), axis.text.y = element_text(face = "bold")) + theme_bw()
#theme(axis.text.y = element_text(size = 14, face = "bold"))

ggsave(paste("Pathways_Barplot_", "SelectTh1Th17", ".pdf", sep = ""), width=18, height=4+nrow(df)*0.25)

###
#############################################################################################################################################################
###

table(immune.combined$reclustered_clusters)

#0   1   2   3   4   5   6 
#341 259 198 180 117 105  94 

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################