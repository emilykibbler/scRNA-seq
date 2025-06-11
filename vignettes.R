## PBMC 1k ----------------
# Following this vignette tutorial:
# https://satijalab.org/seurat/articles/pbmc3k_tutorial

source("libraries.R")
load_packages()

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./data/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k", 
                           min.cells = 3, 
                           min.features = 200)
pbmc


# My demux of similar data: folder output
# See cellranger_intro.sh for details
pbmc.data <- Read10X(data.dir = "./data/raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc.new <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k", 
                           min.cells = 3, 
                           min.features = 200)
pbmc.new
# looks good

# Single file, from my cell ranger run, to demonstrate Read10X_h5
pbmc.data.new <- Read10X_h5("./data/raw_feature_bc_matrix.h5")
pbmc.new <- CreateSeuratObject(counts = pbmc.data.new, 
                           project = "pbmc3k", 
                           min.cells = 3, 
                           min.features = 200)
pbmc.new
# looks good as well

# Following tutorial again
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", 
                           "nCount_RNA", 
                           "percent.mt"), 
        ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Make plots a little larger by having only one legend
ggarrange(plotlist = list(plot1, plot2),
          legend = "bottom",
          common.legend = TRUE)

### Filter, normalize, and scale ---------------

# Cells with >2500 features (unique genes) are likely to have >1 cell in the GEM
# Cells with <200 features are poor quality
# Cells with high fraction of mitochondrial DNA are also low quality--likely dying
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# Plot again, after filtering
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggarrange(plotlist = list(plot1, plot2),
          legend = "bottom",
          common.legend = TRUE)

VlnPlot(pbmc, features = c("nFeature_RNA", 
                           "nCount_RNA", 
                           "percent.mt"), 
        ncol = 3)

# Normalize data
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", # default value
                      scale.factor = 10000) # also default value

# Note that there is another function, SCTransform(),
# which does not assume every cell has the same amount of RNA
# and could be used instead for a custom workflow

pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", # default
                             nfeatures = 2000) # default

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggarrange(plotlist = list(plot1, plot2),
          common.legend = TRUE)

all.genes <- rownames(pbmc)
# Scaling: mean expression is 0, variance across cell is 1
# Eliminates skew coming from highly expressed genes
pbmc <- ScaleData(pbmc, features = all.genes)

### PCA --------------

pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(pbmc)
# The "elbow," point of diminishing returns, happens between the 7th and 12th component
# The tutorial chooses 10 as the cutoff

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

### UMAP ------------

pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap", label = TRUE)

dir.create("./data/output")
saveRDS(pbmc, file = "./data/output/pbmc_tutorial.rds")

### Finding differentially expressed features (cluster biomarkers) -----------------

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc, 
                                ident.1 = 0, 
                                logfc.threshold = 0.25, 
                                test.use = "roc", # finds the "classification power" of any feature
                                only.pos = TRUE)

# Expression probability distributions across clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# Visualizes feature expression on tSNE or PCA plot
# Tutorial recommends these to look at
recommended_features <- c("MS4A1", "GNLY", "CD3E",
             "CD14", "FCER1A", "FCGR3A",
             "LYZ", "PPBP", "CD8A")

FeaturePlot(pbmc, features = recommended_features)

RidgePlot(pbmc, features = recommended_features)

# Compare cells, call by barcode sequence
CellScatter(pbmc, 
            cell1 = #FIXME, 
            cell2 = #FIXME, 
            features = recommended_features)

DotPlot(pbmc, features = recommended_features)

# expression heatmap for given cells and features
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", 
                     "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()

# Markers and corresponding cell types are known in this case
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T",
                     "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()


plot <- DimPlot(pbmc, 
                reduction = "umap", 
                label = TRUE, 
                label.size = 4.5) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  theme(axis.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))
plot

dir.create("./data/output/images")
ggsave(filename = "./data/output/images/pbmc3k_umap.jpg", 
       height = 7, 
       width = 12, 
       plot = plot, 
       quality = 50)
saveRDS(pbmc, file = "./data/output/pbmc3k_final.rds")



