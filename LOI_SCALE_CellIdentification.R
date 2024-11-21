# This Rscript processes single-cell RNA-seq count matrices into a combined dataset for downstream analysis.
#
# Key steps include:
# 1. Conversion of raw count matrices into Seurat objects with filtering for high-quality cells.
# 2. Normalization, identification of highly variable features, and scaling of data for each sample.
# 3. Merging of individual samples into a single dataset and correction of batch effects using Harmony.
# 4. Dimensionality reduction (PCA, UMAP, t-SNE) and clustering of cells to identify subpopulations.
# 5. Integration of marker genes with variable features to refine cell type identification using scSorter.
# 6. Fine-tuning UMAP parameters to optimize cluster granularity and visualize cell types.



# Set working directory
setwd("path/to/your/data")

# Load required packages
packages <- c("RColorBrewer", "Seurat", "SeuratDisk", "SingleCellExperiment", 
              "scater", "scran", "harmony", "scSorter", "dplyr")
lapply(packages, library, character.only = TRUE)

# Set seed for reproducibility
set.seed(6666)

# Function to preprocess individual samples
preprocess_sample <- function(count_matrix, tech_label) {
  # Create Seurat object from count matrix
  seurat_obj <- CreateSeuratObject(counts = count_matrix, project = tech_label)
  
  # Filter cells based on feature and count thresholds
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA < 20000 & nCount_RNA > 200 & nFeature_RNA > 200)
  
  # Normalize, find variable features, and scale data
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  
  # Add technical label
  seurat_obj[["tech"]] <- tech_label
  return(seurat_obj)
}

# Load count matrices and preprocess samples
sample_files <- list.files("path/to/counts", pattern = "*.txt", full.names = TRUE) # Adjust path for public use
sample_labels <- gsub(".txt", "", basename(sample_files))

# Apply preprocessing to all samples
samples <- lapply(seq_along(sample_files), function(i) {
  counts <- read.table(sample_files[i], header = TRUE, row.names = 1)
  preprocess_sample(count_matrix = counts, tech_label = sample_labels[i])
})

# Merge all samples into a combined dataset
combined <- Reduce(function(x, y) merge(x, y), samples)

# Save the combined dataset as an H5Seurat object
#SaveH5Seurat(object = combined, filename = "breast_all.h5seurat")
# Reload the combined dataset
#combined <- LoadH5Seurat("breast_all.h5seurat")

# Calculate mitochondrial percentage
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# Plot feature relationships
FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter combined dataset for high-quality cells
filtered_combined <- subset(combined, nCount_RNA > 200)

# Normalize and scale data
filtered_combined <- NormalizeData(filtered_combined, normalization.method = "LogNormalize")
filtered_combined <- FindVariableFeatures(filtered_combined, selection.method = "vst")
filtered_combined <- ScaleData(filtered_combined, model.use = "negbinom")

# Run PCA and visualize
filtered_combined <- RunPCA(filtered_combined)
ElbowPlot(filtered_combined, ndims = 20)

# Dimensionality reduction
filtered_combined <- RunUMAP(filtered_combined, dims = 1:20, verbose = FALSE)
filtered_combined <- RunTSNE(filtered_combined, dims = 1:20, verbose = FALSE)

# Batch correction with Harmony
combined_harmony <- RunHarmony(filtered_combined, group.by.vars = 'tech', plot_convergence = FALSE)

# Dimensionality reduction using Harmony embeddings
combined_harmony <- combined_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  RunTSNE(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# Visualize UMAP before and after Harmony
before <- DimPlot(filtered_combined, reduction = 'umap', group.by = 'tech') + ggtitle('Before Batch Effect Correction')
after <- DimPlot(combined_harmony, reduction = 'umap', group.by = 'tech') + ggtitle('After Batch Effect Correction')
before | after

# Marker genes and scSorter
anno <- read.csv("path/to/markers.csv", sep = ";") # Replace path
expr <- combined_harmony@assays$RNA@data
topgenes <- head(combined_harmony@assays$RNA@var.features, 2000)
topgene_filter <- rowSums(expr[topgenes, ] != 0) > ncol(expr) * 0.1
topgenes <- topgenes[topgene_filter]

picked_genes <- unique(c(anno$Marker, topgenes))
expr_subset <- expr[rownames(expr) %in% picked_genes, ]
rts <- scSorter(expr_subset, anno)

# Add cell type information
combined_harmony$cell_subtypes <- rts$Pred_Type
DimPlot(combined_harmony, group.by = 'cell_subtypes') + ggtitle("Cell Subtypes")

# Fine-tune UMAP parameters
combined_harmony <- RunUMAP(
  combined_harmony,
  reduction = "pca", 
  dims = 1:10,
  umap.param = list(n_neighbors = 20, min_dist = 0.5, metric = "cosine")
)
DimPlot(combined_harmony, group.by = 'seurat_clusters', label = TRUE)
