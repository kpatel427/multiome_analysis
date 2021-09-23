# Script to perform integration
# Goals: To align same cell types across conditions.
# setwd("/Volumes/target_nbl_ngs/KP/singleCellProjects/multiomeProject_MatkarS")

library(Seurat)
library(ggplot2)
library(Rcpp)
library(openxlsx)
library(here)
library(dplyr)
library(stringr)
library(tibble)
library(markdown)
library(XML)
library(RCurl)
library(AnnotationHub) # BiocManager::install("AnnotationHub")
library(gridExtra)
library(cowplot)
set.seed(1234)


# read in data -------------------
split_seurat <- readRDS('data/split_seurat.rds')

# In order to perform any downstream analysis on these T cells (e.g. DE analysis between control/treatment, ligand-receptor analysis, etc.) 
# it is generally important that the cells of the same cell types are present in the same clusters.

# Oftentimes, when clustering cells from multiple conditions there are condition-specific clusters and integration can help ensure the same cell types cluster together.

# The goal of integration is to ensure that the cell types of one condition/dataset align with 
# the same celltypes of the other conditions/datasets (e.g. control macrophages align with stimulated macrophages).

# If cell types are present in one dataset, but not the other, then the cells will still appear as a separate sample-specific cluster.


# testing out if cells cluster by condition --------------------
# using seurat_phase obj from previous step
seurat_phase <- RunPCA(seurat_phase, verbose = FALSE)
seurat_phase <- RunUMAP(seurat_phase, dims = 1:30, verbose = FALSE)

seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30, verbose = FALSE)
seurat_phase <- FindClusters(seurat_phase, verbose = FALSE)


a1 <- DimPlot(seurat_phase, label = TRUE)
a2 <- DimPlot(seurat_phase, 
        label = FALSE,
        group.by = 'sample')
a3 <- DimPlot(seurat_phase, 
        label = FALSE,
        group.by = 'seq_folder')


plot <- plot_grid(a1, a2, a3, labels=c("Cluster all cells", "Cluster by Source", "Cluster by Samples"), ncol = 2, nrow = 2)

ggsave(plot, filename = 'figures/clusterCells_by_condition_PDX_beforeIntegration.pdf', width = 20, height = 10)


# integration makes sense because cells cluster together by condition
# Condition-specific clustering of the cells indicates that we need to integrate the cells across conditions to ensure that cells of the same cell type cluster together.
# we will "integrate" or "harmonize" the groups to overlay cells that are similar or have a "common set of biological features" between groups. 

# Select the most variable features to use for integration ---------------
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# run PCA on each object in the list, which is required for running the alternative reciprocal PCA workflow --------------------
split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- ScaleData(x, features = integ_features, verbose = FALSE)
  x <- RunPCA(x, features = integ_features, verbose = FALSE)
})

# Prepare the SCTransform object for integration -------------------------
split_seurat <- PrepSCTIntegration(object.list = split_seurat,
                                   anchor.features = integ_features)



# Perform CCA, find the best buddies or anchors and filter incorrect anchors --------------------
# CCA is prohibitively computationally expensive, hence using rpca reduction for better runtimes
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        reduction = 'rpca',
                                        anchor.features = integ_features)

# Integrate across samples ------------------------------------
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Save integrated seurat object ---------------------------
#saveRDS(seurat_integrated, "results/integrated_seurat.rds")
#seurat_integrated <- readRDS("results/integrated_seurat.rds")

# UMAP visualization - to visualize integrated data -------------------
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA (We want to see a good overlap of both conditions by PCA)
a1 <- PCAPlot(seurat_integrated,
        split.by = "sample",
        group.by = 'sample')  

a2 <- PCAPlot(seurat_integrated,
        split.by = "seq_folder",
        group.by = 'seq_folder')  

# We can see with the PCA mapping that we have a good overlay of both conditions by PCA.

ggsave(a1, filename = 'figures/PCAplot_by_condition_postIntegration.pdf', width = 10, height = 10)
ggsave(a2, filename = 'figures/PCAplot_by_sample_postIntegration.pdf', width = 10, height = 10)


# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP     
d1 <- DimPlot(seurat_integrated, group.by = 'sample')                             
ggsave(d1, filename = 'figures/UMAPplot_groupby_condition_postIntegration.pdf', width = 10, height = 10)

d2 <- DimPlot(seurat_integrated, group.by = 'sample', split.by = 'sample')                             
ggsave(d2, filename = 'figures/UMAPplot_splitby_condition_postIntegration.pdf', width = 10, height = 10)

plot <- plot_grid(d1,d2, labels=c("Grouped by Samples", "Split by Samples"), ncol = 2)
ggsave(plot, filename = 'figures/UMAPplot_by_condition_postIntegration.pdf', width = 20, height = 10)


# Save integrated Seurat Object
saveRDS(seurat_integrated, "results/integrated_seurat.rds")




