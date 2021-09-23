# Utilize methods for evaluating the selection of PCs to use for clustering
# Perform clustering of cells based on significant PCs
# setwd("~/KP/singleCellProjects/multiomeProject")

library(Seurat)
library(ggplot2)
library(harmony) # 1.0
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
library(purrr)
set.seed(1234)

# Now that we have our high quality cells integrated, we want to know the different cell types present within our population of cells.
# Determining how many PCs to include in the clustering step is
# therefore important to ensure that we are capturing the majority of the variation, or cell types, present in our dataset.

seurat_integrated <- readRDS('results/integrated_seurat.rds')

# Elbow Plot --------------------------------------------------
ElbowPlot(object = seurat_integrated,
          ndims = 43)

# quantitative approach to choose PCs
# We can calculate where the principal components start to elbow by taking the larger value of:
#
# 1. The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
# 2. The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs

# Create a dataframe with values
plot_df <- data.frame(pct = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


# Cluster the cells ---------------------------------------------
# The resolution is an important argument that sets the "granularity" of the downstream clustering and will need to be optimized for every individual experiment.
# For datasets of 3,000 - 5,000 cells, the resolution set between 0.4-1.4 generally yields good clustering.
# Increased resolution values lead to a greater number of clusters, which is often required for larger datasets.


# Determine the K-nearest neighbor graph
# Using first 40 PCs to cluster
seurat_integrated <- FindNeighbors(object = seurat_integrated,
                                   dims = 1:40)

# Determine the clusters for various resolutions
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))


# ....Explore resolutions ---------------
seurat_integrated@meta.data %>%
  View()


# ....Assign identity of clusters ----------------------
# We will start with a resolution of 0.8 by assigning the identity of the clusters using the Idents() function.
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
umap_0.8 <- DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6) + ggtitle('Resolution = 0.8')

grid.arrange(umap_0.4, umap_0.6, umap_0.8, umap_1, umap_1.4)

# Going ahead with 0.8 resolution
ggsave(umap_0.8, filename = 'figures/UMAPplot_clustering_integrated_object.pdf', width = 10, height = 10)



# Clustering QC ---------------------------------------------
# Goals:
# Evaluate whether clustering artifacts are present
# Determine the quality of clustering with PCA and UMAP plots and understand when to re-cluster
# Assess known cell type markers to hypothesize cell type identities of clusters


# ....Segregation of clusters by sample -----------------
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
# n_cells <- FetchData(seurat_integrated,
#                      vars = c("ident", "orig.ident")) %>%
#   dplyr::count(ident, orig.ident) %>%
#   tidyr::spread(ident, n)
#
# # View table
# View(n_cells)

# UMAP of cells in each cluster by sample
qc1 <- DimPlot(seurat_integrated,
        label = TRUE,
        split.by = "sample")  + NoLegend()

ggsave(qc1, filename = 'figures/UMAP_QC1_postClustering_cluster_by_sample.pdf', width = 10, height = 10)

# Generally, we expect to see the majority of the cell type clusters to be present in all conditions;
# however, depending on the experiment we might expect to see some condition-specific cell types present.
# These clusters look pretty similar between conditions, which is good since we expected similar cell types to be present in both PDX and Parental tumors


# ....Segregation of clusters by cell cycle phase -----------------

# Explore whether clusters segregate by cell cycle phase
qc2 <- DimPlot(seurat_integrated,
        label = TRUE,
        split.by = "Phase")  + NoLegend()

# !!!
ggsave(qc2, filename = 'figures/UMAP_QC2_postClustering_cluster_by_cellCyclePhase.pdf', width = 10, height = 10)


# ....Segregation of clusters by various sources of uninteresting variation -----------------

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

qc3 <- FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

ggsave(qc3, filename = 'figures/FeaturePlot_QC3_postClustering.pdf', width = 10, height = 10)


# ....Exploration of the PCs driving the different clusters -----------------
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated,
                     vars = columns)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated,
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data,
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc),
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE,
                         low = "grey90",
                         high = "blue")  +
    geom_text(data=umap_label,
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>%
  plot_grid(plotlist = .)

# Examine PCA results
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)




# ....Exploration of known cell type markers -----------------
aa1 <- DimPlot(object = seurat_integrated,
        reduction = "umap",
        label = TRUE)

# The SCTransform normalization was performed only on the 3000 most variable genes, so many of our genes of interest may not be present in this data.
# Hence need to normalize data of genes not present in 3000 most variable genes
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)


FeaturePlot(seurat_integrated,
            reduction = "umap",
            features =c('HBB', 'HBA2'),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
# unable to detect most commonly known markers
