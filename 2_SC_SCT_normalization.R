# Script to normalize and perform PCA
  # |__Understand normalizing counts is necessary for accurate comparison between cells
  # |__Understand how similarity in cellular gene expression between cells can be evaluated by Principal Components Analysis (PCA)
  # |__ Executing the normalization, variance estimation, and identification of the most variable genes for each sample

# Goals:
# To accurately normalize and scale the gene expression values to account for differences in sequencing depth and overdispersed count values.
# To identify the most variant genes likely to be indicative of the different cell types present.

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
set.seed(1234)

#load('data/seurat_filtered.RData')
str(filtered_seurat)

# Exploring Sources of unwanted variation
# 1. Normalize the counts -------------------------
seurat_phase <- NormalizeData(filtered_seurat)
str(seurat_phase)


# 2. Evaluating effects of cell cycle --------------
load('~/KP/supporting_files/singleCell/cell_cycle_markers.RData')

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Score cells for cell cycle
# Assign each cell a score based on its expression og G2/M and S phase markers. The function calculates cell cycle phase scores based on canonical markers
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells
View(seurat_phase@meta.data)


# Next we would like to determine whether cell cycle is a major source of variation in our dataset using PCA
# To perform PCA: Find most variable features > scale data > runPCA
# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase,
                                     selection.method = "vst",
                                     nfeatures = 2000,
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)


# Perform PCA to evaulate similarities/differences between cell cycle phase
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
cc1 <- DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")

cc2 <- DimPlot(seurat_phase,
               reduction = "pca",
               group.by= "Phase",
               split.by = "Phase")
cc3 <- DimPlot(seurat_phase,
               reduction = "pca",
               group.by= "Phase",
               split.by = "sample")
plot <- plot_grid(cc1, cc2, cc3, labels=c("Grouped all cells together", "Split by Cell Cycle Phase", "Split by Source"), ncol = 2, nrow = 2)


ggsave(plot, filename = 'figures/dimplot_cellCycle_before_regressing.pdf', width = 20, height = 10)

# Since cells separate entirely by phase, we need to regress out cell cycle scores
# regressing it out (Takes a long time!)
# seurat_phase <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_phase))
# regressed_cc_seurat_phase <- seurat_phase

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
# regressed_cc_seurat_phase <- RunPCA(regressed_cc_seurat_phase, features = VariableFeatures(regressed_cc_seurat_phase), nfeatures.print = 10)


# Plot the PCA colored by cell cycle phase post-regression
# plot2 <- DimPlot(regressed_cc_seurat_phase,
#         reduction = "pca",
#         group.by= "Phase")
# ggsave(plot2, filename = 'figures/dimplot_cellCycle_post_regressing.pdf', width = 10, height = 10)
#
# saveRDS(regressed_cc_seurat_phase, "data/regressed_cc_seurat_phase.rds")


# 3. Evaluating effects of mitochondrial gene expression --------------
# Next, we want to evaluate expression of mitochondrial genes' expression
# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio,
                                     breaks=c(-Inf, 0.003867, 0.011655, 0.013933, Inf),
                                     labels=c("Low","Medium","Medium high", "High"))



# Plot the PCA colored by mitoFr
# Plot the PCA colored by cell cycle phase
m1 <- DimPlot(seurat_phase,
               reduction = "pca",
               group.by= "mitoFr")

m2 <- DimPlot(seurat_phase,
               reduction = "pca",
               group.by= "mitoFr",
               split.by = "mitoFr")
m3 <- DimPlot(seurat_phase,
               reduction = "pca",
               group.by= "mitoFr",
               split.by = "sample")
plot <- plot_grid(m1, m2, m3, labels=c("Grouped all cells together", "Split by Mito Fractions", "Split by Source"), ncol = 2, nrow = 2)


ggsave(plot, filename = 'figures/dimplot_mitoFr_before_regressing.pdf', width = 20, height = 10)

# since they show different pattern of scatter, mitoRatio will be required to be regressed out

# 4. SCTransform --------------
# SCTransform automatically accounts for cellular sequencing depth by regressing out sequencing depth (nUMIs)
# Split seurat object by condition to regress cell cycle scoring, mitoRatio and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "seq_folder")
options(future.globals.maxSize = 4000 * 1024^2)


for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("S.Score","G2M.Score","mitoRatio"))
}


# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat.rds")


# Split seurat object by condition to perform cell cycle scoring and SCT on conditions (PDX & Parental)

# Important NOTE: By default, after normalizing, adjusting the variance, and regressing out uninteresting sources of variation,
# SCTransform will rank the genes by residual variance and output the 3000 most variant genes. If the dataset has larger cell numbers,
# then it may be beneficial to adjust this parameter higher using the variable.features.n argument.

# split_seurat <- SplitObject(seurat_phase, split.by = "sample")
# options(future.globals.maxSize = 4000 * 1024^2)
#
#
# for (i in 1:length(split_seurat)) {
#   split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("S.Score","G2M.Score"))
# }
#
# # finding variable features after SCTransform
# # split_seurat$PDX@assays$SCT@var.features[1:20]
# # split_seurat$Parental@assays$SCT@var.features[1:20]
#
# # Save the split seurat object
# saveRDS(split_seurat, "data/split_seurat.rds")


#' #' 4a. SCTransform for regressed_cc_seurat_phase (using seurat obj which is manually regressed for cc genes)--------------
#' regressed_cc_seurat_phase <- readRDS("data/regressed_cc_seurat_phase.rds")
#' # split by conditions (PDX & Parental)
#' split_regressed_cc_seurat_phase <- SplitObject(regressed_cc_seurat_phase, split.by = "sample")
#' options(future.globals.maxSize = 4000 * 1024^2)
#'
#'
#' for (i in 1:length(split_regressed_cc_seurat_phase)) {
#'   split_regressed_cc_seurat_phase[[i]] <- SCTransform(split_regressed_cc_seurat_phase[[i]])
#' }
#'
#' # finding variable features after SCTransform
#' #split_regressed_cc_seurat_phase$PDX@assays$SCT@var.features[1:20]
#' #split_regressed_cc_seurat_phase$Parental@assays$SCT@var.features[1:20]
#'
#' # Save the split seurat object
#' saveRDS(split_regressed_cc_seurat_phase, "data/split_regressed_cc_seurat_phase.rds")
#'
