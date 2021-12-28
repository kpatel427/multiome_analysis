# script to run doubletFinder to detect doublets in each tumor sample
# setwd("/mnt/isilon/maris_lab/target_nbl_ngs/KP/singleCellProjects/multiomeProject_MatkarS")
# steps:
# - run doublet finder on individual samples (not merged object)
# - filter out low quality cells for each sample
# - preprocess using standard workflow steps
# - run doubletfinder

library(Seurat)
library(ggplot2)
library(Rcpp)
library(DoubletFinder) # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(dplyr)
library(Matrix)
library(fields)
library(KernSmooth)
library(ROCR)
library(parallel)
library(stringr)

set.seed(1234)

# specifying outs directory
outs <- "/mnt/isilon/maris_lab/target_nbl_ngs/Mosse_SMatkar_Multiome/YMosse_SMatkar_Multiome_08172021_release/"

# read in count data and create Seurat objects --------------------------------------------------------------------------------------

# Create a Seurat object for each sample
for(file in c("FelixLRX1_Multiome", "FelixLRX2_Multiome", "FelixParental1_Multiome", "FelixParental2_Multiome")){
  seurat_data <- Read10X(data.dir = paste0(outs,file,"/outs/filtered_feature_bc_matrix"))
  counts <- seurat_data$`Gene Expression`
  seurat_obj <- CreateSeuratObject(counts = counts, 
                                   min.features = 500, 
                                   min.cells = 10,
                                   project = file)
  assign(file, seurat_obj)
  
  
  # filtering low quality cells -----------
  
  # Compute percent mito ratio
  
  mitoRatio <- PercentageFeatureSet(object = get(file), pattern = "^MT-")
  names(mitoRatio) <- 'mitoRatio'
  assign(file,AddMetaData(get(file), mitoRatio))

  
  # filtering
  assign(file, subset(x = get(file), subset = nFeature_RNA > 500
                      & nFeature_RNA < 4000 
                      & nCount_RNA < 16000 
                      & mitoRatio < 2))
  
  
  # preprocess standard workflow ---------------
  assign(file, NormalizeData(get(file)))
  assign(file, FindVariableFeatures(get(file), selection.method = "vst", nfeatures = 2000))
  assign(file, ScaleData(get(file)))
  assign(file, RunPCA(get(file)))
  assign(file, RunUMAP(get(file), dims = 1:20))
  assign(file, FindNeighbors(object = get(file)))
  assign(file, FindClusters(object = get(file)))
  
}


# doubletFinder function ----------------
# Questions to ask Liron-
# What pN and pK value were used?


detectDoublet <- function(seurat.obj){
  
  print("Finding pK...")
  # without ground truth
  seurat.obj.list <- paramSweep_v3(seurat.obj, PCs = 1:20, sct = FALSE)
  sweep.seurat.obj <- summarizeSweep(seurat.obj.list, GT = FALSE)
  bcmvn_seurat.obj <- find.pK(sweep.seurat.obj)
  
  pK <- bcmvn_seurat.obj %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  annotations <- get(file)@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(get(file)@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies 
  print("Run doubletFinder...")
  #assign(seurat.obj,doubletFinder_v3(seurat.obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE))
  seurat.doublets <- doubletFinder_v3(seurat.obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  
  # doublet1 <- DimPlot(FelixLRX1_Multiome, reduction = 'umap', group.by = 'DF.classifications_0.25_0.09_591')
  # doublet2 <- DimPlot(FelixLRX1_Multiome, reduction = 'umap', group.by = 'DF.classifications_0.25_0.09_528')
  # doublet_plots <- gridExtra::grid.arrange(doublet1, doublet2, ncol = 2)
  # ggsave(doublet_plots, filename = 'figures/doubletfinder_plot1.pdf', width = 10, height = 10)
  
  
  
  # create doublet groupings and visualize results
  #DF.class <- names(get(seurat.obj)@meta.data) %>% str_subset("DF.classifications")
  #pANN <- names(get(seurat.obj)@meta.data) %>% str_subset("pANN")
  
  # p1 <- ggplot(bcmvn_felixlrx1, aes(x=pK, y=BCmetric)) +
  #   geom_bar(stat = "identity") + 
  #   ggtitle(paste0("pKmax=",pK)) +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # p2 <- DimPlot(FelixLRX1_Multiome, group.by = DF.class)
  # p3 <- FeaturePlot(FelixLRX1_Multiome, features = pANN)
  # 
  # plots <- gridExtra::grid.arrange(p1,p2,p3, ncol = 2)
  # ggsave(plots, filename = 'figures/doubletfinder_plots.pdf', width = 10, height = 10)
  # 
  
  return(seurat.doublets)
  

}


# run doublet finder --------------

FelixLRX1_Multiome.doublet <- detectDoublet(FelixLRX1_Multiome)
table(FelixLRX1_Multiome.doublet@meta.data$DF.classifications_0.25_0.11_454)
# Doublet Singlet 
# 454    7428 
FelixLRX2_Multiome.doublet <- detectDoublet(FelixLRX2_Multiome)
table(FelixLRX2_Multiome.doublet@meta.data$DF.classifications_0.25_0.12_454)
# Doublet Singlet 
# 454    7333 
FelixParental1_Multiome.doublet <- detectDoublet(FelixParental1_Multiome)
table(FelixParental1_Multiome.doublet@meta.data$DF.classifications_0.25_0.23_454)
# Doublet Singlet 
# 454    7784 
FelixParental2_Multiome.doublet <- detectDoublet(FelixParental2_Multiome)
table(FelixParental2_Multiome.doublet@meta.data$DF.classifications_0.25_0.15_454)
# Doublet Singlet 
# 454    5595 

# filter out doublets -------------------
# before
FelixLRX1_Multiome
# remove doublets
Idents(FelixLRX1_Multiome.doublet) <- FelixLRX1_Multiome.doublet@meta.data$DF.classifications_0.25_0.11_454
FelixLRX1_Multiome.filtered <- subset(FelixLRX1_Multiome.doublet, idents = "Singlet")


# before
FelixLRX2_Multiome
# remove doublets
Idents(FelixLRX2_Multiome.doublet) <- FelixLRX2_Multiome.doublet@meta.data$DF.classifications_0.25_0.12_454
FelixLRX2_Multiome.filtered <- subset(FelixLRX2_Multiome.doublet, idents = "Singlet")



# before
FelixParental1_Multiome
# remove doublets
Idents(FelixParental1_Multiome.doublet) <- FelixParental1_Multiome.doublet@meta.data$DF.classifications_0.25_0.23_454
FelixParental1_Multiome.filtered <- subset(FelixParental1_Multiome.doublet, idents = "Singlet")


# before
FelixParental2_Multiome
# remove doublets
Idents(FelixParental2_Multiome.doublet) <- FelixParental2_Multiome.doublet@meta.data$DF.classifications_0.25_0.15_454
FelixParental2_Multiome.filtered <- subset(FelixParental2_Multiome.doublet, idents = "Singlet")


# merge objects ---------------------------
merged_seurat <- merge(x = FelixLRX1_Multiome.filtered, 
                       y = c(FelixLRX2_Multiome.filtered, FelixParental1_Multiome.filtered, FelixParental2_Multiome.filtered), 
                       add.cell.id = c("LRX1", "LRX2", "parent1", "parent2"))

# Create .RData object to load at any time
save(merged_seurat, file=paste0("data/",Sys.Date(),"_merged_filtered_for_doublets_seurat.RData"))

merged_seurat <- NormalizeData(object = merged_seurat)
merged_seurat <- FindVariableFeatures(object = merged_seurat)
merged_seurat <- ScaleData(object = merged_seurat)
merged_seurat <- RunPCA(object = merged_seurat)
merged_seurat <- FindNeighbors(object = merged_seurat)
merged_seurat <- FindClusters(object = merged_seurat, resolution = 0.8)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:20)

d1 <- DimPlot(object = merged_seurat, reduction = "umap", group.by = 'seurat_clusters')
d2 <- DimPlot(object = merged_seurat, reduction = "umap", group.by = 'orig.ident')
plots <- gridExtra::grid.arrange(d1,d2, ncol = 2)
ggsave(plots, filename = paste0("figures/",Sys.Date(),"_filteration_changes_doublets_removed_umap.pdf"), width = 20, height = 10)






