# Script to preprocess and QC snRNASeq data
# setwd("/Volumes/target_nbl_ngs/KP/singleCellProjects/multiomeProject_MatkarS")

library(Seurat)
library(ggplot2)
library(harmony) # 1.0
library(Rcpp)
library(DoubletFinder) # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(openxlsx)
library(here)
library(dplyr)
library(stringr)
library(tibble)
library(markdown)
set.seed(1234)

# specifying outs directory
outs <- "/Volumes/target_nbl_ngs/Mosse_SMatkar_Multiome/YMosse_SMatkar_Multiome_08172021_release/"

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
}

# Check the metadata in the new Seurat objects
head(FelixLRX1_Multiome@meta.data)
head(FelixParental1_Multiome@meta.data)

# Create a merged Seurat object
# This will make it easier to run the QC steps for both sample groups together and enable us to easily compare the data quality for all the samples.
merged_seurat <- merge(x = FelixLRX1_Multiome, 
                       y = c(FelixLRX2_Multiome, FelixParental1_Multiome, FelixParental2_Multiome), 
                       add.cell.id = c("LRX1", "LRX2", "parent1", "parent2"))


# let's look at metadata for merged_seurat
tail(merged_seurat@meta.data)
View(merged_seurat@meta.data)


# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)


# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^LRX"))] <- "PDX"
metadata$sample[which(str_detect(metadata$cells, "^parent"))] <- "Parental"

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")


# Assessing Quality metrics --------------------------------------------------------------------------------------

# to generate a QC report 
#install.packages("tinytex")
#tinytex::install_tinytex()

# Run generateReport.R to generate a HTML QC report


# 1. cell counts
# Visualize the number of cell counts per sample
p1 <- metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

p2 <- metadata %>% 
  ggplot(aes(x=sample, fill=seq_folder)) + 
  geom_bar(position = 'dodge') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

p1 + p2
# We see higher number (over 15,000) of cells (Parental Samples have even higher number of cells compared to PDX), 
# which is quite a bit more than the 12-13,000 expected. It is clear that we likely have some junk 'cells' present.


# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500, that is the low end of what we expect. 
# If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)



# genes detected per cell
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500)
# Majority of the cells appears to have a high gene count

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NGenes: Parental and PDX")

# UMIs vs. genes detected
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point(alpha = 0.4) + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Cells that are poor quality are likely to have low genes and UMIs per cell, and correspond to the data points in the bottom left quadrant of the plot. 
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs.


# mitochondrial counts ratio
# This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells. 
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
# Majority of cells have way under 20% of mitochondrial ratio (we have used a cutoff of 0.2)


# complexity
# We can evaluate each cell in terms of how complex the RNA species are by using a measure called the novelty score. 
# The novelty score is computed by taking the ratio of nGenes over nUMI.
# we expect the novelty score to be above 0.80 for good quality cells.
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
# Almost all our cells seem to have a good complexity.

preQC_meta <- metadata

# Filtering --------------------------------------
# ...Cell-level Filtering ---------

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, subset = nGene > 500
                          & nGene < 4000 
                          & nUMI < 16000 
                          & mitoRatio < 0.2)


# ...Gene-level Filtering ---------
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene (to keep genes expressed in 10 or more cells)
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)


# Create .RData object to load at any time
save(filtered_seurat, file="data/seurat_filtered.RData")
#load("data/seurat_filtered.RData")

postQC_meta <- filtered_seurat@meta.data


















