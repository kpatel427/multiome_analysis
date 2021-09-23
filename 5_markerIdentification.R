# Marker Identification
# Understand how to determine markers of individual clusters
# Understand the iterative processes of clustering and marker identification

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
library(XML)
library(RCurl)
library(AnnotationHub) # BiocManager::install("AnnotationHub")
library(multtest) # BiocManager::install('multtest')
library(metap) # install.packages('metap')
library(reactome.db)
library(org.Hs.eg.db)
library(fgsea)
set.seed(1234)

# Goals:
# To determine the gene markers for each of the clusters
# To identify cell types of each cluster using markers
# To determine whether there's a need to re-cluster based on cell type markers, perhaps clusters need to be merged or split

# Process:
# 1. First identify conserved markers for each cluster grouped by condition
# 2. Do the clusters corresponding to the same cell types have biologically meaningful differences?


# Identification of conserved markers in all conditions -------------------------
# Since we have samples representing different conditions in our dataset, our best option is to find conserved markers.

# For performing differential expression after integration, we switch back to the original data
DefaultAssay(seurat_integrated) <- "RNA"

# testing out on one cluster
markers <- FindConservedMarkers(seurat_integrated, ident.1 = 0, grouping.var = "sample", verbose = FALSE)
head(markers)

# running on all clusters (check total number of clusters by = Idents(seurat_integrated))
load('~/KP/supporting_files/singleCell/human_gene_annotations.RData')

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample") %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# function call
# Iterate function across all clusters
conserved_markers <- map_dfr(c(0:19), get_conserved)
write.table(conserved_markers, file = 'results/conserved_marks_each_cluster.txt', col.names = T, row.names = F, sep = '\t', quote = F)

# Extract top 5 markers per cluster
top5 <- conserved_markers %>%
  mutate(avg_fc = (Parental_avg_log2FC + PDX_avg_log2FC) /2) %>%
  group_by(cluster_id) %>%
  top_n(n = 5,
        wt = avg_fc)

# Visualize top 5 markers per cluster
View(top5)


# Plot interesting marker gene expression for cluster 20
# FeaturePlot(object = seurat_integrated,
#             features = c("ARL6IP1", "CENPF", "CENPE", "MKI67", "ASPM"),
#             order = TRUE,
#             min.cutoff = 'q10',
#             label = TRUE,
#             repel = TRUE)

# Plot top 5 markers for all clusters
for(i in 0:19){
  print(i)
  f1 <- FeaturePlot(object = seurat_integrated,
              features = top5[top5$cluster_id == i, "gene"] %>%
                pull(gene),
              split.by = 'sample',
              label = TRUE)


  v1 <- VlnPlot(object = seurat_integrated,
          features = top5[top5$cluster_id == i, "gene"] %>%
            pull(gene))



  #ggsave(f1, filename = paste0('figures/featurePlot_top5_conservedmarkers_cluster_',i,'.pdf'), width = 15, height = 10)
  ggsave(v1, filename = paste0('figures/VlnPlot_top5_conservedmarkers_cluster_',i,'.pdf'), width = 15, height = 10)
}

# generate a dotplot with top 5 markers from each cluster
markers.to.plot <- c(unique(conserved_markers %>%
                    mutate(avg_fc = (Parental_avg_log2FC + PDX_avg_log2FC) /2) %>%
                    group_by(cluster_id) %>%
                    top_n(n = 1,
                          wt = avg_fc) %>%
                    dplyr::select(gene)))
markers.to.plot <- unique(markers.to.plot$gene)
dp1 <- DotPlot(seurat_integrated, features = markers.to.plot, cols = c("blue", "red"), split.by = "sample") +
  RotatedAxis() + ggtitle('Top Conserved Marker DE in each cluster')



# Identification of all markers for each cluster -------------------------

# Find markers for every cluster compared to all remaining cells, report only the positive ones
combined_markers <- FindAllMarkers(object = seurat_integrated,
                                   only.pos = TRUE,
                                   logfc.threshold = log(2))

View(combined_markers)

# Combine markers with gene descriptions
ann_comb_markers <- inner_join(x = combined_markers,
                               y = annotations[, c("gene_name", "description")],
                               by = c("gene" = "gene_name")) %>%
                    unique()

# Rearrange the columns to be more intuitive
ann_comb_markers <- ann_comb_markers[ , c(6, 7, 2:4, 1, 5,8)]

# Order the rows by p-adjusted values
ann_comb_markers <- ann_comb_markers %>%
  dplyr::arrange(cluster, p_val_adj)


# generate a dotplot with top 5 markers from each cluster
comb.markers.to.plot <- c(unique(ann_comb_markers %>%
                              group_by(cluster) %>%
                              top_n(n = 1,
                                    wt = avg_log2FC) %>%
                              dplyr::select(gene)))
comb.markers.to.plot <- unique(comb.markers.to.plot$gene)
dp2 <- DotPlot(seurat_integrated, features = comb.markers.to.plot, cols = c("blue", "red"), split.by = "sample") +
  RotatedAxis() + ggtitle('Top Marker DE in each cluster')


dp <- grid.arrange(dp1, dp2)
ggsave(dp, filename = 'figures/DotPlot_markers_DE_each_cluster.pdf', width = 20, height = 15)


# Identifying gene markers for each cluster --------------------------------------
# Using the FindMarkers() function to determine the genes that are differentially expressed between two specific clusters.
# should it be compared against same clusters in each condition or different clusters across conditions?


# testing it out on one cluster

seurat_integrated$celltype.sample <- paste(Idents(seurat_integrated), seurat_integrated$sample, sep = "_")
seurat_integrated$celltype <- Idents(seurat_integrated)
Idents(seurat_integrated) <- "celltype.sample"
DE.markers.allClusters <- FindMarkers(seurat_integrated, ident.1 = "0_PDX", ident.2 = "0_Parental", verbose = FALSE, logfc.threshold = log(2), min.pct = 0.25)
head(DE.markers.allClusters, n = 15)

plots <- VlnPlot(seurat_integrated, features = c("ALCAM", "KAZN", "PCDH7"), split.by = "sample", group.by = 'sample',
        pt.size = 0, combine = FALSE)
grid.arrange(grobs=plots)





# in a loop for all clusters
seurat_integrated$celltype.sample <- paste(Idents(seurat_integrated), seurat_integrated$sample, sep = "_")
seurat_integrated$celltype <- Idents(seurat_integrated)
Idents(seurat_integrated) <- "celltype.sample"
DE.allClusters <- data.frame()
for (i in 0:19){
  print(i)
  try({
    ident1 <- paste0(i,"_PDX")
    ident2 <- paste0(i,"_Parental")
    condition.diffgenes <- FindMarkers(seurat_integrated, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold = log(2))
    condition.diffgenes$cluster <- i
    DE.allClusters <- rbind(DE.allClusters, condition.diffgenes)
  })
}

DE.allClusters <- DE.allClusters %>%
  rownames_to_column(var = 'gene')
write.table(DE.allClusters, file = 'results/DE.features.allClusters.txt', col.names = T, row.names = F, sep = '\t', quote = F)


# performing GSEA using conserved markers ------------------------
conserved_markers <- read.delim('results/conserved_marks_each_cluster.txt', header = T)
# filter out genes withy non-significant P-vals (p-adjusted < 0.05)
conserved_markers <- conserved_markers[conserved_markers$Parental_p_val_adj < 0.05 & conserved_markers$PDX_p_val_adj < 0.05,]

# ...convert gene symbol to EntrezID -----------
hs <- org.Hs.eg.db
my.symbols <- c(as.character(conserved_markers$gene))
df <- AnnotationDbi::select(hs,
                            keys = my.symbols,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")


# merge df with df.res
conserved_markers_gsea <- merge(conserved_markers, df, by.x = 'gene', by.y = 'SYMBOL')
conserved_markers_gsea <- na.omit(conserved_markers_gsea)
# subset
conserved_markers_gsea <- conserved_markers_gsea[,c(1:4,7:9,12,15,16)]

# generate average log fold changes in both Parental and PDX
conserved_markers_gsea$avg_log2FC <- conserved_markers_gsea$Parental_avg_log2FC + conserved_markers_gsea$PDX_avg_log2FC / 2


# for each cluster, rank genes
out.gsea <- conserved_markers_gsea %>%
  dplyr::select(1,2,9:11) %>%
  arrange(cluster_id) %>%
  group_by(cluster_id) %>%
  group_map(~performGSEA(.x$avg_log2FC, .x$ENTREZID, .y$cluster_id))


performGSEA <- function(lfc, entrzid, cluster){
  print(cluster)
  # create ranks
  ranks <- lfc
  names(ranks) <- entrzid
  head(ranks)

  # plot ranked fold changes
  #barplot(sort(ranks, decreasing = T))

  # load pathways
  my_pathways <- reactomePathways(names(ranks))

  out <- fgsea(pathways = my_pathways,
        stats = ranks,
        minSize=15,
        maxSize=500,
        nperm=100000)

  # plotting gsea table
  topPathwaysUp <- out[ES > 0][head(order(pval), n=20), pathway]
  topPathwaysDown <- out[ES < 0][head(order(pval), n=20), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  pdf(paste0("figures/GSEA_cluster_",cluster,"_table_genes_top20_upDown_regulated_pathways.pdf"), width = 20, height = 10)
  plot.new()
  plotGseaTable(my_pathways[topPathways], ranks, out,
                gseaParam = 0.5)
  dev.off()

  return(out)
  print("done!")
}

# spot check! - checks out
# test <- conserved_markers_gsea %>%
#   group_by(cluster_id) %>%
#   filter(cluster_id == 0)
#
# test <- test[,c(1,2,9:11)]
# test <- test[!duplicated(test),]
# ranks <- test$avg_log2FC
# names(ranks) <- test$ENTREZID
#
# # load pathways
# my_pathways <- reactomePathways(names(ranks))
#
# test.gsea <- fgsea(pathways = my_pathways,
#       stats = ranks,
#       minSize=15,
#       maxSize=500,
#       nperm=100000)
