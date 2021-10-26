# script to perform trajectory analysis
# setwd("~/KP/singleCellProjects/multiomeProject_MatkarS")
library(tidyverse)
library(monocle)
library(Seurat)
library(devtools)
library(roxygen2)
# install_github('hb-gitified/cellrangerRkit', user = 'kpatel427',
#                 auth_token = '' )
library(cellrangerRkit)
#BiocManager::install("slingshot")
library(slingshot)
library(RColorBrewer)
library(scales)
library(pals)
library(grDevices)
library(ggbeeswarm)
library(tradeSeq)
library(pheatmap)
library(reactome.db)
library(org.Hs.eg.db)
library(fgsea)




# reading in split_seurat after SCTransform ------------------
split_seurat <- readRDS("data/split_seurat.rds")

FelixLRX1_Multiome <- split_seurat$FelixLRX1_Multiome
FelixLRX2_Multiome <- split_seurat$FelixLRX2_Multiome
FelixParental1_Multiome <- split_seurat$FelixParental1_Multiome
FelixParental2_Multiome <- split_seurat$FelixParental2_Multiome

# merging seurat object post SCTransormation into one object
merged_seurat <- merge(x = FelixLRX1_Multiome, 
                       y = c(FelixLRX2_Multiome, FelixParental1_Multiome, FelixParental2_Multiome), 
                       add.cell.id = c("LRX1", "LRX2", "parent1", "parent2"))

saveRDS(merged_seurat, file = paste0('data/',Sys.Date(),'_beforeIntegration.rds'))

View(merged_seurat@meta.data)


# just a quick check!
# getting variable features from SCT assay
VariableFeatures(merged_seurat[["SCT"]]) <- rownames(merged_seurat[["SCT"]]@scale.data)

merged_seurat <- RunPCA(merged_seurat, verbose = FALSE)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:30, verbose = FALSE)

merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30, verbose = FALSE)
merged_seurat <- FindClusters(merged_seurat, verbose = FALSE)

merged_seurat$sample <- gsub('PDX','LRX', merged_seurat$sample)

a1 <- DimPlot(merged_seurat, label = TRUE)
a2 <- DimPlot(merged_seurat, 
              label = FALSE,
              group.by = 'sample')
a3 <- DimPlot(merged_seurat, 
              label = FALSE,
              group.by = 'seq_folder')

plot <- plot_grid(a1, a2, a3, labels=c("Cluster all cells", "Cluster by Source", "Cluster by Samples"), ncol = 2, nrow = 2)

ggsave(plot, filename = paste0('figures/',Sys.Date(),'_postSCTransform_preIntegration_clustering.pdf'), width = 15, height = 10)
saveRDS(merged_seurat, file = paste0('data/',Sys.Date(), '_postSCTransform_beforeIntegration.rds'))

# reading in seurat objects pre and post SCTransform ---------------
post.sct <- readRDS('~/KP/singleCellProjects/multiomeProject_MatkarS/data/2021-10-11_postSCTransform_beforeIntegration.rds')
pre.sct <- readRDS('~/KP/singleCellProjects/multiomeProject_MatkarS/data/2021-10-11_preSCTransform_beforeIntegration.rds')


# getting overlapping ADRN and MES genes
adrn.mes.sig <- read.delim('~/KP/singleCellProjects/PRAME_ARDN_MES/ADRN_MES_signature_genes_list.txt', header = F)
mesench <- data.frame('V1'=adrn.mes.sig[which(adrn.mes.sig[,2] == "MES"),1])
adrenal <- data.frame('V1'=adrn.mes.sig[which(adrn.mes.sig[,2] == "ADRN"),1])
mesench$group <- 'MES'
adrenal$group <- 'ADRN'

ADRN.overlap <- unique(intersect(adrenal$V1, VariableFeatures(pre.sct)))
MES.overlap <- unique(intersect(mesench$V1, VariableFeatures(pre.sct)))



# starting with pre-SCT data
class(pre.sct)
head(rownames(pre.sct))
head(colnames(pre.sct))

# visualizations
t <- ggplot(pre.sct@meta.data, aes(nCount_RNA, nFeature_RNA, color = seq_folder)) +
  geom_point(size = 0.5) +
  scale_color_brewer(type = "qual", palette = "Set2", name = "cell type") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  # Make points larger in legend
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(x = "Total UMI counts", y = "Number of genes detected")

ggsave(t, filename = 'figures/test.pdf', width = 15, height = 10)

# looking at PCs
e <- ElbowPlot(pre.sct, ndims = 70)
ggsave(e, filename = 'figures/elbowplot.pdf', width = 15, height = 10)

# looking at reduced dimensions and clustering
d <- DimPlot(pre.sct, reduction = "umap",
             group.by = "seq_folder", pt.size = 0.5, label = TRUE, repel = TRUE) +
  scale_color_brewer(type = "qual", palette = "Set2")
ggsave(d, filename = 'figures/dimplot.pdf', width = 15, height = 10)


# Slingshot -----------------------------------
# convert seurat to singleCellExperiment obj
pre.sct.sc <- as.SingleCellExperiment(pre.sct)
post.sct.sc <- as.SingleCellExperiment(post.sct)


# running slingshot ------------------------------
# ....using seurat object to calculate lineages -------------------
pre.sct.ss <- slingshot(Embeddings(pre.sct, "umap"), clusterLabels = pre.sct$seurat_clusters, 
                        start.clus = 4, stretch = 0)

post.sct.ss <- slingshot(Embeddings(post.sct, "umap"), clusterLabels = post.sct$seurat_clusters, 
                        start.clus = 10, stretch = 0, omega = TRUE)

# to test whether we have disjoint trajectories - no disjoint trajectories detected
# to give a different start cluster - 7: Still does not change the lineage structure
pre.sct.ss.test <- slingshot(Embeddings(pre.sct, "umap"), clusterLabels = pre.sct$seurat_clusters, 
                        start.clus = 7, stretch = 0, omega = TRUE)

# using singleCellExperiment object
#sce <- slingshot(pre.sct.sc, clusterLabels = pre.sct.sc$seurat_clusters, reducedDim = 'PCA')

# assign colors to cells
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

# pre.SCT
cell_colors <- cell_pal(pre.sct$seq_folder, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(pre.sct$seurat_clusters, hue_pal())

# post.SCT
cell_colors <- cell_pal(post.sct$seq_folder, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(post.sct$seurat_clusters, hue_pal())

# ...plotting trajectory ------------------
# for seurat object
plot(Embeddings(pre.sct, "umap"), col = cell_colors, pch = 16, cex = 0.5)


pdf(file = 'figures/trajectoryAnalysis_slingshot_trajectories_clusters.pdf', width = 10, height = 10)
plot(Embeddings(pre.sct, "umap"), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(SlingshotDataSet(pre.sct.ss), lwd=2, col='black', type = "lineages", show.constraints = TRUE)
dev.off()


# to get pseudotime for each cell across lineages
slingPseudotime(post.sct.ss, na = TRUE)




# ....using SCE object -------------------
#sce <- slingshot(pre.sct.sc, clusterLabels = pre.sct.sc$seurat_clusters, reducedDim = 'umap', approx_points = 200, start.clus = 4) # step was run on test.R and sbatched using qsub-test.R script
load('data/slingobj.RData')
pre.sce <- sce

load('data/slingobj_postSCT.RData')
post.sce <- sce

summary(sce$slingPseudotime_1)
str(sce)

# ....to visualize inferred lineage for dta with points colored by pseudotime ---------------
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot.new()
pdf('figures/trajectoryAnalysis_slingshot_SCEobj_UMAP.pdf', width = 10, height = 10)
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black', type = 'lineages', show.constraints = TRUE)
dev.off()

# ....plotting pseudotime variable (one per lineage) for cells in all clusters ---------------
time <- as.data.frame(as.matrix(sce@colData$slingPseudotime_3))
cluster.info <- as.data.frame(as.matrix(sce@colData$seurat_clusters))
samples <- as.data.frame(as.matrix(sce@colData$seq_folder))

slingshot_df <- cbind(slingPseudotime_1 = time,
                      seurat_clusters = cluster.info,
                      seq_folder = samples)

names(slingshot_df) <- c('slingPseudotime_3', 'seurat_clusters', 'seq_folder')



plot.new()
pdf('figures/slingshot_orderedby_pseudotime3_postSCT.pdf', width = 10, height = 10)
ggplot(slingshot_df, aes(x = slingPseudotime_3, y = seurat_clusters, 
                         colour = seq_folder)) +
  geom_quasirandom(groupOnX = FALSE) + 
  theme_classic() +
  xlab("Third Slingshot pseudotime") + 
  ylab("cluster number") +
  ggtitle("Cells ordered by Slingshot pseudotime")
dev.off()






# ....identifying temporally dynamic genes --------------------------------
# finding genes that change their expression over the course of development.

#library(tradeSeq)

# lin <- getLineages(pre.sct.sc, clusterLabels = pre.sct.sc$seurat_clusters, start.clus = 4, reducedDim = "UMAP")
# crv <- getCurves(lin)
# a2 <- plotGeneCount(curve = crv, clusters = pre.sct.sc$seurat_clusters)
# ggsave(a2, filename = 'figures/trajectoryCurves_tradeSeq_preSCT_UMAP.pdf', width = 10, height = 10)
# 
# # determining number of knots
# icMat <- evaluateK(counts = counts(sce), k=3:10, nGenes = 200,
#                    pseudotime = slingPseudotime(crv, na = FALSE),
#                    cellWeights = slingCurveWeights(crv))
# 
# ggsave(icMat, filename = 'figures/PseudotimeAnalysis_determineKnots_pre.SCT_UMAP.pdf', width = 10, height = 10)
# 
# # 1. fit negative binomial GAM --------
# # since fitGAM runtimes are long, run it on most variable features from seurat object
# hvg <- VariableFeatures(pre.sct)
# length(rownames(sce)[rownames(sce) %in% hvg])
# 
# # to be fed into genes parameter
# idx <- which(rownames(sce) %in% hvg)
# 
# # following lines of code were run on a separate script on a scheduler (test2.R)
# # fitting GAM on top variable genes from seurat object
# # ordering rownames in sce as per variableFeatures in SCE
# sce <- sce[order(match(rownames(sce),hvg)),] 
# 
# sce <- fitGAM(sce, nknots = 6, genes = 1:10)
# 
# # 2. test for dynamic expression -----------
# ATres <- associationTest(sce, lineages = TRUE)

# steps above are run in separate scripts - test2.R, test3.R and test4.R


# ____1. Pre.SCT per lineage ------------------
load('data/associateTestResult_perLineages_k6_preSCT_UMAP.RData') # pre.sct, per lineage

# We want to inspect two lineages i.e. lineage 1 and 2
# pre.sct.ss@metadata$lineages
# pre.sct.ss is a "PseudotimeOrdering" object i.e. it has lineage information

lineage1.DE.genes <- rownames(ATres)[which(p.adjust(ATres$pvalue_1, "fdr") <= 0.05)]

length(unique(intersect(lineage1.DE.genes, ADRN.overlap)))
length(unique(intersect(lineage1.DE.genes, MES.overlap)))

# only parental cluster
lineage6.DE.genes <- rownames(ATres)[which(p.adjust(ATres$pvalue_6, "fdr") <= 0.05)]

length(unique(intersect(lineage6.DE.genes, ADRN.overlap)))
length(unique(intersect(lineage6.DE.genes, MES.overlap)))



# ____2. Pre.SCT global ------------------
load('data/associateTestResult_global_k6_preSCT_UMAP.RData') # pre.sct global

global.DE.genes <- rownames(ATres)[which(p.adjust(ATres$pvalue, "fdr") <= 0.05)]

length(unique(intersect(global.DE.genes, ADRN.overlap)))
length(unique(intersect(global.DE.genes, MES.overlap)))


# ____3. Post.SCT per lineage ------------------
load('data/associateTestResult_perLineages_k6_postSCT_UMAP.RData') # post.sct, per lineage

# We want to inspect two lineages i.e. lineage 1 and 2
# pre.sct.ss@metadata$lineages
# pre.sct.ss is a "PseudotimeOrdering" object i.e. it has lineage information

lineage1.DE.genes <- rownames(ATres)[which(p.adjust(ATres$pvalue_1, "fdr") <= 0.05)]

length(unique(intersect(lineage1.DE.genes, ADRN.overlap)))
length(unique(intersect(lineage1.DE.genes, MES.overlap)))

# only parental cluster
lineage6.DE.genes <- rownames(ATres)[which(p.adjust(ATres$pvalue_6, "fdr") <= 0.05)]

length(unique(intersect(lineage6.DE.genes, ADRN.overlap)))
length(unique(intersect(lineage6.DE.genes, MES.overlap)))


# ____4. Post.SCT global ------------------
load('data/associateTestResult_global_k6_postSCT_UMAP.RData') # post.sct global

global.DE.genes <- rownames(ATres)[which(p.adjust(ATres$pvalue, "fdr") <= 0.05)]

length(unique(intersect(global.DE.genes, ADRN.overlap)))
length(unique(intersect(global.DE.genes, MES.overlap)))


# 3. plotting a heatmap to visualize expression of DE genes across lineages ------------
#topgenes <- rownames(ATres[order(ATres$pvalue_1), ])[1:200]
topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:200]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$logcounts[topgenes, pst.ord]


#write.table(topgenes, file = 'data/topgenes_assoTest_preSCT.txt', col.names = F, row.names = F, quote = F)


col_anno <- colData(sce)[,c(1,12)]
col_anno.ordered <- col_anno[colnames(heatdata),]
any(rownames(col_anno.ordered) == colnames(heatdata))

col_anno.ordered <- as.data.frame(col_anno.ordered)

# setting annotation colors
n <- 18
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, n)

cluster_cols <- setNames(col, unique(col_anno.ordered$seurat_clusters))

anno_colors <- list(
  orig.ident = c('parent2' = '#6D9886', 'parent1' = 'royalblue', 'LRX1' = '#FF7777', 'LRX2' = '#90AACB'))

# plot heatmap
p <- pheatmap(heatdata,
              annotation_col = col_anno.ordered,
              cluster_cols = F,
              show_colnames = F,
              fontsize_row = 3,
              annotation_colors = anno_colors)
#ggsave(p, filename = paste0('figures/',Sys.Date(),'_heatmap_pseudotime1_200genes_preSCT.pdf'), width = 10, height = 10)
ggsave(p, filename = paste0('figures/',Sys.Date(),'_heatmap_global_pseudotime_200genes_preSCT.pdf'), width = 10, height = 10)


rownames(heatdata[p$tree_row[["order"]],])
gene.clusters <- c("IQCJ-SCHIP1","FGF13","KCNIP4","PCDH7","AL163541.1","XYLT1","FMN1","CREB5","PTPRG","PTPRM","PTCHD1-AS","PRKN","ROBO2","MT-ND3","CCSER1","CLSTN2","AC073050.1","SLC8A1","SGCZ","ROBO1","MAST4","CHGB","ZNF804A","AC069277.1")

intersect(rownames(heatdata[p$tree_row[["order"]],]),mesench$V1)
intersect(rownames(heatdata[p$tree_row[["order"]],]),adrenal$V1)


wnt.signaling <- c("APC","RHOA","ATF4","AXIN1","BCL9","BIRC5","BTRC","CBY1","GCFC2","CAMK2A","CAMK2B","CAMK2D","CAMK2G","CCND1","CDH1","CDH10","CDH11","CDH12","CDH13","CDH15","CDH16","CDH17","CDH18","CDH19","CDH2","CDH20","CDH22","CDH23","CDH24","CDH26","CDH3","CDH4","CDH5","CDH6","CDH7","CDH8","CDH9","CEBPB","CER1","CHUK","CLDN1","CREB1","CREB3","CREB3L4","CREBBP","CSNK1A1","CSNK1D","CSNK1E","CSNK1G1","CSNK1G2","CSNK1G3","CTNNA1","CTNNB1","DAAM1","DAAM2","DKK1","DKK2","DKK4","DVL1","ENC1","FOSL1","FRZB","FZD1","FZD10","FZD2","FZD3","FZD4","FZD5","FZD6","FZD7","FZD8","FZD9","GAST","GNA11","GNA12","GNA13","GNA14","GNA15","GNAI1","GNAI2","GNAI3","GNAL","GNAO1","GNAS","GNAT1","GNAT2","GNAZ","GNB1","GNB2","GNB3","GNB4","GNB5","GNG10","GNG11","GNG12","GNG13","GNG2","GNG3","GNG4","GNG5","GNG7","GNG8","GNGT1","GNGT2","GSK3B","HNF4A","JUN","KREMEN1","KREMEN2","LEF1","LRP5","LRP6","MAP2K4","MAP2K7","MAP3K1","MAP3K10","MAP3K11","MAP3K12","MAP3K13","MAP3K14","MAP3K2","MAP3K3","MAP3K4","MAP3K5","MAP3K6","MAP3K7","MAP3K8","MAP3K9","MAPK10","MAPK8","MAPK9","MMP7","MYC","NFAT5","NFATC1","NFATC2","NFATC3","NFATC4","NFE2L1","NRCAM","PSMD6","DCHS1","PLAUR","PLCB1","PLCB2","PLCB3","PLCB4","PLCD1","PLCD3","PLCD4","PLCE1","PLCG1","PLCG2","PLCZ1","PPARD","PPP3CA","PPP3CB","PPP3CC","PPP3R1","PPP3R2","PRKCA","PRKCB","PRKCD","PRKCE","PRKCG","PRKCH","PRKCI","PRKD1","PRKD3","PRKCQ","PRKCZ","PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7","PSMB1","PSMB10","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7","PSMB8","PSMB9","PSMC1","PSMC2","PSMC3","PSMC4","PSMC5","PSMC6","PSMD1","PSMD2","PSMD3","PSMD4","PSMD5","PSMD7","PTGS2","PYGO1","PYGO2","RAC1","REST","ROCK2","SMARCA4","HNF1A","TCF12","TCF15","TCF19","HNF1B","TCF20","TCF21","TCF23","TCF3","TCF4","TCF7","TCF7L2","ZEB1","TEAD1","TEAD4","TFAM","TP53","UBB","UBC","UBD","UBE2D2","VIM","WIF1","WNT1","WNT10A","WNT10B","WNT11","WNT9A","WNT9B","WNT16","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","ZNF354A")
wnt <- intersect(topgenes, wnt.signaling)

# canon.wnt.signaling <- c("CDH15","CDH13","CDH12","CDH11","CDH1","DKK4","DKK3","DKK2","DKK1","LEF1","TCF4","TCF25","TCF21","CDH16","CDH17","CDH18","CDH19","CDH2","CDH24","CDH26","CDH3","CDH4","CDH5","CDH6","CDH8","CTNNB1","CSNK1A1","GSK3B","ANAPC10","ANAPC11","ANAPC13","ANAPC15","ANAPC16","ANAPC2","ANAPC5","WNT1","WNT10A","WNT11","WNT16","WNT2","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT9A")
# unique(intersect(canon.wnt.signaling, global.topgenes))
# 
# 
# wnt.hedgehog.notch <- c("AFP","ACVR1","ALPL","RUNX3","RIPK4","APC","ARRB1","ARRB2","AXIN1","AXIN2","BAG6","BCL9","BCL9L","TGFBI","BMP4","BMP6","BMP7","BMPR2","BRD2","CACYBP","PROM1","CDC73","CDCP1","CDX2","CELSR2","CSNK1A1","CSNK1D","CSNK1E","CCN2","CTNNB1","CTNND1","CUX1","CYLD","DDX5","DIXDC1","DKK1","DKK2","DLK1","DLL1","DLL3","DLL4","DVL2","DVL3","EMD","ENG","FSCN1","FOXD3","FOXK2","FOXP1","FOXP2","FOXP3","FUT4","FZD5","FZD6","GDF15","GLI1","GLI2","GREM1","GSK3A","GSK3B","HES1","ID2","ID3","IGF2BP1","ITCH","JAG1","JAG2","WWC1","KIF3A","KLF4","LATS1","LEF1","LEFTY1","LFNG","LIN28A","LIN28B","LRP5","LRP6","MAML1","MAML2","MARK3","MESD","MIB1","AMHR2","MOB1A","MOB1B","STK4","MSX1","MSI1","MYB","NANOG","NDRG1","NDRG4","NEDD4L","CSPG4","NCSTN","NKD1","NKD2","NKX2-2","NKX2-5","NLK","NOTCH1","NOTCH2","NOTCH3","NOTCH4","NR6A1","NUMB","POU5F1","PARN","PAX2","PAX3","PAX5","PAX9","PPM1A","PSEN1","PSEN2","PSENEN","PTCH1","PTCH2","PTPN14","RBPJ","RING1","RUVBL1","ZFYVE9","SAV1","SCRIB","SERPINE1","SFRP1","SHH","SIX1","SKIL","SMAD1","SMAD2","SMAD3","SMAD4","SMAD5","SMAD9","SMURF1","SMURF2","SNAI1","SNAI2","SOX2","SP1","SPARC","SUFU","TAB1","TAB2","ADAM17","MAP3K7","WWTR1","TBK1","TCF12","TCF7L1","TCF7","TCF7L2","ZEB1","TDGF1","TEAD1","TEAD2","TEAD3","TEAD4","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR3","TRIM33","TLE1","TLE2","TLE3","TLE4","TRIB2","USP15","UTF1","VANGL1","WBP2","WIF1","WNT3A","WNT5A","WNT5B","AMER1","YAP1")
# unique(intersect(wnt.hedgehog.notch, lineage.topgenes))
# 
# 
# non.canon.wnt <- c("DAG1","TK1","PFN1","PFN2","PFN4","PXN","ROCK2","CDC42","PPP3CA","NR2C2","NFATC1","NFATC3","TCF21","TCF25","TCF4","LEF1","WNT1","WNT10A","WNT11","WNT16","WNT2","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT9A","FZD1","FZD10","FZD2","FZD4","FZD5","FZD9","MAP2K3","MAP2K6","PRKCE","PRKCA","PRKCB","PRKCD","PRKCG","PRKCH","PRKCI","PRKCZ","CSNK1A1","RHOA","RAC1","CAMK2A","CAMK2B","CAMK2G","MAPK14","MAPK1","MAPK6","MAPK7","CTNNB1","NLK","ATF2")
# unique(intersect(non.canon.wnt, lineage.topgenes))

#hippo <- c("PATJ","TPTEP2-CSNK1E","APC2","YAP1","YWHAQ","RASSF1","FZD10","FRMD6","WTIP","CSNK1D","CSNK1E","CCN2","CTNNA1","CTNNA2","CTNNB1","GDF7","AMOT","RASSF6","DLG1","AFP","DLG2","DLG3","DLG4","DVL1","DVL2","DVL3","FGF1","WWC1","FBXW11","CRB1","SCRIB","FZD2","WWTR1","LATS2","AMH","BBC3","GLI2","CRB2","CTNNA3","GSK3B","APC","BIRC2","BIRC3","BIRC5","ID1","ID2","BMP8A","ITGB2","AREG","GDF6","LLGL2","LLGL1","SMAD1","SMAD2","SMAD3","SMAD4","SMAD7","MYC","NF2","SERPINE1","PARD6A","LEF1","WNT16","WNT4","PPP1CA","PPP1CB","PPP1CC","PPP2CA","PPP2CB","PPP2R1A","PPP2R1B","PPP2R2A","PPP2R2B","PPP2R2C","MOB1A","PRKCI","PPP2R2D","PRKCZ","PARD3","CCND1","ACTB","SAV1","PALS1","BMP2","BMP4","BMP5","BMP6","BMP7","BMP8B","BMPR1A","BMPR1B","BMPR2","SNAI2","SOX2","STK3","TCF7","TCF7L2","TEAD1","TEAD4","TEAD3","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","ACTG1","TP53BP2","TP73","WNT1","WNT2","WNT3","WNT5A","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT10B","WNT11","WNT2B","WNT9A","WNT9B","YWHAB","YWHAE","YWHAG","YWHAH","YWHAZ","FZD5","FZD3","FRMD1","WNT10A","WNT5B","GDF5","AXIN1","AXIN2","FZD1","FZD4","FZD6","FZD7","FZD8","FZD9","TCF7L1","PARD6G","PARD6B","TEAD2","AJUBA","NKD1","NKD2","CCND2","BTRC","CCND3","WNT3A","LIMD1","LATS1","DLG5","MOB1B","CDH1")
# unique(intersect(hippo, global.topgenes))

robo.receptor.signaling <- c("ROBO1","SLIT2","GPC1","RAC1","ABL1","ABL2","CAP2","CAP1","CLASP2","CLASP1","ROBO2","ROBO3","CDC42","NCK1","PAK1","PAK2","PAK3","SOS1","SOS2","NCK2","SRGAP3","PAK4","SRGAP2","PAK6","PAK5","SRGAP1","ARHGAP39","BUB1B-PAK6","PFN1","PFN2","VASP","EVL","ENAH")
robo <- intersect(topgenes, robo.receptor.signaling)


signaling.by.slit <- c("ABL1","ABL2","ACTR2","ACTR3","ARHGEF1","ARHGEF10","ARHGEF11","ARHGEF12","ARHGEF15","ARHGEF16","ARHGEF2","ARHGEF3","ARHGEF4","ARHGEF5","ARHGEF6","ARHGEF7","ARHGEF9","CDC42","CLASP1","CLASP2","CXCL12","CXCR4","DCC","ENAH","GNA11","GNA12","GNA13","GNA14","GNA15","GNAS","GNB1","GNB2","GNB3","GNB4","GNB5","GNG2","GNG3","GNG4","MCF2L","NET1","ROBO1","SLIT1","SLIT2","SLIT3","SRGAP1","SRGAP2","SRGAP3","WASL")
slit <- intersect(topgenes, signaling.by.slit)


EMT <- c("ARHGEF18","CGN","F11R","FKBP1A","MIR6869","PARD3","PARD6A","PRKCZ","RHOA","RPS27A","SMURF1","TGFB1","TGFBR1","TGFBR2","UBA52","UBB","UBC")
intersect(topgenes, EMT)


# 4. extracting genes of interest ----------
# reshape data
# converting sparse matrix to data frame
summ <- summary(heatdata)

df.heatdata <- data.frame(gene = rownames(heatdata)[summ$i],
           cells = colnames(heatdata)[summ$j],
           logcounts = summ$x)

# merge df.heatdata with col_anno.ordered
df.heatdata <- merge(df.heatdata, col_anno.ordered, by.x = 'cells', by.y = 'row.names')

# comparing mean expression of ADRN/MES genes between Parental and LRX cells
unique(intersect(topgenes, ADRN.overlap))
unique(intersect(topgenes, MES.overlap))

#genes.of.interest <- rbind(adrenal, mesench)
genes.of.interest <- unique(c(wnt, robo, slit, 'SOX6', 'MCC', 'BCL2'))
genes.of.interest <- as.data.frame(genes.of.interest)
names(genes.of.interest)[1] <- 'V1'

df.heatdata <- merge(df.heatdata, genes.of.interest, by.x = 'gene', by.y = 'V1')


 
# plotting expression 
df.heatdata$seurat_clusters <- factor(df.heatdata$seurat_clusters, levels=c(13,5,6,3,10,8,0,18,15,9,19,7,17,16,11,4,2,14,1,12))
a1 <- ggplot(df.heatdata, aes(seurat_clusters, logcounts, fill = group)) +
  geom_boxplot() +
  facet_wrap(. ~ gene)
#ggsave(a1, filename = 'figures/boxplot_pseudotime1_top200genes_preSCT.pdf', width = 15, height = 10)
#ggsave(a1, filename = 'figures/boxplot_global_pseudotime_top200genes_postSCT.pdf', width = 15, height = 10)


df.heatdata$group <- ifelse(grepl('LRX', df.heatdata$orig.ident), 'LRX', 'Parental')
a2 <- ggplot(df.heatdata, aes(seurat_clusters, logcounts)) +
  geom_boxplot() +
  facet_wrap(. ~ gene)
#ggsave(a2, filename = 'figures/boxplot_global_GOI_pseudotime_top200genes_preSCT.pdf', width = 15, height = 10)
#ggsave(a2, filename = 'figures/boxplot_GOI_pseudotime1_top200genes_preSCT.pdf', width = 15, height = 10)

ggsave(a2, filename = paste0('figures/',Sys.Date(),'_boxplot_test_pseudotime1_top200genes_preSCT.pdf'), width = 15, height = 10)



p1 <- ggplot(df.heatdata, aes(reorder(gene, logcounts), logcounts, fill = orig.ident)) +
  geom_boxplot() +
  facet_wrap(. ~ group) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(p1, filename = 'figures/boxplot_pseudotime2_top100genes_preSCT.pdf', width = 15, height = 10)


p2 <- ggplot(df.heatdata, aes(logcounts, fill = orig.ident)) +
  geom_density(alpha = 0.4) +
  facet_wrap(. ~ group) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(p2, filename = 'figures/densityplot_pseudotime1_top50genes_preSCT.pdf', width = 15, height = 10)



# 5. GSEA on DE genes for lineages -------------------------------------------
# top 100,200 and 500 genes give very generalized hits

# hs <- org.Hs.eg.db
# 
# my.symbols <- c(as.character(topgenes))
# df <- AnnotationDbi::select(hs, 
#                             keys = my.symbols,
#                             columns = c("ENTREZID", "SYMBOL"),
#                             keytype = "SYMBOL")
# 
# 
# # getting entrezsymbol for top genes
# top.genes.df <- as.data.frame(topgenes)
# top.genes.gsea <- merge(top.genes.df, df, by.x = 'topgenes', by.y = 'SYMBOL')
# top.genes.gsea <- na.omit(top.genes.gsea)
# 
# # setting ranks using wald stat for lineage 1
# ranks <- ATres[rownames(ATres) %in% top.genes.gsea$topgenes, 'waldStat_1']
# names(ranks) <- top.genes.gsea$ENTREZID
# head(ranks)
# 
# # load pathways
# my_pathways <- reactomePathways(names(ranks))
# 
# eaRes <- fgsea(pathways = my_pathways, stats = ranks, nperm = 100000, minSize = 15)
# out <- eaRes




#' trying monocle! -------------

# # creating feature data
# gene_annotation <- as.data.frame(rownames(pre.sct@reductions[["pca"]]@feature.loadings),
#                                  row.names = rownames(pre.sct@reductions[["pca"]]@feature.loadings))
# colnames(gene_annotation) <- "gene_short_name"
# 
# # creating pheno data
# cell_metadata <- as.data.frame(pre.sct@assays[["RNA"]]@counts@Dimnames[[2]],
#                                row.names = pre.sct@assays[["RNA"]]@counts@Dimnames[[2]])
# colnames(cell_metadata) <- "barcode"
# 
# # creating a CellDataSet for Monocle
# New_matrix <- pre.sct@assays[["RNA"]]@counts
# New_matrix <- New_matrix[rownames(pre.sct@reductions[["pca"]]@feature.loadings), ]
# expression_matrix <- New_matrix
# 
# 
# cds_from_seurat <- newCellDataSet(expression_matrix,
#                                   phenoData = new("AnnotatedDataFrame", data = cell_metadata),
#                                   featureData = new("AnnotatedDataFrame", data = gene_annotation),
#                                   expressionFamily=negbinomial.size())
# 
# # exploring slot names
# slotNames(cds_from_seurat)
# 
# 
# 
# # Estimating size factors and dispersions (required)
# cds_from_seurat <- estimateSizeFactors(cds_from_seurat)
# cds_from_seurat <- estimateDispersions(cds_from_seurat)
# 
# 
# # tally number of cells expressing a gene and # of genes expressed among all cells (recommended step for all monocle analysis)
# cds_from_seurat <- detectGenes(cds_from_seurat, min_expr = 0.1)
# print(head(fData(cds_from_seurat), n = 20))
# 
# summary(fData(cds_from_seurat)$num_cells_expressed)
# 
# expressed_genes <- row.names(subset(fData(cds_from_seurat),
#                                     num_cells_expressed >= 10))
# 
# 
# 
# # to find number of genes expressed per cell
# head(pData(cds_from_seurat))
# 
# 
# 
# summary(pData(cds_from_seurat)$num_genes_expressed)
# 
# 
# # standardise to Z-distribution
# x <- pData(cds_from_seurat)$num_genes_expressed
# x_1 <- (x - mean(x)) / sd(x)
# summary(x_1)
# 
# library(ggplot2)
# # I like the default theme of cowplot
# library(cowplot)
# 
# df <- data.frame(x = x_1)
# g <- ggplot(df, aes(x)) +
#   geom_histogram(bins = 50) +
#   geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
# ggsave(g, filename = 'figures/QC_monocle_num_of_genes.pdf', width = 10, height = 10)
# 
# 
# # clustering cells without marker genes
# disp_table <- dispersionTable(cds_from_seurat)
# head(disp_table)
# 
# # select genes with mean expression >= 0.1 to use in clustering step
# table(disp_table$mean_expression>=0.1)
# 
# unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
# 
# cds_from_seurat <- setOrderingFilter(cds_from_seurat, unsup_clustering_genes$gene_id)
# p <- plot_ordering_genes(cds_from_seurat)
# # The plot_ordering_genes() function plots mean expression against the empirical dispersion and highlights the set of genes (as black dots) that will be used for clustering.
# ggsave(p, filename = 'figures/QC_monocle_custering.pdf', width = 10, height = 10)
# 
# 
# # Constructing Single Cell Trajectories
# # The trajectory analysis consists of three stages:
# # 1. Choose genes that define progress
# # 2. Reduce the dimensionality of the data
# # 3. Order cells in pseudotime
# 
# 
# # just using expressed_genes (genes_expressed >= 10 cells)
# cds_from_seurat_subset <- cds_from_seurat[expressed_genes,]
# cds_from_seurat_subset
# 
# # ordering - based on genes that differ between clusters
# cds_from_seurat_subset <- detectGenes(cds_from_seurat, min_expr = 0.1)
# fData(cds_from_seurat_subset)$use_for_ordering <- fData(cds_from_seurat_subset)$num_cells_expressed > 0.05 * ncol(cds_from_seurat_subset)
# head(fData(cds_from_seurat_subset))
# 
# # how many genes are used?
# table(fData(cds_from_seurat_subset)$use_for_ordering)
# 
# # FALSE  TRUE 
# # 973  1027 
# 
# 
# p1 <- plot_pc_variance_explained(cds_from_seurat_subset, return_all = FALSE)
# ggsave(p1, filename = 'figures/trajectoryAnalysis_monocle_PCA_preSCT.pdf', width = 10, height = 10)
# 
# 
# # reduce Dimensions
# cds_from_seurat_subset <- reduceDimension(cds_from_seurat_subset,
#                                  max_components = 2,
#                                  norm_method = 'log',
#                                  num_dim = 20,
#                                  reduction_method = 'tSNE',
#                                  verbose = TRUE)
# 
# cds_from_seurat_subset <- clusterCells(cds_from_seurat_subset, verbose = FALSE)
# 
# p2 <- plot_rho_delta(cds_from_seurat_subset, rho_threshold = 2, delta_threshold = 10)
# ggsave(p2, filename = 'figures/trajectoryAnalysis_monocle_clustering_rhoDelta_preSCT.pdf', width = 10, height = 10)
# 
# # use rho = 2 and delta = 10 to cluster the cells again.
# cds_from_seurat_subset <- clusterCells(cds_from_seurat_subset,
#                               rho_threshold = 2,
#                               delta_threshold = 10,
#                               skip_rho_sigma = T,
#                               verbose = FALSE)
# 
# table(pData(cds_from_seurat_subset)$Cluster)
# 
# p3 <- plot_cell_clusters(cds_from_seurat_subset)
# ggsave(p3, filename = 'figures/trajectoryAnalysis_monocle_clustering_preSCT.pdf', width = 10, height = 10)
# 
# 
# # performing the differential gene expression analysis across all cell clusters.
# # to use the top 1,000 most significantly differentially expressed genes as the set of ordering genes and perform the dimension reduction and the trajectory analysis (using the orderCells() function).
# 
# clustering_DEG_genes <- differentialGeneTest(cds_from_seurat_subset,
#                                              fullModelFormulaStr = '~Cluster',
#                                              cores = 8)
# 
# dim(clustering_DEG_genes)
# 
# clustering_DEG_genes %>% arrange(qval) %>% head()
# 
# 
# my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
# cds_from_seurat_subset <- setOrderingFilter(cds_from_seurat_subset, ordering_genes = my_ordering_genes)
# cds_from_seurat_subset <- reduceDimension(cds_from_seurat_subset, method = 'DDRTree')
# cds_from_seurat_subset <- orderCells(cds_from_seurat_subset)
# 
# p4 <- plot_cell_trajectory(cds_from_seurat_subset, color_by = "Cluster")
# ggsave(p4, filename = 'figures/trajectoryAnalysis_monocle_trajectory_preSCT.pdf', width = 10, height = 10)
# 
# 
# # Finding Genes that Change as a Function of Pseudotime
# # Once we have a trajectory, we can use differentialGeneTest() to find genes 
# # that have an expression pattern that varies according to pseudotime.
# 
# # pseudotime is now a column in the phenotypic data as well as the cell state
# head(pData(cds_from_seurat_subset))



