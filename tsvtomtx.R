library(tidyverse)
library(DropletUtils)

umi <- read.delim('~/matrix.tsv', header = T)
genes <- read.delim('~/genes.tsv', header = T)

genes.subset <- data.frame(genes = genes$gene)
umi <- cbind(umi, genes.subset)
umi <- umi %>%
  column_to_rownames(var = 'genes')


# convert tsv to sparse matrix 
umi.sparse <- Matrix(as.matrix(umi), sparse = TRUE)

write10xCounts(
  path = '~/out',
  x = umi.sparse,
  overwrite = FALSE)
