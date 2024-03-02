library(Seurat)
library(escape)


WORKDIR <- 'home/rstudio/levenberg/'
PLOTSDIR <- '~/levenberg/plots/'
use.name <- 'Full'

fibroblasts <- readRDS(paste0(WORKDIR, use.name, '.new.fibroblasts.rds'))

gene.sets <- getGeneSets(library = "C2", subcategory = 'CP:WIKIPATHWAYS')

ES <- enrichIt(obj = fibroblasts[["RNA"]]$counts,
               gene.sets = gene.sets,
               groups = 1500, cores = 40,
               min.size = 5, #method = 'UCell',
               ssGSEA.norm = T)
ES$cluster <- fibroblasts$seurat_clusters

saveRDS(ES, paste0(WORKDIR, use.name, 'wiki.fibroblasts.rds'))