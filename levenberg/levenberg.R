library(Seurat)
library(dittoSeq)
library(SingleR)
library(celldex)
library(dplyr)
library(escape)
library(plyr)
library(cowplot)

sc1 = Read10X_h5('~/Documents/levenberg/filtered_feature_bc_matrix_SL111_old.h5')
sc2 = Read10X_h5('~/Documents/levenberg/filtered_feature_bc_matrix_SL112_old.h5')
sc3 = Read10X_h5('~/Documents/levenberg/filtered_feature_bc_matrix_SL161.h5')
sc4 = Read10X_h5('~/Documents/levenberg/filtered_feature_bc_matrix_SL162.h5')
sc5 = Read10X_h5('~/Documents/levenberg/filtered_feature_bc_matrix_SL111.h5')
sc6 = Read10X_h5('~/Documents/levenberg/filtered_feature_bc_matrix_SL112.h5')

N = 10000
raw = cbind(sc1[,sample(ncol(sc1),N)],sc2[,sample(ncol(sc2),N)],sc3[,sample(ncol(sc3),N)],sc4[,sample(ncol(sc4),N)],sc5[,sample(ncol(sc5),N)],sc6[,sample(ncol(sc6),N)])

colnames(raw) = make.unique(colnames(raw))
#source = c(rep('Pl.My.Pr.Hu',ncol(sc1)),rep('My.Pr.Hu',ncol(sc2)),rep('My.Hu',ncol(sc3)),rep('Pl.My.Hu',ncol(sc4)))
source = c(rep('Pl.My.Pr.Hu.1',N),rep('My.Pr.Hu.1',N),rep('My.Hu',N),rep('Pl.My.Hu',N),rep('Pl.My.Pr.Hu.2',N),rep('My.Pr.Hu.2',N))

createBasicSC = function(data,source) {
  sc <- CreateSeuratObject(data, project = "muscle",min.cells = 10,min.features = 200)
  names(source) = colnames(sc)
  sc$source = source[colnames(sc)]
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  #VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  sc <- SCTransform(sc, vars.to.regress = "percent.mt", verbose = FALSE)
  sc <- RunPCA(sc)
  #ElbowPlot(sc)
  sc <- RunUMAP(sc, reduction = "pca", dims = 1:30)
  sc <- FindNeighbors(sc, dims = 1:30)
  sc <- FindClusters(sc, resolution = 0.5)
  sc
}
#sc$source = plyr::mapvalues(sc$source, from=c('S1','S2'), to = c('WithPlacenta','NoPlacenta'))

sc.all = createBasicSC(raw,source)
sc.only3 = createBasicSC(sc3[,sample(ncol(sc3),N)],rep('My.Hu',N))

sc = sc.all
use.name = 'All2'
#dittoDimPlot(sc, reduction.use = "umap", var = "nFeature_RNA",size = 0.5)+dittoDimPlot(sc, reduction.use = "umap", var = "nCount_RNA",size = 0.5)+dittoDimPlot(sc, reduction.use = "umap", var = "percent.mt",size = 0.5)

p0 = dittoDimPlot(sc, reduction.use = "umap", "source",size = 0.5)+ggtitle('By sample')
ggsave(paste0('~/Documents/levenberg/',use.name,'.plot0.png'),width=8,height=7)


p1 = dittoDimPlot(sc, reduction.use = "umap", "source",size = 0.5,split.by = 'source')+ggtitle('By sample')
p1
ggsave(paste0('~/Documents/levenberg/',use.name,'.plot1.png'),width=8,height=7)

new.cluster.ids <- c('Pericytes',
                     'Myoblasts',
                     'Pericytes',
                     'Myoblasts',
                     'Myoblasts',
                     'Placenta',
                     'Myoblasts',
                     'Myoblasts',
                     'Per_Myo',
                     'Pericytes',
                     'Myoblasts',
                     'Myoblasts',
                     'Per_Myo',
                     'Per_Myo',
                     'HUVEC')

#new.cluster.ids <- c('Pericyte','Myoblast','Myoblast','Pericyte','Myoblast','Myoblast','Placenta','Pericyte','Myoblast','Pericyte','Myoblast','Myoblast','Myoblast','Myoblast','HUVEC')
new.cluster.ids2 = new.cluster.ids
new.cluster.ids = make.unique(new.cluster.ids)
names(new.cluster.ids) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids)
sc$idents= sc@active.ident


p2 = dittoDimPlot(sc, reduction.use = "umap", "idents",size = 0.5,split.by = 'source')+ggtitle('By clusters')
p2
ggsave(paste0('~/Documents/levenberg/',use.name,'.plot2.png'),width=7,height=7)

p2.1 = dittoDimPlot(sc, reduction.use = "umap", "idents",size = 0.5,split.by = 'idents')+ggtitle('By clusters')
p2.1
ggsave(paste0('~/Documents/levenberg/',use.name,'.plot2.1.png'),width=7,height=7)


p2.0 = dittoDimPlot(sc, reduction.use = "umap", "idents",size = 0.5)+ggtitle('By clusters')
ggsave(paste0('~/Documents/levenberg/',use.name,'.plot2.0.png'),width=7,height=7)

names(new.cluster.ids2) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids2)
sc$idents2 = sc@active.ident

p2.3 = dittoDimPlot(sc, reduction.use = "umap", "idents2",size = 0.5)+ggtitle('By clusters')
ggsave(paste0('~/Documents/levenberg/',use.name,'.plot2.3.png'),width=7,height=7)


p3 = dittoBarPlot(sc, "source", group.by = "idents")+ggtitle('')
p3
ggsave(paste0('~/Documents/levenberg/',use.name,'.bars1.png'),width=7,height=7)

p3.1 = dittoBarPlot(sc, "source", group.by = "idents",scale='count')+ggtitle('')
p3.1
ggsave(paste0('~/Documents/levenberg/',use.name,'.bars2.png'),width=7,height=7)

p3.2 = dittoBarPlot(sc, "idents2", group.by = "source",scale='count')+ggtitle('')
p3.2
ggsave(paste0('~/Documents/levenberg/',use.name,'.bars3.png'),width=7,height=7)

p3.3 = dittoBarPlot(sc, "source", group.by = "idents",scale='count')+ggtitle('')
p3.3
ggsave(paste0('~/Documents/levenberg/',use.name,'.bars4.png'),width=7,height=7)


tbl = table(sc$source,sc$idents)
corrplot(t(t(tbl)/colSums(tbl)),is.corr=F,method = 'square', addCoef.col = 'black')


tbl = table(sc$source,sc$idents)
corrplot(t(t(tbl)/rowSums(tbl)),is.corr=F,method = 'square', addCoef.col = 'black')
#sc <- FindClusters(sc, resolution = 80)

cowplot::plot_grid(plot_grid(p1,p2,labels = c('A', 'B')),plot_grid(p3,labels = c('C')),nrow=2)
ggsave('~/Documents/levenberg/plot1.png',width=12,height=10)
sc = readRDS('~/Documents/levenberg/levengberg_new.rds')
singler_ref <- celldex::BlueprintEncodeData()

sc <- FindClusters(sc, resolution = 80)
singler_results <- SingleR(test = GetAssayData(sc, assay = 'RNA', slot = 'data'),
                           ref = singler_ref,
                           labels = singler_ref@colData@listData$label.main,
                           clusters = sc@meta.data$SCT_snn_res.80
)


x <- data.frame(row.names = rownames(singler_results), labels=singler_results$labels)
sc <- AddMetaData(sc, metadata=x[sc$SCT_snn_res.80, 1], col.name = "singler")
saveRDS(sc,'~/Documents/levenberg/levengberg_new.rds')

dittoDimPlot(sc, reduction.use = "umap", "singler", main = "By cell type")

sc@active.ident= sc$idents
markers.cluster <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25,min.diff.pct = 0.1, logfc.threshold = 0.25)
top5 <- markers.cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p4 = dittoDotPlot(sc, vars = unique(top5$gene), group.by = "idents",max.color='blue',min.color='white')+ylab('Clusters')
p4
ggsave(paste0('~/Documents/levenberg/',use.name,'.plo4.png'),width=12,height=6)

plot_grid(plot_grid(p1,p2,labels = c('A', 'B')),plot_grid(p3,labels = c('C')),plot_grid(p4,labels = c('D')),nrow=3)
ggsave(paste0('~/Documents/levenberg/',use.name,'.figure1.png'),width=12,height=14)


pericytes = subset(sc,idents = levels(sc)[grepl('Per',levels(sc))])
myocytes = subset(sc,idents = levels(sc)[grepl('Myo',levels(sc))])

markers.cluster <- FindAllMarkers(pericytes, only.pos = TRUE, min.pct = 0.25,min.diff.pct = 0.1, logfc.threshold = 0.25)
top5 <- markers.cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p4.per = dittoDotPlot(pericytes, vars = unique(top5$gene), group.by = "idents",max.color='blue',min.color='white')+ylab('Clusters')
p4.per

markers.cluster <- FindAllMarkers(myocytes, only.pos = TRUE, min.pct = 0.25,min.diff.pct = 0.1, logfc.threshold = 0.25)
top5 <- markers.cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p4.myo = dittoDotPlot(myocytes, vars = unique(top5$gene), group.by = "idents",max.color='blue',min.color='white')+ylab('Clusters')
p4.myo
ggsave(paste0('~/Documents/levenberg/',use.name,'.markers.group.png'),width=12,height=8)


plot_grid(p4.per,p4.myo,labels = c('A', 'B'),nrow=2)


genes = c('KRT7','ALPP','HLA-G','ACTA2','VIM','FAP','COL1A1','COL1A2','MCAM','CSPG4','PDGFRB','PECAM1','CDH5','VWF')
dittoDotPlot(sc, vars = genes[genes %in% rownames(sc)], group.by = "idents",max.color='blue',min.color='white')+ylab('Clusters')


genes = c('MYOD1','MYF5','TNFRSF12A','TTN','MUSK','FLNC','FGFR4','STK40','CXCL1','PAX7','ACTA2','PDGFRB','COLEC11','','MCAM','HIGD1B','ZIC1','ANGPT1','TBX18','CSPG4','DES','PTH1R')
dittoDotPlot(sc, vars = genes[genes %in% rownames(sc)], group.by = "idents",max.color='blue',min.color='white')+ylab('Clusters')

genes = c('DES',rownames(sc)[grepl('MYH',rownames(sc))],rownames(sc)[grepl('MYOG',rownames(sc))])
dittoDotPlot(sc, vars = genes[genes %in% rownames(sc)], group.by = "idents",max.color='blue',min.color='white')+ylab('Clusters')


genes = c('KLF5','PHX2B','SOX11','WNT3A',"PTGDS","MYPN" )
dittoDotPlot(sc, vars = genes[genes %in% rownames(sc)], group.by = "SCT_snn_res.0.5",max.color='blue',min.color='white')+ylab('Clusters')


library(gridExtra)
library(grid)
library(ggrepel)

p = list()
markers=list()

createCorrelationPlot = function(sc,celltype) {
  sc2 = subset(sc,idents2==celltype)
  sc2$idents = factor(sc2$idents)
  
  avg = data.frame(S1 = log1p(rowMeans(sc2@assays$RNA@data[,sc2$source=='My.Pr.Hu'])),
                   S2 = log1p(rowMeans(sc2@assays$RNA@data[,sc2$source=='Pl.My.Pr.Hu'])),
                   gene = rownames(sc2@assays$RNA@data),row.names=rownames(sc2@assays$RNA@data))
  
  markers <- FindMarkers(sc2,group.by = 'source', ident.1 = "My.Pr.Hu", ident.2 = "Pl.My.Pr.Hu", verbose = FALSE,min.diff.pct=0.1)
  avg$label = avg$gene
  avg$label[!(avg$label %in% rownames(markers)[markers$p_val_adj< 0.0001 & abs(markers$avg_log2FC)>0.3])] = ''
  p = ggplot(avg, aes(S1,S2,label=label)) + geom_point()+theme_classic()+xlab('No Placenta')+ylab('With Placenta')+ggtitle(celltype)+
    geom_text_repel(label=avg[,'label'],box.padding = 0.5, max.overlaps = Inf)+geom_point(data = avg[avg$label != "",], color = "red")
  list(p,markers)
}

p1 = createCorrelationPlot(sc,'Myoblasts')
p2 = createCorrelationPlot(sc,'Pericytes')
p3 = createCorrelationPlot(sc,'HUVEC')
cowplot::plot_grid(p1[[1]],p2[[1]])





##
for (i in 1:length(levels(sc2$idents))) {
  
  cl <- subset(sc2, idents= levels(sc2$idents)[i])
  avg = data.frame(S1 = log1p(rowMeans(cl@assays$RNA@data[,cl$source=='My.Pr.Hu'])),
                   S2 = log1p(rowMeans(cl@assays$RNA@data[,cl$source=='Pl.My.Pr.Hu'])),
                   gene = rownames(cl@assays$RNA@data),row.names=rownames(cl@assays$RNA@data))
  Idents(cl) = cl$source
  
  avg$label = avg$gene
  avg$label[abs(avg[,1]-avg[,2]) < 0.5] = ''
  
  p[[i]] = ggplot(avg, aes(S1,S2,label=label)) + geom_point()+theme_classic()+xlab('No Placenta')+ylab('With Placenta')+ggtitle(levels(sc@active.ident)[i])+
    geom_text_repel(label=avg[,'label'],box.padding = 0.5, max.overlaps = Inf)+geom_point(data = avg[avg$label != "",], color = "red")
  markers[[i]] <- FindMarkers(cl, ident.1 = "My.Pr.Hu", ident.2 = "Pl.My.Pr.Hu", verbose = FALSE,min.diff.pct=0.2)
  markers[[i]]$Cluster = i
  markers[[i]]$Gene = rownames(markers[[i]])
  markers[[i]]$rank = 1:nrow(markers[[i]])
}
dev.off()

markers = do.call(rbind,markers)


res <- marrangeGrob(p, nrow = 3, ncol = 3)
ggexport(res, filename = "~/Documents/levenberg/scatters_small.png")

res <- marrangeGrob(p[1:9], nrow = 3, ncol = 3)
ggexport(res, filename = "~/Documents/levenberg/scatters1.png")
res <- marrangeGrob(p[10:13], nrow = 2, ncol = 3)
ggexport(res, filename = "~/Documents/levenberg/scatters2.png")

dittoPlot(sc,var='IGFBP2',plots='vlnplot',group.by = 'seurat_clusters',color.by = 'source')

dittoPlot(sc,var='VEGFA',plots='ridgeplot',group.by = 'source')


dittoDotPlot(sc, vars = rownames(sc)[grepl('MYO',rownames(sc))], group.by = "SCT_snn_res.0.2",max.color='blue',min.color='white')+ylab('Clusters')


##

sc2 = sc[,!(sc$SCT_snn_res.0.2 %in% c(4,6))]

sc2 <- SCTransform(sc2, vars.to.regress = "percent.mt", verbose = FALSE)

sc2 <- RunPCA(sc2)
ElbowPlot(sc2)

sc2 <- RunUMAP(sc2, reduction = "pca", dims = 1:15)
p1 = dittoDimPlot(sc2, reduction.use = "umap", "source",size = 0.5)
p1

sc2 <- FindNeighbors(sc2, dims = 1:15)
sc2 <- FindClusters(sc2, resolution = 0.2)
p2 = dittoDimPlot(sc2, reduction.use = "umap", "SCT_snn_res.0.2",size = 0.5)

p1+p2



###

expr = data.table::fread('~/Documents/levenberg/GSE143704_DeMicheli_HumanMuscleAtlas_normalizeddata.txt',header=T,sep='\t')
expr <- as.matrix(expr,rownames=1)
md = read.table('~/Documents/levenberg/GSE143704_DeMicheli_HumanMuscleAtlas_metadata.txt',header=T,row.names = 1,sep='\t')



sc <- FindClusters(sc, resolution = 80)

singler_results <- SingleR(test = GetAssayData(sc, assay = 'RNA', slot = 'data'),
                           ref = expr,
                           labels = md$cell_annotation
                           # ,clusters = sc@meta.data$SCT_snn_res.80
)
x <- data.frame(row.names = rownames(singler_results), labels=singler_results$labels)
sc <- AddMetaData(sc, metadata=x[sc$SCT_snn_res.80, 1], col.name = "HumanMuscleAtlas")

p <- dittoDimPlot(sc, reduction.use = "umap", "HumanMuscleAtlas", main = "HumanMuscleAtlas")

tbl = table(sc$HumanMuscleAtlas,sc$idents)
tbl = t(t(tbl)/colSums(tbl))
corrplot(tbl,is.corr = F)

p1 <- dittoDimPlot(sc, reduction.use = "umap", "idents", main = "Clusters")

p2 <- dittoDimPlot(sc, reduction.use = "umap", "HumanMuscleAtlas", main = "HumanMuscleAtlas")

dittoDimPlot(sc,reduction.use = "umap", "idents",split.by = 'source', main = "Clusters")



###
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- ReadVelocity("~/Documents/levenberg/SL_111.loom")

sl <- as.Seurat(x = ldat)
sl <- SCTransform(object = sl, assay = "spliced")
sl <- RunPCA(object = sl, verbose = FALSE)
sl <- FindNeighbors(object = sl, dims = 1:20)
sl <- RunUMAP(object = sl, dims = 1:20)
sl <- FindClusters(object = sl,resolution = 0.2)


sl <- RunVelocity(object = sl, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sl)))
names(x = ident.colors) <- levels(x = sl)
cell.colors <- ident.colors[Idents(object = sl)]
names(x = cell.colors) <- colnames(x = sl)
show.velocity.on.embedding.cor(emb = Embeddings(object = sl, reduction = "umap"), vel = Tool(object = sl, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

dittoDimPlot(sl, reduction.use = "umap", "seurat_clusters")

markers.cluster <- FindAllMarkers(sl, only.pos = TRUE, min.pct = 0.25,min.diff.pct = 0.1, logfc.threshold = 0.25)
top5 <- markers.cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
dittoDotPlot(sl, vars = unique(top5$gene), group.by = "SCT_snn_res.0.2",max.color='blue',min.color='white')+ylab('Clusters')

saveRDS(sl,file='~/Documents/levenberg/sl_withplacenta_velo_seurat.rds')

expr = data.table::fread('~/Documents/levenberg/GSE143704_DeMicheli_HumanMuscleAtlas_normalizeddata.txt',header=T,sep='\t')
expr <- as.matrix(expr,rownames=1)
md = read.table('~/Documents/levenberg/GSE143704_DeMicheli_HumanMuscleAtlas_metadata.txt',header=T,row.names = 1,sep='\t')

singler_results <- SingleR(test = GetAssayData(sl, assay = 'RNA', slot = 'data'),
                           ref = expr,
                           labels = md$cell_annotation
                           #,clusters = sc@meta.data$SCT_snn_res.80
)
x <- data.frame(row.names = rownames(singler_results), labels=singler_results$labels)
sc <- AddMetaData(sc, metadata=x[sc$SCT_snn_res.80, 1], col.name = "HumanMuscleAtlas")




#### new analysis
sc = readRDS('~/Documents/levenberg/levengberg_new.rds')
sc <- FindClusters(sc, resolution = 1)
tbl = table(sc$source,sc$seurat_clusters)
dittoDimPlot(sc, reduction.use = "umap", "seurat_clusters",size = 0.5,do.label = T)+ggtitle('By clusters')

corrplot(t(ceiling(100*t(tbl)/colSums(tbl))),is.corr=F,method = 'square', addCoef.col = 'black')
new.cluster.ids <- c('0'='Myoblasts',
                     '1'='Myoblasts',
                     '2'='Pericytes',
                     '3'='Myoblasts',
                     '4'='Pericytes',
                     '5'='Placenta',
                     '6'='Myoblasts',
                     '7'='Pericytes',
                     '8'='Myoblasts',
                     '9'='Pericytes',
                     '10'='Myoblasts',
                     '11'='Myoblasts',
                     '12'='Pericytes',
                     '13'='Myoblasts',
                     '14'='Pericytes',
                     '15'='Pericytes',
                     '16'='Myoblasts',
                     '17'='Pericytes',
                     '18'='Pericytes',
                     '19'='HUVEC',
                     '20'='Pericytes'
)

#new.cluster.ids <- c('Pericyte','Myoblast','Myoblast','Pericyte','Myoblast','Myoblast','Placenta','Pericyte','Myoblast','Pericyte','Myoblast','Myoblast','Myoblast','Myoblast','HUVEC')
new.cluster.ids2 = new.cluster.ids
new.cluster.ids = make.unique(new.cluster.ids)
names(new.cluster.ids) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids)
sc$idents= sc@active.ident
names(new.cluster.ids2) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids2)
sc$idents2 = sc@active.ident

runSeurat = function(sc) {
  sc <- SCTransform(sc, vars.to.regress = "percent.mt", verbose = FALSE)
  sc <- RunPCA(sc)
  #ElbowPlot(sc)
  sc <- RunUMAP(sc, reduction = "pca", dims = 1:30)
  sc <- FindNeighbors(sc, dims = 1:30)
  sc <- FindClusters(sc, resolution = 0.5)
  sc
}

#sc = runSeurat(subset(sc,source %in% c('My.Pr.Hu', 'Pl.My.Pr.Hu')))

library(monocle3)
library(SeuratWrappers)
library(patchwork)

my = runSeurat(subset(sc,idents2=='Myoblasts' & source=='Pl.My.Hu'))

cds <- as.cell_data_set(my)
cds <- cluster_cells(cds, resolution=1e-2)

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
plot_cells(cds,
           color_cells_by = "idents",
           label_groups_by_cluster=F,
           label_leaves=F,
           label_branch_points=F,cell_size = 1,
           label_roots = F,label_cell_groups = F)

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 4]))
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")

cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)

rowData(cds)$gene_short_name <- row.names(rowData(cds))

head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)

deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)

##
library(SingleCellExperiment)
library(destiny)
my = runSeurat(subset(sc,idents2=='Myoblasts'))
Idents(my) = 'source'
my2 = runSeurat(subset(x = my, downsample = 1000))

sce <- as.SingleCellExperiment(my2)
dm <- DiffusionMap(sce, verbose = TRUE)
library(ggplot2)
cellLabels <- sce$idents
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  DC4 = eigenvectors(dm)[, 4],
                  Samples = cellLabels)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Samples)) +
  geom_point()  + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

sce$pseud_dm1 <- rank(eigenvectors(dm)[,1])      # rank cells by their dpt dm1
sce$pseud_dm2 <- rank(eigenvectors(dm)[,2])      # rank cells by their dpt dm2
sce$pseud_dm1R <- rank(-eigenvectors(dm)[,1])    # rank cells by their dpt dm1 reverse order
sce$pseud_dm2R <- rank(-eigenvectors(dm)[,2])    # rank cells by their dpt dm2 reverse order

SortedDM1 <- data.frame(DM1Sort = as.data.frame(colData(sce))$pseud_dm1,
                        Samples = as.data.frame(colData(sce))$idents)
SortedDM2 <- data.frame(DM2Sort = as.data.frame(colData(sce))$pseud_dm2,
                        Samples = as.data.frame(colData(sce))$idents)
SortedDM1R <- data.frame(DM1SortR = as.data.frame(colData(sce))$pseud_dm1R,
                         Samples = as.data.frame(colData(sce))$idents)
SortedDM2R <- data.frame(DM2SortR = as.data.frame(colData(sce))$pseud_dm2R,
                         Samples = as.data.frame(colData(sce))$idents)

ggplot(SortedDM1, aes(x=SortedDM1[,1], y=Samples,color=Samples)) +
  geom_jitter() + xlab("Diffusion component 1 (DC1)") + ylab("Samples") +
  ggtitle("Cells ordered by DC1")
ggplot(SortedDM2, aes(x=SortedDM2[,1], y=Samples,color=Samples)) +
  geom_jitter() + xlab("Diffusion component 2 (DC2)") + ylab("Samples") +
  ggtitle("Cells ordered by DC2")

ggplot(SortedDM1R, aes(x=SortedDM1R[,1], y=Samples,color=Samples)) +
  geom_jitter() + xlab("Minus Diffusion component 1 (DC1)") + ylab("Samples") +
  ggtitle("Cells ordered by reversed DC1")
ggplot(SortedDM2R, aes(x=SortedDM2R[,1], y=Samples,color=Samples)) +
  geom_jitter() + xlab("Minus Diffusion component 2 (DC2)") + ylab("Samples") +
  ggtitle("Cells ordered by reversed DC2")

##

library(slingshot)

sce <- as.SingleCellExperiment(my2)

geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)
pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')
plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)
cl1 <- mclust::Mclust(rd1)$classification
colData(sce)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[as.factor(as.character(sce$idents))], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

Idents(my) = my$idents
markers.cluster <- FindAllMarkers(my, only.pos = TRUE, min.pct = 0.25,min.diff.pct = 0.1, logfc.threshold = 0.25)
top5 <- markers.cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p4.myo = dittoDotPlot(my, vars = unique(top5$gene), group.by = "idents",max.color='blue',min.color='white')+ylab('Clusters')
p4.myo
