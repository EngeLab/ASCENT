# setwd("/wrk/data/martin/all/bcell_trajectory/peds-bm-atlas")
# source('scDataAnalysis_Utilities.R')
# remotes::install_github("carmonalab/UCell")
setwd("/wrk/resources/dntr2")
library(tidyverse)
library(Seurat)
# library(UCell)

## build the normal trajectory (in process_healthydonor_final.R) ####
seurat.rna <- readRDS("/wrk/data/martin/all/bcell_trajectory/peds-bm-atlas/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC.rds")
DimPlot(seurat.rna,group.by="Ctype")

# ## revovery all the previous transformation done for thess healthy donor cells ####
# vegs = VariableFeatures(seurat.rna)
# feature.load = seurat.rna@reductions$pca@feature.loadings
# # get the coefficient from ScaleData for each gene
# coef.lm = centers.vegs = sd.vegs = list()
# data0 = seurat.rna@assays$RNA@data
# mdata0 = seurat.rna@meta.data
# for(gene0 in vegs){
#   expr0 = data0[gene0, ]
#   perc.mito0 = mdata0[, 'perc.mito']
#   numi0 = mdata0[, 'nCount_RNA']
#   model0 = lm(expr0 ~ perc.mito0 + numi0)
#   residual0 = residuals(model0)
#   coef.lm[[gene0]] = coef(model0)
#   
#   centers.vegs[[gene0]] = mean(residual0)
#   sd.vegs[[gene0]] = sd(residual0)
# }
# coef.lm = do.call('rbind', coef.lm)
# centers.vegs = do.call('c', centers.vegs)
# sd.vegs = do.call('c', sd.vegs)
# 
# ## save coefficients for future use
# save(coef.lm, centers.vegs, sd.vegs, file = 
#        'coef_hd.RData')
# 
# 
# ## recovery the umap model
# set.seed(42) 
# npc = 20
# umap_model = uwot::umap(as.matrix(seurat.rna@reductions$pca@cell.embeddings[, 1:npc]),
#                         n_neighbors = 30, min_dist = 0.3,
#                         ret_model = T, metric = 'cosine')

health_ann = seurat.rna@meta.data
health_ann = subset(health_ann, select = c('sample', 'nCount_RNA', 'nFeature_RNA', 'perc.mito',  'Ctype'))
health_ann = data.table::data.table(health_ann, keep.rownames = T)
#setkey(health_ann, rn) # Fucks with ordering.
Idents(object = seurat.rna) <-  "Ctype"
# seurat.rna <- ScaleData(seurat.rna, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat.rna))
seurat.b.rna <- seurat.rna[,health_ann$Ctype %in% c("HSPC","LMPP","CLP","Pre-pro-B","Pro-B","Pre-B","Immature-B","Mature-B")]
DimPlot(seurat.rna, group.by='Ctype')
a <- DimPlot(seurat.b.rna, group.by='Ctype')
b <- FeaturePlot(seurat.b.rna, "pseudotime")
plot_grid(a,b)

# seurat.b.rna <- ScaleData(seurat.b.rna, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat.b.rna))

# S/M phase regressed out in b-cell lineage only
# a <- FindAllMarkers(seurat.b.rna, only.pos=TRUE, test.use='LR', slot="scale.data")
# a.default <- FindAllMarkers(seurat.b.rna, only.pos=TRUE, slot="scale.data")

seurat.b.rna.2 <- seurat.rna[,health_ann$Ctype %in% c("HSPC","LMPP","CLP","Pre-pro-B","Pro-B","Pre-B","Immature-B","Mature-B")]

# S/M phase regressed out in whole data set
# a.2 <- FindAllMarkers(seurat.b.rna.2, only.pos=TRUE, test.use='LR', slot="scale.data")
# a.2.default <- FindAllMarkers(seurat.b.rna.2, only.pos=TRUE, slot="scale.data")

# Without regressing S/M out.
# a.2.noreg <- FindAllMarkers(seurat.b.rna.2, only.pos=TRUE, slot="data")
a.2.noreg <- FindAllMarkers(seurat.b.rna.2, only.pos=TRUE, slot="data", logfc.threshold=0.5)

# marker.genes <- a.2 %>% filter(pct.2 < 0.8) %>% select(gene) %>% unlist() %>% unique()
# 
# a.adj.maxpct2 <- a.2 %>% filter(pct.2 < 0.8 & p_val_adj < 1e-10) %>% group_by(cluster) %>% slice_min(n=400, order_by=p_val_adj) %>% ungroup
# a.s <- split(a.adj.maxpct2$gene, a.adj.maxpct2$cluster)
# 
# a.def.adj.maxpct2 <- a.2.default %>% filter(pct.2 < 0.8 & p_val_adj < 1e-10) %>% group_by(cluster) %>% slice_min(n=400, order_by=p_val_adj) %>% ungroup
# a.def.s <- split(a.def.adj.maxpct2$gene, a.def.adj.maxpct2$cluster)

a.noreg.adj.maxpct2 <- a.2.noreg %>% filter(pct.2 < 0.8 & p_val_adj < 1e-10) %>% group_by(cluster) %>% slice_min(n=400, order_by=p_val_adj) %>% ungroup
table(a.noreg.adj.maxpct2$cluster)
a.noreg.s <- split(a.noreg.adj.maxpct2$gene, a.noreg.adj.maxpct2$cluster)


# all.rna <- readRDS('/wrk/data/martin/all/bcell_trajectory/all/rna_seurat.Rds')
# print(load("aonghus/rna_seurat_wAnnot211101.Rda"))
all.rna <- readRDS("results/ALL-220128/rna_seurat.Rds")
DimPlot(all.rna)

# aa <- unique(unlist(a.s)[unlist(a.s) %in% rownames(all.rna)])

# all.rna <- MyAddModuleScore(all.rna, features=a.s, pool=aa, assay='RNA', nbin=25)
# all.rna <- AddModuleScore(all.rna, features=a.s, assay='RNA', nbin=25)
all.rna <- AddModuleScore(all.rna, features=a.noreg.s, assay='RNA', nbin=25)
colnames(all.rna@meta.data)[(ncol(all.rna@meta.data)-length(a.noreg.s)+1):ncol(all.rna@meta.data)] <- paste0(names(a.noreg.s),"_score")
# all.rna <- AddModuleScore_UCell(all.rna, features=a.noreg.s, assay='RNA')
seurat.b.rna.2 <- AddModuleScore(seurat.b.rna.2, features=a.noreg.s, assay='RNA', nbin=25)
colnames(seurat.b.rna.2@meta.data)[(ncol(seurat.b.rna.2@meta.data)-length(a.noreg.s)+1):ncol(seurat.b.rna.2@meta.data)] <- paste0(names(a.noreg.s),"_score")


T.NK <- all.rna@reductions$umap@cell.embeddings[,1] > 3 & all.rna@reductions$umap@cell.embeddings[,2] < -6

score_names <- make.names(names(a.noreg.s))
colnames(all.rna@meta.data)[40:47] <- score_names

FeaturePlot(object = all.rna, features = c("HSPC","LMPP","CLP","Pre.pro.B","Pro.B","Pre.B","Immature.B","Mature.B"))
# FeaturePlot(object = all.rna, features = c("HSPC","LMPP","CLP","Pre.pro.B","Pro.B","Pre.B","Immature.B","Mature.B"), order=T)
# FeaturePlot(object = all.rna, features = c("HSPC","LMPP","CLP","Pre.pro.B","Pro.B","Pre.B","Immature.B","Mature.B"), min.cutoff="q1",max.cutoff="q99", order=T)
# FeaturePlot(object = all.rna, features = c("HSPC","LMPP","CLP","Pre.pro.B","Pro.B","Pre.B","Immature.B","Mature.B"),order = T,min.cutoff= "q10",max.cutoff = "q90")

FeaturePlot(object = all.rna, features = 'HSPC')
FeaturePlot(object = all.rna, features = 'LMPP')
FeaturePlot(object = all.rna, features = 'CLP')
FeaturePlot(object = all.rna, features = 'Pre-pro-B')
FeaturePlot(object = all.rna, features = 'Pre-B')
FeaturePlot(object = all.rna, features = 'Pro-B')
FeaturePlot(object = all.rna, features = 'Mature-B')
FeaturePlot(object = all.rna, features = 'Immature-B')

#Dimensionality reduction
my.svd <- svd(all.rna[[]][!T.NK,40:47])
plot(my.svd$u[,1:2])

my.dist <- dist(all.rna[[]][!T.NK,40:47])
#fit <- cmdscale(my.dist, eig=TRUE, k=2)
fit <- cmdscale(my.dist, eig=TRUE, k=1)
my.dist <- dist(all.rna[[]][!T.NK,48:55])
fit2 <- cmdscale(my.dist, eig=TRUE, k=1)
my.dist <- dist(all.rna[[]][!T.NK,paste0(names(a.s), ".noreg")])
fit.noreg <- cmdscale(my.dist, eig=TRUE, k=1)

my.dist <- dist(seurat.b.rna.2[[]][,43:50])
fit.ref <- cmdscale(my.dist, eig=TRUE, k=1)

library(caTools)

plotme("HSPC.noreg", obj=fit.noreg, pcol='red', lcol='darkred')
plotme("CLP.noreg", obj=fit.noreg, pcol='red', lcol='darkred')
plotme("LMPP.noreg", obj=fit.noreg, pcol='red', lcol='darkred')
plotme("Mature-B.noreg", obj=fit.noreg, pcol='red', lcol='darkred')
plotme("Immature-B.noreg", obj=fit.noreg, pcol='red', lcol='darkred')
plotme("Pre-B.noreg", obj=fit.noreg, pcol='red', lcol='darkred')
plotme("Pre-pro-B.noreg", obj=fit.noreg, pcol='red', lcol='darkred')
plotme("Pro-B.noreg", obj=fit.noreg, pcol='red', lcol='darkred')
plotme("HSPC", obj=fit, pcol='red', lcol='darkred')
plotme("Immature-B", pcol='orange', lcol='yellow', add.to.plot=T)
plotme("Pre-B", pcol='darkgreen', lcol='green', add.to.plot=T, obj=fit2)
legend("topright", legend=c('HSPC', 'Pre-B', 'Immature-B'), fill=c('red', 'green', 'yellow'))
dev2bitmap("Module_scores_by_diffscale.png", type="png16m", res=200)
plotme("Mature-B", pcol='darkblue', lcol='blue', add.to.plot=T)

plotme("LMPP", pcol='green', lcol='darkgreen', add.to.plot=T, obj=fit2)
plotme("Pre-B", obj=fit)
plotme("Pre-B", pcol='darkgreen', lcol='green', add.to.plot=T, obj=fit2)
plotme("Pro-B", pcol='blue', lcol='green', add.to.plot=T, obj=fit2)
plotme("Pre-pro-B", pcol='blue', lcol='green', add.to.plot=T)
plotme("Pre-pro-B")
plotme("Immature-B", pcol='orange', lcol='yellow', add.to.plot=T)
plotme("Mature-B", pcol='darkblue', lcol='blue', add.to.plot=T)
plotme("Mature-B", )
plotme <- function(name, pcol='black', lcol='red', add.to.plot=F, obj=fit, ...) {
  if(add.to.plot) {
    points(all.rna[[]][!T.NK,][order(obj$points),name], pch=19, cex=0.3, col=pcol)
  } else {
    plot(all.rna[[]][!T.NK,][order(obj$points),name], pch=19, cex=0.3, col=pcol, ...)
  }
  lines(runmean(all.rna[[]][!T.NK,][order(obj$points),name], k=100), lwd=3, col=lcol)
}

annot <- read.table("/wrk/data/martin/all/bcell_trajectory/all/all_metadata_long-211101.tsv", sep="\t", header=T)

rownames(annot) <- annot$cell_id
annot.all.rna <- annot[colnames(all.rna@assays$RNA),]

library(beeswarm)
boxplot(split(fit$points, paste(annot.all.rna$patient_id, annot.all.rna$timepoint)[!T.NK]))
beeswarm(split(fit.noreg$points, paste(annot.all.rna$patient_id, annot.all.rna$timepoint)[!T.NK]))
boxplot(split(fit$points, annot$cell_phase))
boxplot(split(fit$points, annot$seurat_clusters))

o <- !T.NK & !grepl('MRD', annot.all.rna$patient_id)
o2 <- !grepl('MRD', annot.all.rna$patient_id[!T.NK])
o <- !T.NK & grepl('MRD', annot.all.rna$patient_id)
o2 <- grepl('MRD', annot.all.rna$patient_id[!T.NK])
beeswarm(split(fit.noreg$points[o2], paste(annot.all.rna$patient_id, annot.all.rna$timepoint)[o]))
png(filename="all_bees_MRD.png", width=10000, height=600, res=50)
beeswarm(split(fit.noreg$points[o2], paste(annot.all.rna$dna_class, annot.all.rna$timepoint)[o]))
dev.off()

# Add Module scores to original object
seurat.b.rna <- AddModuleScore(seurat.b.rna, features=a.s, assay='RNA', nbin=25)
seurat.b.rna <- AddModuleScore(seurat.b.rna, features=a.def.s, assay='RNA', nbin=25)
seurat.b.rna <- AddModuleScore(seurat.b.rna, features=a.noreg.s,  assay='RNA', nbin=25)

colnames(seurat.b.rna@meta.data)[43:50] <- paste0(names(a.s))
colnames(seurat.b.rna@meta.data)[59:66] <- paste0(names(a.s), ".def")
colnames(seurat.b.rna@meta.data)[67:74] <- paste0(names(a.s), ".noreg")
colnames(seurat.b.rna@meta.data)[75:82] <- paste0(names(a.s), ".noreg2")
colnames(seurat.rna@meta.data)[43:50] <- paste0(names(a.s))
DimPlot(seurat.b.rna, group.by='Ctype')

FeaturePlot(object = seurat.b.rna, features = paste0(names(a.s), ".2"))
FeaturePlot(object = seurat.b.rna, features = paste0(names(a.s), ".def"))
FeaturePlot(object = seurat.b.rna, features = paste0(names(a.s), ".noreg"))
FeaturePlot(object = seurat.b.rna, features = paste0(names(a.s), ".noreg2"))
FeaturePlot(object = seurat.b.rna, features = 'HSPC')
FeaturePlot(object = seurat.b.rna, features = 'LMPP')
FeaturePlot(object = seurat.b.rna, features = 'CLP')
FeaturePlot(object = seurat.b.rna, features = 'Pre-pro-B')
FeaturePlot(object = seurat.b.rna, features = 'Pre-B')
FeaturePlot(object = seurat.b.rna, features = 'Pro-B')
FeaturePlot(object = seurat.b.rna, features = 'Mature-B')
FeaturePlot(object = seurat.b.rna, features = 'Immature-B')

DimPlot(seurat.rna, group.by='Ctype')
FeaturePlot(object = seurat.rna, features = names(a.s))

FeaturePlot(object = seurat.b.rna, features = paste0(names(a.s), '.noreg2'))

my.dist.s <- dist(seurat.b.rna[[]][,43:50])
my.dist.s <- dist(seurat.b.rna[[]][,paste0(names(a.s), '.noreg2')])
fit.noreg <- cmdscale(my.dist.s, eig=TRUE, k=1)
beeswarm(split(fit.noreg$points, seurat.b.rna@meta.data[,'Ctype_updated']))
