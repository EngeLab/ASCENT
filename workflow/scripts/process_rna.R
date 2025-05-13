library(tidyverse)
library(Seurat)

#Since seurat updated accidentally 
options(Seurat.object.assay.version = 'v3')

# Convenience functions 
dna_to_cellid <- function(x){ sub("([A-Z]{3}[0-9]{5})(D?)_([A-Z][0-9]{2}).*","\\1.\\3", x) }
rna_to_cellid <- function(x){ sub("([A-Z]{3}[0-9]{5})(R?).([A-Z][0-9]{2}).*","\\1.\\3", x) }
cellid_to_plateid <- function(x){ sub("([A-Z]{3}[0-9]{5}).*","\\1", x) }
norm.log.counts <- function(counts) {
  norm.fact <- colSums(counts)
  counts.norm <- t( apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
  counts.log <- log2(counts.norm)
  rownames(counts.log) <- rownames(counts)
  counts.log
}
get.cutoff.lognorm <- function(my.counts.log, quantile.cut=0.001, gene.name=NULL) {
  cl.act <- my.counts.log[gene.name,]
  cl.act.m <- median(cl.act)
  cl.act.sd <- sqrt(sum((cl.act[cl.act > cl.act.m] - cl.act.m)^2)/(sum(cl.act  > cl.act.m)-1))
  my.cut <- qnorm(p=quantile.cut, mean=cl.act.m, sd=cl.act.sd)
  my.cut
}

# Filtering parameters
mincount <- as.integer(snakemake@params[["min_count"]]) # Set very low for immune cells
mingenes <- as.integer(snakemake@params[["min_genes"]]) # Only drops ~30 extra cells after count + ACTB filter
minquant <- as.numeric(snakemake@params[["quantile_cut"]])
genequant <- as.character(snakemake@params[["quantile_gene"]])
perpl_rna <- as.integer(snakemake@params[["rna_perplexity"]])

# Metadata: without filters it has 5-cell bulk and 0-cell controls (NA)
meta <- read_tsv(snakemake@input[["meta_wsort"]])
dim(meta)
meta_wfacs <- read_tsv(snakemake@input[["meta_wfacs"]])
dim(meta_wfacs)

# Counts
counts <- data.table::fread(snakemake@input[["counts"]], sep="\t", stringsAsFactors=F, header=T, data.table=F)
rownames(counts) <- counts[[1]]
counts[[1]] <- NULL
colnames(counts) <- rna_to_cellid(colnames(counts)) # Convert to cell_id
dim(counts)

# Check cycle gene ids
if(any(grepl("^ENSG", rownames(counts)))){
  gene_info <- read_tsv(snakemake@params[["gene_info"]], skip = 1, col_names = c("gene_id", "gene_name", "gene_type"))
  # Ensembl gene ids: convert seurat cc genes
  cc_genes <- lapply(cc.genes.updated.2019, function(x){
    gene_info$gene_id[match(x, gene_info$gene_name)]
  })
} else {
  cc_genes <- cc.genes.updated.2019
}


# TODO: Keep only single cells (drop 5-cell bulks and non-sorted wells)
if("Sorted Count" %in% colnames(meta)){
  meta.f <- filter(meta, `Sorted Count`==1 | is.na(`Sorted Count`) | (patient_id %in% c("MRD_ALL35","ALL3","ALL5") & !grepl("A01$|P24$",cell_id)) )
  meta_wfacs.f <- filter(meta_wfacs, `Sorted Count`==1 | is.na(`Sorted Count`) | (patient_id %in% c("MRD_ALL35","ALL3","ALL5") & !grepl("A01$|P24$",cell_id)) )
} else {
  meta.f <- meta
}
dim(meta.f)

include_cells <- colnames(counts) %in% meta.f$cell_id
ercc.idx <- grepl("^ERCC-", rownames(counts)) # ERCC spike-ins
counts.f <- as.matrix(counts[!grepl("^__|_",rownames(counts)) & !ercc.idx, include_cells]) # Throws HTseq rows, ERCC and underscored genes (2)
counts.ercc <- as.matrix(counts[ercc.idx,include_cells])
ercc.frac <- colSums(counts.ercc)/(colSums(counts.f)+colSums(counts.ercc))

# Get mt genes if available
if(!is.null(snakemake@params[["mt_genes"]])){
  mtgenes <- readLines(snakemake@params[["mt_genes"]])  
  mt.idx <- rownames(counts) %in% mtgenes # MT-genes
  mt.frac <- colSums(counts[mt.idx,])/colSums(counts)
} else {
  mt.frac <- NULL
}

counts.f[1:5,1:5]
dim(counts.f)

# Metadata for Seurat obj
meta_seurat <- data.frame(cell_id=colnames(counts.f)) %>% 
  mutate(ercc_frac=ercc.frac[match(cell_id, names(ercc.frac))],
         mt_frac=ifelse(!is.null(mt.frac), 
                                 mt.frac[match(cell_id, names(mt.frac))], 
                                 NA)
  ) %>% left_join(meta, by="cell_id")

rownames(meta_seurat) <- meta_seurat$cell_id
colnames(meta_seurat) <- make.names(colnames(meta_seurat))

# Load Seurat obj
filtered.counts <- counts.f %>% 
  subset(rowSums(.)>0) %>% # Drop genes (rows) with 0 counts in all samples
  subset(select=colSums(.)>mincount) %>% # Drop cells (cols) with total read count < minimum count
  subset(select=colSums(. > 0)>mingenes) %>% # Drop cells with genes below minimum gene (Default 500)
  as.matrix()
dim(filtered.counts)
n_countsfilt <- ncol(counts.f)-ncol(filtered.counts)
countsfilter <- setdiff(colnames(counts.f), colnames(filtered.counts))
counts.log1 <- norm.log.counts(filtered.counts) # Log2 of counts per million mapped reads (log2 CPM)
# Apply 2nd filter, low quantile gene (eg ACTB) expression
actb_cutoff <- get.cutoff.lognorm(counts.log1, quantile.cut=minquant, gene.name=genequant)
good.cells <- counts.log1[genequant,] > actb_cutoff
final.counts <- filtered.counts[,good.cells]
cat("^Dropped ",n_countsfilt," (low counts) + ",sum(!good.cells)," (low ", genequant, ") of ", ncol(counts.f)," cells (",round((n_countsfilt+sum(!good.cells))/ncol(counts.f)*100,1),"%)\n",sep="")

# QC plots from filtering
colsumhist <- colSums(counts.f) %>% enframe("cell_id","colsum") %>% mutate(plate_id=cellid_to_plateid(cell_id))
phist1 <- colsumhist %>% ggplot(aes(log10(colsum+1))) + ggridges::geom_density_ridges(aes(y=plate_id)) + 
  geom_vline(xintercept=log10(mincount), color="red") + labs(x="log10(read counts)", y="Plate ID")
genehist <- colSums(counts.f > 0) %>% enframe("cell_id","genes") %>% mutate(plate_id=cellid_to_plateid(cell_id))
phist2 <- genehist %>% ggplot(aes(genes)) + ggridges::geom_density_ridges(aes(y=plate_id)) + 
  geom_vline(xintercept=mingenes, color="red") + labs(x="genes detected", y="Plate ID")
actbhist <- counts.log1[genequant,] %>% enframe("cell_id","ACTB") %>% mutate(plate_id=cellid_to_plateid(cell_id))
phist3 <- actbhist %>% ggplot(aes(ACTB)) + ggridges::geom_density_ridges(aes(y=plate_id)) + 
  geom_vline(xintercept=actb_cutoff, color="red") + labs(x="log2(ACTB+1)",y="Plate ID")
pqc <- cowplot::plot_grid(phist1,phist2,phist3,ncol=3)
pdf(snakemake@output[["plot_qc"]],height=10,width=6)
plot(pqc)
dev.off()

pcadims <- ifelse(ncol(final.counts) < 100, 10, 30)
perpl <- ifelse(ncol(final.counts) < 100, floor(ncol(final.counts)/3), 50)
clustres <- ifelse(ncol(final.counts) < 100, 0.5, 1.5)
neighb <- ifelse(ncol(final.counts) < 30, 5, 30L)

dat <- CreateSeuratObject(final.counts, names.field=1, names.delim="\\.", assay="RNA")
dat <- AddMetaData(dat, meta_seurat)
dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 1e4, verbose=F)
dat <- FindVariableFeatures(dat, selection.method="dispersion", nfeatures = 2000, verbose=F) # Seurat feature finder (vst, dispersion)
dat <- ScaleData(dat, verbose=F, rownames(dat))
dat <- RunPCA(dat, npcs=pcadims, verbose=F)
dat <- FindNeighbors(dat, dims=1:pcadims, verbose=F)
dat <- FindClusters(dat, resolution = clustres, algorithm=1, n.iter=1000, n.start=100, verbose=F)
dat <- RunTSNE(dat, n.components=2, dims=1:pcadims,perplexity=perpl, verbose=F)
dat <- RunUMAP(dat, n.components=2, dims=1:pcadims,verbose=F,n.neighbors = neighb)

# Add cell cycle scoring
cc.cutoff <- snakemake@params[["cc_expr_cutoff"]]
sgenes <- cc_genes$s.genes[cc_genes$s.genes %in% rownames(dat)]
mgenes <- cc_genes$g2m.genes[cc_genes$g2m.genes %in% rownames(dat)]
cc.df <- data.frame(mscore=log2(colSums(as.matrix(dat@assays$RNA@data[mgenes,]))), sscore=log2(colSums(as.matrix(dat@assays$RNA@data)[sgenes,])))
cc.df.cl <- cc.df %>% 
  rownames_to_column("cell_id") %>% 
  mutate(cell_phase=case_when(mscore < cc.cutoff & sscore < cc.cutoff ~ "G1",
                              mscore > sscore ~ "G2M",
                              TRUE ~ "S"))
cells.phase <- select(cc.df.cl, cell_id, cell_phase)
write_tsv(cells.phase, snakemake@output[["rna_phase"]])
# Add to Seurat object
dat$cell_phase <- cc.df.cl$cell_phase
dat$m_score <- cc.df.cl$mscore
dat$s_score <- cc.df.cl$sscore


pdf(snakemake@output[["dimplot"]], width=8, height=8)
DimPlot(dat, reduction="tsne", group.by="patient_id", label=T)
DimPlot(dat, reduction="umap", group.by="patient_id", label=T)
if("timepoint" %in% colnames(dat@meta.data) && !all(is.na(dat$timepoint))) DimPlot(dat, reduction="umap", group.by="timepoint")
DimPlot(dat, reduction="tsne", group.by="plate_id")
DimPlot(dat, reduction="umap", group.by="plate_id")
DimPlot(dat, reduction="tsne", group.by="cell_phase") + scale_color_manual(values=c("G2M"="#757de8", "S"="#d73027", "G1"="#cfcfcf"))
DimPlot(dat, reduction="umap", group.by="cell_phase") + scale_color_manual(values=c("G2M"="#757de8", "S"="#d73027", "G1"="#cfcfcf"))
dev.off()

write_rds(dat, snakemake@output[["seurat_obj"]])
