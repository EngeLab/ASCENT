# Clustering and heatmap of LLR breakpoints in DNTR-seq
library(tidyverse)
library(uwot)
library(dbscan)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(GenomicRanges)
library(copynumber)
library(ape)
library(ggtree)
library(scales)
options(scipen = 999)

# Functions
source(snakemake@params[["src_general"]])
manhattan.dist <- function(factor, segs, na.rm=F) {
  if(na.rm) segs <- segs[!is.na(segs)]
  sum(abs((segs*factor)-round(segs*factor)), na.rm = na.rm)
}

# Subset for cells of interest
cells.all <- colnames(data.table::fread(snakemake@input[["copynumber"]], header = T, nrows = 1))
cells.idx <- which(cells.all %in% snakemake@params[["cells"]])

cn <- data.table::fread(snakemake@input[["copynumber"]], select = cells.idx) %>% as.matrix
counts <- data.table::fread(snakemake@input[["counts"]], select = cells.idx) %>% as.matrix

cat(snakemake@wildcards[["patient_id"]],"n=",ncol(counts),"with",nrow(counts),"bins (before filtering)\n")
dim(counts)

# Sanity check
if (ncol(cn) != ncol(counts)) stop("Copy number and raw count matrix do not match. Exiting...")
if (!all(colnames(cn) == colnames(counts))) stop("Discrepant cell names between cn and count matrix. Exiting...")

if(!is.null(snakemake@input[["rna_phase"]])){
  cells.phase <- read_tsv(snakemake@input[["rna_phase"]]) %>% 
    filter(cell_id %in% dna_to_cellid(colnames(cn))) %>%  # Keep only cells also in the copy number matrix == snakemake@params[["cells"]]
    mutate(dna_library_id=colnames(cn)[match(cell_id,dna_to_cellid(colnames(cn)))])
  cells.replicating <- cells.phase$cell_id[cells.phase$cell_phase %in% c("S", "G2M", "G2/M", "M", "G2")] 
} else {
  cells.phase <- tibble(cell_id = NA, cell_phase = NA, dna_library_id = NA)
  cells.replicating <- c()
}

bc <- colSums(counts)

# Load and filter bins by mappability and GC content
bins <- read_tsv(snakemake@input[["bins"]], col_names = c("chr", "start", "end")) %>%
  mutate(
    chr = factor(chr, levels = unique(chr)),
    bin_unfilt = seq_along(chr),
    chr_short = sub("chr", "", chr),
    chr_int = as.integer(case_when(chr_short == "X" ~ "23", chr_short == "Y" ~ "24", TRUE ~ chr_short))
  )
bins$end_cum <- cumsum(as.numeric(bins$end)-as.numeric(bins$start))
bins.global <- bins %>% mutate(end_cum = cumsum(as.numeric(end) - as.numeric(start)))
map <- as.numeric(readLines(snakemake@input[["map"]]))
gc <- as.numeric(readLines(snakemake@input[["gc"]]))
if(!is.null(snakemake@input[["badbins"]])) badbins <- read.table(snakemake@input[["badbins"]],header=T,stringsAsFactors=F) else badbins <- NULL
good.bins <- gc > snakemake@params[["gc_min"]] & map > snakemake@params[["map_min"]] & !bins$bin_unfilt %in% badbins$bin_unfilt
bins.f <- bins[good.bins, ] %>% mutate(bin=1:nrow(.))
gc.f <- gc[good.bins]

dim(cn)
cn[1:5, 1:5]
cells_per_plate <-table(cellid_to_plateid(colnames(cn)), useNA = "always")
cells_per_plate

# Filter out cells with low total reads (and rare cells with negative CNs because of scaling errors)
good.cells <- bc > snakemake@params[["counts_min"]] & !apply(cn, 2, function(x) { any(x < 0) })
good.cells.ids <- names(bc)[which(good.cells)]
cat("Good cells by bincount:\n")
table(good.cells)
cells_per_plate_qc <- table(cellid_to_plateid(colnames(cn[, good.cells])), useNA = "always")
if (length(cells_per_plate_qc) == length(cells_per_plate)) {
  cat("^^Good cells as fractions per plate:\n")
  round(cells_per_plate_qc / cells_per_plate, 3)
} else {
  cat("^^Good cells per plate:\n")
  cells_per_plate_qc
}

# NEW: Keep only non-replicating cells (according to RNA): any cells not classified as S or G2M by Seurat
good.cells2 <- good.cells & !dna_to_cellid(names(bc)) %in% cells.replicating
good.cells2.ids <- names(bc)[which(good.cells2)]
cells.replicating.pass <- dna_to_cellid(names(bc)[good.cells & dna_to_cellid(names(bc)) %in% cells.replicating])
cat(
  "Number of replicating cells:",
  length(cells.replicating.pass),
  "(total",length(cells.replicating),")\n"
)
cat("Good cells by bincount and non-cycling phase:", sum(good.cells2), "\n")

cn.f <- cn[, good.cells2]
dim(cn.f)

cells_per_plate_qc2 <- table(cellid_to_plateid(colnames(cn.f)), useNA = "always")
if (length(cells_per_plate_qc2) == length(cells_per_plate)) {
  cat("^^Final set as fractions per plate:\n")
  round(cells_per_plate_qc2 / cells_per_plate, 3)
} else {
  cat("^^Final set, cells per plate:\n")
  cells_per_plate_qc2
}

### If too few cells per group/patient -- do not cluster anything, just plot and add all to one class (=0)
# Run UMAP on multiple components (since its done for classification, not visualization primarily)
# NOTE: No truncation of copy number is done before classification

# Run UMAP
if(snakemake@params[["umap_log"]]) cn.umap <- log(cn.f + 1) else cn.umap <- cn.f
set.seed(42)
start_time <- Sys.time()
cat("Clustering", snakemake@wildcards[["patient_id"]],"(dim:",dim(cn.umap),"; log",snakemake@params[["umap_log"]],") with neigbours=",snakemake@params[["umap_neighbors"]],"\n")
umap <- umap(t(cn.umap), n_threads = snakemake@threads, n_sgd_threads = 0, 
             metric=snakemake@params[["umap_metric"]], 
             n_neighbors=snakemake@params[["umap_neighbors"]], 
             n_components=snakemake@params[["umap_n_components"]], 
             min_dist=snakemake@params[["umap_mindist"]], 
             spread=snakemake@params[["umap_spread"]],
             pca=NULL)
end_time <- Sys.time()
cat("Time (min):", difftime(end_time, start_time, units="mins"),"\n")

# Classify UMAP embeddings with HBDSCAN (https://umap-learn.readthedocs.io/en/latest/clustering.html#umap-enhanced-clustering)
hc <- hdbscan(umap, minPts = snakemake@params[["min_cells"]])

# Add any replicating cells back
if (length(cells.replicating.pass) > 0) {
  rep_cells.df <- cells.phase %>% 
    filter(cell_id %in% cells.replicating.pass) %>% 
    select(dna_library_id, dna_class_byPatient=cell_phase)
} else {
  rep_cells.df <- tibble()
}

# Collect clone assignments and add back in replicating cells (held out for classification)
cl.df <- data.frame(
  dna_library_id = colnames(cn.umap),
  dna_class_byPatient = as.character(hc$cluster),
  dna_class_prob = hc$membership_prob,
  dna_class_outlierScore = hc$outlier_scores
) %>% 
  bind_rows(rep_cells.df) %>% # Add in replicating cells
  mutate(
    patient_id = snakemake@wildcards[["patient_id"]],
    dna_class = paste(patient_id, dna_class_byPatient, sep = "_"),
    dna_class = factor(dna_class, levels = unique(str_sort(dna_class, numeric = T)))
  )
table(cl.df$patient_id, cl.df$dna_class_byPatient)

# Make clone-level profiles
clones.ids <- split(cl.df$dna_library_id, cl.df$dna_class)
cn.clones <- lapply(clones.ids, function(i) cn[,i])
cn.clones.profiles <- lapply(cn.clones, function(x) {
  if(is.matrix(x)) Rfast::rowMedians(x) else x # If single cell use vector
})
cn.clones.mtx <- do.call(cbind, cn.clones.profiles)
cn.clones.mtx.smooth <- smooth_bins_mean(cn.clones.mtx, 50)[["mat"]] # x/n bin smoothing by mean
cn.clones.hclust <- hclust(dist(t(cn.clones.mtx.smooth), method="manhattan"))

# Divergence from median profile (per bin)
# Fixed cutoff (fraction of bins divergent from median clone profile: 
diff.cutoff = snakemake@params[["max_clone_diff_frac"]]
diff.clone.cutoff = snakemake@params[["max_clone_excl_frac"]] #solrun should this be max_clone_excl_frac?
tree.dist.cutoff <- snakemake@params[["tree_dist_max"]]

cn.clones.diff <- setNames(lapply(names(cn.clones), function(i){
  if(is.matrix(cn.clones[[i]])) {
    colSums(cn.clones[[i]] != cn.clones.profiles[[i]])
  }
}), nm=names(cn.clones))
diff <- setNames(unlist(cn.clones.diff,use.names=F) / nrow(cn),
                 nm=unlist(lapply(cn.clones.diff,names),use.names=F))
diff.frac.excl <- unlist(lapply(cn.clones.diff, function(x) sum(x/nrow(cn) > diff.cutoff)/length(x)))
# Find bad clones
bad.clones <- names(which(diff.frac.excl > diff.clone.cutoff))

# Tree with all cells. Smooth cn integer matrix by window of bins (w) for efficient tree calculation
cn.smooth.l <- smooth_bins_mean(cn[,good.cells.ids], bins=bins.f, w=50) # with S/G2M but drop bad QC
cn.smooth <- cn.smooth.l$mat
dim(cn.smooth)

# Create distance matrix
dist.mtx <- make_tree(cn.smooth, return_dist=T, cores=8, dist_method = "manhattan")
qdist <- apply(dist.mtx, 1, quantile, 0.01, na.rm=T) # distance to nearest 1%?

# Full class table:
# Move bad cells/clones to undetermined (class 0)
# a) individual cells above divergence or tree cutoff 
# b) whole clones with high excl fraction to class 0 (undetermined)
cl.all <- cl.df %>%
  mutate(tree_dist=qdist[dna_library_id],
         clone_diff_frac=diff[dna_library_id],
         qc_clone_diff=clone_diff_frac < diff.cutoff,
         qc_whole_clone=! dna_class %in% bad.clones,
         qc_tree_dist=tree_dist < tree.dist.cutoff,
         qc_all=qc_clone_diff & qc_whole_clone & qc_tree_dist,
         cell_phase=cells.phase$cell_phase[match(dna_library_id,cells.phase$dna_library_id)], # Add RNA phase for all available
         dna_class_revised=case_when(
           cell_phase %in% c("S","G2M") ~ as.character(dna_class),
           ! qc_all ~ paste(patient_id, "0", sep="_"),
           TRUE ~ as.character(dna_class)),
  ) %>% 
  mutate(
    dna_class=factor(dna_class, levels=unique(str_sort(dna_class, numeric=T))),
    dna_class_revised=factor(dna_class_revised, levels=unique(str_sort(dna_class_revised, numeric=T)))
  )

# Get clone colors
clone_colors <- get_clones_col(c(paste0(snakemake@wildcards[["patient_id"]],"_",0),
                                 levels(cl.all$dna_class))) # Always add class 0
# Get metadata to annotate clones
m <- read_tsv(snakemake@input[["meta_per_cell"]]) %>% filter(dna_library_id %in% colnames(cn))
cl.all.m <- cl.all %>% left_join(m)

### 2. Run multipcf segmentation
# Get normal cell counts for better bin-level normalization
# TODO: Make optional
cn.norm.all <- data.table::fread(snakemake@input[["normal_bincounts"]],data.table=F)[good.bins,]
normals <- read_tsv(snakemake@params[["normals"]], col_names = c("id","sex","read_length")) %>% filter(id %in% colnames(cn.norm.all)) # Filter for correct read length = matching ids in cn matrix of normals
autochr <- bins.f$chr %in% paste0("chr",1:22)
xchr <- bins.f$chr %in% "chrX"
ychr <- bins.f$chr %in% "chrY"
norm.fact.x <- c(rep(2, sum(autochr)+sum(xchr)), rep(NA,sum(ychr)))
if(snakemake@params[["exclude_y"]]){
  norm.fact.y <- c(rep(2, sum(autochr)), rep(1,sum(xchr)), rep(NA,sum(ychr)))
} else {
  norm.fact.y <- c(rep(2, sum(autochr)), rep(1,sum(xchr) + sum(ychr)))
}
cn.norm.x <- rowSums(cn.norm.all[,normals$id[which(normals$sex=="female")]])
cn.norm.y <- rowSums(cn.norm.all[,normals$id[which(normals$sex=="male")]])
norm2.x <- cn.norm.x/norm.fact.x
norm2.y <- cn.norm.y/norm.fact.y
norm3.x <- norm2.x/mean(norm2.x,na.rm=T)
norm3.y <- norm2.y/mean(norm2.y,na.rm=T)
norm.comb <- apply(cbind(norm3.x,norm3.y),1,mean,na.rm=T)
norm.factors <- norm.comb/mean(norm.comb,na.rm=T)
if(snakemake@params[["exclude_y"]]) norm.factors[ychr] <- 1 # No Y chr normalization
# Debug, plot norm factors
# plot(1:length(norm.factors), norm.factors,cex=.1, col=factor(bins.f$chr))

# Get chr arm info
cyto <- read.table(snakemake@params[["cytoband_file"]],col.names=c("chr","start","end","cytoband","chromatin"), sep="\t")
arms <- cyto %>% mutate(arm=substr(cytoband,1,1)) %>% 
  filter(chr %in% unique(bins.f$chr)) %>% 
  mutate(chr=factor(chr, levels=paste0("chr",c(1:22,"X","Y")))) %>%  
  group_by(chr, arm) %>% 
  summarize(start=min(start),end=max(end))

# With GRanges
arms.gr <- GRanges(seqnames=arms$chr, ranges=IRanges(arms$start, arms$end), arm=arms$arm)
bins.gr <- GRanges(seqnames=bins.f$chr, ranges=IRanges(bins.f$start, bins.f$end))
gr.intersect <- findOverlaps(bins.gr, arms.gr)
bins.f$arm <- arms.gr$arm[subjectHits(gr.intersect)]

# Debug plot arms
# plot(bins.f$chr_int, bins.f$end, col=factor(bins.f$arm),cex=0.2, pch=21)
# Can be done with new dplyr versions (>1.1.0) but can mess up conda environments
# by <- join_by(chr, within(x$start, x$end, y$start, y$end))
# bin.arms <- inner_join(bins.f, arms, by)

# Prepare data
# rowSum(counts x cells) x clones
clones.all <- split(cl.all$dna_library_id, cl.all$dna_class_revised)
clones.l <- clones.all[!grepl("_S$|_G2M$|_0$",names(clones.all))]

d <- tibble(clone=names(clones.l),
            ids=map(clone,function(x) clones.l[[x]]),
            counts=map(ids, ~counts[good.bins, .x]),
            counts.m=map(counts, ~rowSums(.x)),
            norm=map(counts.m, function(.x) (.x+1)/mean(.x+1)), 
            norm.adj=map(norm, function(.x) .x/norm.factors),
            norm.gc=map(norm.adj,~lowess.gc(gc.f,.x)),
            lrr=map(norm.gc, ~log2(.x)))

lrr.mat <- as.matrix(do.call(cbind, d$lrr))
colnames(lrr.mat) <- d$clone
mpcf.input <- cbind(cbind("chr"=bins.f$chr, "pos"=bins.f$start),lrr.mat)
dim(mpcf.input)
mpcf.input[1:5,1:ncol(mpcf.input)]

cat("Starting multipcf for sample",snakemake@wildcards[["patient_id"]],"...")
res.mpcf <- copynumber::multipcf(as.data.frame(mpcf.input), pos.unit = "bp", arms = bins.f$arm, gamma = snakemake@params[["gamma"]], normalize = FALSE, fast = F, verbose = T, return.est = F)
cat("done.\n")

s <- seq(snakemake@params[["min_scale"]], snakemake@params[["max_scale"]], length.out = 1000) # Reasonable ploidies/scaling factors
d2 <- d %>% 
  mutate(mpcf.segcounts=lapply(res.mpcf[, 6:ncol(res.mpcf)], function(x) 2^rep.int(x, res.mpcf$n.probes)),
         mpcf.scale_tests=map(mpcf.segcounts, ~sapply(s, function(x) manhattan.dist(x, .x, na.rm=TRUE))),
         mpcf.scale_factor=map(mpcf.scale_tests, ~s[which.min(.x)]),
         mpcf.cn=map2(mpcf.segcounts, mpcf.scale_factor, ~round(.x*.y)),
         mpcf.segspruned=pmap(list(mpcf.cn, norm.adj, mpcf.scale_factor), function(cn, n, s) prune.neighbors(cn, n, s, thresh1=snakemake@params[["prune1"]], thresh2=snakemake@params[["prune2"]], exclude_high_cn=8)),
         mpcf.cnpruned=map(mpcf.segspruned, ~rep(.x$cpy, .x$size))
  )

df.mpcf <- lapply(d2$clone, function(clone){
  i <- which(d2$clone==clone)
  tibble(bin=1:nrow(bins.f), bin_unfilt=bins.f$bin_unfilt, bp=bins.f$end_cum, chr=bins.f$chr, bp_chr_start=bins.f$start, bp_chr_end=bins.f$end,
         raw=d2$norm[[i]]*d2$mpcf.scale_factor[[i]], # Normalized
         gc=d2$norm.gc[[i]]*d2$mpcf.scale_factor[[i]], # GC + normal-panel corrected
         seg=d2$mpcf.segcounts[[i]]*d2$mpcf.scale_factor[[i]],
         cn=d2$mpcf.cn[[i]],
         cn.pruned=d2$mpcf.cnpruned[[i]])
})
names(df.mpcf) <- d2$clone

df.chr <- bins.f %>% group_by(chr) %>% summarize(line=max(end_cum),label_pos=round(mean(end_cum))) %>% mutate(chr_short=sub("chr","",chr))

# Segment matrix by clone after multipcf
clone.mtx <- do.call(cbind, d2$mpcf.cnpruned)
all.segs <- res.mpcf[,1:5]
# Loop over unique states and find matches in whole cn matrix
cn.states.all <- apply(unique(clone.mtx), 1, function(v){
  # Note: need to transpose matrix to match vector
  which(colSums(v == t(clone.mtx)) == ncol(clone.mtx))
})
cn.seg.list <- lapply(cn.states.all, function(x){
  x. <- c(0,cumsum( diff(x)-1 ))
  l <- lapply( split(x, x.), range )
  data.frame( matrix(unlist(l), nrow=length(l), byrow=T) )
})
all.segs.pruned <- bind_rows(cn.seg.list) %>% arrange(X1,X2)
colnames(all.segs.pruned) <- c("start","end")

#Output the segments 
write_tsv(all.segs.pruned, snakemake@output[["segments"]])


clone_labels <- lapply(seq_along(d2$clone), function(i){
  paste0("Clone ",d2$clone[[i]], " (n=",length(d2$ids[[i]])," cells; sf=",round(d2$mpcf.scale_factor[[i]],2),")")
}) %>% unlist %>% setNames(., nm=d2$clone)

max.cn.plot <- max(sapply(df.mpcf, function(x) x$cn.pruned))+1
# p.clones <- plot_clone(df.mpcf, labels=clone_labels, df.chr = df.chr, show_detail = T, cn.trunc = max.cn.plot)
# ggsave(snakemake@output[["clone_profiles"]], p.clones, height=length(clone_labels)*2)

## New grid with multipcf segments: identify identical clones
# Full cn bins with new clones: to be used for consensus profiles of S/G2M/0 class clones that cannot go in the multipcf segmentation
clones.ids.rev <- split(cl.all$dna_library_id, cl.all$dna_class_revised)
cn.clones.rev <- lapply(clones.ids.rev, function(i) cn[,i])
cn.clones.profiles.rev <- lapply(cn.clones.rev, function(x) {
  if(is.matrix(x)) Rfast::rowMedians(x) else x # If single cell use vector
})
cn.clones.mtx.rev <- do.call(cbind, cn.clones.profiles.rev)

# Pruned segs to cn grid
clones.mtx.full <- segs_to_mtx(all.segs.pruned,do.call(cbind, d2$mpcf.cnpruned))
# segs.idx <- rowSums(clones.mtx.full==2)!=ncol(clones.mtx.full) # All non diploid
segs.idx <- colSums(apply(clones.mtx.full, 1, function(x) x[[1]]==x))!=ncol(clones.mtx.full) # All diverging segments
clones.mtx <- clones.mtx.full[segs.idx,]
sm.clones <- grep("_S$|_G2M$|_0$",colnames(cn.clones.mtx.rev),value = T)
clones.mtx.sg2m <- segs_to_mtx(all.segs.pruned[segs.idx,], round(cn.clones.mtx.rev[,sm.clones]))
clones.mtx.all <- cbind(clones.mtx, clones.mtx.sg2m)
if(length(sm.clones)==1) colnames(clones.mtx.all)[ncol(clones.mtx.all)] <- sm.clones
clones.mtx.hclust <- hclust(dist(t(clones.mtx.all),method="manhattan"))

### Plot consensus grid profiles
seg.grid.labs <- all.segs.pruned[segs.idx,] %>% 
  mutate(start_bp=bins.f$start[start],
         start_bp_cum=bins.f$end_cum[start],
         end_bp=bins.f$end[end],
         end_bp_cum=bins.f$end_cum[end],
         chr=bins.f$chr[start],
         size=end_bp-start_bp,
         label=paste0(chr,":",round(start_bp/1e6,1),"-",round(end_bp/1e6,1)),
         bin_idx=1:nrow(.))
clones.df <- as.data.frame(clones.mtx.all) %>% mutate(bin_idx=1:nrow(.)) %>% pivot_longer(-bin_idx, names_to="clone", values_to="cn") %>%
  mutate(clone=factor(clone, levels=clones.mtx.hclust$labels[clones.mtx.hclust$order]),
         excluded=grepl("*_S$|*_G2M$|*_0$",clone))
p.grid <- clones.df %>%
  ggplot(aes(x=as.factor(bin_idx), y=clone, fill=as.factor(cn))) + geom_tile() +
  geom_point(data=tibble(clone=levels(clones.df$clone), excluded=grepl("*_S$|*_G2M$|*_0$",clone)), aes(x=-Inf, y=clone, color=clone), size=5, pch=15, show.legend = F,inherit.aes = F) +
  scale_color_manual(values=clone_colors) +
  scale_fill_manual(values=cn.colors[as.character(c(0:max(clones.df$cn)))]) + 
  theme_dntr(axis_ticks = T, axis_labels = F) +
  scale_x_discrete(expand=c(0.01,0.01), labels=seg.grid.labs$label) +
  facet_grid(excluded~.,scales="free", space="free_y") +
  theme(strip.text=element_blank(), legend.position=c(0,0), 
        axis.text.y=element_text(face="bold"), axis.text.x=element_text(angle = 45, hjust=1, vjust=1, size=8)) +
  guides(fill=guide_colorsteps(title.position="top", title=paste0("copy number [0-",max(clones.df$cn),"]"), 
                               title.hjust=0, title.theme=element_text(size=8), frame.colour = "black", 
                               label=F, barwidth=5,barheight=0.7, show.limits = T)) +
  labs(subtitle="Profile per common segment. All clones")
if(any(!is.na(cl.all.m$timepoint))){
  times <- unique(cl.all.m$timepoint)
  p.hist <- cl.all.m %>% 
    filter(dna_class_revised %in% clones.df$clone) %>% 
    group_by(dna_class_revised, timepoint) %>% summarize(n=n()) %>% 
    mutate(clone=factor(dna_class_revised, levels=clones.mtx.hclust$labels[clones.mtx.hclust$order]),
           excluded=grepl("*_S$|*_G2M$|*_0$",clone)) %>% 
    ggplot(aes(x=clone, y=n, fill=timepoint)) + geom_bar(stat="identity", show.legend = T, position="stack") +
    scale_fill_manual(values=timepoint_cols[times]) + guides(fill=guide_legend(title = NULL))
} else {
  p.hist <- cl.all %>% 
    group_by(dna_class_revised) %>% summarize(n=n()) %>% 
    mutate(clone=factor(dna_class_revised, levels=clones.mtx.hclust$labels[clones.mtx.hclust$order]),
           excluded=grepl("*_S$|*_G2M$|*_0$",clone)) %>% 
    ggplot(aes(x=clone, y=n)) + geom_bar(stat="identity", show.legend = F)  
}
p.hist <- p.hist +
  scale_y_continuous(expand=c(0.01,0.01), breaks=pretty_breaks(3)) +
  facet_grid(excluded~.,scales="free", space="free_y") + coord_flip() +
  theme_dntr(axis_labels = F, axis_ticks=T) +
  theme(strip.text=element_blank(), axis.text.y=element_blank(), plot.margin = unit(c(0,0,0,0),"pt"))
p.consensus_profiles <- p.grid + p.hist + plot_layout(widths=c(6,1), guides="collect") & theme(legend.position="bottom")

# Segment plots
clone_labels.all = setNames(paste0(names(clones.ids.rev)," n=",sapply(clones.ids.rev,length)), nm=names(clones.ids.rev))
segs.chr <- unique(bins.f$chr[bins.f$bin %in% all.segs.pruned[segs.idx,"start"]])
segs.cuts <- all.segs.pruned[segs.idx,] %>% mutate(chr=bins.f$chr[start])
plist.chr <- lapply(segs.chr, function(chr){
  cuts <- seg.grid.labs[seg.grid.labs$chr %in% chr,]
  p <- plot_clone(df.mpcf, labels=clone_labels.all,
                  df.chr = df.chr,
                  region=chr,
                  show_detail = T,
                  cn.trunc = max.cn.plot) 
  p + geom_vline(xintercept=unique(c(cuts$start_bp_cum,cuts$end_bp_cum)),linetype="solid",color="red", alpha=0.5)
})

pdf(snakemake@output[["consensus_clone_profiles"]])
p.consensus_profiles
plist.chr
dev.off()

### Merge identical clones
dup_clones.idx <- duplicated(t(clones.mtx)) | duplicated(t(clones.mtx), fromLast = TRUE)
dup_clones <- names(which(dup_clones.idx))
dup_clones.df <- data.frame(old=colnames(clones.mtx),
                            merged=colnames(clones.mtx))
dup_clones.df$merged[dup_clones.idx] <- str_sort(dup_clones, numeric = T)[1]
clones.mtx.dedup <- clones.mtx[,unique(dup_clones.df$merged)]
cl.all.merged <- cl.all %>% 
  mutate(dna_class_merge=case_when(
    grepl("_S$|_G2M$|_0$",dna_class_revised) ~ dna_class_revised,
    TRUE ~ dup_clones.df$merged[match(dna_class_revised, dup_clones.df$old)]),
    dna_class_merge=factor(dna_class_merge, levels=str_sort(unique(dna_class_merge),numeric=T))
  )
cl.all.merged.m <- cl.all.merged %>% left_join(m)
# Debug, tabulation
# table(cl.all.merged$dna_class_revised, cl.all.merged$dna_class_merge)

# Save final clone set
write_tsv(cl.all.merged, snakemake@output[["clones"]])

# TODO: Recreate pseudobulk profiles with merged profiles
# Start from d (new lists + raw counts -> all.segs.pruned)
# Pass through pruning again optionally

### Plots: Trees
# Make unfiltered and filtered trees
keep.ids <- filter(cl.all.merged, qc_all==T) %>% pull(dna_library_id)
filtered.ids <- setdiff(cl.all.merged$dna_library_id, keep.ids)
tree.ids <- list(all=colnames(cn.smooth),
                 filt=keep.ids,
                 filt_norep=intersect(good.cells2.ids,keep.ids))

# Use smooth cn matrix for tree creation
trees <- lapply(names(tree.ids), function(i){
  idx <- tree.ids[[i]]
  make_tree(cn.smooth[,idx], cores=snakemake@threads, dist_method="manhattan", tree_method=snakemake@params[["tree_method"]])
})

# Color by clones
clones.list <- lapply(tree.ids, function(i){
  cl <- cl.all.merged %>% filter(dna_library_id %in% i)
  split(cl$dna_library_id, cl$dna_class)
})
trees.clones <- lapply(seq_along(trees), function(i){
  tree <- trees[[i]]
  tree <- groupOTU(tree, clones.list[[i]], group_name = "clone")
  clones <- unique(names(clones.list[[i]]))
  p <- ggtree(tree, aes(color=droplevels(clone)), ladderize = F) +
    layout_dendrogram() +
    scale_color_manual(values = clone_colors[clones], name = names(clone_colors[clones])) +
    geom_treescale(x = 10) +
    guides(color=guide_legend(override.aes=list(size=3), nrow=2, title.position = "top", title = "Clones"),
           shape=guide_legend(title.position="top", title="Excluded")) +
    theme(legend.position="bottom")
  if(i==1) p <- p + geom_tippoint(aes(shape=label %in% filtered.ids), color="red") + scale_shape_manual(values=c("TRUE"=21))
  return(p)
})

if(any(!is.na(m$timepoint))){
  timepoints.list <- lapply(tree.ids, function(i){
    cl <- cl.all.merged.m %>% filter(dna_library_id %in% i)
    split(cl$dna_library_id, cl$timepoint)
  })
  # timepoint_cols # from general.R
  trees.time <- lapply(seq_along(trees), function(i){
    tree <- trees[[i]]
    tree <- groupOTU(tree, timepoints.list[[i]], group_name = "timepoint")
    times <- unique(names(timepoints.list[[i]]))
    p.tree <- ggtree(tree, aes(color=droplevels(timepoint)), ladderize = F) +
      layout_dendrogram() +
      scale_color_manual(values = timepoint_cols[times]) +
      geom_treescale(fontsize=4, linesize=1, offset=20) +
      guides(color=guide_legend(override.aes=list(size=3), nrow=3, title.position = "top", title = "Timepoints"),
             shape=guide_legend(title.position="top", title="Excluded")) +
      theme(legend.position="bottom")
    if(i==1) p.tree <- p.tree + geom_tippoint(aes(shape=label %in% filtered.ids), color="red") + scale_shape_manual(values=c("TRUE"=21))
    return(p.tree)
  })
}

# Select 1+ clones/timepoints
# sel <- c("_S","_G2M","_0",NA)
# sel.idx <- grepl(paste(sel,collapse="$|"),names(clone_colors))
# highlight_colors <- clone_colors[sel.idx]
# 
# trees.clones[[1]] + scale_color_manual(values = highlight_colors, name = names(highlight_colors), na.value="grey")
# trees.clones[[2]] + scale_color_manual(values = highlight_colors, name = names(highlight_colors), na.value="grey")

### Plots: QC measures per clone
cn.clones.df <- as.data.frame(round(cn.clones.mtx.smooth)) %>% mutate(bin=1:nrow(.)) %>% pivot_longer(-bin, names_to="clone", values_to="cn") %>%
  mutate(clone=factor(clone, levels=cn.clones.hclust$labels[cn.clones.hclust$order]),
         excluded=grepl("*_S$|*_G2M$|*_0$",clone))
p.grid <- cn.clones.df %>%
  ggplot(aes(x=bin, y=clone, fill=as.factor(cn))) + geom_tile() +
  geom_point(data=tibble(clone=levels(cn.clones.df$clone), excluded=grepl("*_S$|*_G2M$|*_0$",clone)), aes(x=-Inf, y=clone, color=clone), size=5, pch=15, show.legend = F,inherit.aes = F) +
  scale_color_manual(values=clone_colors) +
  scale_fill_manual(values=cn.colors[as.character(c(0:max(cn.clones.df$cn)))]) + 
  theme_dntr(axis_ticks = T, axis_labels = F) +
  scale_x_continuous(expand=c(0.01,0.01)) +
  facet_grid(excluded~.,scales="free", space="free_y") +
  theme(strip.text=element_blank(), legend.position=c(0,0), axis.text.y=element_text(face="bold")) +
  guides(fill=guide_colorsteps(title.position="top", title=paste0("copy number [0-",max(cn.clones.df$cn),"]"), title.hjust=0, title.theme=element_text(size=8), frame.colour = "black", 
                               label=F, barwidth=5,barheight=0.7, show.limits = T)) +
  labs(subtitle="Median profile per smoothed bin. All clones")
if(any(!is.na(cl.all.merged.m$timepoint))){
  times <- unique(cl.all.merged.m$timepoint)
  p.hist <- cl.all.merged.m %>% 
    group_by(dna_class, timepoint) %>% summarize(n=n()) %>% 
    mutate(dna_class=factor(dna_class, levels=cn.clones.hclust$labels[cn.clones.hclust$order]),
           excluded=grepl("*_S$|*_G2M$|*_0$",dna_class)) %>% 
    ggplot(aes(x=dna_class, y=n, fill=timepoint)) + geom_bar(stat="identity", show.legend = T, position="stack") +
    scale_fill_manual(values=timepoint_cols[times]) + guides(fill=guide_legend(title = NULL))
} else {
  p.hist <- cl.all.merged %>% 
    group_by(dna_class) %>% summarize(n=n()) %>% 
    mutate(dna_class=factor(dna_class, levels=cn.clones.hclust$labels[cn.clones.hclust$order]),
           excluded=grepl("*_S$|*_G2M$|*_0$",dna_class)) %>% 
    ggplot(aes(x=dna_class, y=n)) + geom_bar(stat="identity", show.legend = F)  
}
p.hist <- p.hist +
  scale_y_continuous(expand=c(0.01,0.01)) +
  facet_grid(excluded~.,scales="free", space="free_y") + coord_flip() +
  theme_dntr(axis_labels = F, axis_ticks=T) +
  theme(strip.text=element_blank(), axis.text.y=element_blank(), plot.margin = unit(c(0,0,0,0),"pt"))
p.avg_profiles <- p.grid + p.hist + plot_layout(widths=c(6,1), guides="collect") & theme(legend.position="bottom")

# Save plots. Plot distributions and bad clones
pdf(snakemake@output[["qc_clones"]],width=12,height=10)
layout(matrix(c(1,1,3,3,2,2,4,4,5,5,6,7), nrow=4, ncol=3,byrow=F))
boxplot(split(cl.all.merged$clone_diff_frac, cl.all.merged$dna_class),las=2,cex=.2, main="Fraction of divergent bins per cell (from clone median)", xlab="clone")
hist(cl.all.merged$clone_diff_frac,breaks=100,main=paste0("Divergent bins per cell (from clone median; n=",nrow(cl.all.merged),")"))
abline(v=diff.cutoff,col="red")
text(x=diff.cutoff, y=50, label=paste0("n=",sum(!cl.all.merged$qc_clone_diff & cl.all.merged$dna_class_byPatient!=0, na.rm=T)),col="red", adj=-0.1)
boxplot(split(log2(cl.all.merged$tree_dist), cl.all.merged$dna_class), las=2,cex=.2, main="Distance to nearest 1% of single-cell tree (log2)")
hist(log2(cl.all.merged$tree_dist),breaks=100,main="Distance to nearest 1% of single-cell tree (log2)")
abline(v=log2(tree.dist.cutoff),col="red")
barplot(diff.frac.excl,las=2,main="Fraction of clone excluded (bin divergence)",col=ifelse(grepl("*_S$|*_G2M$|*_0$", names(diff.frac.excl)), "red", "grey"))
abline(h=diff.clone.cutoff,col="red")
barplot(table(cl.all.merged$dna_class),las=2, main="cells per clone (original)")
barplot(table(cl.all.merged$dna_class_revised),las=2, main="cells per clone (revised)")

p1 <- ggplot(cl.all.merged, aes(log2(tree_dist), clone_diff_frac, color=cell_phase)) + geom_point(pch=19, size=0.8, alpha=1) + 
  geom_point(data=cl.all.merged[!cl.all.merged$qc_all & !cl.all.merged$dna_class_byPatient %in% c("0","S","G2M"),], color="red", size=1.5, pch=21, alpha=0.3) + 
  facet_wrap(~dna_class) + geom_vline(xintercept=log2(tree.dist.cutoff),color="grey",linetype="dotted") + geom_hline(yintercept=diff.cutoff,color="grey",linetype="dotted") +
  labs(x="L1 dist to nearest 1% of cells in tree, log2", y="Diff from clone (frac of bins)") +
  scale_color_manual(values=c(phase_cols,"NA"="black")) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  theme_dntr(axis_ticks = T, legend="bottom") +
  guides(color=guide_legend(override.aes=list(size=3,alpha=1), title=NULL))
cl.crossfreq <- table("merged"=cl.all.merged$dna_class_merge, "old"=cl.all.merged$dna_class)
# Scaled by column (eg fraction of old classes)
p2 <- ggplot(as.data.frame(scale(cl.crossfreq)), aes(old, merged, fill=Freq)) + geom_tile(show.legend=F) + scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + labs(x="Original clone", y="Merged/final clone")

(p1 + p2 + p.avg_profiles + plot_layout(widths=c(1,1,2), guides="collect") & theme(legend.position="bottom")) / trees.clones[[1]] + plot_layout(heights = c(1,2))


trees.clones[[1]] / trees.clones[[2]] / trees.clones[[3]] + 
  plot_annotation(tag_levels="a", 
                  subtitle=paste0("Trees a) raw (n=",length(unique(trees.clones[[1]]$data$label)),
                                  ") b) excl outliers (n=",length(unique(trees.clones[[2]]$data$label)),
                                  ") c) excl outliers + replicating (n=",length(unique(trees.clones[[3]]$data$label)),")")) + 
  plot_layout(guides="collect") & theme(legend.position="bottom")

if(any(!is.na(m$timepoint))){
  trees.time[[1]] / trees.time[[2]] / trees.time[[3]] + plot_layout(guides="collect") & theme(legend.position="bottom")
}

dev.off()

### Final cn matrix by clone
clones.ids.final <- split(cl.all.merged$dna_library_id, cl.all.merged$dna_class_merge)
cn.clones.final <- lapply(clones.ids.final, function(i) cn[,i])
cn.clones.profiles.final <- lapply(cn.clones.final, function(x) {
  if(is.matrix(x)) Rfast::rowMedians(x) else x # If single cell use vector
})
cn.clones.mtx.final <- do.call(cbind, cn.clones.profiles.final)

### Plots: Heatmaps
cn.smooth10 <- smooth_bins_mean(cn[,good.cells.ids], bins=bins.f, w=snakemake@params[["heatmap_smooth_factor"]]) # with S/G2M but drop bad QC
mtx <- t(cn.smooth10$mat)
bins.m <- cn.smooth10$bins

cn.distL1 <- amap::Dist(mtx, method="manhattan", nbproc = 8) # Note! Transposed
cn.hclust <- hclust(cn.distL1, method="complete")
hc_order <- cn.hclust$labels[cn.hclust$order]

# Metadata for plotting
qc <- read_tsv(snakemake@input[["dna_qc"]]) %>% filter(dna_library_id %in% colnames(cn))
meta.us <- cl.all.merged.m %>% 
  left_join(qc,by="dna_library_id") %>% 
  mutate(dna_class_byPatient=factor(dna_class_byPatient, levels=c(str_sort(unique(dna_class_byPatient),numeric=T)[-1],0)),
         clone_revised_short=sub(".*_(.*)$","\\1",dna_class_revised),
         clone_revised_short=factor(clone_revised_short,levels=c(str_sort(unique(clone_revised_short),numeric=T)[-1],0)),
         clone_merge_short=sub(".*_(.*)$","\\1",dna_class_merge),
         clone_merge_short=factor(clone_merge_short,levels=c(str_sort(unique(clone_merge_short),numeric=T)[-1],0))
  )
chr_labels <- factor(bins.m$chr_short,levels=c(1:22,"X","Y"))

# Addon: reorder clones by average "damage" as calc by absolute distance from 2 copies (diploid)
mean_clones_dmg <- sort(colSums(abs(cn.clones.mtx.final - 2)))
mean_clones_clust.all <- sub(".*_(.*)$","\\1", names(mean_clones_dmg))
mean_clones_clust <- c(mean_clones_clust.all[! mean_clones_clust.all %in% c("G2M","S","0") & mean_clones_clust.all %in% levels(meta.us$clone_merge_short)],"S","G2M","0")
meta.us$clone_merge_short <- factor(meta.us$clone_merge_short, levels=mean_clones_clust)

meta <- meta.us[match(rownames(mtx), meta.us$dna_library_id),]

#Check correct order
all(rownames(mtx)==meta$dna_library_id)

# Colors and legends
col_fun = circlize::colorRamp2(c(0, 2, 4, 5, 20, 50), c("blue","white","red", "brown", "green", "yellow"))
col_fun_dups = circlize::colorRamp2(c(0,0.1,0.4), c("white","white","black"))
col_fun_prob = circlize::colorRamp2(c(0,1), c("white","darkgreen"))
rna_qc_colors <- setNames(c("black","white"),nm=c("PASS","FAIL"))

lgd_list = list(
  Legend(title="copy number", col_fun = col_fun, break_dist = 1, title_position = "lefttop-rot")
  # Legend(title="clone", labels=names(clone_colors), legend_gp=gpar(fill=clone_colors), title_position = "lefttop-rot")
  # Legend(title="rna qc", labels=names(rna_qc_colors), legend_gp=gpar(fill=rna_qc_colors), title_position = "lefttop-rot")
  # Legend(title="dups", col_fun = col_fun_dups, title_position = "lefttop-rot")
)

# Annotations
ha <- rowAnnotation(clone=meta$dna_class,
                    clone_rev=meta$dna_class_revised,
                    clone_merge=meta$dna_class_merge,
                    time=meta$timepoint,
                    cloneDiff=meta$clone_diff_frac,
                    rna_qc=ifelse(is.na(meta$cell_phase),"FAIL","PASS"),
                    dups=meta$bam_dup_frac,
                    classProb=meta$dna_class_prob,
                    # scp=anno_barplot(meta$scp,bar_width = 1, border = F, baseline="min"),
                    reads=anno_barplot(log10(meta$bam_read_pairs), bar_width=1, border=F, baseline = "min"),
                    col=list(clone=clone_colors, clone_rev=clone_colors, clone_merge=clone_colors, time=timepoint_cols,
                             rna_qc=rna_qc_colors, dups=col_fun_dups, classProb=col_fun_prob),
                    show_legend=T, annotation_name_gp = gpar(fontsize=8))

png(snakemake@output[["heatmap"]], width=2000,height=1400,units="px",res=150)
ht = Heatmap(mtx, name="cn", col=col_fun,
             heatmap_legend_param=list(title="cn",col_fun=col_fun,break_dist=1),
             cluster_columns = F,
             row_split = meta$clone_merge_short,
             
             # Opt 1: outside clustering
             row_order = hc_order, # Named vector matching rownames in mtx
             
             # Opt 2: inside clustering (+/- reordering of non-visible dendrogram)
             # cluster_row_slices = F, # Keep slice order = factor order of clones
             # cluster_rows= T, # Cluster within each clone, to L1 dist/complete
             # clustering_distance_rows = "manhattan",
             # clustering_method_rows="complete",
             # row_dend_reorder = F,
             # show_row_dend = F,
             
             column_split = chr_labels,
             row_gap = unit(2, "mm"),
             column_gap = unit(0, "mm"),
             border=T,
             row_title_gp=gpar(fontsize=8),
             column_title_gp=gpar(fontsize=8),
             column_title_side = "bottom",
             show_row_names = F,
             show_column_names = F,
             right_annotation = ha,
             use_raster = F,
             show_heatmap_legend=F)
ht <- draw(ht, annotation_legend_list=packLegend(list=lgd_list))
dev.off()

png(snakemake@output[["heatmap_plate"]], width=2000,height=1400,units="px",res=150)
ht.plate = Heatmap(mtx, name="cn", col=col_fun,
                   heatmap_legend_param=list(title="cn",col_fun=col_fun,break_dist=1),
                   cluster_columns = F,
                   row_split = meta$plate_id,
                   
                   # Opt 1: Set ordering from outside clustering
                   row_order = hc_order, # Named vector matching rownames in mtx
                   # cluster_row_slices = T, # Cluster row slices by mean profiles?
                   
                   # # Opt 2: inside clustering (+/- reordering of non-visible dendrogram)
                   # cluster_rows= T, # Cluster within each clone, to L1 dist/complete
                   # clustering_distance_rows = "manhattan",
                   # clustering_method_rows="complete",
                   # row_dend_reorder = F,
                   # show_row_dend = F,
                   
                   column_split = chr_labels,
                   row_gap = unit(2, "mm"),
                   column_gap = unit(0, "mm"),
                   border=T,
                   row_title_gp=gpar(fontsize=8),
                   column_title_gp=gpar(fontsize=8),
                   column_title_side = "bottom",
                   show_row_names = F,
                   show_column_names = F,
                   right_annotation = ha,
                   use_raster = F,
                   show_heatmap_legend=F)
ht.plate <- draw(ht.plate, annotation_legend_list=packLegend(list=lgd_list))
dev.off()

# Plot unfiltered heatmap (split by raw clone + filtered out, showing new rev clone too)
# All cells from unfiltered cn matrix
cn.all.smooth.l <- smooth_bins_mean(cn, bins=bins.f, w=snakemake@params[["heatmap_smooth_factor"]])
mtx.raw <- t(cn.all.smooth.l$mat)

# All NA cells can not get distance computed
excl.idx <- rowSums(is.na(mtx.raw)) == ncol(mtx.raw)
mtx.dist <- amap::Dist(mtx.raw[!excl.idx,], method="manhattan",nbproc=8) # Note! Transposed
mtx.hc <- hclust(mtx.dist, method="complete")

# Get order by cell ids and add NA cells at the end of vector
hc_order.raw <- mtx.hc$labels[mtx.hc$order]
hc_order.all <- c(hc_order.raw, setdiff(rownames(mtx.raw),hc_order.raw))

mtx.bins <- cn.all.smooth.l$bins
mtx.chr_labels <- factor(mtx.bins$chr_short,levels=c(1:22,"X","Y"))
mtx.meta.us <- m %>% left_join(qc, by="dna_library_id") %>% 
  mutate(
    clone=as.character(cl.all.merged$dna_class)[match(dna_library_id,cl.all.merged$dna_library_id)],
    clone_rev=as.character(cl.all.merged$dna_class_revised)[match(dna_library_id,cl.all.merged$dna_library_id)],
    clone_rev=ifelse(grepl("*_0$",clone_rev),NA,clone_rev),
    clone_merge=as.character(cl.all.merged$dna_class_merge)[match(dna_library_id,cl.all.merged$dna_library_id)],
    clone_merge=ifelse(grepl("*_0$",clone_rev),NA,clone_merge),
    dna_qc=ifelse(dna_library_id %in% good.cells.ids, "PASS", "FAIL"),
    clone_qc=ifelse(dna_library_id %in% names(which(cl.all.merged$qc_all)), "PASS", "FAIL"),
    cell_phase=cells.phase$cell_phase[match(dna_library_id,cells.phase$dna_library_id)]) # Add RNA phase for all available)
mtx.meta <- mtx.meta.us[match(rownames(mtx.raw), mtx.meta.us$dna_library_id),] # Order by cn matrix to plot

ha2 <- rowAnnotation(dna_qc=mtx.meta$dna_qc,
                     clone_qc=mtx.meta$clone_qc,
                     clone=mtx.meta$clone,
                     clone_rev=mtx.meta$clone_rev,
                     clone_merge=mtx.meta$clone_merge,
                     time=mtx.meta$timepoint,
                     # cloneDiff=mtx.meta$clone_diff_frac,
                     rna_qc=ifelse(is.na(mtx.meta$cell_phase),"FAIL","PASS"),
                     dups=mtx.meta$bam_dup_frac,
                     # scp=anno_barplot(mtx.meta$scp,bar_width = 1, border = F, baseline=0),
                     reads=anno_barplot(log10(mtx.meta$bam_read_pairs+1), bar_width=1, border=F, baseline = "min"),
                     col=list(rna_qc=rna_qc_colors, dna_qc=rna_qc_colors, clone_qc=rna_qc_colors,
                              dups=col_fun_dups, clone=clone_colors, clone_rev=clone_colors, clone_merge=clone_colors,
                              time=timepoint_cols),
                     show_legend=T, annotation_name_gp = gpar(fontsize=8))

png(snakemake@output[["heatmap_raw"]], width=2000,height=1400,units="px",res=150)
ht.raw = Heatmap(mtx.raw, name="cn", col=col_fun,
                 heatmap_legend_param=list(title="cn",col_fun=col_fun,break_dist=1),
                 cluster_columns = F,
                 row_split = mtx.meta$dna_qc,
                 row_order = hc_order.all, # Named vector matching rownames in mtx.raw
                 column_split = mtx.chr_labels,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(0, "mm"),
                 border=T,
                 row_title_gp=gpar(fontsize=8),
                 column_title_gp=gpar(fontsize=8),
                 column_title_side = "bottom",
                 show_row_names = F,
                 show_column_names = F,
                 right_annotation = ha2,
                 use_raster = F,
                 show_heatmap_legend=F)
ht.raw <- draw(ht.raw, annotation_legend_list=packLegend(list=lgd_list))
dev.off()


to_save <- c("bins", "good.bins", "gc", 
             "counts", "cn",
             "good.cells", "good.cells.ids", "good.cells2", "good.cells2.ids", 
             "cells.phase", "cells.replicating.pass", 
             "qc", "meta", "cn.umap", "umap", "cl.all",
             "all.segs.pruned","d2","clone.mtx","df.mpcf")
cat("^Saving Rda\n"); save(list = to_save, file=paste0(snakemake@output[["rda"]]))
