library(tidyverse)
library(GenomicRanges)
library(copynumber)
library(furrr)
options(scipen = 999)
try(RhpcBLASctl::blas_set_num_threads(snakemake@threads), silent = TRUE)
# Functions
dna_to_cellid <- function(x){ sub("([A-Z]{3}[0-9]{5})(D?)_([A-Z][0-9]{2}).*","\\1.\\3", x) }
manhattan.dist <- function(factor, segs) sum(abs((segs*factor)-round(segs*factor)))
safe_future_map <- function(...) {
  RhpcBLASctl::blas_set_num_threads(1)  # Set to 1 for each worker
  future_map(...)
}
safe_future_map2 <- function(...) {
  RhpcBLASctl::blas_set_num_threads(1)  # Set to 1 for each worker
  future_map2(...)
}
gc.map.correct <- function(bin.counts, gc, map, gc.coefs) {
  predq <- function(qest, gc, scale=1) {
    (qest[1]+gc*qest[2]+gc^2*qest[3])/scale
  }
  gc.scale <- predq(gc.coefs, 0.45, 1)
  weight <- map*predq(gc.coefs, gc, gc.scale)
  bin.counts/weight
}
lowess.gc = function(x, y) {
  low = lowess(x, log(y), f=0.05, iter=10);
  z = approx(low$x, low$y, x)
  return(exp(log(y) - z$y))
}

pid <- snakemake@wildcards[["patient_id"]]

# Parameters
counts_min <- snakemake@params[["counts_min"]]
binsize <- as.integer(snakemake@wildcards[["binsize"]]) 
cell.ids <- snakemake@params[["cells"]]

# Load bins
bins <- read_tsv(snakemake@input[["goodbins"]], skip=1, col_names = c("chr","start","end","idx","bin"), show_col_types = F) %>% 
  mutate(
    chr_short=sub("chr","",chr),
    chr_int=as.integer(case_when(chr_short == "X" ~ "23", chr_short == "Y" ~ "24", TRUE ~ chr_short))
  )
bins$chr <- factor(bins$chr, levels = unique(bins$chr))

# Get chr arm info
cyto <- read.table(snakemake@params[["cytoband_file"]],col.names=c("chr","start","end","cytoband","chromatin"), sep="\t")
arms <- cyto %>% mutate(arm=substr(cytoband,1,1)) %>% 
  filter(chr %in% unique(bins$chr)) %>% 
  mutate(chr=factor(chr, levels=paste0("chr",c(1:22,"X","Y")))) %>%  
  group_by(chr, arm) %>% 
  summarize(start=min(start),end=max(end))

# With GRanges
arms.gr <- GRanges(seqnames=arms$chr, ranges=IRanges(arms$start, arms$end), arm=arms$arm)
bins.gr <- GRanges(seqnames=bins$chr, ranges=IRanges(bins$start, bins$end))
gr.intersect <- findOverlaps(bins.gr, arms.gr)
bins$arm <- arms.gr$arm[subjectHits(gr.intersect)]

# GC normalization of a) normal panel of cells and b) cells of interest in patient
# Load gc and mapability tracks, filtered for good bins
gc <- as.numeric(readLines(snakemake@input[["gc"]]))[bins$idx] 
map <- as.numeric(readLines(snakemake@input[["map"]]))[bins$idx]

# All cells from patient and QC
cells.all <- colnames(data.table::fread(snakemake@input[["counts"]], header = T, nrows = 1))
cells.idx <- which(cells.all %in% cell.ids)
counts <- as.matrix(data.table::fread(snakemake@input[["counts"]], select = cells.idx))

cat(pid,"n=",ncol(counts),"with",nrow(counts),"bins (before filtering)\n")
dim(counts)

# 1. Filter for good cells
bc <- colSums(counts)
good.cells <- bc > counts_min
good.cells.ids <- names(bc)[which(good.cells)]

# 2. Filter out replicating cells (by RNA)
if(!is.null(snakemake@input[["rna_phase"]])){
  cells.phase <- read_tsv(snakemake@input[["rna_phase"]]) %>%
    filter(cell_id %in% dna_to_cellid(colnames(counts))) %>%  # Keep only cells also in the copy number matrix == snakemake@params[["cells"]]
    mutate(dna_library_id=colnames(counts)[match(cell_id,dna_to_cellid(colnames(counts)))])
  cells.replicating <- cells.phase$cell_id[cells.phase$cell_phase %in% c("S", "G2M", "G2/M", "M", "G2")] 
} else {
  cells.phase <- tibble(cell_id = NA, cell_phase = NA, dna_library_id = NA)
  cells.replicating <- c()
}
good.cells2 <- good.cells & !dna_to_cellid(names(bc)) %in% cells.replicating
good.cells2.ids <- names(bc)[which(good.cells2)]
counts.f <- counts[bins$idx,good.cells2.ids]

# Scale factors from normal panel: subset normal counts for 1) good bins 2) NOT cells in this particular patient
if(!is.null(snakemake@input[["normal_bincounts"]])){
  n.counts.all <- data.table::fread(snakemake@input[["normal_bincounts"]],data.table=F)
  exclude_self <- !colnames(n.counts.all) %in% colnames(counts.f)
  n.counts <- n.counts.all[bins$idx, exclude_self]
  # n.meta <- read_tsv(snakemake@config[["dna"]][["normal_cells"]],col_names = c("id","sex","read_length")) %>% 
  #   filter(id %in% colnames(n.counts))
  # n.meta <- n.meta[match(colnames(n.counts), n.meta$id),] # Ensure same order
  autochr <- bins$chr %in% paste0("chr",1:22)
  xchr <- bins$chr %in% "chrX"
  ychr <- bins$chr %in% "chrY"
  # Normalize autosomal and X chr separately
  n.auto.counts <- rowSums(n.counts[autochr,])
  n.x.counts <- rowSums(n.counts[xchr,])
  # FT + lowess
  n.auto.ft <- sqrt(n.auto.counts) + sqrt(n.auto.counts+1)
  n.auto.ft.lowess <- lowess.gc(gc[autochr], n.auto.ft)
  n.x.ft <- sqrt(n.x.counts) + sqrt(n.x.counts+1)
  n.x.ft.lowess <- lowess.gc(gc[xchr], n.x.ft)
  # Gc + map 
  n.auto.gc.coefs <- summary(lm(n.auto.counts ~ gc[autochr] + I(gc[autochr]^2)))$coefficients[,1]
  n.auto.gcmap <- gc.map.correct(n.auto.counts, gc[autochr], map[autochr], n.auto.gc.coefs)
  n.auto.gcmap.mc <- n.auto.gcmap/median(n.auto.gcmap)
  n.x.gc.coefs <- summary(lm(n.x.counts ~ gc[xchr] + I(gc[xchr]^2)))$coefficients[,1]
  n.x.gcmap <- gc.map.correct(n.x.counts, gc[xchr], map[xchr], n.x.gc.coefs)
  n.x.gcmap.mc <- n.x.gcmap/median(n.x.gcmap)
  # Collect factors
  n_factors.ft <- c(n.auto.ft.lowess, n.x.ft.lowess, rep(1, sum(ychr))) # Set chrY factor to 1 (=no correction)
  n_factors.gcmap <- c(n.auto.gcmap.mc, n.x.gcmap.mc, rep(1, sum(ychr))) # Set chrY factor to 1 (=no correction)
} else {
  # No scaling
  n_factors.ft <- rep(1, nrow(bins))
  n_factors.gcmap <- rep(1, nrow(bins))
}

# NEW 241028: lowess.gc + norm depth
ft <- apply(counts.f, 2, function(x) sqrt(x) + sqrt(x + 1))
ft.lowess <- apply(ft, 2, function(x) lowess.gc(gc, x))
ft.lowess.nc <- ft.lowess/n_factors.ft

gcmap.mc.nc <- apply(counts.f, 2, function(x){
  gc.coefs <- summary(lm(x ~ gc + I(gc^2)))$coefficients[,1]
  r <- gc.map.correct(x, gc, map, gc.coefs)
  r.mc <- r/mean(r)
  r.mc/n_factors.gcmap
})
#Exclude cells that have bins with negative values after GC normalisation
indices <- which(gcmap.mc.nc < (0), arr.ind = TRUE)
exclude.cells<-colnames(gcmap.mc.nc)[unique(indices[,2])]
cat(paste0("Cells excluded because of GC filtering:", exclude.cells, "...\n"))
gc_norm_f<-gcmap.mc.nc[,!colnames(gcmap.mc.nc)%in%exclude.cells]


mpcf.input.norm <- cbind(cbind("chr"=bins$chr, "pos"=bins$start), ft.lowess.nc)
# gamma <- as.numeric(sub("^0", "0.", snakemake@wildcards[["gamma"]]))
gamma <- as.numeric(snakemake@wildcards[["gamma"]])
mpcf_fast <- ifelse(is.null(snakemake@params[["mpcf_exact"]]), T, !snakemake@params[["mpcf_exact"]])
# Parallelize by chromosome
plan("multicore", workers = snakemake@threads)
chr.idx <- split(seq_len(nrow(mpcf.input.norm)), bins$chr)
input.split <- lapply(chr.idx, function(i) mpcf.input.norm[i, , drop = FALSE])
bins.split <- lapply(chr.idx, function(i) bins[i,])

bins_dropped <- list()

#mpcf does not run when a chr arm only contains one bin - this happens in some of the resolutions, like 40kb and 100kb for example (chr21p arm)
for (i in seq_along(bins.split)) {
  df <- bins.split[[i]]
  arm_counts <- table(df$arm)
  drop_arms <- names(arm_counts[arm_counts == 1])
  
  bins_dropped[[i]] <- df %>% filter(arm %in% drop_arms)
  bins.split[[i]] <- df %>% filter(!(arm %in% drop_arms))
}

input.split <- Map(function(bins_df) {
  mpcf.input.norm[which(bins$chr %in% bins_df$chr & bins$start %in% bins_df$start), , drop = FALSE]
}, bins.split)

cat(paste0("Running joint segmentation in ",pid," at gamma=",gamma," with fast=",mpcf_fast," (n=",sum(good.cells2)," filtered cells)...\n"))
res.split <- safe_future_map2(input.split, bins.split, function(d, b)
  copynumber::multipcf(as.data.frame(d), pos.unit="bp", arms=b$arm, gamma=gamma, normalize=F, fast=mpcf_fast, verbose=F, return.est=F)
)

# Combine back to joint res object
res <- do.call(rbind, res.split)

# Add back dropped bins (if any)
if (length(bins_dropped) > 0) {
  dropped_all <- do.call(rbind, bins_dropped)
  input_dropped <- mpcf.input.norm[with(dropped_all, which(bins$chr %in% chr & bins$start %in% start)), , drop = FALSE]
  
  res_dropped <- data.frame(
    chrom = sub("^chr", "", dropped_all$chr),
    arm       = dropped_all$arm,
    start.pos = dropped_all$start,
    end.pos   = dropped_all$end,
    n.probes  = 1,
    input_dropped[, -c(1:2), drop = FALSE],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  res2 <- rbind(res, res_dropped)
  res2$chrom <- factor(res2$chrom, levels = unique(res2$chrom)) 
  res<-res2[order(res2$chrom, res2$start.pos), ]
}
write_tsv(res, snakemake@output[["mpcf"]])

# Recalculate segment means: get segment means from matrix of gc-normalized counts (instead of using multipcf output)
all.segs <- res[,1:5] %>% 
  mutate(
    start = c(1, head(cumsum(n.probes) + 1, -1)),
    end = cumsum(n.probes),
    size = (end - start) + 1
  )
# Make list of segment indices
seg.means <- map(1:nrow(all.segs), function(i) {
  seg <- all.segs[i,]
  r <- gc_norm_f[seg$start:seg$end, ]
  if(is.null(nrow(r)) || nrow(r)==1) {
    return(r)
  } else { 
    return(colMeans(r))
  }
}) %>% bind_rows %>% as.matrix
segcounts.mtx <- seg.means[rep(1:nrow(seg.means), times = all.segs$n.probes), ]
# segcounts.mtx.n<-apply(segcounts.mtx, 2, function(x) x/mean(x)) # mean centered upstream instead
segcounts <- setNames(split(segcounts.mtx, col(segcounts.mtx)),nm=colnames(segcounts.mtx))

s <- seq(snakemake@params[["min_scale"]], snakemake@params[["max_scale"]], length.out = 100)

scale.tests <- safe_future_map(segcounts, function(segs_vec) sapply(s, function(x) manhattan.dist(x, segs_vec)))
scale.factors <- map(scale.tests, ~s[which.min(.x)])
q1 <- quantile(unlist(scale.factors), 0.25)
q3 <- quantile(unlist(scale.factors), 0.75)
iqr <-q3 - q1
outliers <- which(unlist(scale.factors) < (q1 - 1.5*iqr) | 
                    unlist(scale.factors) > (q3 + 1.5*iqr))
write.table(t(as.data.frame(scale.factors)), snakemake@output[["scalefactor"]], quote=FALSE, col.names = FALSE)


pdf(snakemake@output[["scaleplot"]])
par(mfrow=c(1,2))
plot(s, scale.tests[[1]], type="l", ylim=range(unlist(scale.tests)), col=rgb(0,0,0,0.1))
for(i in 2:length(scale.tests)) {
  lines(s, scale.tests[[i]], col=rgb(0,0,0,0.1))
}
# Highlight outliers in red
for(i in outliers) {
  lines(s, scale.tests[[i]], col=rgb(1,0,0,0.5), lwd=0.5)
}
hist(unlist(scale.factors),breaks=50, main=NULL, xlab="scale factor")
dev.off()

cn.mpcf <- map2(segcounts, scale.factors, ~round(.x*.y))
cn.mpcf.mtx <- do.call(cbind, cn.mpcf)
cn.mpcf.mtx[1:5,1:10]
write_tsv(as.data.frame(cn.mpcf.mtx), snakemake@output[["cn"]])
write_tsv(as.data.frame(all.segs), snakemake@output[["sc_segs"]])
