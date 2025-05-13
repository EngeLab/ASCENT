# Collect, QC and re-call segments at clone level
library(tidyverse)
library(uwot)
library(dbscan)
library(ComplexHeatmap)
library(circlize)
library(furrr)
library(patchwork)
source(snakemake@params[["src_general"]])
options(scipen=999)
try(RhpcBLASctl::blas_set_num_threads(snakemake@threads), silent = TRUE)

# Parameters/setup
gamma <- as.numeric(snakemake@wildcards[["gamma"]])
min_bins <- snakemake@params[["bin_filter"]]
min_clone_size <- snakemake@params[["min_cells"]]
k_threshold <- snakemake@params[["kmeans_threshold"]]
k_max <- snakemake@params[["kmeans_max"]]
dim_reduction_method <- snakemake@wildcards[["method"]]
bin_diff_max <- snakemake@params[["bin_diff_max"]]
bin_diff_zlim <- snakemake@params[["bin_diff_zlim"]]
clone_frac_diff <- snakemake@params[["clone_frac_diff"]]

# Debug
cat(snakemake@wildcards[["patient_id"]],"gamma:",gamma,
    "min bins per segment:",min_bins, "min clone size:",min_clone_size,"k-means threshold:",k_threshold, "k-means max:",k_max,"\n")

# Run on one instance
res <- read_tsv(snakemake@input[["mpcf"]])
res$filter <- res$n.probes < min_bins
res$start <- c(1, cumsum(res$n.probes)[-length(res$n.probes)] + 1)
res$end <- cumsum(res$n.probes)
res.f <- res[!res$filter,]
sum(res$filter)

# Load correct bins
bins <- read_tsv(snakemake@input[["bins"]], skip=1, col_names = c("chr","start","end","idx","bin")) %>%
  mutate(
    chr_short=sub("chr","",chr),
    chr_int=as.integer(case_when(chr_short == "X" ~ "23", chr_short == "Y" ~ "24", TRUE ~ chr_short)),
    chr=factor(chr, levels=unique(chr))
  )

# New functions
zscale <- function(x) {
  # If the whole vector is a single value, return 0 instead of NaN
  if(all(x == x[1])) return(setNames(rep(0,length(x)),nm=names(x)))
  (x - mean(x)) / sd(x)
}

manhattan.dist <- function(factor, segs) {
  sum(abs((segs*factor)-round(segs*factor)))
}

get_clones_states <- function(x, min_cells=5, plot=F, bins=NULL){
  uniq_mtx <- unique(t(x))
  # Vectorized counting
  row_strings <- apply(t(x), 1, paste, collapse="")
  uniq_strings <- apply(uniq_mtx, 1, paste, collapse="")
  counts <- tabulate(match(row_strings, uniq_strings))
  
  # Assign clone numbers (0 for < min_cells, sequential for others) 
  clone_ids <- numeric(length(counts)) # initialize with all 0
  valid_clone <- counts >= min_cells
  clone_ids[valid_clone] <- seq_len(sum(valid_clone))
  cell_clones <- clone_ids[match(row_strings, uniq_strings)]
  names(cell_clones) <- colnames(x)
  
  if(plot){
    smooth <- smooth_bins_subset(t(uniq_mtx), w=50, bins=bins.f)
    states_smooth <- smooth$mat
    colnames(states_smooth) <- clone_ids
    
    chr_labels <- factor(smooth$bins$chr_short,levels=c(1:22,"X","Y"))
    col_fun = circlize::colorRamp2(c(0, 2, 4, 5, 20, 50), c("blue","white","red", "brown", "green", "yellow"))
    
    row_annot <- rowAnnotation(
      "cells" = anno_barplot(
        counts,
        width = unit(4, "cm"),
        border = FALSE,
        gp = gpar(fill = "darkgrey"),
        axis_param = list(
          side = "top",
          labels_rot = 0
        )
      ),
      annotation_name_rot = 90  # Rotate the "cells" title
    )
    
    ht <- ComplexHeatmap::Heatmap(t(states_smooth), col=col_fun,
                                  name = "cn",
                                  show_row_names = F,
                                  show_column_names =F,
                                  cluster_rows=F, 
                                  cluster_columns = F,
                                  row_split = factor(clone_ids, levels=c(1:max(cell_clones), 0)),
                                  column_split = chr_labels,
                                  show_row_dend = F,
                                  show_column_dend = F, 
                                  row_gap = unit(2, "mm"),
                                  column_gap = unit(0, "mm"),
                                  border=T,
                                  right_annotation = row_annot)
    draw(ht)
  }
  
  return(list(clones=cell_clones, counts=counts, states=uniq_mtx))
}

get_clones_umap <- function(x, min_cells=NULL, metric="manhattan", cluster_dims=30, threads=1, n_nb=20, min_dist=0.1, spread=1.2, seed=42, init="spectral"){
  if(!is.null(seed)) set.seed(seed)
  
  # Find and handle duplicates
  xt <- t(x)
  dups <- duplicated(xt)
  xunique <- xt[!dups,]
  n <- nrow(xunique)
  cat("umap n=",n," (unique sets)\n")
  if(n < n_nb | n < cluster_dims) {
    n_nb <- n-1
    cluster_dims <- n-1
    init <- "random"
  }
  
  # Run UMAP on unique points
  umap_result <- umap(xunique,
                      nn_method = "annoy",
                      n_threads = threads, 
                      n_sgd_threads = 0, 
                      metric=metric, 
                      n_neighbors=n_nb, 
                      n_components=cluster_dims, 
                      min_dist=min_dist, 
                      spread=spread,
                      init=init,
                      pca=NULL)
  
  # Create initial results for unique points
  rall_unique <- data.frame(
    cell=rownames(xunique),
    dim1=umap_result[,1],
    dim2=umap_result[,2]
  )
  
  # Run HDBSCAN on unique points
  hc_unique <- hdbscan(umap_result[,1:cluster_dims], minPts = min_cells)
  rall_unique$clone <- hc_unique$cluster
  
  # Handle duplicates if they exist
  if(sum(dups) > 0){
    xdup <- xt[which(dups),]
    # matches <- xdup %*% t(xunique)
    matches <- tcrossprod(xdup, xunique)
    if(is.matrix(xdup) && nrow(xdup)>1){
      xdup_mag <- rowSums(xdup^2)  
      # Compare each row against its own magnitude
      match_idx <- sapply(1:nrow(matches), function(i) {
        which(matches[i,] == xdup_mag[i])[1]
      })
      # Create results for duplicates
      r_dups <- data.frame(
        cell=rownames(xdup),
        dim1=umap_result[match_idx,1],
        dim2=umap_result[match_idx,2],
        clone=hc_unique$cluster[match_idx]  # Assign same cluster as matching point
      )
    } else {
      xdup_mag <- sum(xdup^2)
      match_idx <- which(matches == xdup_mag)
      r_dups <- data.frame(
        cell=names(which(dups)),
        dim1=umap_result[match_idx,1],
        dim2=umap_result[match_idx,2],
        clone=hc_unique$cluster[match_idx]  # Assign same cluster as matching point
      )
    }
    
    # Combine unique and duplicate results
    rall <- rbind(rall_unique, r_dups)
  } else {
    rall <- rall_unique
  }
  
  # Reorder results to match input data
  rall.o <- rall[match(colnames(x), rall$cell),]
  
  # Calculate cluster centers for non-noise points
  clone_centers <- rall.o %>%
    filter(clone != 0) %>%
    group_by(clone) %>%
    summarize(
      center_x = mean(dim1),
      center_y = mean(dim2)
    )
  
  # Create plot
  plot(rall.o[,c("dim1","dim2")],
       main=paste0("n=",nrow(rall.o),", ",
                   length(unique(rall.o$clone[rall.o$clone != 0]))," clones"),
       sub=paste0("(dims=",paste0(cluster_dims),"; neighb=", n_nb, 
                  "; spread=",spread,"; dist=", min_dist,")"),
       col=factor(rall.o$clone), 
       cex=.8, pch=19)
  
  # Add clone numbers at cluster centers
  if(nrow(clone_centers) > 0) {
    text(clone_centers$center_x, 
         clone_centers$center_y, 
         labels=clone_centers$clone,
         font=2,
         cex=1)
  }
  
  # Mark duplicates with red circles
  if(sum(dups) > 0) {
    points(rall.o[which(dups), c("dim1","dim2")], col="red", pch=21, cex=1.2)
    message("umap: found ", sum(dups), " duplicate points out of ", length(dups))
  }
  
  return(list(
    cluster=rall.o$clone, 
    coordinates=rall.o
  ))
}

get_clones_dist <- function(x, min_cells=NULL, threads=1, metric="manhattan", plot=F){
  require(dbscan)
  dist <- amap::Dist(t(x), method=metric, nbproc = threads)
  hc <- hdbscan(dist, minPts = min_cells)
  if(plot){
    ht <- ComplexHeatmap::Heatmap(as.matrix(dist),
                                  name = "L1 dist",
                                  show_row_names = F,
                                  show_column_names = F, 
                                  row_split = hc$cluster, 
                                  column_split= hc$cluster,
                                  show_row_dend = F, 
                                  show_column_dend = F, 
                                  column_title = "L1 dist")
    draw(ht)
  }
  return(list(cluster=hc$cluster, dist=dist))
}

mtx <- as.matrix(read_tsv(snakemake@input[["cn"]]))
if(nrow(mtx) != nrow(bins)) stop("Unequal length of bin and mtx")

#Adjust copy numbers according to logodds ratios 
logodds<-read_tsv(snakemake@input[["logodds"]])

mtx_adjusted <- sweep(mtx, 2, logodds$multiplication[match(colnames(mtx), logodds$dna_library_id)] %>% replace_na(1), `*`)


# Filter out small segments
f.idx <- rep(res$filter, res$n.probes)
sum(f.idx)
rownames(mtx_adjusted) <- 1:nrow(mtx_adjusted) # Add original bin index as rowname. Allows subsetting AFTER filtering segments by original bin indices
mtx.f <- mtx_adjusted[!f.idx, ]
bins.f <- bins[!f.idx,]
nrow(mtx.f) == nrow(bins.f) # TRUE

#Use winsorised counts for clustering: 
mtx.w <- mtx.f
mtx.w[mtx.w >= 8] <- 8
unique_states <- nrow(unique(t(mtx.w)))

# Take first bin in each segment (should be identical copy state across segment in a single cell)
segs.mtx <- mtx.w[rownames(mtx.w) %in% res.f$start,] 
# Sanity check
# rr <- apply(res.f, 1, function(r){
#   binnames <- r[["start"]]:r[["end"]]
#   mm <- mtx.f[rownames(mtx.f) %in% binnames,]
#   all(apply(mm, 2, function(x) length(unique(x))) == 1)  
# })

# 241024: Run all dim reductions, but use "primary" and "secondary" for plotting
pdf(snakemake@output[["qc_dimred"]],width=12, height=12)
if(dim_reduction_method=="umap"){
  if(unique_states < snakemake@params[["umap_unique_threshold"]]){
    warning(paste0("Too low complexity for UMAP (n=",unique_states,"). Getting clones from unique copy states"))
    r <- get_clones_states(mtx.w, min_cells=min_clone_size, plot=T, bins=bins.f)
    clones.umap <- list(cluster=r$clones)
  } else {
    cat("Calculating UMAP...")
    clones.umap <- get_clones_umap(mtx.w, min_cells=min_clone_size, 
                                   threads=snakemake@threads, 
                                   metric=snakemake@params[["umap_metric"]],
                                   cluster_dims=as.integer(snakemake@params[["umap_cluster_dims"]]), 
                                   n_nb=as.integer(snakemake@params[["umap_neighbors"]]), 
                                   min_dist=as.numeric(snakemake@params[["umap_min_dist"]]), 
                                   spread=as.numeric(snakemake@params[["umap_spread"]]), 
                                   seed=as.integer(snakemake@params[["umap_seed"]]))
    cat("done.\n")  
  }
  clone.df <- tibble(
    dna_library_id=colnames(mtx.w),
    clone=clones.umap$cluster
  )
}
if(dim_reduction_method=="L1"){
  cat("Manhattan dist...")
  mtx.smooth <- smooth_bins_subset(mtx.w, w=10, bins=bins.f)[["mat"]]
  clones.dist <- get_clones_dist(mtx.smooth, min_cells=min_clone_size, threads=snakemake@threads, metric="manhattan", plot=T)
  cat("done.\n")
  clone.df <- tibble(
    dna_library_id=colnames(mtx.w),
    clone=clones.dist$cluster
  )
}
if(dim_reduction_method=="lv"){
  cat(paste0("Running SNN + louvain in",snakemake@wildcards[["patient_id"]],"\n"))
  snn <- Seurat::FindNeighbors(clones.dist$dist)
  snn.lv <- Seurat::FindClusters(snn[["snn"]], algorithm=1, resolution=1, group.singletons=F)
  clones.snn <- setNames(as.numeric(snn.lv[[1]]),nm=colnames(snn[["snn"]])) # Avoid 0 by turning factor numeric
  clone.df <- tibble(
    dna_library_id=colnames(mtx.w),
    clone=clones.snn
  )
}
dev.off()
clone.df <- clone.df %>% mutate(across(everything(), as.character))
print(table(clone.df$clone))

### QC Clones
# CN matrix x clone x bins (median profile)
# Drop class 0 = unclassified by hdbscan, unless the only clone left
clone.ids <- split(clone.df$dna_library_id, clone.df$clone)
if(length(names(clone.ids))>1 && "0" %in% names(clone.ids)) {
  clone.ids <- clone.ids[names(clone.ids) != "0"]
}
cn.clones <- lapply(clone.ids, function(i) mtx.w[,i])
seg.clones <- lapply(clone.ids, function(i) segs.mtx[,i])
# DEBUG: Inspect individual clones
# Heatmap(t(seg.clones[[4]]),cluster_rows = T, cluster_columns = F, col = col_fun, show_row_names = F, show_column_names = F)

# !! Rounding median non-integers (from even-number cell vectors with ~even splits). Stay with integer copy numbers for QC purposes, not in clone-level segmentation
clone.mtx.list <- lapply(cn.clones, function(x) {
  if(is.matrix(x)) round(matrixStats::rowMedians(x)) else median(x)
})
clone.mtx <- do.call(cbind, clone.mtx.list)

# CN matrix x clone x segment (median profile)
# 1. Original segs to clones
cn.clones.seg <- lapply(cn.clones, function(m){
  apply(res.f, 1, function(x){
    o.idx <- x[["start"]]:x[["end"]] # Original bin index
    r <- m[rownames(m) %in% o.idx,]
    round(median(r))
  })
})
clone.mtx.seg <- do.call(cbind, cn.clones.seg)

### Correlation of clones
cor_clust_plot <- function(x, label="bins", order=NULL){
  if(ncol(x)<=1) {
    p <- ggplot() + annotate("text", x=0.5, y=0.5, label=paste0("Too few clones (n=",ncol(x),")")) + theme_void()
    return(p)
  }
  diff <- outer(1:ncol(x), 1:ncol(x), 
                Vectorize(function(i, j) sum(x[,i]!=x[,j])))
  colnames(diff) <- colnames(x)
  rownames(diff) <- colnames(x)
  hc <- hclust(dist(diff,method="manhattan"),method="complete")
  df <- as.data.frame(diff) %>% rownames_to_column("clone") %>% 
    pivot_longer(-clone,names_to = "comparison", values_to = label)
  if(!is.null(order)){
    df <- df %>% mutate(clone=factor(clone, levels=names(order)),
                        comparison=factor(comparison, levels=names(order)))
  } else {
    df <- df %>% mutate(clone=factor(clone, levels=hc$order),
                        comparison=factor(comparison, levels=hc$order))
  }
  
  # Plot correlation grid
  df %>%
    ggplot(aes(x=clone, y=comparison, fill=.data[[label]]/nrow(x))) + geom_tile(show.legend=T) +
    geom_text(aes(label=.data[[label]]), color="black", fontface="bold", show.legend = F, size=3.1) +
    geom_text(aes(label=.data[[label]], color=.data[[label]]), fontface="bold", show.legend = F, size=3) +
    scale_fill_viridis_c() + 
    scale_color_viridis_c(direction = -1) +
    theme_dntr(axis_ticks = T, axis_labels = F) +
    scale_x_discrete(expand=c(0.01,0.01)) +
    scale_y_discrete(expand=c(0.01,0.01))
  # theme(strip.text=element_blank(), 
  # legend.direction = "vertical", 
  # legend.position = "right")
}
sort_from_diploid <- function(x){
  # Get closest to diploid
  mean_clones_dmg <- sort(colSums(abs(x - 2)))
  diploid <- names(mean_clones_dmg)[which.min(mean_clones_dmg)]
  diff <- colSums(x[,diploid] != x)
  o <- sort(diff)
  return(o)
}
# Plot 1: Bin + segment level
p1 <- cor_clust_plot(clone.mtx, label="bins", order=sort_from_diploid(clone.mtx)) + 
  cor_clust_plot(clone.mtx.seg, label="segs", order=sort_from_diploid(clone.mtx.seg)) +
  plot_annotation("Clone vs clone difference, sorted from diploid (lower left).\nBins (left) and segments (right)")

# Filter and qc segments and divergent cells
# For each clone (single cells), compare internal consistency per segment
cnsegs.qc <- lapply(cn.clones, function(m){
  rsegs <- apply(res.f, 1, function(x){
    o.idx <- x[["start"]]:x[["end"]] # Original bin index
    r <- m[rownames(m) %in% o.idx,]
    # Fraction of bins in segment that do not match clone median
    # Will always be 0 or 1? Because mpcf should return same value for whole segment
    #!!Note: Rounded median to avoid non-integer cn counts (for QC purposes)
    colSums(round(median(r)) != r) / nrow(r)
  })
  # By cell: fraction of segments that do not match median. Per cell. No weight by segment length! 
  cells <- rowSums(rsegs) / sum(!res$filter)
  # By segment: Fraction of cells that do not match median
  segs <- colSums(rsegs) / ncol(m)
  list("cells"=cells, "segs"=segs)
})
# cnsegs.qc.bins <- lapply(cnsegs.qc, function(x){
#   sapply(1:length(x), function(i){
#     rep(x[[i]],res.f$n.probes[[i]])
#   }) %>% unlist
# })

# Plots
cnsegs.qc.segs <- do.call(rbind, cnsegs.qc %>% map("segs"))
rownames(cnsegs.qc.segs) <- names(cnsegs.qc)
cnsegs.qc.cells <- cnsegs.qc %>% map("cells")

# Calculate cell-level distance to median of clone at bin-level (summarized as fraction of bins discordant per CELL)
cn.clones.diff <- setNames(lapply(names(cn.clones), function(i){
  # clone.mtx.list -- vectors of each clone, filtered bin length
  # cn.clones -- matrices of cells x bins per clone
  if(is.matrix(cn.clones[[i]])) {
    colSums(cn.clones[[i]] != clone.mtx.list[[i]]) / nrow(mtx.w)
  }
}), nm=names(cn.clones))
diff.frac.bins <- setNames(unlist(cn.clones.diff,use.names=F),
                           nm=unlist(lapply(cn.clones.diff,names),use.names=F))
diff.frac.segs <- setNames(unlist(cnsegs.qc.cells,use.names=F),
                           nm=unlist(lapply(cnsegs.qc.cells,names),use.names=F))
# Z-scale by clone (x-mean)/sd
zscale.bins <- lapply(cn.clones.diff, zscale)
zscale.segs <- lapply(cnsegs.qc.cells, zscale)
diff.frac.bins.zscale <- setNames(unlist(zscale.bins,use.names=F),
                                  nm=unlist(lapply(zscale.bins,names)))
diff.frac.segs.zscale <- setNames(unlist(zscale.segs,use.names=F),
                                  nm=unlist(lapply(zscale.segs,names),use.names=F))

# Remove outliers based differential bin fraction (z-scaled and max 20% different)
if(is.null(bin_diff_max) | is.null(bin_diff_zlim)) stop("No bin difference filters set. Stopping")
diff_filter_z <- abs(diff.frac.bins.zscale) > bin_diff_zlim
diff_filter_abs <- diff.frac.bins > bin_diff_max
diff_filter <- diff_filter_z | diff_filter_abs
cat(sum(diff_filter), "of", length(diff_filter),"clone-level outliers dropped:\n",
    paste0("zscale to median >",bin_diff_zlim," n=",sum(diff_filter_z), 
           " and absolute diff to clone median >",bin_diff_max*100,"% n=",sum(diff_filter_abs)),"\n")
diff.good.cells <- names(which(!diff_filter))
diff.bad.cells <- names(which(diff_filter))
# Drop empty clones after diff filter
keep.idx <- sapply(seg.clones, function(x) sum(!colnames(x) %in% diff.bad.cells) >= min_clone_size)
cn.clones_f <- lapply(cn.clones[keep.idx], function(x) {
  x[, colnames(x) %in% diff.good.cells, drop = FALSE]
})
seg.clones_f <- lapply(seg.clones[keep.idx], function(x) {
  x[, colnames(x) %in% diff.good.cells, drop = FALSE]
})
# Add filtered clone id
clone.df$clone_filtered <- ifelse(clone.df$dna_library_id %in% diff.bad.cells, 0, clone.df$clone)
# Recalculate bin level diff fraction in filtered clones
cnsegs_f.qc <- lapply(cn.clones_f, function(m){
  rsegs <- apply(res.f, 1, function(x){
    o.idx <- x[["start"]]:x[["end"]] # Original bin index
    r <- m[rownames(m) %in% o.idx,]
    # Fraction of bins in segment that do not match clone median
    # Will always be 0 or 1? Because mpcf should return same value for whole segment
    #!!Note: Rounded median to avoid non-integer cn counts (for QC purposes)
    colSums(round(median(r)) != r) / nrow(r)
  })
  # By cell: fraction of segments that do not match median. Per cell. No weight by segment length! 
  cells <- rowSums(rsegs) / sum(!res$filter)
  # By segment: Fraction of cells that do not match median
  segs <- colSums(rsegs) / ncol(m)
  list("cells"=cells, "segs"=segs)
})
cnsegs_f.qc.cells <- cnsegs_f.qc %>% map("cells")

# Plot qc before and after diff bin filter
pdf(snakemake@output[["qc_corplot"]])
plot(p1)
par(mfrow=c(2,1),mar=c(2,4,2,2))
o <- order(sapply(cnsegs.qc.cells, mean))
boxplot(cnsegs.qc.cells[o],xlab="clone",main="Fraction of bins different from clone (before/after filtering)",las=2)
boxplot(cnsegs_f.qc.cells[o],xlab="clone",las=2)
dev.off()

# Subclustering by k-means
kneedle <- function(x, y, S = NULL, i=NULL, plot=T) {
  # Normalize x and y to [0, 1]
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  # Calculate differences from the line connecting the first and last points
  diffs <- (1 - x_norm) * y_norm[1] + x_norm * y_norm[length(y_norm)] - y_norm
  
  # 2. Use k at maximum difference
  knee_index <- which.max(diffs)
  
  if(plot) {
    plot(x, diffs)
    # points(x[close_points],diffs[close_points],col="grey",pch=19,xlab="k")
    points(x[knee_index],diffs[knee_index],col="orange",pch=19,xlab="k")
    points(x[which.max(diffs)],diffs[which.max(diffs)],col="red",pch=19,xlab="k")
    abline(h=S,col="grey")
  }
  if(!is.null(i)) title(paste0("clone ",i, "; k=",x[knee_index],"; threshold: ",diffs[which.max(diffs)] > S))
  
  # Apply sensitivity threshold
  # Enough that the maximum diff has surpassed threshold
  if (diffs[which.max(diffs)] > S) {
    return(x[knee_index])
  } else {
    return(NULL)
  }
}
pdf(snakemake@output[["qc_kmplot"]])
k.list <- lapply(sort(as.integer(names(cn.clones_f))), function(i){
  i <- as.character(i)
  n_cells <- ncol(cn.clones_f[[i]])
  n_states <- nrow(unique(t(cn.clones_f[[i]])))
  k_max <- min(k_max, n_cells - 1, n_states)
  cat(snakemake@wildcards[["patient_id"]], ", gamma",gamma,": clone ",i,"/",length(names(cn.clones_f)),"k_max:",k_max, "n_cells:",n_cells,"n_states:",n_states,"\n")
  # Smooth to 10-bin windows before k-means (for speed)
  x_smooth <- smooth_bins_subset(cn.clones_f[[i]], w=10)[["mat"]]
  wss <- sapply(1:k_max, function(k) {
    kmeans(t(x_smooth), centers=k, nstart=10)$tot.withinss
  })
  par(mfrow=c(1,2))
  plot(1:k_max, wss, main=paste0("clone ",i), type="b", xlab="k", ylab="total within cluster sum of sq")
  if(k_max>1) {
    k <- kneedle(1:k_max, wss, S=k_threshold, i=i, plot=T)
  } else {
    k <- NULL
  }
  if(!is.null(k)){
    cat("PASS threshold. Split at k=",k,"\n")
    kres <- kmeans(t(cn.clones_f[[i]]), centers=k, nstart=10)
    return(kres$cluster)
  } else {
    cat("FAIL threshold, no split\n")
    # If no k passes threshold, all in k_clone = 1
    return(setNames(rep(1,n_cells),nm=colnames(cn.clones_f[[i]])))
  }
})
dev.off()
k.vec <- unlist(k.list)

### Heatmaps
dna_qc <- read_tsv(snakemake@input[["dna_qc"]])
meta_patient <- read_tsv(snakemake@input[["metadata_patient"]])
if(!"timepoint" %in% colnames(meta_patient)) meta_patient$timepoint <- NA
if(!is.null(snakemake@input[["rna_phase"]])){
  rna_phase <- read_tsv(snakemake@input[["rna_phase"]])
} else {
  rna_phase <- tibble(cell_id = NA, cell_phase = NA, dna_library_id = NA)
}

# Collect metadata, ensure same order of cells as mtx
meta <- clone.df[match(colnames(mtx), clone.df$dna_library_id),] %>%
  left_join(dna_qc,by="dna_library_id") %>%
  mutate(rna_qc=ifelse(dna_to_cellid(dna_library_id) %in% rna_phase$cell_id, "PASS", "FAIL"),
         diff_bins=diff.frac.bins[dna_library_id],
         diff_segs=diff.frac.segs[dna_library_id],
         diff_bins_z=diff.frac.bins.zscale[dna_library_id],
         diff_segs_z=diff.frac.segs.zscale[dna_library_id],
         kmeans=unname(k.vec[dna_library_id]),
         kmeans=ifelse(is.na(kmeans),0,kmeans),
         timepoint=meta_patient$timepoint[match(dna_library_id, meta_patient$dna_library_id)],
  ) %>% 
  group_by(clone_filtered, kmeans) %>% 
  mutate(n_clone_comb=n(),
         kmeans_rev=ifelse(n_clone_comb < min_clone_size, 0, kmeans),
         clone_comb=paste(clone_filtered,kmeans_rev,sep="_"))

# Create new clone level cn matrix
clone.mtx.list_f <- lapply(cn.clones_f, function(x) {
  if(is.matrix(x)) round(matrixStats::rowMedians(x)) else median(x)
})
clone.mtx_f <- do.call(cbind, clone.mtx.list_f)

# Reorder by average distance from diploid (per CLONE)
mean_clones_dmg <- sort(colSums(abs(clone.mtx_f - 2)))
mean_clones_clust.all <- sub(".*_(.*)$","\\1", names(mean_clones_dmg))
mean_clones_clust <- c(mean_clones_clust.all[! mean_clones_clust.all %in% c("G2M","S","0")],"S","G2M","0")
meta$clone_dmgOrder <- factor(meta$clone_filtered, levels=mean_clones_clust)

# Plot colors
col_fun = circlize::colorRamp2(c(0, 2, 4, 5, 20, 50), c("blue","white","red", "brown", "green", "yellow"))
col_fun_dups = circlize::colorRamp2(c(0,0.1,0.4), c("white","white","black"))
col_fun_prob = circlize::colorRamp2(c(0,1), c("white","darkgreen"))
rna_qc_colors <- setNames(c("black","white"),nm=c("PASS","FAIL"))
diff_filter_colors <- setNames(c("black","white"),nm=c("TRUE","FALSE"))

lgd_list = list(
  Legend(title="copy number", col_fun = col_fun, break_dist = 1, title_position = "lefttop-rot")
)

### Heatmap inputs (smoothed)
smooth <- smooth_bins_subset(mtx_adjusted, bins=bins, w=10) # TODO: Add S/G2M but drop bad QC
ht_mtx <- smooth$mat
chr_labels <- factor(smooth$bins$chr_short,levels=c(1:22,"X","Y"))

# Genome level
bins.all <- read_tsv(snakemake@input[["all_bins"]],col_names=c("chr","start","end")) %>% 
  mutate(idx=seq_along(chr),
         bin=seq_along(chr),
         chr_short=sub("chr","",chr),
         chr_int=as.integer(case_when(chr_short == "X" ~ "23", chr_short == "Y" ~ "24", TRUE ~ chr_short)),
         chr=factor(chr, levels=unique(chr)))
mtx.e <- matrix(
  NA,
  nrow = max(bins.all$idx),
  ncol = ncol(mtx),
  dimnames = list(
    1:max(bins.all$idx),
    colnames(mtx)
  )
)
mtx.e[bins$idx, ] <- mtx_adjusted
smooth <- smooth_bins_subset(mtx.e, bins.all, w=10)
ht_mtx <- smooth$mat
chr_labels <- factor(smooth$bins$chr_short,levels=c(1:22,"X","Y"))

# Order within clones by global L1
smooth.distL1 <- amap::Dist(t(smooth$mat), method="manhattan", nbproc = 8) # Note! Transposed
smooth.hclust <- hclust(smooth.distL1, method="complete")
smooth.hc_order <- smooth.hclust$labels[smooth.hclust$order]
smooth.hc_order.vec <- setNames(seq_along(smooth.hc_order),nm=smooth.hc_order)

# Set final order of cells: clone dmg (as factor) > kmeans dmg (sort) > by L1 dist (sort)
# Reorder subclones (k-means) by average distance from diploid
meta.km.filt <- meta %>% filter(kmeans_rev!=0)
clone.ids.km <- split(meta.km.filt$dna_library_id, meta.km.filt$clone_comb)
clones.km.list <- lapply(clone.ids.km, function(i) round(matrixStats::rowMedians(mtx.f[,i])))
clones.km.mtx <- do.call(cbind, clones.km.list)
mean_clones_dmg.km <- colSums(abs(clones.km.mtx - 2))
meta$kmeans_diff <- mean_clones_dmg.km[meta$clone_comb]
meta$L1_order <- smooth.hc_order.vec[meta$dna_library_id]
cells.final_order <- meta %>% arrange(clone_filtered,kmeans_diff,L1_order) %>% pull(dna_library_id)

# Final clones and metadata
# Filter: clones where x% of cells are x% different from median clone profile
clone_level.bin_filter <- names(which(sapply(cnsegs_f.qc.cells, function(x) sum(x > bin_diff_max)/length(x)) > clone_frac_diff))

# Revised clone assignment (after bin filters and kmeans subclustering)
meta <- meta %>% ungroup() %>% 
  mutate(clone_rev=case_when(clone_filtered=="0" ~ "0",
                             kmeans_rev=="0" ~ "0",
                             clone_filtered %in% clone_level.bin_filter ~ "0",
                             diff_bins > bin_diff_max ~ "0",
                             .default=clone_comb))

# Make median profiles of revised clones - merge identical
clone.ids.rev <- split(meta$dna_library_id, meta$clone_rev)
if(length(names(clone.ids.rev))>1 && "0" %in% names(clone.ids.rev)){
  clone.ids.rev <- clone.ids.rev[names(clone.ids.rev) != "0"]
}
clone.rev.mtx <- do.call(cbind, lapply(clone.ids.rev, function(i) round(matrixStats::rowMedians(mtx.w[,i]))))

# Merge duplicates, if any
dup_idx <- duplicated(t(clone.rev.mtx))
if(any(dup_idx)){
  # Get mapping of duplicates to first occurrence
  clone_map <- setNames(
    apply(t(clone.rev.mtx), 1, function(x) {
      colnames(clone.rev.mtx)[which(apply(t(clone.rev.mtx), 1, function(y) all(x == y)))[1]]
    }),
    nm=colnames(clone.rev.mtx)
  )
  # Show which clones were merged
  dups <- split(names(clone_map), clone_map)
  message("Merged clones (mapping to first occurrence):")
  print(dups[lengths(dups) > 1])
  meta$clone_merged <- ifelse(meta$clone_rev==0, 0, clone_map[meta$clone_rev])
} else {
  meta$clone_merged <- meta$clone_rev
}

# Set final clone ids
meta <- meta %>% 
  group_by(clone_merged) %>%
  mutate(n_merged=n(),
         clone_final=case_when(n_merged < min_clone_size ~ "0",
                               .default=clone_merged))

# Reorder final clones by distance from diploid
order_clones_from_diploid <- function(clone.list, mtx, excluded=c("0","S","G2M")){
  clone.mtx <- do.call(cbind, lapply(clone.list, function(i) {
    if(length(i)>1) round(matrixStats::rowMedians(mtx[,i])) else median(mtx[,i])
  }))
  s <- names(sort(colSums(abs(clone.mtx - 2))))
  # Move special clones/excluded last in order
  reorder <- c(setdiff(s, excluded), intersect(s, excluded))
  return(reorder)
}
clone.ids.final <- split(meta$dna_library_id, meta$clone_final)
clone_dmgOrder.final <- order_clones_from_diploid(clone.ids.final, mtx.w)
meta$clone_final_dmgOrder <- factor(meta$clone_final, levels=clone_dmgOrder.final)

# Build annotation
clone_colors <- get_clones_col(meta$clone)
kmeans_all <- unique(c(meta$kmeans, meta$kmeans_rev))
clone_colors_km <- get_clones_col(kmeans_all[!is.na(kmeans_all)])
clone_colors_rev <- get_clones_col(meta$clone_rev)
clone_colors_final <- get_clones_col(meta$clone_final)

ha <- rowAnnotation(clone=meta$clone,
                    clone_kmeans=meta$kmeans,
                    clone_rev=meta$clone_rev,
                    clone_final=meta$clone_final,
                    timepoint=meta$timepoint,
                    diff_bins=if(diff(range(meta$diff_bins, na.rm = T)) > 0) meta$diff_bins else rep_len(F,nrow(meta)),
                    diff_bins_z=if(diff(range(meta$diff_bins_z, na.rm = T)) > 0) meta$diff_bins_z else rep_len(F,nrow(meta)),
                    diff_filter=abs(meta$diff_bins_z) > bin_diff_zlim | meta$diff_bins > bin_diff_max,
                    rna_qc=meta$rna_qc,
                    reads=anno_barplot(log10(meta$bam_read_pairs), bar_width=1, border=F, baseline = "min"),
                    col=list(clone=clone_colors,
                             clone_kmeans=clone_colors_km,
                             clone_rev=clone_colors_rev,
                             clone_final=clone_colors_final,
                             rna_qc=rna_qc_colors,
                             timepoint=timepoint_cols,
                             diff_filter=diff_filter_colors),
                    show_legend=T, annotation_name_gp = gpar(fontsize=8), na_col = "white")


# 1st heatmap: on primary assigned clones
png(snakemake@output[["heatmap"]], width=2400,height=1400,units="px",res=150)
ht = Heatmap(t(ht_mtx), name="cn", col=col_fun,
             heatmap_legend_param=list(title="cn",col_fun=col_fun,break_dist=1),
             cluster_columns = F,
             row_split = meta$clone_dmgOrder,
             # Opt 1: outside clustering
             row_order = smooth.hc_order, # Named vector matching rownames in mtx
             cluster_row_slices = F,
             cluster_rows= F, # Cluster within each clone, to L1 dist/complete
             row_dend_reorder = F,
             show_row_dend = F,
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
             show_heatmap_legend=F,
             na_col="#cccccc")
ht <- draw(ht, annotation_legend_list=packLegend(list=lgd_list),
           column_title = paste0(snakemake@wildcards[["patient_id"]], " (n=",ncol(mtx),"); gamma=",gamma,
                                 " (",nrow(res.f),"/",nrow(res)," segs); ", dim_reduction_method,
                                 "; min_clone_size=",min_clone_size,"; min_segment=",min_bins,"; k_sens=",k_threshold,
                                 "; frags per bin: ",round(mean(meta$bam_read_pairs)/nrow(mtx.f),1),
                                 " (",as.numeric(snakemake@wildcards[["binsize"]])/1e3,"kb bins)"),
           column_title_side = "top")
dev.off()

# 2nd heatmap: on revised subclustered clones
png(snakemake@output[["heatmap_revised"]], width=2400,height=1400,units="px",res=150)
ht = Heatmap(t(ht_mtx), name="cn", col=col_fun,
             heatmap_legend_param=list(title="cn",col_fun=col_fun,break_dist=1),
             cluster_columns = F,
             row_split = meta$clone_final_dmgOrder,
             # Opt 1: outside clustering
             row_order = cells.final_order, # 
             cluster_row_slices = F,
             cluster_rows= F,
             row_dend_reorder = F,
             show_row_dend = F,
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
             show_heatmap_legend=F,
             na_col="#cccccc")
ht <- draw(ht, annotation_legend_list=packLegend(list=lgd_list),
           column_title = paste0(snakemake@wildcards[["patient_id"]], " (n=",ncol(mtx),"); gamma=",gamma,
                                 " (",nrow(res.f),"/",nrow(res)," segs); ", dim_reduction_method,
                                 "; min_clone_size=",min_clone_size,"; min_segment=",min_bins,"; k_sens=",k_threshold,
                                 "; frags per bin: ",round(mean(meta$bam_read_pairs)/nrow(mtx.f),1),
                                 " (",as.numeric(snakemake@wildcards[["binsize"]])/1e3,"kb bins)"),
           column_title_side = "top")
dev.off()

# Save clones
clone.df.final <- meta %>% ungroup %>% select(dna_library_id, clone, clone_comb, clone_filtered, clone_rev, clone_final, diff_bins, diff_bins_z, kmeans, kmeans_rev, kmeans_diff)
write_tsv(clone.df.final, snakemake@output[["clones"]])