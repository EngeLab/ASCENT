### NOTE! Some functions are duplicated from general.R

# Basic
dna_to_cellid <- function(x){ 
  sub("([A-Z]{3}[0-9]{5})(D?)_([A-Z][0-9]{2}).*","\\1.\\3", x) 
}
rna_to_cellid <- function(x){ 
  sub("([A-Z]{3}[0-9]{5})(R?).([A-Z][0-9]{2}).*","\\1.\\3", x) 
}
filter_vec <- function(x, filter.idx){
  x[filter.idx] <- NA
  return(x)
}
zscale <- function(x) {
  # If the whole vector is a single value, return 0 instead of NaN
  if(all(x == x[1])) return(setNames(rep(0,length(x)),nm=names(x)))
  (x - mean(x)) / sd(x)
}
expand_gaps_vec <- function(x, length.out, idx){
  if(!is.vector(x)) stop("Not a vector")
  if(length(x)!=length(idx)) stop("Length of vector and idx should be equal")
  r <- rep(NA, length.out)
  r[idx] <- x
  return(r)
}
lowess.gc = function(x, y) {
  low = lowess(x, log(y), f=0.05, iter=10);
  z = approx(low$x, low$y, x)
  return(exp(log(y) - z$y))
}
gc.map.correct <- function(bin.counts, gc, map, gc.coefs) {
  predq <- function(qest, gc, scale=1) {
    (qest[1]+gc*qest[2]+gc^2*qest[3])/scale
  }
  gc.scale <- predq(gc.coefs, 0.45, 1)
  weight <- map*predq(gc.coefs, gc, gc.scale)
  bin.counts/weight
}
manhattan.dist <- function(factor, segs, na.rm=F) {
  if(na.rm) segs <- segs[!is.na(segs)]
  sum(abs((segs*factor)-round(segs*factor)))
}
make_cn_colorscale <- function(max_cn) {
  if (max_cn <= 5) {
    return(circlize::colorRamp2(c(0, 2, 4, 5), c("blue","white","red", "brown")))
  } else if(max_cn <= 20) {
    return(circlize::colorRamp2(c(0, 2, 4, 5, 20), c("blue","white","red", "brown", "green")))
  } else {
    return(circlize::colorRamp2(c(0, 2, 4, 5, 20, 50), c("blue","white","red", "brown", "green", "yellow")))
  }
}
get_divergent_segs <- function(mtx){
  apply(mtx, 1, function(x) length(unique(x)) > 1)
}
get_identical_segs <- function(mtx){
  apply(mtx, 1, function(x) length(unique(x))==1)
}

phase_cols <- c("G2M"="purple","S"="darkorange","G1"="#cfcfcf")

get_clones_col <- function(x){
  n_clones = length(unique(x))
  if(n_clones == 1){
    colors <- "#f2f3f4" # Grey
  }
  else if (n_clones <= 21){
    colors <- Polychrome::kelly.colors(n=n_clones+1)[-2] # Drop white
  } else if(n_clones <= 26){
    colors <- Polychrome::alphabet.colors(n_clones)
  } else if(n_clones <= 36){
    colors <- Polychrome::palette36.colors(n_clones)
  } else {
    colors <- Polychrome::createPalette(n_clones,  c("#ff0000", "#00ff00", "#0000ff"))
  }
  names(colors) <- as.character(sort(unique(x)))
  colors[grepl("*G2M$|*S$", names(colors))] <- c("purple","darkorange")
  if(any(grepl("*_0$",names(colors)))) colors[grepl("*_0$",names(colors))] <- "black"
  return(colors)
}

smooth_bins_subset <- function(m, w, bins=NULL){
  sub <- m[seq(1, nrow(m), w), ]
  # Debug: Expand back out: rep(sub,each=10)
  # Faster version for mean across bin windows
  if(!is.null(bins)){
    bins.sub <- bins[seq(1, nrow(bins), w),]
    rownames(sub) <- bins.sub$id
    return(list(mat=sub, bins=bins.sub))
  } else {
    rownames(sub) <- 1:nrow(sub)
    return(list(mat=sub, bins=NULL))
  }
}

smooth_vector <- function(vector, bin_size) {
  n <- length(vector)
  smoothed <- sapply(seq(1, n, by = bin_size), function(i) {
    any(vector[i:min(i + bin_size - 1, n)])
  })
  return(smoothed)
}

timepoint_cols <- c("diagnosis"="#fc8d59",
                    "day 2"="#fee08b",
                    "day 3"="#d9ef8b",
                    "day 5"="#91cf60",
                    "day 15"="#1a9850",
                    "day 29"="#136b39",
                    "relapse"="#d73027",
                    "relapse 1"="#d73027",
                    "relapse 1d29"="#d77c27",
                    "relapse 2"="#002984",
                    "relapse 2d6"="#4f0084",
                    "relapse 2d13"="#4f0084",
                    "relapse 2d10"="#840077")

load_and_join <- function(annotation_files, join_by=NULL) {
  data_list <- map(annotation_files, function(x) read_tsv(x, show_col_types = FALSE))
  joined_data <- data_list[[1]]
  if (length(data_list) > 1) {
    for (i in 2:length(data_list)) {
      joined_data <- joined_data %>%
        left_join(data_list[[i]], by = join_by)
    }
  }
  return(joined_data)
}

theme_dntr <- function (font_size = 14, font_family = "Helvetica", line_size = 0.5, axis_labels=T, axis_ticks=F, legend="bottom") {
  half_line <- font_size/2
  small_rel <- 0.857
  small_size <- small_rel * font_size
  th <- theme_grey(base_size = font_size, base_family = font_family) %+replace% 
    theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black", linewidth=1), 
          legend.justification = "top", legend.background = element_blank(), legend.key = element_blank(),
          legend.key.size = unit(1, "lines"), strip.background = element_blank(), strip.text=element_text(hjust=0,size=12),
          rect = element_rect(fill = "transparent", color = "black", linewidth = 1, linetype = "solid"), axis.line = element_blank(),
          text = element_text(family = font_family, face = "plain", color = "black", size = font_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9, margin = margin(), debug = FALSE),
          # axis.text.x = element_text(color="black", margin = margin(t = small_size/4), vjust = 1), 
          # axis.text.y = element_text(color="black", margin = margin(r = small_size/4), hjust = 1), 
          axis.title.x = element_text(size=ifelse(axis_labels,NA,0), margin = margin(t = small_size/4, b = small_size/4)),
          axis.title.y = element_text(size=ifelse(axis_labels,NA,0), angle = 90, margin = margin(r = small_size/4, l = small_size/4)), 
          axis.ticks = element_line(linewidth=ifelse(axis_ticks,0.5,0)),
          axis.text = element_text(size=ifelse(axis_ticks,NA,0), margin = margin(t = small_size/4, b = small_size/4)))
  if(legend=="bottom"){
    th <- th %+replace%
      theme(legend.position="bottom", legend.direction = "horizontal")
  } else if(!is.null(legend)) {
    th <- th %+replace%
      theme(legend.position=legend)
  }
  return(th)
}

# Done
load_bins <- function(bins_file, gc_file, map_file, cytoband_file){
  # bins file: chr, start, end
  # gc and map file: one value per line/bin
  # cytoband_file: chr, start, end, cytoband as first columns
  bins <- read_tsv(bins_file,
                   col_names=c("chr","start","end"), 
                   col_select=1:3, show_col_types = F)
  cytobands <- read_tsv(cytoband_file,
                        col_names=c("chr","start","end","cytoband"),
                        col_select=1:4, show_col_types = F)
  gc <- as.numeric(readLines(gc_file))
  map <- as.numeric(readLines(map_file))
  
  stopifnot(
    "Unequal length of bins and gc/mappability" = 
      nrow(bins) == length(gc) && nrow(bins) == length(map)
  )
  
  arms <- cytobands %>% mutate(arm=substr(cytoband,1,1)) %>% 
    filter(chr %in% unique(bins$chr)) %>% 
    group_by(chr, arm) %>% 
    summarize(start=min(start),end=max(end),.groups = "drop")
  arms.gr <- GRanges(seqnames=arms$chr, ranges=IRanges(arms$start, arms$end), arm=arms$arm)
  bins.gr <- GRanges(seqnames=bins$chr, ranges=IRanges(bins$start, bins$end))
  gr.intersect <- findOverlaps(bins.gr, arms.gr, type = "within")
  bins$arm <- arms.gr$arm[subjectHits(gr.intersect)]
  
  bins$chr <- factor(bins$chr, levels=paste0("chr",c(1:22,"X","Y"))) 
  bins$gc <- gc
  bins$map <- map
  bins$id <- 1:nrow(bins)
  
  
  return(bins)
}

# Calculate scale factors from panel of normal diploid cells sequenced similary
#   normal_counts: matrix Normal cell panel bincounts
#   bins: data.frame Table of bins and positions, including gc and mappability fraction per bin
#   good_bins_idx: vector Index of "good" bin positions in bins. Or logical vector
#   gc: vector If gc counts not in bins, provide separately 
#   map: vector If mapability fraction not in bins, provide separately
#   exclude: vector Ids of cells to exclude from normal panel (for self exclusion)
#   plot: logical Plot results
# Returns:
#   named list ft and gcmap for respective normalization factors
calc_normal_scale_factors <- function(normal_counts, bins, good_bins_idx=NULL, gc=NULL, map=NULL, exclude=NULL, plot = FALSE) {
  # Read in data
  if(!is.matrix(normal_counts)){
    if(is.character(normal_counts) && file.exists(normal_counts)){
      normal_counts <- as.matrix(data.table::fread(normal_counts))
    } else {
      stop("normal_counts is not a matrix and cannot be read from file")
    }
  }
  
  # Exclude defined cells
  if(!is.null(exclude)) {
    cell_idx <- !colnames(normal_counts) %in% exclude
    if(sum(!cell_idx)>0) message(paste0("Excluded ",sum(!cell_idx)," cells (self-exclusion)"))
  } else {
    cell_idx <- T
  }
  # Subset for good bins if provided
  if(!is.null(good_bins_idx)){
    counts <- normal_counts[good_bins_idx, cell_idx]
    bins <- bins[good_bins_idx, ]
  } else {
    counts <- normal_counts[, cell_idx]
  }
  
  # Get mappability and gc from bins data if not provided
  if(is.null(gc)) gc <- bins$gc
  if(is.null(map)) map <- bins$map
  
  autochr <- bins$chr %in% paste0("chr", 1:22)
  xchr <- bins$chr %in% "chrX"
  ychr <- bins$chr %in% "chrY"
  
  # Calculate summed counts for autosomes and X chromosome
  auto.counts <- rowSums(counts[autochr, ])
  x.counts <- rowSums(counts[xchr, ])
  
  # Freeman-Tukey + lowess normalization
  auto.ft <- sqrt(auto.counts) + sqrt(auto.counts + 1)
  auto.ft.lowess <- lowess.gc(gc[autochr], auto.ft)
  x.ft <- sqrt(x.counts) + sqrt(x.counts + 1)
  x.ft.lowess <- lowess.gc(gc[xchr], x.ft)
  
  # GC + mappability normalization
  auto.gc.coefs <- summary(lm(auto.counts ~ gc[autochr] + I(gc[autochr]^2)))$coefficients[, 1]
  auto.gcmap <- gc.map.correct(auto.counts, gc[autochr], map[autochr], auto.gc.coefs)
  auto.gcmap.mc <- auto.gcmap/median(auto.gcmap)
  
  x.gc.coefs <- summary(lm(x.counts ~ gc[xchr] + I(gc[xchr]^2)))$coefficients[, 1]
  x.gcmap <- gc.map.correct(x.counts, gc[xchr], map[xchr], x.gc.coefs)
  x.gcmap.mc <- x.gcmap/median(x.gcmap)
  
  # Combine factors for all chromosomes (autosomes, X, and Y)
  factors.ft <- c(auto.ft.lowess, x.ft.lowess, rep(1, sum(ychr)))
  factors.gcmap <- c(auto.gcmap.mc, x.gcmap.mc, rep(1, sum(ychr)))
  
  if(plot){
    # Plot for debug
    colors <- Polychrome::alphabet.colors(length(unique(bins$chr)))
    plot(factors.gcmap, cex=.1, col=adjustcolor(colors[factor(bins$chr)], alpha=0.1), main="gcmap scaling factors from normal cell panel")
    lines(caTools::runmean(factors.gcmap, k=50), col="black")
    abline(h=mean(factors.gcmap),col="red")
    plot(factors.ft, cex=.1, col=adjustcolor(colors[factor(bins$chr)], alpha=0.1), main="ft+lowess scaling factors from normal cell panel")
    lines(caTools::runmean(factors.ft, k=50), col="black")
    abline(h=mean(factors.ft),col="red")
  }
  
  return(list("ft"=factors.ft, "gcmap"=factors.gcmap))
}

create_pseudobulk_analysis <- function(counts_matrix, bins_info, good_bins, cell_metadata, normal_counts = NULL,
                                       params = list(min_bins = 10, min_cells = 8, gamma = 1, recall_mad_filter = 7, scale_range = c(1.5, 4.5))
) {
  # Input validation
  stopifnot(
    "Bins info must match counts" = 
      nrow(counts_matrix) == nrow(bins_info),
    "Required columns missing in bins" = 
      all(c("chr", "start", "end", "gc", "map", "arm") %in% colnames(bins_info)),
    "Good bins index is incorrect length" = 
      length(good_bins) > 0 && length(good_bins) < nrow(bins_info),
    "Required columns missing in metadata" = 
      all(c("dna_library_id", "clone") %in% colnames(cell_metadata)),
    "Cell with clones missing in counts data" = 
      all(cell_metadata$dna_library_id %in% colnames(counts_matrix))
  )
  
  cat(sprintf("Processing %d cells with clone information (%d total)\n", sum(!is.na(cell_metadata$clone)), ncol(counts_matrix)))
  cat(sprintf("Read in %d good bins (%d total)\n", length(good_bins), nrow(bins_info)))
  
  # Exclude clone 0 for clone-level processing
  clone_keep_idx <- !is.na(cell_metadata$clone) & cell_metadata$clone!="0"
  clone_order <- gtools::mixedsort(unique(cell_metadata$clone[clone_keep_idx]))
  clone_counts <- table(factor(cell_metadata$clone[clone_keep_idx], levels=clone_order))
  cat(sprintf("Loaded %d clones (%d cells unclassified/class 0):\n", length(clone_counts), sum(cell_metadata$clone %in% "0")))
  for(i in seq_along(clone_counts)) {
    cat(sprintf("  Clone %s: %d cells\n", names(clone_counts)[i], clone_counts[i]))
  }
  
  cnv_data <- list(
    bins = list(
      all = bins_info %>% mutate(
        id = row_number(),
        chr = factor(chr, levels=paste0("chr", c(1:22, "X", "Y"))),
        chr_int = as.integer(case_when(
          chr == "chrX" ~ "23", 
          chr == "chrY" ~ "24", 
          TRUE ~ sub("chr", "", chr))),
        end_cum = cumsum(as.numeric(end)-as.numeric(start)),
        good = id %in% good_bins
      ),
      good = NULL
    ),
    
    # Cells: ordered by counts_matrix
    cells = tibble(
      dna_library_id = colnames(counts_matrix),
      raw_counts = map(seq_len(ncol(counts_matrix)), ~counts_matrix[, .x]),
      clone_initial = cell_metadata$clone[match(colnames(counts_matrix), cell_metadata$dna_library_id)],
      !!!cell_metadata[match(colnames(counts_matrix), cell_metadata$dna_library_id), 
                       !names(cell_metadata) %in% c("dna_library_id", "clone")]
    ),
    
    # Clone level data
    clones = tibble(
      clone_id = clone_order,
      revision = "initial",
      parent_clone = NA_character_,
      cells = map(clone_id, ~cell_metadata$dna_library_id[cell_metadata$clone %in% .x]),
      # Clone-level count summaries
      raw_counts = list(NULL),
      ft_lowess = list(NULL),
      ft_lowess_normal = list(NULL),
      gcmap = list(NULL),
      gcmap_normal = list(NULL),
      segment_means = list(NULL),  
      segment_bins = list(NULL),  
      scale_tests = list(NULL),
      scale_factor = list(NULL),
      cn = list(NULL),    # Clone-specific copy numbers
      residuals = list(NULL)
    ),
    
    # Segment level data after processing: [initial, filtered, ...]
    # Updates d$clones$segment_means and (if not null) cn_integers based on latest iteration, since its always clone based
    segments = list(),
    
    parameters = params
    
    # Optional: state tracking
  )
  
  # Set good bins
  cnv_data[["bins"]][["good"]] <- cnv_data[["bins"]][["all"]][good_bins,]
  
  # Add normal panel scaling if provided
  if(!is.null(normal_counts)) {
    cat("Calculating normal cell scale factors\n")
    nf <- calc_normal_scale_factors(normal_counts=normal_counts, 
                                    bins=cnv_data[["bins"]][["all"]], 
                                    good_bins_idx=cnv_data[["bins"]][["good"]]$id, 
                                    exclude=colnames(counts_matrix),
                                    plot=F)
    cnv_data[["bins"]][["good"]]$normal_factors_ft <- nf$ft
    cnv_data[["bins"]][["good"]]$normal_factors_gcmap <- nf$gcmap
  }
  
  return(cnv_data)
}

normalize_counts <- function(data, methods = c("ft_lowess", "gcmap"), subset=NULL) {
  # Usage: data <- normalize_counts(data, methods = c("ft_lowess", "gcmap"))
  methods <- match.arg(methods, several.ok = TRUE)
  bins <- data[["bins"]][["good"]]
  
  cat(paste0("\n\nNormalizing raw counts with methods: ", paste(methods,collapse=", "),"\n"))
  cat("Normal panel scaling:",paste(!is.null(bins$normal_factors_ft),!is.null(bins$normal_factors_gcmap),collapse=", "),"\n")
  
  run_idx = TRUE
  if(!is.null(subset)){
    run_idx <- data[["clones"]][["clone_id"]] %in% subset
    if(sum(run_idx)==0) stop("Subset not matching any indices in clones table")
    cat("Using subset of clones:",paste(subset,collapse=", "),"\n")
  }
  
  data[["clones"]][run_idx,] <- data[["clones"]][run_idx,] %>%
    mutate(
      raw_counts = map(cells, function(cells) {
        cell_counts <- data[["cells"]][["raw_counts"]][data[["cells"]][["dna_library_id"]] %in% cells]
        Reduce(`+`, cell_counts)[bins$id] # Sum bins per clone
      })
    )
  
  if("ft_lowess" %in% methods) {
    data[["clones"]][run_idx,] <- data[["clones"]][run_idx,] %>%
      mutate(
        ft_lowess = map(raw_counts, function(counts) {
          ft <- sqrt(counts) + sqrt(counts + 1)
          lowess.gc(bins$gc, ft)
        })
      )
    
    if(!is.null(bins$normal_factors_ft)) {
      data[["clones"]][run_idx,] <- data[["clones"]][run_idx,] %>%
        mutate(
          ft_lowess_normal = map2(ft_lowess, list(bins$normal_factors_ft), 
                                  function(ft, norm) ft/norm)
        )
    }
  }
  
  if("gcmap" %in% methods) {
    data[["clones"]][run_idx,] <- data[["clones"]][run_idx,] %>%
      mutate(
        gcmap = map(raw_counts, function(counts) {
          gc.coefs <- summary(lm(counts ~ gc + I(gc^2), 
                                 data = data.frame(counts = counts, gc = bins$gc)))$coefficients[,1]
          gcmap <- gc.map.correct(counts, bins$gc, bins$map, gc.coefs)
          gcmap/median(gcmap)
        })
      )
    
    if(!is.null(bins$normal_factors_gcmap)) {
      data[["clones"]][run_idx,] <- data[["clones"]][run_idx,] %>%
        mutate(
          gcmap_normal = map2(gcmap, list(bins$normal_factors_gcmap), 
                              function(gc, norm) gc/norm)
        )
    }
  }
  
  return(data)
}

call_segments <- function(data, 
                          gamma = 1,
                          norm_segments = "ft_lowess_normal",
                          norm_ratio = "gcmap_normal",
                          segs_slot = "initial",
                          weights = NULL,
                          verbose = TRUE) {
  stopifnot(
    "Empty/missing data in norm_segments slot"=
      all(!sapply(data[["clones"]][[norm_segments]], is.null)),
    "Empty/missing data in norm_ratio slot"=
      all(!sapply(data[["clones"]][[norm_ratio]], is.null))
  )
  
  bins <- data[["bins"]][["good"]]
  norm.mtx <- do.call(cbind, data[["clones"]][[norm_segments]])
  colnames(norm.mtx) <- data[["clones"]][["clone_id"]]
  
  mpcf.input <- data.frame(
    chr = bins$chr_int,
    pos = bins$start,
    norm.mtx,
    check.names = F
  )
  
  # Optional Y input for segment ratio calculation
  if(!is.null(norm_ratio)) {
    normY.mtx <- do.call(cbind, data[["clones"]][[norm_ratio]])
    colnames(normY.mtx) <- data[["clones"]][["clone_id"]]
    mpcf.inputY <- data.frame(
      chr = bins$chr_int,
      pos = bins$start,
      normY.mtx,
      check.names=F
    )
  } else {
    mpcf.inputY <- NULL
  }
  
  # Default cell-based weights: none
  if(is.null(weights)) weights <- 1
  
  # Run segmentation
  res <- copynumber::multipcf(
    data = mpcf.input,
    Y = mpcf.inputY,
    pos.unit = "bp",
    arms = bins$arm,
    gamma = gamma,
    normalize = FALSE,
    fast = FALSE,
    verbose = TRUE,
    return.est = TRUE,
    w = weights
  )
  
  # Process results
  segs <- res$segments[,1:5] %>%
    mutate(
      seg_idx = 1:n(),
      start = c(1, cumsum(n.probes)[-length(n.probes)] + 1),
      end = cumsum(n.probes),
      chr = factor(case_when(
        chrom %in% 1:22 ~ paste0("chr", chrom),
        chrom == 23 ~ "chrX",
        chrom == 24 ~ "chrY"
      ), levels=paste0("chr",c(1:22,"X","Y"))),
      filtered = FALSE,
      merged_from = NA,
      high_residual = FALSE
    ) %>% select(chr, everything(), -chrom)
  
  data[["segments"]][[segs_slot]] <- segs
  
  # Update data structure
  # NOTE: Could just store segments, but would have to keep track of segment boundaries/sizes at all times
  data[["clones"]][["segment_means"]] <- as.list(as.data.frame(res$segments[,-c(1:5)]))
  data[["clones"]][["segment_bins"]] <- as.list(as.data.frame(res$estimates[,-c(1:2)]))
  
  return(data)
}

merge_small_segments <- function(data, current="initial", revision="merged", min_bins_filter, boundary_filter=NULL, update_clones=T) {
  if(is.null(data[["segments"]][[current]])) stop("Segment slot empty")
  merged <- data$segments[[current]]
  merged$merged_from <- as.character(merged$merged_from)
  
  is_boundary_segment <- function(idx, df) {
    chr_arm_segments <- which(df$chr == df$chr[idx] & df$arm == df$arm[idx])
    return(idx == min(chr_arm_segments) || idx == max(chr_arm_segments))
  }
  
  find_segments_to_merge <- function(df, min_bins_filter, boundary_filter) {
    if(!is.null(boundary_filter)) {
      boundary_idx <- which(sapply(1:nrow(df), function(i) {
        !df$filtered[i] && is_boundary_segment(i, df) && df$n.probes[i] < boundary_filter
      }))
      if(length(boundary_idx) > 0) return(boundary_idx[1])
    }
    small_idx <- which(!df$filtered & df$n.probes < min_bins_filter)
    if(length(small_idx) > 0) return(small_idx[which.min(df$n.probes[small_idx])])
    return(numeric(0))
  }
  
  while(length(i <- find_segments_to_merge(merged, min_bins_filter, boundary_filter)) > 0) {
    cat("Processing segment:", i, "chr", merged$chr[i], merged$arm[i],
        "n.probes:", merged$n.probes[i],
        if(!is.null(boundary_filter) && is_boundary_segment(i, merged)) " (boundary segment)" else "", "\n")
    
    left_neighbor <- if(i > 1 &&
                        merged$chr[i-1] == merged$chr[i] &&
                        merged$arm[i-1] == merged$arm[i]) i-1 else NA
    right_neighbor <- if(i < nrow(merged) &&
                         merged$chr[i] == merged$chr[i+1] &&
                         merged$arm[i] == merged$arm[i+1]) i+1 else NA
    
    cat("-- Left neighbor:", left_neighbor, "Right neighbor:", right_neighbor, "\n")
    
    left_size <- if(!is.na(left_neighbor)) merged$n.probes[left_neighbor] else 0
    right_size <- if(!is.na(right_neighbor)) merged$n.probes[right_neighbor] else 0
    cat("-- Left size:", left_size, "Right size:", right_size, "\n")
    
    current_merged <- if(!is.na(merged$merged_from[i])) strsplit(merged$merged_from[i], ",")[[1]] else character(0)
    
    if(left_size >= right_size && !is.na(left_neighbor)) {
      cat("-- Merging with left neighbor\n")
      merged$end.pos[left_neighbor] <- merged$end.pos[i]
      merged$n.probes[left_neighbor] <- merged$n.probes[left_neighbor] + merged$n.probes[i]
      if(is.na(merged$merged_from[left_neighbor])) {
        merged$merged_from[left_neighbor] <- merged$seg_idx[i]
      } else {
        merged$merged_from[left_neighbor] <- paste(c(merged$merged_from[left_neighbor], merged$seg_idx[i], current_merged), collapse=",")
      }
      # merged$filtered[i] <- TRUE
      merged <- merged[-i,]
    } else if(!is.na(right_neighbor)) {
      cat("-- Merging with right neighbor\n")
      merged$start.pos[right_neighbor] <- merged$start.pos[i]
      merged$n.probes[right_neighbor] <- merged$n.probes[right_neighbor] + merged$n.probes[i]
      if(is.na(merged$merged_from[right_neighbor])) {
        merged$merged_from[right_neighbor] <- merged$seg_idx[i]
      } else {
        merged$merged_from[right_neighbor] <- paste(c(merged$merged_from[right_neighbor], merged$seg_idx[i], current_merged), collapse=",")
      }
      # merged$filtered[i] <- TRUE
      merged <- merged[-i,]
    } else {
      # NOTE: Make segment values NA, goes to clones slot (below)
      cat(paste0("-- No valid neighbors to merge with! Setting values to NA\n"))
      merged$filtered[i] <- TRUE
    }
  }
  
  # Add new segments and update clone metrics
  data[["segments"]][[revision]] <- merged[order(as.numeric(row.names(merged))), ]
  # Update clones, but run cn integer calculation separately
  if(update_clones) data <- update_clone_metrics(data, norm_slot = "gcmap_normal", segs_slot = revision, calc_cn=FALSE)
  
  return(data)
}

calc_segment_means <- function(norm_data, segs) {
  # Calculate mean per segment using the normalized values
  vapply(1:nrow(segs), function(j) {
    bins_in_seg <- segs$start[j]:segs$end[j]
    if("filtered" %in% colnames(segs) && segs$filtered[j]) {
      NA_real_
    } else {
      mean(norm_data[bins_in_seg]) 
    }
  }, numeric(1))
}

mask_high_residuals <- function(data, segs_slot = NULL, new_slot = "masked", max_residual=0.3, clone_filter_fraction=0.3, update_clones=T){
  # Default to latest iteration of segments
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])] 
  if(!segs_slot %in% names(data[["segments"]])) stop("Specified segment slot not found")
  cat(sprintf("Using segments: %s\n", segs_slot))
  
  residuals <- sapply(data$clones$residuals, function(x) x > max_residual)
  high_residual_vals <- rowMeans(residuals, na.rm=TRUE) > clone_filter_fraction
  
  clone_names <- data[["clones"]][["clone_id"]]
  high_residual_clones <- apply(residuals, 1, function(x) {
    valid_high <- which(!is.na(x) & x)  # Get indices where valid and TRUE
    paste(clone_names[valid_high], collapse=",")  # Use clone names for the high residuals
  })
  
  cat(sprintf("Found %d segments with high residuals", sum(high_residual_vals,na.rm = T)))
  
  # If no masking - return original data
  if(sum(high_residual_vals,na.rm = T)==0) return(data)
  
  # Otherwise print to new slot
  if(new_slot %in% names(data[["segments"]])){
    new_slot <- ifelse(!grepl("\\d+$", new_slot), 
                       paste0(new_slot,2), 
                       paste0(new_slot,as.numeric(gsub(".*?(\\d+)$", "\\1", new_slot))+1))
  }
  cat(sprintf("Writing segments to: %s\n", new_slot))
  
  data[["segments"]][[new_slot]] <- data[["segments"]][[segs_slot]] %>% 
    mutate(
      high_residual=high_residual_vals,
      high_residual_clones=high_residual_clones,
      filtered=ifelse(filtered, filtered, high_residual)
    )
  
  if(update_clones)  data <- update_clone_metrics(data, norm_slot = "gcmap_normal", segs_slot=new_slot, calc_cn=T)
  
  return(data)
}

mask_segment <- function(data, mask=NULL, segs_slot = NULL, update_clones=T){
  # Masks specific segments by original index, and returns them to the same slot
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])] 
  if(!segs_slot %in% names(data[["segments"]])) stop("Specified segment slot not found")
  if(! length(mask)>0 & ! all(mask %in% data[["segments"]][[segs_slot]][["seg_idx"]])) stop("Incorrect segments defined for masking")
  cat(sprintf("\nWill mask (index): %s\n", mask))
  cat(sprintf("\nUsing segments: %s\n", segs_slot))
  
  data[["segments"]][[segs_slot]] <- data[["segments"]][[segs_slot]] %>% 
    mutate(
      filtered=ifelse(seg_idx %in% mask, TRUE, filtered)
    )
  if(update_clones)  data <- update_clone_metrics(data, norm_slot = "gcmap_normal", segs_slot=segs_slot, calc_cn=T)
  
  return(data)
}

update_clone_metrics <- function(data, norm_slot="gcmap_normal", segs_slot=NULL, clone_slot=NULL, calc_cn=TRUE, subset=NULL, ...){
  # Default to latest added segments
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])]
  stopifnot("Not recognized normalization" = norm_slot %in% names(data[["clones"]]),
            "No segments to process" = !is.null(data[["segments"]]), 
            "Unrecognized segment slot" = segs_slot %in% names(data[["segments"]]))
  
  run_idx = TRUE # Default run all clones
  update_msg = "\n\n> Updating clone metrics (all clones)\n"
  if(!is.null(subset)){
    run_idx <- data[["clones"]][["clone_id"]] %in% subset
    if(sum(run_idx)==0) stop("Subset not found in clone ids\n")
    update_msg = paste0("\n\n> Updating subset of clones (",paste(subset,collapse=", "),")")
  }
  
  # TODO: Break out to separate function (update_clone_composition)
  if(!is.null(clone_slot)) {
    cat("> Updating clone assignments from slot:",clone_slot,"\n")
    current <- data[["clones"]] %>% select(current_clone=clone_id, cell=cells) %>% unnest(cell)
    cross_tab <- data[["cells"]] %>% 
      select(cell=dna_library_id, new_clone=!!clone_slot) %>% 
      left_join(current, by="cell")
    
    # Cross tabulate and check for clones that changed: 
    # 1. across rows/current: excluded cell from original clone
    # 2. across columns/new: added from NA originally
    clone_table <- table(current=cross_tab$current_clone, new=cross_tab$new_clone, useNA="always")
    print(clone_table)
    changed_clones <- unique(
      c(names(which(apply(clone_table[!is.na(rownames(clone_table)), ], 1, function(x) sum(x > 0) > 1))),
        names(which(apply(clone_table[, !colnames(clone_table) %in% c(NA,"0")], 2, function(x) sum(x > 0) > 1))))
    )
    cat("Clones with changed composition:\n")
    print(changed_clones)
    
    for(clone_id in changed_clones) {
      cat(">>> Updating",clone_id,"\n")
      new_cells <- cross_tab %>%
        filter(new_clone == clone_id) %>%
        pull(cell)
      
      idx <- which(data[["clones"]]$clone_id == clone_id)
      current_cells <- data[["clones"]]$cells[[idx]]
      
      cat(sprintf("%d -> %d",length(current_cells), length(new_cells)),"\n")
      dropped <- setdiff(current_cells, new_cells)
      added <- setdiff(new_cells, current_cells)
      
      if(length(dropped) > 0) cat(sprintf("Dropped (%d): %s", length(dropped), paste(dropped,collapse=", ")), "\n")
      if(length(added) > 0) cat(sprintf("Added (%d): %s", length(added), paste(added,collapse=", ")), "\n")
      cat("\n")
      
      data[["clones"]]$cells[[idx]] <- new_cells
      data[["clones"]]$revision[[idx]] <- sub("clone_","",clone_slot)
      data[["clones"]]$cn[[idx]] <- NA
      
      # Update summed clone counts and normalizations
      data <- normalize_counts(data, methods = c("ft_lowess", "gcmap"), subset=clone_id)
      cat("\n\n")
    }
  }
  
  data# Second run or re-calculating cn
  if(calc_cn && all(sapply(data[["clones"]][["cn"]][run_idx], is.null))) cn_version="\nCalculating" else cn_version="\nRe-calculating"
  cat(update_msg,
      "\nSegment slot:", segs_slot,
      "\nNormalized counts:", norm_slot,
      cn_version,"scale factor/cn integers:", calc_cn,"\n")
  
  segment_means_list <- lapply(data[["clones"]][[norm_slot]][run_idx], function(x) calc_segment_means(x, data[["segments"]][[segs_slot]]))
  data[["clones"]][["segment_means"]][run_idx] <- segment_means_list
  data[["clones"]][["segment_bins"]][run_idx] <- lapply(segment_means_list, function(x) rep(x, data[["segments"]][[segs_slot]][["n.probes"]])) 
  if(calc_cn) data <- calc_cn_integers(data, subset=subset, ...)
  
  return(data)
}

calc_cn_integers <- function(data, subset=NULL){
  cat("\nCalculating scale factors and copy number integers\n")
  
  run_idx = TRUE # Default run all clones
  if(!is.null(subset)) {
    run_idx <- data[["clones"]][["clone_id"]] %in% subset
    if(sum(run_idx)==0) stop("Subset does not match any clones")
    cat("Using subset of clones:",paste(subset,collapse=", "),"\n")
  }
  
  data[["clones"]][run_idx,] <- data[["clones"]][run_idx,] %>%
    mutate(
      scale_factor = map(cells, ~{
        data$cells %>%
          filter(dna_library_id %in% .x) %>%
          pull(correct_scalefactor) %>%
          na.omit() %>%
          { if (length(.) == 0) NA_real_ else median(.) }
      }),
      cn = map2(segment_bins, scale_factor, ~round(.x * .y)),
      residuals = map2(segment_means, scale_factor, ~abs(.x * .y - round(.x * .y)))
    )
  return(data)
}

group_by_arm <- function(regions, segments, max_gap=1) {
  # Create identifier for each chromosome arm
  group_ids <- sapply(regions, function(r) {
    chr <- segments$chr[r]
    if(chr %in% c(23, 24, "chrX", "chrY")) {
      return(paste0(chr))  # Just use chromosome for X and Y
    } else {
      return(paste0(chr, "_", segments$arm[r]))  # Use chromosome_arm for others
    }
  })
  
  # Split into list by chromosome/arm
  groups <- split(regions, group_ids)
  
  # Further split each group by gaps if needed
  final_groups <- unlist(lapply(groups, function(arm_regions) {
    if(length(arm_regions) <= 1) return(list(arm_regions))
    
    arm_regions <- sort(arm_regions)
    sub_groups <- list()
    current_group <- arm_regions[1]
    
    for(i in 2:length(arm_regions)) {
      if(arm_regions[i] - tail(current_group, 1) <= max_gap + 1) {
        current_group <- c(current_group, arm_regions[i])
      } else {
        sub_groups <- c(sub_groups, list(current_group))
        current_group <- arm_regions[i]
      }
    }
    sub_groups <- c(sub_groups, list(current_group))
    return(sub_groups)
  }), recursive=FALSE)
  
  return(final_groups)
}

#TODO(VZ): Do not sub-split if less than min_cells/min_clone_size
split_mixed_clones <- function(data, clones=NULL, segs_slot=NULL, cells_slot="gcmap_normal", 
                               residual_threshold = 0.3, min_cells = 10,
                               improvement_threshold = 0.8,
                               total_improvement_threshold = 1,
                               plot=TRUE, verbose=TRUE, update_clones=TRUE) {
  stopifnot("cells_slot should be gcmap or gcmap_normal"=
              cells_slot %in% c("gcmap","gcmap_normal"))
  # Unless specific slot needed, default to latest
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])]
  segs <- data[["segments"]][[segs_slot]]
  
  # Initialize assignments data frame with all cells: from current clone assignments data$clone
  assignments <- tibble(
    cell = unlist(data$clones$cells),
    original_clone = rep(data$clones$clone_id, sapply(data$clones$cells, length)),
    current_clone = original_clone,
  )
  
  cat(sprintf("\n^ Looking for admixed clones \nUsing segments in: %s \nUsing normalized counts in: %s\n\n",segs_slot,cells_slot))
  
  # Store split results separately
  split_results <- list()
  
  clone_indices <- if (!is.null(clones)) {
    which(data$clones$clone_id %in% clones)
  } else {
    seq_along(data$clones$clone_id)
  }
  
  for(i in clone_indices) {
    clone_id <- data$clones$clone_id[i]
    n_cells <- length(data$clones$cells[[i]])
    cat("\nProcessing clone:", clone_id, "with", n_cells, "cells\n")
    
    if(clone_id=="0"){
      cat("Skipping - clone 0\n")
      next
    } else if(n_cells < min_cells) {
      cat(paste0("Skipping ",clone_id," (too few cells)\n"))
      next
    }
    
    # Calculate initial residuals
    residuals <- data$clones$residuals[[i]]
    high_res_regions <- which(residuals > residual_threshold)
    
    if(length(high_res_regions) == 0) {
      cat("No high residual regions found\n")
      next
    }
    
    # Group by chromosome arm
    region_groups <- group_by_arm(high_res_regions, segs)
    
    # Calculate score for each group
    group_scores <- sapply(region_groups, function(regions) {
      sum(residuals[regions]^2 * log(segs$n.probes[regions]))
    })
    
    # Order groups by score
    ordered_groups <- region_groups[order(group_scores, decreasing=TRUE)]
    ordered_scores <- sort(group_scores, decreasing=TRUE)
    
    # Get cells in this clone
    clone_cell_idx <- which(assignments$original_clone == clone_id)
    
    # Initialize current nodes for this clone
    current_nodes <- list(list(
      name = clone_id,
      cells = assignments$cell[clone_cell_idx],
      # Index to pull cell-level info from data[["cells"]]
      cells_idx = match(assignments$cell[clone_cell_idx], data$cells$dna_library_id), 
      profile = data$clones$segment_means[[i]] * data$clones$scale_factor[[i]]
    ))
    
    cat(sprintf("Evaluating %d high residual regions\n", length(ordered_groups)))
    # Process each region group in order
    for(g in seq_along(ordered_groups)) {
      regions <- ordered_groups[[g]]
      if(verbose){
        cat("\n>",g,"of",length(ordered_groups),
            "\nTesting group", g, "chromosome", unique(segs$chr[regions]),
            "arm", unique(segs$arm[regions]),
            "\nregions:", paste(regions, collapse=","),
            "\nscore:", ordered_scores[g], "\n")
      }
      
      # Try splitting each current node
      new_nodes <- list()
      
      for(node_idx in seq_along(current_nodes)) {
        current_node <- current_nodes[[node_idx]]
        
        if(length(current_node$cells) < min_cells) {
          new_nodes[[length(new_nodes) + 1]] <- current_node
          next
        }
        
        # Get data for these regions using relative indices
        na.idx <- which(is.na(residuals))
        
        # Create cell-level normalizations if needed (gcmap or gcmap_normal)
        if(is.null(data[["cells"]][[cells_slot]]) | any(sapply(data$cells[[cells_slot]][current_node$cells_idx], is.null))){
          bins <- data[["bins"]][["good"]]
          # If slot doesn't exist, create it with NULL - otherwise proceed to "fill in" 
          if(is.null(data[["cells"]][[cells_slot]])) data[["cells"]][[cells_slot]] <- list(NULL)
          data[["cells"]][[cells_slot]][current_node$cells_idx] <- lapply(data[["cells"]][["raw_counts"]][current_node$cells_idx], function(x){
            x <- x[bins$id] # subset for good bins
            gc.coefs <- summary(lm(x ~ bins$gc + I(bins$gc^2)))$coefficients[,1]
            r <- gc.map.correct(x, bins$gc, bins$map, gc.coefs)
            rm <- r/mean(r) # mean instead of median at single-cell level
            if(cells_slot == "gcmap_normal") rm <- rm/bins$normal_factors_gcmap
            return(rm)
          })
        }
        
        # Summarize over current segments and add to data$clones$segment_means_byCell
        # NOTE: Safer to recalculate based on cells needed in this loop, rather than store
        gcmap_cells <- data[["cells"]][[cells_slot]][current_node$cells_idx]
        gcmap_segs <- lapply(gcmap_cells, function(x){
          r <- sapply(1:nrow(segs), function(i) {
            idx <- segs[["start"]][i]:segs[["end"]][i]
            mean(x[idx])
          })
          r[na.idx] <- NA
          r
        })
        gcmap_segs.mtx <- do.call(cbind, gcmap_segs)
        colnames(gcmap_segs.mtx) <- current_node$cells
        
        region_counts <- t(gcmap_segs.mtx[regions,, drop = FALSE]) # transpose for kmeans
        
        current_residuals <- abs(current_node$profile - round(current_node$profile))
        region_ssr <- ifelse(length(regions)>1,
                             sum(current_residuals[regions]^2, na.rm=TRUE),
                             sum(current_residuals[regions]^2, na.rm=TRUE))
        total_ssr <- sum(current_residuals^2, na.rm=TRUE)
        
        # Try k=2 clustering
        set.seed(42)
        km <- kmeans(region_counts, centers=2)
        
        # Calculate new profiles
        subclone_profiles <- list()
        for(k in 1:2) {
          cluster_cells <- km$cluster == k
          if(sum(cluster_cells) > 0) {
            subclone_profiles[[k]] <- rowMeans(gcmap_segs.mtx[,cluster_cells, drop=FALSE] *
                                                 data$clones$scale_factor[[i]])
          }
        }
        
        # Calculate improvements
        new_region_residuals <- sapply(subclone_profiles, function(profile) {
          abs(profile[regions] - round(profile[regions]))
        })
        new_region_ssrs <- ifelse(length(regions)>1,
                                  colSums(new_region_residuals^2, na.rm=TRUE),
                                  sum(new_region_residuals^2, na.rm=TRUE))
        
        new_total_residuals <- sapply(subclone_profiles, function(profile) {
          abs(profile - round(profile))
        })
        new_total_ssrs <- ifelse(length(regions)>1,
                                 colSums(new_total_residuals^2, na.rm=TRUE),
                                 sum(new_total_residuals^2, na.rm=TRUE))
        
        region_improvement <- mean(new_region_ssrs) / region_ssr
        total_improvement <- mean(new_total_ssrs) / total_ssr
        
        accept_split <- region_improvement < improvement_threshold && total_improvement < total_improvement_threshold
        
        # Create position labels with better formatting
        position_label <- if(length(regions) > 1) {
          paste(sapply(regions, function(i) {
            sprintf("seg_idx=%d, %d%s:%s-%s (%d bins)",
                    segs$seg_idx[i],
                    segs$chr[i],
                    segs$arm[i],
                    format(segs$start.pos[i], scientific = FALSE, big.mark = ","),
                    format(segs$end.pos[i], scientific = FALSE, big.mark = ","),
                    segs$n.probes[i])
          }), collapse = "\n")
        } else {
          sprintf("seg_idx=%d, %d%s:%s-%s (%d bins)",
                  segs$seg_idx[regions],
                  segs$chr[regions],
                  segs$arm[regions],
                  format(segs$start.pos[regions], scientific = FALSE, big.mark = ","),
                  format(segs$end.pos[regions], scientific = FALSE, big.mark = ","),
                  segs$n.probes[regions])
        }
        
        # Calculate percentage changes (now will be positive when residuals increase)
        region_pct_change <- (mean(new_region_ssrs) - region_ssr)/region_ssr * 100
        total_pct_change <- (mean(new_total_ssrs) - total_ssr)/total_ssr * 100
        
        if(plot){
          # Collect eval info
          eval_info <- sprintf("region SSR: %.2f (%.2f) %+.0f%%\ntotal SSR: %.2f (%.2f) %+.0f%%",
                               mean(new_region_ssrs), region_ssr, region_pct_change,
                               mean(new_total_ssrs), total_ssr, total_pct_change)
          
          
          # Set up the plot with adjusted margins for title space
          par(mar = c(4, 4, 6, 2) + 0.1)  # Increase top margin
          
          # Create main plot
          plot(subclone_profiles[[1]], type="l", col="red",
               main="",
               sub="",
               # ylim=c(0, max(sapply(subclone_profiles, function(x) max(x,na.rm=T)))*1.2),
               ylim=c(0, max(sapply(subclone_profiles, function(x) max(x,na.rm=T)))*1.5),
               xlab="", ylab="Copy number estimate")
          
          # Add split lines
          lines(subclone_profiles[[2]], col="blue")
          lines(current_node$profile, col="black", lty=2)
          
          # Add highlighted region
          rect(xleft = min(regions)-0.5, xright = max(regions)+0.5,
               ybottom = par("usr")[3], ytop = par("usr")[4],
               col = rgb(0.9, 0.9, 0.2, 0.2), border = NA)
          
          # Add legend
          legend("topright",
                 legend=c("Original", paste0("Split 1 (",sum(km$cluster==1),")"),
                          paste0("Split 2 (",sum(km$cluster==2),")")),
                 col=c("black", "red", "blue"),
                 lty=c(2,1,1),
                 cex = 0.8,
                 bg = c(0,0,0,alpha=0.3),  # transparent background
                 box.lwd = 0.5,           # thinner box
                 xjust = 1,               # right justify
                 yjust = 1,               # top justify
                 y.intersp = 0.8) 
          
          # Add evaluation info at top left
          legend("topleft", legend = eval_info, bty = "n", cex = 0.8)
          
          # Add titles using mtext
          mtext(sprintf("Clone %s, node %s: %s",
                        clone_id, current_node$name,
                        ifelse(accept_split, "PASS", "FAIL")),
                side = 3, line = 4, cex = 1.2, font = 2)
          
          mtext(sprintf("Split evaluation at region(s) %s",
                        paste(regions, collapse="-")),
                side = 3, line = 3, cex = 1)
          
          # Add position label(s) at bottom, adjusting line position based on number of regions
          if(length(regions) > 1) {
            mtext(position_label, side = 1, line = 3, cex = 0.8)
          } else {
            mtext(position_label, side = 1, line = 2.5, cex = 0.8)
          }
          
          # Reset par to original settings
          par(mar = c(5, 4, 4, 2) + 0.1)
        }
        
        if(verbose){
          # Print evaluation
          cat("\nRegional evaluation:")
          cat("\nOriginal SSR:", region_ssr)
          cat("\nKm clusters:")
          print(table(km$cluster))
          cat("\nNew clone SSRs:", paste(new_region_ssrs, collapse=", "))
          cat("\nRegion improvement:", region_improvement)
          cat("\n\nTotal evaluation:")
          cat("\nOriginal SSR:", total_ssr)
          cat("\nNew clone SSRs:", paste(new_total_ssrs, collapse=", "))
          cat("\nTotal improvement:", total_improvement, "\n")
        }
        
        if(accept_split) {
          cat("\nSplit ACCEPTED\n")
          
          # Update assignments for this split
          for(k in 1:2) {
            cluster_cells_idx <- km$cluster == k
            cluster_cells <- current_node$cells[cluster_cells_idx]
            new_clone_name <- paste0(current_node$name, "_", k)
            match_idx <- match(cluster_cells, assignments$cell)
            assignments$current_clone[match_idx] <- new_clone_name
            
            # Add new node
            new_nodes[[length(new_nodes) + 1]] <- list(
              cells = cluster_cells,
              cells_idx = current_node$cells_idx[cluster_cells_idx],
              profile = subclone_profiles[[k]],
              name = new_clone_name
            )
          }
          
          # Store split results
          split_results[[paste0(clone_id, "_group", g, "_node", current_node$name)]] <- list(
            regions = regions, 
            cluster_assignments = km$cluster,
            cell_ids = current_node$cells,
            region_ssrs = new_region_ssrs,
            total_ssrs = new_total_ssrs,
            region_improvement,
            total_improvement,
            subclone_profiles = subclone_profiles,
            parent = current_node$name
          )
          cat("\nAfter storing split for clone", clone_id, 
              "\nKey:", paste0(clone_id, "_group", g, "_node", current_node$name),
              "\nCurrent structure of split_results:\n")
          str(split_results,max.level = 2)
          
        } else {
          cat("\nSplit REJECTED - Reason: ",
              if(region_improvement >= improvement_threshold) "Insufficient regional improvement"
              else "Total SSR not improved on average", "\n")
          # Keep original node if split rejected
          new_nodes[[length(new_nodes) + 1]] <- current_node
        }
      }
      
      # Update current nodes for next region
      current_nodes <- new_nodes
    }
    cat("Finished clone", clone_id,"\n")
    cat(length(split_results),"\n")
    cat(names(split_results),"\n")
  }
  
  # Return both the split results and the final assignments
  if(length(split_results)>0){
    cat(sprintf("\n>>> Finished: %d splits made. Updated clone assignments.\n",length(split_results)))
    diffs <- assignments %>% filter(original_clone!=current_clone)
    diffs_clone <- assignments %>% filter(original_clone %in% unique(diffs$original_clone))
    print(table(diffs_clone[,2:3]))
    
    # Return new clone assignments 
    # Return class 0 or new clone set
    data[["cells"]][["clone_split"]] <- ifelse(
      data$cells$clone_initial == "0", "0",
      assignments$current_clone[match(data$cells$dna_library_id, assignments$cell)]
    )
    # data[["cells"]][["clone_split"]] <- assignments$current_clone[match(data$cells$dna_library_id, assignments$cell)]
    
    # Update clone info
    if(update_clones){
      data <- update_split_clones(data, split_results)
      all_new_clones <- as.vector(sapply(split_results, function(x) paste0(x$parent, "_", c(1,2))))
      data <- normalize_counts(data, methods="gcmap", subset=all_new_clones)
      data <- update_clone_metrics(data, subset=all_new_clones, calc_cn=TRUE)
    }
    
    # Return split information to new subclones. Should append if multiple
    data[["clones"]][["splits"]] <- vector("list", nrow(data[["clones"]]))
    for (split_name in names(split_results)) {
      sr <- split_results[[split_name]]
      parent_clone <- sr$parent
      new_clone_ids <- paste0(parent_clone, c("_1", "_2"))
      for(k in 1:2){
        subclone_idx <- which(data[["clones"]]$clone_id == new_clone_ids[k])
        subclone_split_info <- list(
          regions = sr$regions,
          region_ssrs = sr$region_ssrs,
          total_ssrs = sr$total_ssrs,
          region_improvement = sr$region_improvement,
          total_improvement = sr$total_improvement,
          subclone_profiles = sr$subclone_profiles[[k]],
          parent = parent_clone
        )
        data[["clones"]][["splits"]][[subclone_idx]] <- subclone_split_info
      }
    }
    
  } else {
    cat("No splits passed threshold.\n")
  }
  
  return(data)
}

# Function to update clones tibble with split results
update_split_clones <- function(data, split_results) {
  new_clones_list <- lapply(names(split_results), function(split_name) {
    split <- split_results[[split_name]]
    parent_id <- split$parent
    
    # Get cells for each subclone from cluster assignments
    cells_1 <- names(split$cluster_assignments[split$cluster_assignments == 1])
    cells_2 <- names(split$cluster_assignments[split$cluster_assignments == 2])
    new_clone_names <- paste(parent_id, c(1,2), sep="_") # Because k always k=2 in this setup
    
    # Get raw counts for new clones from parent clone
    counts_1 <- data[["cells"]][["raw_counts"]][match(cells_1, data[["cells"]][["dna_library_id"]])]
    counts_2 <- data[["cells"]][["raw_counts"]][match(cells_2, data[["cells"]][["dna_library_id"]])]
    
    # Create tibble with two rows (one for each split)
    tibble(
      clone_id = new_clone_names,
      revision = "split",
      parent_clone = parent_id,
      cells = list(cells_1, cells_2)
    )
  })
  new_clones_tbl <- bind_rows(new_clones_list)
  
  # Remove parent clones and add new split clones
  data[["clones"]] <- data[["clones"]] %>%
    filter(!clone_id %in% unique(new_clones_tbl$parent_clone)) %>%
    bind_rows(new_clones_tbl)
  
  return(data)
}

sort_from_diploid_root <- function(x, tree_method="fastme.ols", na.rm=TRUE) {
  if(na.rm) {
    keep <- !rowSums(is.na(x)) > 0
    if(sum(!keep)>0) message(paste0("Removed ",sum(!keep)," rows with NAs"))
    x <- x[keep,]
  }
  # Add fake diploid root
  xroot <- cbind(x,2)
  colnames(xroot) <- c(colnames(x),"root")
  dist_mat <- dist(t(xroot), method="manhattan")
  if(tree_method=="fastme.ols"){
    tree <- fastme.ols(dist_mat) # Create and root tree    
  } else if(tree_method=="nj"){
    tree <- nj(dist_mat) # Create and root tree
  } else if(tree_method=="fastme.bal"){
    tree <- fastme.bal(dist_mat) # Create and root tree    
  } else {
    stop("tree_method does not match any known function")
  }
  
  tree <- root.phylo(tree, outgroup = which(tree$tip.label == "root"), resolve.root = TRUE)
  tree <- drop.tip(tree, tip = "root")
  n_tips <- length(tree$tip.label)
  edge_order <- tree$edge[,2]
  tip_indices <- edge_order[edge_order <= n_tips]
  ordered_labels <- rev(tree$tip.label[tip_indices])
  return(ordered_labels)
}

get_clone_types <- function(data, segs_slot=NULL, diploid=NULL, bin_level=FALSE){
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])]
  
  # Create matrix at bin or segment level
  if(bin_level){
    # Expand to include gaps
    cn_bins.l <- lapply(data$clones$cn, function(x) expand_gaps_vec(x, length.out=nrow(data$bins$all), idx=data$bins$good$id))
    m <- do.call(cbind, cn_bins.l)
    colnames(m) <- data$clones$clone_id
  } else {
    cn_segs.l <- lapply(data$clones$cn, function(x) x[data[["segments"]][[segs_slot]]$start])
    m <- do.call(cbind, cn_segs.l)
    colnames(m) <- data$clones$clone_id 
  }
  
  if(is.null(diploid)) {
    diploid <- names(which.min(colMeans(abs(m-2),na.rm = T)))
    cat("Diploid clone not provided. Estimated to be clone",diploid,"\n")
  }
  
  tumor <- !colnames(m) %in% diploid
  seg_clonal <- FALSE
  seg_subclonal <- FALSE
  if(sum(tumor)>1){
    seg_clonal <- get_identical_segs(m[,tumor]) & rowMeans(m[,diploid]!=m[,tumor])>0
    seg_subclonal <- get_divergent_segs(m[,tumor]) 
  }
  clone_type = rep("none", nrow(m))
  clone_type[seg_clonal] = "clonal"
  clone_type[seg_subclonal] = "subclonal"
  return(clone_type)
}

merge_duplicate_clones <- function(data, segs_slot=NULL, clone_slot=NULL){
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])]
  if(is.null(clone_slot)) clone_slot <- names(data[["cells"]])[max(grep("clone_", names(data[["cells"]])))]
  cat("\n\n> Merging duplicate clones",
      "\nSegments slot:",segs_slot,
      "\nClone slot:", clone_slot,"\n")
  
  cn_segs <- lapply(data[["clones"]][["cn"]], function(x) x[data[["segments"]][[segs_slot]][["start"]]])
  m <- do.call(cbind, cn_segs)
  colnames(m) <- data$clones$clone_id
  m <- m[,colnames(m)!="0"]
  
  # Compare full columns and group them
  dups <- duplicated(t(m)) | duplicated(t(m), fromLast=TRUE)
  if(!any(dups)){
    cat("No duplicate clones found\n")
    return(data)
  }
  clone_groups <- split(colnames(m)[dups], apply(m[,dups], 2, paste, collapse=","))
  
  # Create pairs of original and duplicates
  merge_table <- do.call(rbind, lapply(clone_groups, function(x) {
    data.frame(original_clone = x[1], 
               duplicate_clone = x,
               stringsAsFactors = FALSE)
  }))
  rownames(merge_table) <- NULL
  
  if(nrow(merge_table)>0){
    cat(">>> Found duplicate clones (includes self):\n")
    print(merge_table)
    
    # Update merged clone assignments
    new_clones.df <- merge_table %>% mutate(merged_clone=paste0(original_clone,"m"))
    new_clones.names <- unique(new_clones.df$merged_clone)
    merged_from <- split(new_clones.df$duplicate_clone, new_clones.df$merged_clone)
    
    clone_match <- match(data[["cells"]][[clone_slot]], merge_table$duplicate_clone)
    new_clones.vec <- sapply(1:length(clone_match), function(i) {
      ifelse(is.na(clone_match[i]), data[["cells"]][[clone_slot]][i], new_clones.df$merged_clone[clone_match[i]])
    })
    
    
    cat("Updated",sum(data[["cells"]][[clone_slot]] != new_clones.vec, na.rm = T),"cells\n")
    cat("Updated",length(new_clones.names),"clones\n")
    data[["cells"]]$clone_merged <- new_clones.vec 
    
    ### Update merged clone metrics
    # Create tibble with new clones (replaces old ones)
    clone_updates <- tibble(
      clone_id = new_clones.names,
      revision = "merged",
      parent_clone = map_chr(merged_from, ~paste(.x, collapse=";")),
      cells = map(clone_id, function(clone) data[["cells"]][["dna_library_id"]][which(data[["cells"]][["clone_merged"]] == clone)])
    )
    
    # Remove parent clones and add new split clones
    data[["clones"]] <- data[["clones"]] %>%
      filter(!clone_id %in% new_clones.df$duplicate_clone) %>%
      bind_rows(clone_updates)
    
    # Update normalization and segment metrics
    data <- normalize_counts(data, methods="gcmap", subset=new_clones.names)
    data <- update_clone_metrics(data, subset=new_clones.names, calc_cn=TRUE)
    
    return(data)
  }
}
fuzzy_merge_clones <- function(data, segs_slot=NULL, clone_slot=NULL, max_cn=6, max_diff_frac=0.05, max_residual=0.3, ...){
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])]
  if(is.null(clone_slot)) clone_slot <- names(data[["cells"]])[max(grep("clone_", names(data[["cells"]])))]
  cat("\n\n> Merging duplicate clones",
      "\nSegments slot:",segs_slot,
      "\nClone slot:", clone_slot,"\n")
  
  cn_segs <- lapply(data[["clones"]][["cn"]], function(x) x[data[["segments"]][[segs_slot]][["start"]]])
  res_segs <- data[["clones"]][["residuals"]]
  m <- do.call(cbind, cn_segs)
  r <- do.call(cbind, res_segs)
  colnames(m) <- data$clones$clone_id
  m <- m[,colnames(m)!="0"]
  # Mask values over max_cn
  m[m > max_cn] <- NA
  colnames(r)<-data$clones$clone_id
  r <- r[,colnames(r)!="0"]
  
  #If there are less than 10% of the rows that are different - then lets look at the residuals of those rows 
  similarity_check <- function(data, threshold) {
    colnames_data <- colnames(data)
    results <- list()
    
    for (i in 1:(ncol(data) - 1)) {
      for (j in (i + 1):ncol(data)) {
        col1 <- data[, i]
        col2 <- data[, j]
        
        # Compare ignoring NA values
        valid_idx <- !is.na(col1) & !is.na(col2)
        differing_rows <- which(col1 != col2)
        proportion_diff <- length(differing_rows) / sum(valid_idx) # Exclude NAs from the denominator
        
        if (proportion_diff <= threshold) {
          result_key <- paste(colnames_data[i], colnames_data[j], sep = "_vs_")
          results[[result_key]] <- list(
            differing_rows = differing_rows,
            proportion_diff = proportion_diff,
            columns = c(colnames_data[i], colnames_data[j])
          )
        }
      }
    }
    
    return(results)
  }
  
  
  # Run the similarity check
  result <- similarity_check(m, threshold = max_diff_frac)
  
  check_res_values <- function(result, res, threshold = max_residual) {
    valid_pairs <- lapply(result, function(x) {
      # Get column names and differing rows
      cols <- x$columns
      rows <- x$differing_rows
      
      # Extract the values from res for the given columns and rows
      values <- res[rows, cols, drop = FALSE]
      filtered_values <- as_tibble(values) %>%
        filter(rowSums(is.na(.)) < ncol(.))
      higher_values <- apply(filtered_values, 1, max, na.rm = TRUE)
      
      # Check if all values exceed the threshold
      if (all(higher_values > threshold)) {
        return(cols)
      } else {
        return(NULL)
      }
    })
    
    # Filter out NULLs
    valid_pairs <- valid_pairs[!sapply(valid_pairs, is.null)]
    
    # If there are valid pairs, pivot longer
    if (length(valid_pairs) > 0) {
      # Create a dataframe for valid pairs
      valid_df <- do.call(rbind, lapply(names(valid_pairs), function(pair_name) {
        cols <- valid_pairs[[pair_name]]
        data.frame(
          old_clones = cols,
          new_clones = paste0(cols[1], "ma"),
          pair_name = pair_name
        )
      }))
      # Return only the relevant columns
      return(valid_df[, c("old_clones", "new_clones")])
    } else {
      return(data.frame())  # Empty dataframe if no valid pairs
    }
  }
  
  output <- check_res_values(result, r, threshold = max_residual)
  
  if(nrow(output)==0){
    cat("No almost duplicate clones found\n")
    return(data)
  }
  
  if(nrow(output)>0){
    cat(">>> Found duplicate clones (includes self):\n")
    print(output)
    unique_output <- output %>%
      distinct(old_clones, .keep_all = TRUE)
    
    # View result
    output<-unique_output
    
    # Update merged clone assignments
    new_clones.names <- unique(output$new_clones)
    merged_from <- split(output$old_clones, output$new_clones)
    
    clone_match <- match(data[["cells"]][[clone_slot]], output$old_clones)
    new_clones.vec <- sapply(1:length(clone_match), function(i) {
      ifelse(is.na(clone_match[i]), data[["cells"]][[clone_slot]][i], output$new_clones[clone_match[i]])
    })
    
    
    cat("Updated",sum(data[["cells"]][[clone_slot]] != new_clones.vec, na.rm = T),"cells\n")
    cat("Updated",length(new_clones.names),"clones\n")
    data[["cells"]]$clone_merged2 <- new_clones.vec 
    
    ### Update merged clone metrics
    # Create tibble with new clones (replaces old ones)
    clone_updates <- tibble(
      clone_id = new_clones.names,
      revision = "merged2",
      parent_clone = map_chr(merged_from, ~paste(.x, collapse=";")),
      cells = map(clone_id, function(clone) data[["cells"]][["dna_library_id"]][which(data[["cells"]][["clone_merged2"]] == clone)])
    )
    
    # Remove parent clones and add new split clones
    data[["clones"]] <- data[["clones"]] %>%
      filter(!clone_id %in% output$old_clones) %>%
      bind_rows(clone_updates)
    
    # Update normalization and segment metrics
    data <- normalize_counts(data, methods="gcmap", subset=new_clones.names)
    data <- update_clone_metrics(data, subset=new_clones.names, calc_cn=TRUE, ...)
    
    return(data)
  }
}

remove_bad_clones <- function(data, clone_slot=NULL, frequency=0.1, ...){
  if(is.null(clone_slot)) clone_slot <- names(data[["cells"]])[max(grep("clone_", names(data[["cells"]])))]
  if(!is.null(data[["segments"]]$masked)){
    bad_clones <- data[["segments"]]$masked %>%
      filter(high_residual == FALSE) %>%
      select(high_residual_clones) %>%
      pull(high_residual_clones) %>%
      str_split(",") %>%
      unlist() %>%
      table() %>%
      as.data.frame() %>%
      filter(Freq > sum(Freq) *frequency ) %>% pull(1)%>%as.character()
    #Current clone:: 
    clone_to_zero<-data[["clones"]]$clone_id[data[["clones"]]$clone_id%in%bad_clones]
  }else {
    clone_to_zero<-NULL
  }
  if(length(clone_to_zero)==0){
    cat("No bad clones according to current clones \n")
    return(data)
  }
  cat("Bad clones to zero: ")
  print(clone_to_zero)
  
  new_clones.name<-data[["clones"]]$clone_id[!data[["clones"]]$clone_id%in%bad_clones]
  data[["cells"]][["clone_clean"]]<-mgsub(data[["cells"]][[clone_slot]], (clone_to_zero), c(rep("0", length(clone_to_zero))))
  
  #THen clone_updates to data somehow 
  data[["clones"]] <- data[["clones"]] %>%
    filter(clone_id %in% new_clones.name) 
  return(data)
}

refine_segments_from_cn <- function(data, new_slot="refined", update_clones=TRUE) {
  bins <- data[["bins"]][["good"]]
  cn_bins <- do.call(cbind, data[["clones"]][["cn"]])
  
  cat("\n\n>>> Creating minimal segments from change points in current cn-integer profiles\n")
  
  # Create chromosome arm factor levels
  chr_arm <- factor(paste(bins$chr, bins$arm, sep="_"), 
                    levels=paste0("chr", rep(c(1:22, "X", "Y"), each=2), c("_p", "_q")))
  chr_arm.l <- split(1:length(chr_arm), chr_arm)
  
  # Process each chromosome arm to find segments
  refined_segs.l <- lapply(chr_arm.l, function(bin_indices) {
    row_states <- apply(cn_bins[bin_indices,], 1, paste, collapse=",")
    breaks <- which(row_states[-length(row_states)] != row_states[-1])
    
    if(length(breaks) == 0) {
      tibble(
        start = bin_indices[1],
        end = bin_indices[length(bin_indices)]
      )
    } else {
      tibble(
        start = bin_indices[c(1, breaks + 1)],
        end = bin_indices[c(breaks, length(row_states))]
      )
    }
  })
  refined_segs <- do.call(rbind, refined_segs.l)
  
  cat(sprintf("Found %s segments",nrow(refined_segs)),"\n")
  print(sapply(data[["segments"]],nrow))
  initial_segs <- names(data[["segments"]])[1]
  prev_segs <- names(data[["segments"]])[length(names(data[["segments"]]))]
  prev_segs_df <- data[["segments"]][[prev_segs]]
  
  # Find matching original segments
  seg_idx <- sapply(1:nrow(refined_segs), function(i) {
    orig_match <- which(data[["segments"]][[initial_segs]]$start <= refined_segs$start[i] & 
                          data[["segments"]][[initial_segs]]$end >= refined_segs$end[i])
    if(length(orig_match) == 1) orig_match else NA
  })
  
  # Create final segments data frame using bins information
  data[["segments"]][[new_slot]] <- data.frame(
    chr = bins$chr[refined_segs$start],
    arm = bins$arm[refined_segs$start],
    start.pos = bins$start[refined_segs$start],
    end.pos = bins$end[refined_segs$end],
    n.probes = refined_segs$end - refined_segs$start + 1,
    seg_idx = seg_idx,
    start = refined_segs$start,
    end = refined_segs$end,
    filtered = map_lgl(seg_idx, ~ifelse(is.na(.x), FALSE, prev_segs_df$filtered[.x == prev_segs_df$seg_idx])),
    merged_from = NA,
    high_residual = NA,
    high_residual_clones = NA
  )
  
  if(update_clones)  data <- update_clone_metrics(data, norm_slot = "gcmap_normal", segs_slot=new_slot, calc_cn=T)
  
  return(data)
}

recall_cells <- function(data, from="0", mad_cutoff=5, segs_slot=NULL, clone_slot=NULL, plot=TRUE, update_clones=TRUE, ...){
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])]
  if(is.null(clone_slot)) clone_slot <- names(data[["cells"]])[max(grep("clone_", names(data[["cells"]])))]
  
  if(any(!from %in% unique(data[["cells"]][[clone_slot]]))){
    stop(sprintf("From slot (%s) not matching available clones (%s)",clone_slot, paste(unique(data[["cells"]][[clone_slot]]),collapse=", ")))
  } 
  cat(sprintf("\n\n> Clone recall from: %s. Using segs from %s and clones from %s\n", paste(from,collapse=", "), segs_slot, clone_slot))
  
  segs <- data[["segments"]][[segs_slot]]
  seg_idx <- Map(`:`, segs$start, segs$end)
  
  # Calculate for ALL cells
  data$cells <- data$cells %>% 
    mutate(
      raw_counts_segs=map(raw_counts, function(.x){
        r <- .x[data$bins$all$good]
        sapply(seg_idx, function(idx) sum(r[idx]))
      })
    )
  
  data$clones <- data$clones %>% 
    mutate(
      raw_counts_segs = map(raw_counts, ~sapply(seg_idx, function(idx) sum(.x[idx])))
    )
  
  # Create cell and clone matrices
  clone_counts <- do.call(cbind, data$clones$raw_counts_segs)
  colnames(clone_counts) <- data$clones$clone_id
  clone_counts.n <- t(t(clone_counts) / colSums(clone_counts)) # OBS! Have to transpose since R does -first division
  
  cell_counts <- do.call(cbind, data$cells$raw_counts_segs)
  colnames(cell_counts) <- data$cells$dna_library_id[!sapply(data$cells$raw_counts_segs, is.null)]
  cell_counts.n <- t(t(cell_counts) / colSums(cell_counts)) # OBS! Have to transpose since R does -first division
  
  # Mask filtered segments
  na_segs_idx <- which(segs$filtered)
  segs.f <- segs[-na_segs_idx,]
  clone_counts.nf <- clone_counts.n[-na_segs_idx,]
  cell_counts.nf <- cell_counts.n[-na_segs_idx,]
  
  # Print dimensions
  cat(sprintf("Using %d segments from %d clones vs %d cells", nrow(clone_counts.nf), ncol(clone_counts.nf), ncol(cell_counts.nf)),"\n")
  
  # Correlation
  size_factors <- log2(segs.f$n.probes) / segs.f$n.probes
  recall_cor <- Rfast::dista(xnew=t(cell_counts.nf * size_factors), x=t(clone_counts.nf * size_factors), k=0, index=F, type="euclidean", square=T)
  rownames(recall_cor) <- colnames(cell_counts.nf)
  colnames(recall_cor) <- colnames(clone_counts.nf)
  recall_cor_scaled <- rank_to_prob_exp(recall_cor)
  
  # Key part: assign NA cells to unassigned (dna qc fail) or replicating (pass dna, but excluded from clone detection due to S/M phase)
  cell_match_idx <- match(rownames(recall_cor),data$cells$dna_library_id)
  cell_clones <- data[["cells"]][[clone_slot]][cell_match_idx]
  cell_clones_all <- ifelse(is.na(cell_clones), 
                            ifelse(data$cells$cell_phase[cell_match_idx] %in% c("G2M","S"), 
                                   "replicating", 
                                   "unassigned"), 
                            cell_clones)
  
  if(plot){
    clone_colors_col <- get_clones_col(c("0",colnames(recall_cor),"unassigned","replicating"))
    clone_colors_row <- clone_colors_col[cell_clones_all]
    ha_col <- HeatmapAnnotation(clone = colnames(recall_cor), col = list(clone = clone_colors_col))
    ha_row <- rowAnnotation(original_clone = cell_clones_all, col = list(original_clone = clone_colors_row))
    ht <- Heatmap(recall_cor_scaled, name="recall prob",
                  cluster_columns = F,
                  cluster_rows = T,
                  col = colorRamp2(c(0, 1), c("white", "black")),
                  show_column_names = T, 
                  show_row_names = F,
                  top_annotation = ha_col,
                  right_annotation = ha_row) 
    draw(ht)
  }
  
  # Summarize potential recalls
  recall_counts <- table(
    "old" = cell_clones_all,
    "new" = colnames(recall_cor)[unlist(apply(recall_cor, 1, function(x) {
      if(all(is.nan(x))) return(NA)
      which.min(x)
    }))]
  )
  # print(recall_counts)
  
  # Return data to cells df
  l.recall_dist <- asplit(recall_cor, 1)
  l.recall_prob <- asplit(recall_cor_scaled, 1)
  l.recall_dist_min <- sapply(l.recall_dist, min)
  l.recall_clone <- sapply(l.recall_dist, function(x) {
    # Added 250129: If it fails to find a distance, returns NaN
    r <- names(x)[which.min(x)]
    r[length(r)==0] <- NA
    return(r)
  })
  
  # recall_cutoff <- median(l.recall_dist_min) + (mad(l.recall_dist_min) * mad_cutoff)
  # Cutoff based recall from/to confidently assigned clones
  cell_clones_idx.good <- which(!is.na(cell_clones) & ! cell_clones %in% "0")
  recall_cutoff <- median(l.recall_dist_min[cell_clones_idx.good]) + (mad(l.recall_dist_min[cell_clones_idx.good]) * mad_cutoff)
  
  recall_df <- tibble(
    cell=rownames(recall_cor),
    clone_original=cell_clones,
    clone_original_annot=factor(cell_clones_all, 
                                levels=c(setdiff(unique(cell_clones_all), c("0","replicating","unassigned")), 
                                         c("0","replicating","unassigned"))),
    clone_recall=l.recall_clone,
    min_dist=l.recall_dist_min,
    pass_cutoff=min_dist < recall_cutoff,
    display_color = case_when(
      # clone_original_annot == clone_recall & !clone_original_annot %in% c("unassigned", "replicating") ~ "grey",
      clone_original_annot == clone_recall ~ "grey",
      TRUE ~ clone_original_annot
    )
  )
  
  recall_df_stats <- recall_df %>% 
    group_by(clone_original_annot) %>% 
    summarize(n=n(), pass=sum(pass_cutoff), fail=sum(!pass_cutoff))
  
  if(plot){
    p <- recall_df %>% 
      ggplot(aes(log10(min_dist))) + geom_density() + facet_wrap(~clone_original_annot, ncol=1, scales="free_y") +
      geom_vline(xintercept=log10(recall_cutoff), color="red", linetype="dotted") +
      geom_text(data=recall_df_stats, aes(label=paste0(pass,"/",n)), inherit.aes=F, 
                x=max(log10(recall_df$min_dist)), y=Inf, hjust=1, vjust=1.2) +
      theme_dntr() +
      labs(title=paste0(patient_id, " recall: distance to closest matching clone"),
           caption=paste0("cutoff: median + (MAD x ",mad_cutoff,")"))
    plot(p)
    # Dot plot
    matching_df <- recall_df %>% 
      filter(clone_original_annot == clone_recall & !clone_original_annot %in% c("unassigned", "replicating"))
    non_matching_df <- recall_df %>% 
      filter(clone_original_annot != clone_recall | clone_original_annot %in% c("unassigned", "replicating"))
    modified_colors <- clone_colors_col
    modified_colors["correct"] <- "grey70"  # Add grey for correct assignments
    clone_indicators <- data.frame(
      clone = unique(recall_df$clone_recall),
      xstart = as.numeric(factor(unique(recall_df$clone_recall))) - 0.4,
      xend = as.numeric(factor(unique(recall_df$clone_recall))) + 0.4
    )
    p2 <- ggplot(recall_df, aes(x = clone_recall, y = min_dist, color = display_color)) +
      geom_quasirandom(data = matching_df,
                       aes(x = clone_recall, 
                           y = min_dist,
                           color = display_color),
                       alpha = 0.6, size = 0.8, width = 0.4) +
      # Plot non-matching points on top
      geom_quasirandom(data = non_matching_df,
                       aes(x = clone_recall, 
                           y = min_dist,
                           color = display_color),
                       alpha = 0.8,
                       size = 1.2,
                       width = 0.4) +
      geom_segment(data = clone_indicators, aes(x = xstart, xend = xend, y = min(recall_df$min_dist)*0.9, yend = min(recall_df$min_dist*0.9), color = clone),size = 3, show.legend=F) +  # Make the bars thick
      geom_hline(yintercept=recall_cutoff, color="red", linetype="dashed") +
      scale_color_manual(values = modified_colors) +
      scale_y_log10() +  # since distances are very small numbers
      labs(x = "Recalled Clone", y = "log10(min distance)", color = "Original Annotation") +
      theme_dntr(axis_labels = T, axis_ticks = T, legend="right") + labs(x=NULL) +
      theme(plot.margin=margin(b=0), axis.text = element_blank()) +
      guides(color=guide_legend(override.aes=list(size=4), title.position = "top", title = "Original clone"))
    
    # Create proportion data for marginal plot using only passing cells
    proportion_df <- recall_df %>%
      filter(pass_cutoff) %>%
      group_by(clone_recall) %>%
      count(display_color) %>%  # Use display_color to maintain grey for matching
      group_by(clone_recall) %>%
      mutate(prop = n/sum(n))
    
    p_marginal <- ggplot(proportion_df, aes(x = clone_recall, y = prop, fill = display_color)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = modified_colors) +
      theme_dntr(axis_ticks = T, axis_labels = T, legend="none") +
      theme(plot.margin = margin(t = 0)) +
      labs(y=NULL, x="Recalled clone")
    
    plot(p2/p_marginal + plot_layout(heights = c(3, 1)))
  }
  
  match_idx <- match(rownames(recall_cor), data$cells$dna_library_id)
  data$cells$recall_dist <- list(NULL)
  data$cells$recall_prob <- list(NULL)
  data$cells$clone_recall <- NA
  
  data$cells$recall_dist[match_idx] <- l.recall_dist
  data$cells$recall_prob[match_idx] <- l.recall_prob
  recall_df.f <- recall_df %>% 
    mutate(clone_recall=case_when(
      pass_cutoff & clone_original %in% from ~ clone_recall,
      TRUE ~ clone_original
    ))
  data$cells$clone_recall[match_idx] <- recall_df.f$clone_recall
  
  cat("Final recall:\n")
  print(table("original"=data$cells[[clone_slot]], "recall"=data$cells$clone_recall,useNA="always"))
  
  if(update_clones) data <- update_clone_metrics(data, clone_slot=clone_slot, calc_cn=TRUE, ...)
  
  return(data)
}

rank_to_prob_exp <- function(distances, decay_rate = 1) {
  # Get ranks
  ranks <- t(apply(distances, 1, rank))
  
  # Convert ranks to probabilities with exponential decay
  probs <- t(apply(ranks, 1, function(x) {
    p <- exp(-decay_rate * (x - 1))  # Subtract 1 so best rank has no decay
    p/sum(p)  # Normalize to sum to 1
  }))
  return(probs)
}

calc_cell_cn <- function(data, cell_idx=NULL, filter_segs=FALSE, cells_slot="gcmap_normal", segs_slot=NULL, clone_slot=NULL){
  try(RhpcBLASctl::blas_set_num_threads(1))
  try(RhpcBLASctl::omp_set_num_threads(1))
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])]
  if(is.null(clone_slot)) clone_slot <- names(data[["cells"]])[max(grep("clone_", names(data[["cells"]])))]
  stopifnot("cells_slot should be gcmap or gcmap_normal"=
              cells_slot %in% c("gcmap","gcmap_normal"), 
            "segs_slot not found"=
              segs_slot %in% names(data[["segments"]]),
            "clone slot not found"=
              clone_slot %in% names(data[["cells"]]))
  
  segs <- data[["segments"]][[segs_slot]]
  bins <- data[["bins"]][["good"]]
  if((length(cell_idx)==1 & is.logical(cell_idx)) | is.null(cell_idx)) cell_idx <- rep(TRUE, nrow(data[["cells"]]))
  clone_idx <- norm_idx <- cn_idx <- cell_idx
  cat(sprintf("Running %s cells from custom index", sum(cell_idx)),"\n")
  # }
  
  cat(sprintf("Normalizing (%s) %d single cells...\n", cells_slot, sum(norm_idx)))
  data[["cells"]][[cells_slot]][norm_idx] <- mclapply(data[["cells"]][["raw_counts"]][norm_idx], function(x){
    x <- x[bins$id]
    gc.coefs <- summary(lm(x ~ bins$gc + I(bins$gc^2)))$coefficients[,1]
    r <- gc.map.correct(x, bins$gc, bins$map, gc.coefs)
    rm <- r/mean(r)
    if(cells_slot == "gcmap_normal") rm <- rm/bins$normal_factors_gcmap
    return(rm)
  }, mc.cores=threads)
  
  cat(sprintf("Calculating scale factors and cn integers in %d single cells...\n", sum(cn_idx)))
  data[["cells"]][["segment_means"]][cn_idx] <- lapply(data[["cells"]][[cells_slot]][cn_idx], function(x){
    r <- sapply(1:nrow(segs), function(i) {
      idx <- segs[["start"]][i]:segs[["end"]][i]
      mean(x[idx])
    })
    if(filter_segs) r[segs$filtered] <- NA
    r
  })
  data[["cells"]][["segment_bins"]][cn_idx] <- lapply(data[["cells"]][["segment_means"]][cn_idx], function(x) rep(x, segs$n.probes))
  # data[["cells"]][["scale_tests"]][cn_idx] <- mclapply(data[["cells"]][["segment_bins"]][cn_idx], function(x){
  #   sapply(s, function(sf) manhattan.dist(sf, x, na.rm = T))
  # }, mc.cores=threads)
  # data[["cells"]][["scale_factor"]][cn_idx] <- lapply(data[["cells"]][["scale_tests"]][cn_idx], function(x) s[which.min(x)])
  # 
  
  for (i in which(cn_idx)) {
    val <- data[["cells"]]$correct_scalefactor[i]
    if (!is.na(val)) {
      data[["cells"]]$scale_factor[[i]] <- val
    }
  }
  
  data[["cells"]][["cn"]][cn_idx] <- lapply(which(cn_idx), function(i){
    round(data[["cells"]][["segment_bins"]][[i]] * data[["cells"]][["scale_factor"]][[i]])
  })
  return(data)
}



plot_clone_heatmap <- function(data, cn_slot="cn", segs_slot=NULL, bin_level=FALSE, highlight_dups=TRUE, clone_types=NULL, show_chr=FALSE, only_divergent=FALSE, order=TRUE, ...){
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])]
  
  # Create matrix at bin or segment level
  if(bin_level){
    # Expand to include gaps
    cn_bins.l <- lapply(data[["clones"]][[cn_slot]], function(x) expand_gaps_vec(x, length.out=nrow(data$bins$all), idx=data$bins$good$id))
    m <- do.call(cbind, cn_bins.l)
    colnames(m) <- data$clones$clone_id
    if(show_chr) chr_labels <- data$bins$all$chr_int else chr_labels <- NULL
  } else {
    cn_segs.l <- lapply(data[["clones"]][[cn_slot]], function(x) x[data[["segments"]][[segs_slot]]$start])
    m <- do.call(cbind, cn_segs.l)
    colnames(m) <- data$clones$clone_id 
    if(show_chr) chr_labels <- data$segments[[segs_slot]]$chr else chr_labels <- NULL
  }
  
  is_dup = FALSE # Do not color by default
  if(highlight_dups) is_dup <- duplicated(t(m)) | duplicated(t(m), fromLast = TRUE)
  
  if(only_divergent){
    identical_segs <- get_identical_segs(m)
    m <- m[!identical_segs,]
    chr_labels <- chr_labels[!identical_segs]
  } 
  
  column_anno=NULL
  if (!is.null(clone_types)) {
    if (identical(clone_types, TRUE)) { 
      clone_types <- get_clone_types(data, diploid = NULL, bin_level = bin_level)
    } else if (!is.character(clone_types)) { 
      stop("clone_types must be NULL, TRUE, or a character vector") 
    }
    if (only_divergent) clone_types <- clone_types[!identical_segs]
    if (nrow(m) != length(clone_types)) {
      stop(sprintf("Length mismatch: clonality (%d) differs from number of segments (%d)", length(clone_types), nrow(m)))
    }
    
    column_anno = HeatmapAnnotation(
      seg_type = clone_types,
      col = list(seg_type = c("clonal" = "darkred",
                              "subclonal" = "darkblue",
                              "none" = "white")),
      annotation_name_side = "left",
      show_legend = TRUE,
      annotation_legend_param = list(seg_type = list(title = NULL))
    )
  }
  if(all(is.character(order)) & length(order)==length(colnames(m))){
    clone_order <- order
  } else if(order==TRUE){
    clone_order <- sort_from_diploid_root(m, na.rm = T, ...) 
  } else {
    clone_order <- colnames(m)
  }
  
  label_colors <- ifelse(is_dup, "red", "black")
  if (grepl("cn", cn_slot)) {
    col_fun <- make_cn_colorscale(max(m, na.rm = TRUE))  
    
  } else if (grepl("ab_copy", cn_slot)) {
    # Helper to convert R color names to hex
    col2hex <- function(cname) {
      rgb(t(col2rgb(cname)) / 255)
    }
    
    # Get unique genotypes (as strings)
    unique_genotypes <- na.omit(unique(as.character(m)))
    
    # Define color ramps for A/B genotypes
    blue_ramp <- colorRampPalette(c("#b3cde0", "#1E90FF", "#1C66CC"))(2)
    green_ramp <- c("darkgreen", "#a8d08d")  # Custom green for AA and BB
    warm_ramp <- colorRampPalette(c("firebrick", "orangered", "orange", "gold"))(10)
    mixed_ramp <- colorRampPalette(c("darkred", "orangered", "gold", "yellow", "red"))(5)
    
    # Assign colors based on the genotype values
    genotype_colors <- sapply(unique_genotypes, function(gt) {
      len <- nchar(gt)
      if (gt == "AA") {
        return("darkgreen")  # AA gets dark green
      } else if (gt == "BB") {
        return("#a8d08d")  # BB gets light green
      } else if (gt == "AB") {
        return("#f0f0f0")  # AB gets light grey
      } else if (gt == "AAB") {
        return("darkred")  # AAB gets dark red
      } else if (gt == "ABB") {
        return("orangered")  # ABB gets orange-red
      } else if (len < 2) {
        # For genotypes of length 1 (rare case), use blue shades
        blue_ramp[pmin(len + 1, length(blue_ramp))]
      } else if (len == 2) {
        # Pick a nice green for AB and similar
        green_ramp[sample(length(green_ramp), 1)]
      } else if (len == 3) {
        # Mixed genotypes like AAB or ABB
        mixed_ramp[sample(length(mixed_ramp), 1)]
      } else {
        warm_ramp[pmin(len - 2 + 1, length(warm_ramp))]
      }
    }, USE.NAMES = TRUE)
    
    col_fun <- genotype_colors
  } else if (grepl("final_ab", cn_slot)) {
    # Define the color ramps for final_ab values
    blue_ramp <- colorRampPalette(c("#b3cde0", "#1E90FF", "#1C66CC"))(2)
    green_ramp <- c("darkgreen", "#a8d08d")  # Custom green for AA and BB
    warm_ramp <- colorRampPalette(c("firebrick", "orangered", "orange", "gold"))(10)
    mixed_ramp <- colorRampPalette(c("darkred", "orangered", "gold", "yellow", "red"))(5)
    
    # Get unique genotypes (as strings) from the final_ab column
    unique_genotypes <- na.omit(unique(as.character(m)))
    
    # Assign colors based on the values of final_ab
    final_ab_colors <- sapply(unique_genotypes, function(gt) {
      len <- nchar(gt)
      if (grepl("A", gt) | grepl("B", gt)) {
        if (gt == "AA") {
          return("darkgreen")  # AA gets dark green
        } else if (gt == "BB") {
          return("#a8d08d")  # BB gets light green
        } else if (gt == "AB") {
          return("#f0f0f0")  # AB gets light grey
        } else if (gt == "AAB") {
          return("darkred")  # AAB gets dark red
        } else if (gt == "ABB") {
          return("orangered")  # ABB gets orange-red
        } else if (len < 2) {
          # For genotypes of length 1 (rare case), use blue shades
          blue_ramp[pmin(len + 1, length(blue_ramp))]
        } else if (len == 2) {
          # Pick a nice green for AB and similar
          green_ramp[sample(length(green_ramp), 1)]
        } else if (len == 3) {
          # Mixed genotypes like AAB or ABB
          mixed_ramp[sample(length(mixed_ramp), 1)]
        } else {
          warm_ramp[pmin(len - 2 + 1, length(warm_ramp))]
        }
      } else {
        # For non-A/B genotypes (numeric or others), use numeric color scale
        return(make_cn_colorscale(max(m, na.rm = TRUE))(as.numeric(gt)))
      }
    }, USE.NAMES = TRUE)
    
    col_fun <- final_ab_colors
  } else {
    # Default color scale for other cases
    col_fun <- circlize::colorRamp2(
      c(0, 0.2, 0.33, 0.5),
      c("#2c7fb8", "#7fcdbb", "#edf8b1", "white")
    )
  }
  clone_slot <- names(data[["cells"]])[max(grep("clone_", names(data[["cells"]])))]
  title <- paste(paste(names(data$parameters),data$parameters,sep=": "),collapse=", ")
  n_segs <- ifelse(only_divergent, paste(nrow(m), "divergent"), nrow(data[["segments"]][[segs_slot]]))
  subtitle <- paste0("segs_slot: ",segs_slot,
                     ", clones: ",clone_slot,
                     ", ",nrow(data$clones)," clones, ",
                     n_segs, " segments")
  if(sum(is_dup)>0) subtitle <- paste0(subtitle,", ",sum(is_dup)," duplicates")
  ht <- Heatmap(t(m),cluster_columns = F, col=col_fun, name = cn_slot,
                column_split = chr_labels,
                column_gap = unit(0, "mm"),
                column_title_gp=gpar(fontsize=8),
                column_title_side = "bottom",
                row_split = NULL,
                row_gap = unit(0, "mm"),
                border_gp = gpar(col = "black", lwd = 0.5),
                row_order = clone_order,
                cluster_rows=F,
                border=T,
                row_names_gp = gpar(col = label_colors, fontsize=8),
                row_names_side = "left",
                bottom_annotation = column_anno,
                row_title = "clones"
  )
  draw(ht, 
       column_title=paste(title,subtitle,sep="\n"), 
       column_title_side = "top", 
       column_title_gp=gpar(fontsize=10))
}


plot_clone_revisions <- function(data, clone_columns=NULL, show_all = TRUE) {
  if(is.null(clone_columns)) clone_columns <- grep("clone_", names(data), value = TRUE)
  # Get flows between consecutive clone assignments
  flows <- data.frame(source = data[[clone_columns[1]]], 
                      target = data[[clone_columns[2]]])
  
  # Create flow table and convert to data frame
  flow_table <- as.data.frame(table(flows)) %>%
    filter(Freq > 0) %>%
    # Only show changes if show_all is FALSE
    filter(show_all | as.character(source) != as.character(target))
  
  # Summarize target counts for merged endpoints
  target_counts <- flow_table %>%
    group_by(target) %>%
    summarize(total_freq = sum(Freq))
  
  # Create plot
  p <- ggplot(flow_table, aes(x = 1, xend = 2, 
                              y = reorder(source, Freq),
                              yend = target)) +
    geom_segment(aes(size = Freq), color = "steelblue", alpha = 0.6) +
    # Left side labels (sources)
    geom_text(aes(label = paste0(source, " (", Freq, ")")), x = 1, hjust = 1) +
    # Right side labels (targets) - only show once with total count
    geom_text(data = target_counts,
              aes(x = 2, y = target, 
                  label = paste0(target, " (", total_freq, ")")), 
              hjust = -0.1) +
    scale_size_continuous(range = c(0.5, 3)) +
    scale_x_continuous(limits = c(0.5, 2.5),
                       breaks = c(1, 2),
                       labels = c(clone_columns[1], clone_columns[2])) +
    labs(title = "Clone Assignment Changes",
         x = "Clone revisions") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
  
  cat("\nClone revisions:", paste(clone_columns, collapse = " → "), "\n")
  cat("\nFlow summary:\n")
  print(flow_table)
  
  return(p)
}

plot_cell_heatmap <- function(data, filtered=TRUE, clone_slot=NULL, segs_slot=NULL, 
                              smooth_bins=50, annotate=NULL, annotate_colors=NULL, 
                              group_by=NULL, region=NULL){
  if(is.null(clone_slot)) clone_slot <- names(data[["cells"]])[max(grep("clone_", names(data[["cells"]])))]
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])] 
  bins <- data$bins$good
  bins$chr_short <- sub("chr","",bins$chr)
  clone_idx = TRUE
  if(filtered) clone_idx <- !is.na(data[["cells"]][[clone_slot]]) 
  
  cat("\nUsing clone slot:",clone_slot,"\n")
  
  recall_cn <- do.call(cbind, data$cells$cn[clone_idx])
  colnames(recall_cn) <- data$cells$dna_library_id[clone_idx]
  
  # Parse and apply region filtering if specified
  if(!is.null(region)) {
    parsed_region <- parse_region(region)
    if(!is.null(parsed_region)) {
      # Get indices for the specified region
      region_idx <- switch(parsed_region$type,
                           "chromosome" = bins$chr == parsed_region$chr,
                           "arm" = bins$chr == parsed_region$chr & 
                             substr(bins$arm, 1, 1) == parsed_region$arm,
                           "position" = bins$chr == parsed_region$chr & 
                             bins$start >= parsed_region$start & 
                             bins$end <= parsed_region$end
      )
      
      if(!any(region_idx)) {
        stop("No bins found in specified region: ", region)
      }
      
      # Subset the data
      smooth <- smooth_bins_subset(recall_cn[region_idx,], bins=bins[region_idx,], w=smooth_bins)
      ht_mtx <- smooth$mat
      filt_segs.all <- rep(data[["segments"]][[segs_slot]]$filtered, data[["segments"]][[segs_slot]]$n.probes)
      filt_segs <- filt_segs.all[region_idx]
      filt_segs.smooth <- smooth_vector(filt_segs, smooth_bins)
    }
  } else {
    smooth <- smooth_bins_subset(recall_cn, bins=bins, w=smooth_bins)
    ht_mtx <- smooth$mat
    filt_segs <- rep(data[["segments"]][[segs_slot]]$filtered, data[["segments"]][[segs_slot]]$n.probes)
    filt_segs.smooth <- smooth_vector(filt_segs, smooth_bins)
  }
  
  # Build annotation: convert initially NA clones to annotated sets (for plotting)
  meta <- data[["cells"]][clone_idx,]
  meta <- meta %>%
    mutate(
      pass_dna_qc = !is.na(clone_initial),
      across(starts_with("clone_"), 
             ~case_when(
               !is.na(.) ~ .,
               pass_dna_qc == FALSE ~ "fail_qc"
             )))
  
  # Ordering by profile
  non_clones <- c("0","S","G2M","fail_qc")
  clone_list <- split(meta$dna_library_id, meta[[clone_slot]])
  clones_incl <- clone_list[!names(clone_list) %in% non_clones]
  clones_excl <- clone_list[!names(clone_list) %in% names(clones_incl)]
  clone_medians <- lapply(clones_incl, function(i) 
    if(length(i)==1) median(recall_cn[,i]) else Rfast::rowMedians(recall_cn[,i])
  )
  clone_medians.mtx <- do.call(cbind, clone_medians)
  o <- sort_from_diploid_root(clone_medians.mtx)
  meta[[clone_slot]] <- droplevels(factor(meta[[clone_slot]], levels=c(o, non_clones)))
  
  col_fun <- make_cn_colorscale(max(ht_mtx, na.rm=T))
  chr_labels <- factor(smooth$bins$chr_short, levels=c(1:22,"X","Y"))
  clone_vars <- grep("clone_", colnames(meta), value=TRUE)
  clone_colors <- get_clones_col(unique(unlist(lapply(meta[, clone_vars], as.character))))
  names(clone_colors) <- as.character(names(clone_colors))
  recall_data <- NULL
  if("recall_prob" %in% colnames(data$cells)){
    recall_annot <- TRUE
    recall_prob <- do.call(rbind, data$cells$recall_prob[clone_idx])
    recall_colors <- clone_colors[match(colnames(recall_prob), names(clone_colors))] 
    recall_data <- row_anno_barplot(recall_prob, gp = gpar(col=recall_colors, fill = recall_colors))
  }
  
  # Dynamically add annotations
  annot_data <- sapply(clone_vars, function(col) as.character(meta[[col]]), simplify=FALSE)
  annot_colors <- rep(list(clone_colors), length(clone_vars))
  names(annot_colors) <- clone_vars
  if(!is.null(annotate) && annotate %in% colnames(meta)) {
    annot_data[[annotate]] <- meta[[annotate]]
    if(!is.null(annotate_colors)){
      annot_colors[[annotate]] <- annotate_colors
    }
  }
  
  ha <- rowAnnotation(
    df = annot_data,
    col = annot_colors,
    recall = recall_data,
    reads=anno_barplot(log10(meta$bam_read_pairs), bar_width=1, border=F, baseline = "min"),
    show_legend=T,
    annotation_name_gp = gpar(fontsize=8),
    na_col = "white"
  )
  
  column_anno = HeatmapAnnotation(
    masked = filt_segs.smooth,
    col = list(masked = c("TRUE" = "red", "FALSE"=NA)),
    show_annotation_name=FALSE
  )
  
  if(!is.null(group_by)){
    if(!group_by %in% colnames(meta)) stop(paste0("Grouping variable ",group_by, "not in metadata."))
    group_var <- meta[[group_by]]
  } else {
    group_var <- meta[[clone_slot]]
  }
  
  # Add 250129: Empty cells need to be set to NA
  zero_cells <- colSums(is.na(ht_mtx)) == nrow(ht_mtx)
  if(sum(zero_cells) > 0){
    cat(sprintf("%d cells with NA counts excluded: %s", sum(zero_cells), paste(names(which(zero_cells)), collapse=", ")),"\n")
    ht_mtx[is.na(ht_mtx)] <- -1
  }
  
  region_title <- if(!is.null(region)) paste0(", region: ", region) else ""
  
  ht = Heatmap(t(ht_mtx), name="cn", col=col_fun,
               heatmap_legend_param=list(title="cn",col_fun=col_fun, break_dist=1),
               cluster_columns = F,
               row_split = group_var,
               cluster_rows= T,
               cluster_row_slices = F,
               row_dend_reorder = F,
               show_row_dend = T,
               column_split = chr_labels,
               row_gap = unit(2, "mm"),
               column_gap = unit(0, "mm"),
               border=T,
               row_title_rot = 90,
               row_title_gp=gpar(fontsize=8),
               column_title_gp=gpar(fontsize=8),
               column_title_side = "bottom",
               show_row_names = F,
               show_column_names = F,
               right_annotation = ha,
               bottom_annotation = column_anno,
               use_raster = F,
               show_heatmap_legend=F,
               na_col="#cccccc")
  
  ht <- draw(ht)
}

parse_region <- function(region) {
  # Check and standardize the region input
  if (is.null(region)) return(NULL)
  
  if (grepl("^chr[0-9XY]{1,2}[pq]$", region)) {
    # Chromosome arm format: chr12p, chrXp
    chr <- sub("([pq])$", "", region)
    arm <- substr(region, nchar(region), nchar(region))
    return(list(type="arm", chr=chr, arm=arm))
    
  } else if (grepl("^chr[0-9XY]{1,2}:\\d+\\-\\d+[kKmMbB]*$", region)) {
    # Position range format: chr12:0-12Mb, chrX:0-12Mb
    chr <- sub(":.*$", "", region)
    range_str <- sub("^.*:", "", region)
    positions <- strsplit(range_str, "-")[[1]]
    
    convert_position <- function(pos) {
      multiplier <- switch(toupper(substr(pos, nchar(pos), nchar(pos))),
                           "K" = 1e3,
                           "M" = 1e6,
                           "B" = 1e9,
                           1)
      if (multiplier != 1) pos <- substr(pos, 1, nchar(pos)-1)
      as.numeric(pos) * multiplier
    }
    
    start_pos <- convert_position(positions[1])
    end_pos <- convert_position(positions[2])
    
    return(list(type="position", chr=chr, start=start_pos, end=end_pos))
    
  } else if (grepl("^chr[0-9XY]{1,2}$", region)) {
    # Whole chromosome format: chr12, chrX
    return(list(type="chromosome", chr=region))
    
  } else {
    stop("Invalid region format. Use 'chr12', 'chrX', 'chr12p', 'chrXp', or 'chr12:0-12Mb' format")
  }
}

get_segment_highlights <- function(segments, end_cum) {
  seg.grid.labs <- segments %>% 
    mutate(
      start_bp_cum = end_cum[start],
      end_bp_cum = end_cum[end],
      seg_label_pos = (start_bp_cum + end_bp_cum) / 2
    )
  
  cuts <- seg.grid.labs %>%
    select(chr, start_bp_cum, end_bp_cum, seg_idx, filtered) %>% 
    pivot_longer(cols = c(start_bp_cum, end_bp_cum), 
                 names_to = "border", 
                 values_to = "bp") %>%
    arrange(chr, bp) %>%
    group_by(chr) %>%
    mutate(
      drop = case_when(
        lag(filtered) == TRUE ~ TRUE,
        lead(filtered) == TRUE ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(filtered | !drop) %>%
    ungroup() %>%
    arrange(chr, bp)
  
  cuts.labs <- seg.grid.labs %>% 
    select(chr, seg_idx, seg_label_pos, filtered)
  
  list(cuts = cuts, cuts.labs = cuts.labs)
}

plot_clone_detail <- function(data, show_raw=TRUE, show_norm=TRUE, show_segs=TRUE,
                              cn.trunc=4, labels=NULL, region=NULL, clone_ids=NULL,
                              segs_slot=NULL, order=TRUE) {
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])] 
  if(!segs_slot %in% names(data[["segments"]])) stop("Specified segment slot not found")
  
  # Subset clones if specified
  if (!is.null(clone_ids)) {
    clone_idx <- which(data$clones$clone_id %in% clone_ids)
    if (length(clone_idx) == 0) stop("No matching clone IDs found")
  } else {
    clone_idx <- seq_len(nrow(data$clones))
  }
  
  # Create the full plotting dataframe first
  df.list <- lapply(clone_idx, function(i) {
    clone <- data$clones$clone_id[i]
    
    tibble(
      bin = seq_len(nrow(data$bins$all)),
      bp = data$bins$all$end_cum,
      chr = data$bins$all$chr,
      arm = data$bins$all$arm,
      bp_chr_start = data$bins$all$start,
      bp_chr_end = data$bins$all$end,
      raw = expand_gaps_vec(data$clones$raw_counts[[i]], 
                            length.out=nrow(data$bins$all), 
                            idx=data$bins$good$id) * 
        data$clones$scale_factor[[i]],
      gc = expand_gaps_vec(data$clones$gcmap_normal[[i]], 
                           length.out=nrow(data$bins$all), 
                           idx=data$bins$good$id) * 
        data$clones$scale_factor[[i]],
      seg = expand_gaps_vec(data$clones$segment_bins[[i]], 
                            length.out=nrow(data$bins$all), 
                            idx=data$bins$good$id) * 
        data$clones$scale_factor[[i]],
      cn = expand_gaps_vec(data$clones$cn[[i]], 
                           length.out=nrow(data$bins$all), 
                           idx=data$bins$good$id),
      n_cells = length(data$clones$cells[[i]])
    )
  })
  
  names(df.list) <- data$clones$clone_id[clone_idx]
  
  # Combine all data
  df <- bind_rows(df.list, .id="clone")
  if(order){
    cn_segs.l <- lapply(data[["clones"]]$cn, function(x) x[data[["segments"]][[segs_slot]]$start])
    m <- do.call(cbind, cn_segs.l)
    colnames(m) <- data$clones$clone_id
    df$clone <- factor(df$clone, levels=sort_from_diploid_root(m, na.rm = T))
  } else {
    df$clone <- factor(df$clone, levels=names(df.list))
  }
  
  df$bin <- as.integer(df$bin)
  
  # Create chromosome annotation dataframe
  df.chr <- data$bins$all %>%
    group_by(chr) %>%
    summarise(
      line = max(end_cum),
      label_pos = mean(c(min(end_cum), max(end_cum))),
      chr_short = unique(chr)
    ) %>%
    ungroup()
  
  # Convert segment index to all bins
  good_bins_map <- setNames(data$bins$good$id, nm=1:nrow(data$bins$good)) # val = original bin index, name = filtered bin index
  good_end_cum <- data$bins$all$end_cum[good_bins_map]
  seg_highlights <- get_segment_highlights(
    segments=data$segments[[segs_slot]],
    end_cum=good_end_cum
  )
  
  # Filter for region if specified
  if (!is.null(region)) {
    region_info <- parse_region(region)
    
    if (region_info$type == "chromosome") {
      df <- df %>% filter(chr == region_info$chr)
      df.chr <- df.chr %>% filter(chr == region_info$chr)
      seg_highlights$cuts <- seg_highlights$cuts %>% filter(chr == region_info$chr)
      seg_highlights$cuts.labs <- seg_highlights$cuts.labs %>% filter(chr == region_info$chr)
      
    } else if (region_info$type == "arm") {
      df <- df %>% filter(chr == region_info$chr, arm == region_info$arm)
      df.chr <- df.chr %>% filter(chr == region_info$chr)
      seg_highlights$cuts <- seg_highlights$cuts %>% 
        filter(chr == region_info$chr, 
               bp %in% df$bp)  # Only keep cuts within the filtered bins
      seg_highlights$cuts.labs <- seg_highlights$cuts.labs %>% 
        filter(chr == region_info$chr,
               seg_label_pos %in% df$bp)
      
    } else if (region_info$type == "position") {
      df <- df %>% 
        filter(chr == region_info$chr,
               bp_chr_start <= region_info$end,
               bp_chr_end >= region_info$start)
      df.chr <- df.chr %>% filter(chr == region_info$chr)
      seg_highlights$cuts <- seg_highlights$cuts %>% 
        filter(chr == region_info$chr,
               bp >= region_info$start,
               bp <= region_info$end)
      seg_highlights$cuts.labs <- seg_highlights$cuts.labs %>%
        filter(chr == region_info$chr,
               seg_label_pos >= region_info$start,
               seg_label_pos <= region_info$end)
    }
  }
  
  # Set y-axis limits
  if(!is.null(cn.trunc)) {
    ymax <- max(min(cn.trunc + 0.5, max(df$cn, na.rm=TRUE) + 1), 4)
  } else {
    ymax <- max(df$cn, na.rm=TRUE) + 1
  }
  
  # Clone labeller (with counts)
  clone_counts <- sapply(df.list, function(x) unique(x$n_cells))
  clone_labels <- setNames(paste0(names(clone_counts), " (n=",clone_counts,")"), 
                           nm=names(clone_counts))
  
  # Create base plot
  p <- ggplot(df, aes(x=bp)) +
    geom_hline(yintercept=2, color="grey", linetype="solid", linewidth=0.1)
  
  if(show_raw) p <- p + geom_point(aes(y=raw), pch=19, size=0.1, alpha=0.1, color="grey")
  if(show_norm) p <- p + geom_point(aes(y=gc, color=as.factor(cn)), pch=19, size=0.1, alpha=0.5, show.legend=FALSE)
  if(show_segs) {
    #p <- p + geom_point(aes(y=seg), pch=19, size=0.15, col="orange", alpha=0.5) 
    # Add segment highlights
    p <- p + 
      geom_vline(data=seg_highlights$cuts, 
                 aes(xintercept=bp, linetype=filtered), 
                 alpha=0.5, 
                 color="red", 
                 show.legend=FALSE) + 
      # geom_text(data=seg_highlights$cuts.labs, 
      #           aes(x=seg_label_pos, y=0, label=seg_idx), 
      #           vjust=0, 
      #           color="darkgrey", 
      #           size=3) +
      geom_point(data=filter(seg_highlights$cuts.labs, filtered), 
                 aes(x=seg_label_pos, y=0), 
                 pch=4,
                 color="black",
                 size=1.2) +
      scale_linetype_manual(values=c("FALSE"="solid", "TRUE"="dotted"))
  }
  
  p <- p + 
    geom_point(aes(y=cn, color=as.factor(cn)), pch=19, size=0.5, show.legend=FALSE) +
    scale_color_manual(values=setNames(as.character(cut(0:200, 
                                                        breaks=c(-Inf,0,1,2,3,4,20,50,Inf), 
                                                        labels=c("darkblue","blue","grey","red","darkred","brown","green","yellow"))), 
                                       nm=0:200)) +
    scale_x_continuous(expand=c(0,0)) +
    lims(y=c(0, ymax)) +
    theme_dntr(axis_ticks=TRUE) +
    theme(
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank(),
      plot.margin=margin(5,5,0,5)
    ) +
    labs(x=NULL, y="copy number")
  
  if(!is.null(region)){
    p <- p + annotate(
      "text",
      x = -Inf,
      y = -Inf,
      label = region,
      hjust = -0.1,
      vjust = -0.5
    )
  } else {
    p <- p + geom_text(
      data=df.chr,
      aes(x=label_pos, y=0, label=chr_short),
      check_overlap=TRUE,
      hjust=0.5,
      vjust=-0.2
    )
  }
  
  # Add points above cn.trunc
  p <- p + geom_point(
    data=filter(df, cn > cn.trunc),
    aes(y=cn.trunc + 0.5, fill=as.factor(cn)),
    color="black",
    pch=24,
    size=2.5,
    show.legend=FALSE
  )
  
  # Handle faceting
  if(length(unique(df$chr)) > 1 && length(unique(df$chr)) < 24) {
    if(!is.null(labels)) {
      clone.anno.df <- tibble(
        clone = names(labels),
        label = labels,
        chr = unique(df$chr)[1]
      )
      p <- p + geom_text(
        data=clone.anno.df,
        aes(x=-Inf, y=Inf, label=label),
        inherit.aes=FALSE,
        hjust=-0.01,
        vjust=1.5
      )
    }
    p <- p +
      facet_grid(
        factor(clone, levels=levels(df$clone)) ~ chr,
        space="free",
        scales="free",
        labeller=labeller(clone=clone_labels)
      ) +
      theme(
        strip.text.x=element_blank(),
        panel.spacing.x=unit(0, "lines")
      )
  } else {
    p <- p +
      facet_wrap(
        vars(clone),
        ncol=1,
        labeller=labeller(clone=clone_labels),
        scales="free_x"
      ) 
    if(is.null(region)){
      p <- p + geom_vline(xintercept=df.chr$line, color="grey")
    }
  }
  
  plot(p)
}

remove_bad_cells <- function(data, max_diff_bins=NULL, min_seg_size=500, clone_slot=NULL, update_clones=T, segs_slot=NULL, ...){
  if(is.null(clone_slot)) clone_slot <- names(data[["cells"]])[max(grep("clone_", names(data[["cells"]])))]
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])]
  if(! "cn" %in% colnames(data[["cells"]])) stop("Cell-level cn profiles not in data. Run calc_cell_cn() first")
  
  # clone_idx <- !is.na(data[["cells"]][["clone_initial"]]) & data[["cells"]][["clone_initial"]] != "0"
  cat(sprintf("> Removing individual cells with > %s bins difference to clone median. \nClone slot: %s and segs slot: %s", max_diff_bins, clone_slot, segs_slot),"\n")
  clone_idx <- !is.na(data[["cells"]][[clone_slot]]) & data[["cells"]][[clone_slot]] != "0"
  
  clone_list <- split(data[["cells"]][["dna_library_id"]][clone_idx], data[["cells"]][[clone_slot]][clone_idx])
  print(sapply(clone_list, length))
  cells_to_zero <- list()
  segs <- data[["segments"]][[segs_slot]]
  bin_idx <- logical(max(segs$end))
  for(i in 1:nrow(segs)) {
    bin_idx[segs$start[i]:segs$end[i]] <- segs$n.probes[i] > 500
  }
  
  comparison_results <- lapply(names(clone_list), function(clone) {
    cat("Analyzing clone:",clone,"\n")
    lapply(clone_list[[clone]], function(cell){
      # cat(sprintf("Running %s against clone %s", cell, clone),"\n")
      # Compare cell to clone
      cell_cn <- data[["cells"]][["cn"]][[which(data$cells$dna_library_id==cell)]]
      clone_cn <- data[["clones"]][["cn"]][[which(data$clones$clone_id==clone)]]
      diff <- clone_cn[bin_idx] != cell_cn[bin_idx]
      # cat(sprintf("Cell-clone diff %s: %d", cell, sum(diff, na.rm=T)),"\n")
      if(sum(diff, na.rm=T) > max_diff_bins){
        cat(sprintf("Cell %s excluded from clone %s (diff: %d bins)", cell, clone, sum(diff, na.rm=T)),"\n")
        cells_to_zero <- append(cells_to_zero, (cell))
      }
    })
  })
  
  cells_to_move <- unlist(comparison_results) 
  data$cells <- data$cells %>%
    mutate(
      clone_final = if_else(
        .data$dna_library_id %in% cells_to_move, 
        "0",  # Set clone_final to "0" for cells to move
        .data[[clone_slot]]  # Otherwise, retain the value in clone_slot
      )
    )
  # data$cells$clone_final
  print(table("old"=data$cells[[clone_slot]], "new"=data$cells$clone_final))
  data <- update_clone_metrics(data, clone_slot=clone_slot, calc_cn=TRUE, ...)
  return(data)
}

remove_small_clones <- function(data, clone_slot=NULL,  min_size_clone=min_size_clone){
  if(is.null(clone_slot)) clone_slot <- names(data[["cells"]])[max(grep("clone_", names(data[["cells"]])))]
  clone_idx <- !is.na(data[["cells"]][["clone_initial"]])
  clone_list <- split(data[["cells"]][["dna_library_id"]][clone_idx], data[["cells"]][[clone_slot]][clone_idx])
  bad_clones<-names(clone_list)[sapply(clone_list, length) < min_size_clone]
  #Current clone:: 
  clone_to_zero<-data[["clones"]]$clone_id[data[["clones"]]$clone_id%in%bad_clones]
  if(length(clone_to_zero)==0){
    cat("No clone that is too small \n")
    return(data)
  }
  cat("Small clone - moved to zero: ")
  print(clone_to_zero)
  
  new_clones.name<-data[["clones"]]$clone_id[!data[["clones"]]$clone_id%in%bad_clones]
  data[["cells"]][["clone_clean2"]] <- mgsub::mgsub(data[["cells"]][[clone_slot]], (clone_to_zero), c(rep("0", length(clone_to_zero))))
  
  #THen clone_updates to data somehow 
  data[["clones"]] <- data[["clones"]] %>%
    filter(clone_id %in% new_clones.name) 
  return(data)
}

load_parameters <- function(patient_id, yaml_file, defaults = default_params) {
  params <- defaults
  params$patient_id <- patient_id
  
  yaml_data <- yaml::read_yaml(yaml_file)
  if (!patient_id %in% names(yaml_data)) {
    warning(sprintf("Patient %s not found in YAML file. Using default parameters.", patient_id))
    return(params)
  }
  
  # Update with patient specific parameters
  patient_params <- yaml_data[[patient_id]]
  for (param_name in names(patient_params)) {
    params[[param_name]] <- patient_params[[param_name]]
  }
  return(params)
}


plot_clone_heatmap_2 <- function(data, cn_slot="cn", segs_slot=NULL, bin_level=FALSE, highlight_dups=TRUE, clone_types=NULL, show_chr=FALSE, only_divergent=FALSE, order=TRUE, ...){
  if(is.null(segs_slot)) segs_slot <- names(data[["segments"]])[length(data[["segments"]])]
  
  # Create matrix at bin or segment level
  if(bin_level){
    # Expand to include gaps
    cn_bins.l <- lapply(data[["clones"]][[cn_slot]], function(x) expand_gaps_vec(x, length.out=nrow(data$bins$all), idx=data$bins$good$id))
    m <- do.call(cbind, cn_bins.l)
    colnames(m) <- data$clones$clone_id
    if(show_chr) chr_labels <- data$bins$all$chr_int else chr_labels <- NULL
  } else {
    cn_segs.l <- lapply(data[["clones"]][[cn_slot]], function(x) x[data[["segments"]][[segs_slot]]$start])
    m <- do.call(cbind, cn_segs.l)
    colnames(m) <- data$clones$clone_id 
    if(show_chr) chr_labels <- data$segments[[segs_slot]]$chr else chr_labels <- NULL
  }
  
  is_dup = FALSE # Do not color by default
  if(highlight_dups) is_dup <- duplicated(t(m)) | duplicated(t(m), fromLast = TRUE)
  
  if(only_divergent){
    identical_segs <- get_identical_segs(m)
    m <- m[!identical_segs,]
    chr_labels <- chr_labels[!identical_segs]
  } 
  
  #Here we could choose to drop x y 
  # if(only_divergent){ so instead we could remove the segments from x and y  - if chr_labels%in%c("chrX", "chrY") then we can remove those and then update m and chr_labels 
  #   identical_segs <- get_identical_segs(m)
  #   m <- m[!identical_segs,]
  #   chr_labels <- chr_labels[!identical_segs]
  # } 
  
  column_anno=NULL
  if (!is.null(clone_types)) {
    if (identical(clone_types, TRUE)) { 
      clone_types <- get_clone_types(data, diploid = NULL, bin_level = bin_level)
    } else if (!is.character(clone_types)) { 
      stop("clone_types must be NULL, TRUE, or a character vector") 
    }
    if (only_divergent) clone_types <- clone_types[!identical_segs]
    if (nrow(m) != length(clone_types)) {
      stop(sprintf("Length mismatch: clonality (%d) differs from number of segments (%d)", length(clone_types), nrow(m)))
    }
    
    column_anno = HeatmapAnnotation(
      seg_type = clone_types,
      col = list(seg_type = c("clonal" = "darkred",
                              "subclonal" = "darkblue",
                              "none" = "white")),
      annotation_name_side = "left",
      show_legend = TRUE,
      annotation_legend_param = list(seg_type = list(title = NULL))
    )
  }
  if(all(is.character(order)) & length(order)==length(colnames(m))){
    clone_order <- order
  } else if(order==TRUE){
    clone_order <- sort_from_diploid_root(m, na.rm = T, ...) 
  } else {
    clone_order <- colnames(m)
  }
  
  label_colors <- ifelse(is_dup, "red", "black")
  if(grepl("cn", cn_slot)){
    col_fun <- make_cn_colorscale(max(m,na.rm=T))  
  } else if(grepl("baf", cn_slot)) {
    col_fun <- circlize::colorRamp2(
      c(0, 0.2, 0.33, 0.5),
      c("#2c7fb8", "#7fcdbb", "#edf8b1", "white")
    )
  } else if (grepl("allelespecific", cn_slot)){
    col_fun <- c("1"="white", "2"="blue", "3"="red", "4"="darkgreen", "5"="#CCCCCC", "6"="#666666", "7"="black", "8"="#E0E0E0", "9"="grey", na.value="black")  # The colors for each value
    
  }
  
  # scale_color_manual(values = c(
  #   "level1" = "#000000",  # Black
  #   "level2" = "#333333",
  #   "level3" = "#666666",
  #   "level4" = "#999999",
  #   "level5" = "#CCCCCC",
  #   "level6" = "#E0E0E0"   # Light Gray
  # ))
  
  
  clone_slot <- names(data[["cells"]])[max(grep("clone_", names(data[["cells"]])))]
  title <- paste(paste(names(data$parameters),data$parameters,sep=": "),collapse=", ")
  n_segs <- ifelse(only_divergent, paste(nrow(m), "divergent"), nrow(data[["segments"]][[segs_slot]]))
  subtitle <- paste0("segs_slot: ",segs_slot,
                     ", clones: ",clone_slot,
                     ", ",nrow(data$clones)," clones, ",
                     n_segs, " segments")
  if(sum(is_dup)>0) subtitle <- paste0(subtitle,", ",sum(is_dup)," duplicates")
  ht <- Heatmap(t(m),cluster_columns = F, col=col_fun, name = cn_slot,
                column_split = chr_labels,
                column_gap = unit(0, "mm"),
                column_title_gp=gpar(fontsize=8),
                column_title_side = "bottom",
                row_gap = unit(0, "mm"),
                border_gp = gpar(col = "black", lwd = 0.5),
                cluster_rows=F,
                border=T,
                row_names_gp = gpar(col = label_colors, fontsize=8),
                row_names_side = "left",
                bottom_annotation = column_anno,
                row_title = "clones",
                heatmap_legend_param = list(
                  title = "Copy Number",
                  at=c(1,2,3,4,5,6,7,8,9),
                  labels = c("cn_a:1, cn_b:1", "cn_a:1, cn_b:0", "cn_a:2, cn_b:1", "cn_a:2, cn_b:0", "copy 1", "copy 0", "NA", "copy 2"),  # Custom labels for color scale
                  title_gp = gpar(fontsize = 10),  # Legend title font size
                  labels_gp = gpar(fontsize = 8),  # Legend label font size
                  grid_height = unit(5, "mm")  # Height of the legend
                ), na_col="black"
  )
  draw(ht, 
       column_title=paste(title,subtitle,sep="\n"), 
       column_title_side = "top", 
       column_title_gp=gpar(fontsize=10))
}


