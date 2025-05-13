### Naming functions
dna_to_cellid <- function(x){
    sub("([A-Z]{3}[0-9]{5})(D?)_([A-Z][0-9]{2}).*","\\1.\\3", x)
}
rna_to_cellid <- function(x){
    sub("([A-Z]{3}[0-9]{5})(R?).([A-Z][0-9]{2}).*","\\1.\\3", x)
}
cellid_to_plateid <- function(x){
    sub("([A-Z]{3}[0-9]{5}).*","\\1", x)
}
cellid_to_well <- function(x){
    sub("([A-Z]{3}[0-9]{5}).([A-Z][0-9]{2})","\\2", x)
}
cellid_to_bamfile <- function(x, path="/proj/sens2018557/dntr/aligned"){
    r <- sub("([A-Z]{3}[0-9]{5})\\.([A-Z][0-9]{2})","\\1D_\\2.dedup.bam", x)
    if(!is.null(path)) r <- paste(path, r, sep="/")
    return(r)
}
cellid_to_bedfile <- function(x, path="/proj/sens2018557/dntr/aligned"){
    r <- sub("([A-Z]{3}[0-9]{5})\\.([A-Z][0-9]{2})","\\1D_\\2.bed.gz", x)
    if(!is.null(path)) r <- paste(path, r, sep="/")
    return(r)
}
cellid_to_row <- function(x){
    sub("([A-Z]{3}[0-9]{5}).([A-Z])[0-9]{2}","\\2", x)
}
cellid_to_col <- function(x){
    sub("([A-Z]{3}[0-9]{5}).[A-Z]([0-9]{2})","\\2", x)
}

### RNA processing functions
norm.log.counts <- function(counts) {
    norm.fact <- colSums(counts)
    counts.norm <- t( apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
    counts.log <- log2(counts.norm)
    rownames(counts.log) <- rownames(counts)
    counts.log
}

get.cutoff.lognorm <- function(my.counts.log, quantile.cut=0.001, gene.name='ACTB') {
    cl.act <- my.counts.log[gene.name,]
    cl.act.m <- median(cl.act)
    cl.act.sd <- sqrt(sum((cl.act[cl.act > cl.act.m] - cl.act.m)^2)/(sum(cl.act  > cl.act.m)-1))
    my.cut <- qnorm(p=quantile.cut, mean=cl.act.m, sd=cl.act.sd)
    my.cut
}

### Colors and plot themes
theme_dntr <- function (font_size = 14, font_family = "Helvetica", line_size = 0.5, axis_labels=T, axis_ticks=F, legend="bottom") {
    half_line <- font_size/2
    small_rel <- 0.857
    small_size <- small_rel * font_size
    th <- theme_grey(base_size = font_size, base_family = font_family) %+replace% 
        theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black", size=1), 
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

patient_cols <- setNames(c("#e53935","#795548","#757de8","#2196f3","#3AAA80",
                  "#62efff","#1de9b6","#ffc107",
                  "#607d8b","#1a9850","#bbbbbb"),
                  c("ALL1", "ALL2", "ALL3", "ALL4","ALL5",
                    "MRD_ALL71", "MRD_ALL64", "MRD_ALL67", 
                    "MRD_ALL68", "MRD_ALL47", "MRD_ALL35"))

source_cols <- setNames(c("#e53935","#ab000d", # Red
                 "#795548", # Brown
                 "#757de8","#3f51b5", "#002984", # Indigo
                 "#2196f3", "#0069c0", # Blue
                 "#62efff","#00bcd4","#008ba3", # Cyan
                 "#1de9b6","#52c7b8","#009688","#00675b", #Teal
                 "#ffc107","#c79100", #Amber
                 "#607d8b","#34515e", #Blue gray
                 "#cfcfcf"), # Grey
                 nm=c("ALL1_dx", "ALL1_rel", 
                      "ALL2_dx", 
                      "ALL3_dx", "ALL3_rel1", "ALL3_rel2", 
                      "ALL4_dx", "ALL4_rel", 
                      "ALL35_dx", "ALL35_d15", "ALL35_rel", 
                      "ALL47_dx", "ALL47_d2", "ALL47_d3", "ALL47_d15", 
                      "ALL64_dx", "ALL64_d15",
                      "ALL71_dx", "ALL71_d2d5d15", 
                      "other"))

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
                    "relapse 2d10"="#840077",
                    "target"="#136b39",
                    "non_target"="#840077", 
                    "4w_target"="#002984",
                    "0_weeks"="#fee08b", 
                    "4_weeks"="#d77c27", 
                    "prog"="#d73027")

phase_cols <- c("G2M"="purple","S"="darkorange","G1"="#cfcfcf")

# Create copy number colors
cn.colors <- c("darkblue","blue","white","red","darkred",
               grDevices::colorRampPalette(c("darkred","brown","green","yellow"))(46))
names(cn.colors) <- 0:50


### DNA processing functions

# GC normalization by lowess
lowess.gc = function(x, y) {
    low = lowess(x, log(y), f=0.05, iter=10);
    z = approx(low$x, low$y, x)
    return(exp(log(y) - z$y))
}
# L1/Manhattan distance for scale factors
manhattan.dist <- function(factor, segs) sum(abs((segs*factor)-round(segs*factor)))
# euclidean.dist <- function(a, b) sqrt(sum((a - b)^2))

# Prune copy number profiles: applied to single metacell segments after initial segmentation
prune.neighbors <- function(my.segs, norm.bincounts, scale.factor, thresh1=1, thresh2=1, sd.thresh=1, exclude_high_cn=NULL) {
    # Start: find segments, make data frame of them.
    segpos <- which(my.segs[-length(my.segs)] != my.segs[-1])
    int.segs <- data.frame(start=c(1,segpos+1), end=c(segpos, length(my.segs)), cpy=my.segs[c(1,segpos+1)])
    int.segs$size <- int.segs$end-int.segs$start +1
    
    # # Exclude highly amplified regions. Add them back at the end
    # if(!is.null(exclude_high_cn)){
    #   int.segs.high <- int.segs[int.segs$cpy > exclude_high_cn, ]
    #   int.segs <- int.segs[int.segs$cpy < exclude_high_cn, ]
    #   cat("Excluded high CN regions (above cpy",exclude_high_cn,"; n=",nrow(int.segs.high),")\n")
    # }
    
    # Make distribution of adjacent segment scores to determine cutoff (used in round 2)
    winsorize <- function(x, quantile=0.01) {
        minMax <- quantile(x, c(quantile, 1-quantile))
        xret <- x
        xret[xret < minMax[1]] <- minMax[1]
        xret[xret > minMax[2]] <- minMax[2]
        xret
    }
    
    mini <- which.min(int.segs$size)
    minl <- int.segs$size[mini]
    # First round, cut everything that is too small, merge with best neighbor
    while(minl < thresh1 && (is.null(exclude_high_cn) || int.segs$cpy[mini] < exclude_high_cn)) {
        
        cat(".")
        
        m.b <- (int.segs$start[mini]):(int.segs$end[mini])
        m1 <- m2 <- Inf
        if(mini > 1) m1 <- sum((norm.bincounts[m.b]*scale.factor-rep(int.segs$cpy[mini-1]))^2)/length(m.b) #Euclidean, maybe better
        if(mini < nrow(int.segs)) m2 <- sum((norm.bincounts[m.b]*scale.factor-rep(int.segs$cpy[mini+1]))^2)/length(m.b)
        m <- which.min(c(m1,m2))
        if(m1 == m2) {
            # Use the biggest sized fragment if the ploidy is the same on both sides.
            m <- which.max(int.segs$size[c(mini-1,mini+1)]) 
        }
        if(m == 2) {
            mp <- mini+1
            int.segs[mp,'start'] <- int.segs[mini,'start']
            int.segs[mp, 'size'] <- int.segs[mini,'size']+int.segs[mp, 'size']
        } else {
            mp <- mini-1
            int.segs[mp,'end'] <- int.segs[mini,'end']
            int.segs[mp, 'size'] <- int.segs[mini,'size']+int.segs[mp, 'size']
        }
        
        # Prune
        int.segs <- int.segs[-mini,]
        mini <- which.min(int.segs$size)
        minl <- int.segs$size[mini]
    }
    
    msp <- rep(int.segs$cpy, times=int.segs$size)
    segpos2 <- which(msp[-length(msp)] != msp[-1])
    
    # Second round, cut small-ish segments with small diff to neighbor.
    int.segs2 <- data.frame(start=c(1,segpos2+1), end=c(segpos2, length(msp)), cpy=msp[c(1,segpos2+1)])
    int.segs2$size <- int.segs2$end-int.segs2$start+1
    
    mindiffs <- t(as.matrix(as.data.frame(lapply(1:nrow(int.segs2), function(i) {
        l <- int.segs2$size[i]
        m.b <- (int.segs2$start[i]):(int.segs2$end[i])
        m1 <- m2 <- Inf
        if(i > 1) m1 <- sum((norm.bincounts[m.b]*scale.factor-int.segs2$cpy[i-1])^2)/length(m.b) #Euclidean, maybe better
        if(i < nrow(int.segs2)) m2 <- sum((norm.bincounts[m.b]*scale.factor-int.segs2$cpy[i+1])^2)/length(m.b)
        morig <- sum((norm.bincounts[m.b]*scale.factor-rep(int.segs2$cpy[i]))^2)/length(m.b)
        c(morig, min(c(m1,m2)))
    }))))
    
    a <- cbind(int.segs2, mindiffs)
    rownames(a) <- NULL
    a$div <- a[['1']]/a[['2']]
    #    hist(log(a$div), breaks=100)
    # min/max for good cuts. Use median and Winsorize before calculating sd to avoid influence by outliers.
    #    mincut <- median(log(a$div))-sd(log(winsorize(a$div, 0.05))*sd.thresh)
    mincut <- -2
    #    abline(v=mincut)                      
    
    mini <- which.min(int.segs2$size)
    minl <- int.segs2$size[mini]
    
    excluded.segs.ind <- 1
    
    while(minl < thresh2) {
        m.b <- (int.segs2$start[mini]):(int.segs2$end[mini])
        m1 <- m2 <- Inf
        if(mini > 1) m1 <- sum((norm.bincounts[m.b]*scale.factor-int.segs2$cpy[mini-1])^2)/length(m.b) #Euclidean, maybe better
        if(mini < nrow(int.segs2)) m2 <- sum((norm.bincounts[m.b]*scale.factor-int.segs2$cpy[mini+1])^2)/length(m.b)
        morig <- sum((norm.bincounts[m.b]*scale.factor-rep(int.segs2$cpy[mini]))^2)/length(m.b)
        c(morig, min(c(m1,m2)))
        cat(morig, ", ", m1, ", ", m2, "\n")
        div <- log(morig/min(c(m1,m2)))
        if(!is.null(exclude_high_cn) && any(int.segs2$cpy[mini] > exclude_high_cn)){
            cat("Skip high CN segment: ", paste(int.segs2[mini,], collapse=" "), "\n")
            excluded.segs.ind <- excluded.segs.ind+1
        } else if(div < mincut) {
            cat("Do not remove seg with div: ", div, ", and statline: ", paste(int.segs2[mini,], collapse=" "), "\n")
            excluded.segs.ind <- excluded.segs.ind+1
        } else {
            cat("Remove seg with div: ", div, ", and statline: ", paste(int.segs2[mini,], collapse=" "), "\n")
            m <- which.min(c(m1,m2))
            if(m1 == m2) {
                # Use the biggest sized fragment if the ploidy is the same on both sides.
                m <- which.max(int.segs2$size[c(mini-1,mini+1)]) 
            }
            if(m == 2) {
                mp <- mini+1
                int.segs2[mp,'start'] <- int.segs2[mini,'start']
                int.segs2[mp, 'size'] <- int.segs2[mini,'size']+int.segs2[mp, 'size']
            } else {
                mp <- mini-1
                int.segs2[mp,'end'] <- int.segs2[mini,'end']
                int.segs2[mp, 'size'] <- int.segs2[mini,'size']+int.segs2[mp, 'size']
            }
            # Prune
            int.segs2 <- int.segs2[-mini,]
        }
        mini <- order(int.segs2$size)[excluded.segs.ind]
        minl <- int.segs2$size[mini]
    }
    
    msp2 <- rep(int.segs2$cpy, times=int.segs2$size)
    segpos3 <- which(msp2[-length(msp2)] != msp2[-1])
    #    plot(msp2)
    #    abline(v=segpos3)
    
    int.segs3 <- data.frame(start=c(1,segpos3+1), end=c(segpos3, length(msp2)), cpy=msp2[c(1,segpos3+1)])
    int.segs3$size <- int.segs3$end-int.segs3$start+1
    
    int.segs3
}

# Find discordant/disjoint segments across all clones/metacells
# cn (mtx): bins x clones copy number integers
find_disjoint_segs <- function(cn, keep.all=F, min.bin.dist=1){
    # Find all unique copy states across clones 
    cn.uniq <- unique(cn)
    # cn.uniq <- cn.uniq[!rowSums(is.na(cn.uniq))==ncol(cn.uniq),] # Drop NA rows
    # Loop over unique states and find matches in whole cn matrix
    cn.states.all <- apply(cn.uniq, 1, function(v){
        # Note: need to transpose matrix to match vector
        which(colSums(v == t(cn)) == ncol(cn))
    })
    try(if(sum(sapply(cn.states.all, length)) != nrow(cn)) stop("Length of all cn states doesnt match nr of bins"))
    
    # Keep only cn states that are not identical across all clones
    if(keep.all) keep.idx=T else keep.idx <- which(rowSums(cn.uniq[,1]==cn.uniq) != ncol(cn.uniq))
    cn.states <- cn.states.all[keep.idx]
    
    # Make into segments. TODO: Could add a "discordant index" to keep track 
    disc.seg.list <- lapply(cn.states, function(x){
        x. <- c(0,cumsum( diff(x)-1 ))
        l <- lapply( split(x, x.), range )
        data.frame( matrix(unlist(l), nrow=length(l), byrow=T) )
    })
    
    # Join all segments (note: trimming and chr boundary checks in separate fx)
    disc.seg <- bind_rows(disc.seg.list)
    colnames(disc.seg) <- c("start","end")
    # DEBUG: start/end turns non-numeric in some cases?
    if(!all(sapply(disc.seg, is.numeric))) stop("Non-numeric result in disjoint segments")
    disc.seg.filt <- disc.seg %>% filter(end-start >= min.bin.dist)
    o <- order(disc.seg.filt[[1]], disc.seg.filt[[2]])
    return(disc.seg.filt[o,])
}

# If segments cross boundaries (chromosome boundaries): split segments
# segs (df): start, end (in bins)
# bins (df): bin, chr
# bnd (df): chr, start, end (in bins)
split_segs_at_boundaries <- function(segs, bins, bnd, min.bin.dist=1){
    res <- data.frame()
    for(i in 1:nrow(segs)){
        s <- segs[[1]][i] # start
        e <- segs[[2]][i] # end
        chr <- bins$chr[s] # start chr
        if(e > bnd$end[bnd$chr==chr]){ # Does end position cross chr boundary? If so, split into two. 
            # NOTE: Does not cover situation when one continuous segment covers at least 3 chromosomes
            first <- data.frame(chr=chr,
                                start=s,
                                end=bnd$end[bnd$chr==chr])
            second <- data.frame(chr=bins$chr[e],
                                 start=first$end+1,
                                 end=e)
            res <- rbind(res, first, second)
        } else {
            res <- rbind(res, data.frame(chr=chr,start=s,end=e))
        }
    }
    size_select <- (res$end+1)-res$start > min.bin.dist
    if(sum(!size_select)>0) warning(paste0("Filtered out ",sum(!size_select)," segment(s) after chr split:",
                                             paste(res[!size_select,],collapse=","),sep="\n"))
    return(res[size_select,])
}

# Convert segments to clone x bin matrix with CN integers
# segs (df): start, end of segments (in bins)
# cn (mtx): bins x clone matrix with cn integers
segs_to_mtx <- function(segs, cn){
  if(!is.matrix(cn)){
    cn.vec <- sapply(1:nrow(segs), function(i){
      # Single clone, subset vector
      median(cn[segs$start[[i]]:segs$end[[i]]], na.rm=T)
    })
    return(cn.vec)
  } else {
    cn.list <- lapply(1:nrow(segs), function(i){
      idx <- segs$start[[i]]:segs$end[[i]]
      if(length(idx)>1){
        apply(cn[idx,], 2, function(x) median(x, na.rm=T)) 
      } else {
        # If segment = 1 bin
        cn[idx,]
      }
    })
    cn.mat <- as.matrix(bind_rows(cn.list))
    return(cn.mat)
  }
}

# Get ape tree object from CN matrix: clones x bins with CN integers
# tree_init: ape::nj or phangorn::upgma
# excl_clone: exclude clones/metacell, by default clone "0" which is outliers in hdbscan
tree_from_mtx <- function(mtx, root_cn=2, dist_method="manhattan", tree_init="nj", excl_clone="0"){
    require(ape)
    dist <- dist(rbind(t(mtx), "diploid"=root_cn), method=dist_method)
    if(!is.null(excl_clone)) excl <- c("diploid",excl_clone) else excl <- "diploid"
    if(tree_init=="nj"){
        tr <- nj(dist) %>% root(outgroup = "diploid") %>% drop.tip(excl)
    } else if(tree_init=="upgma"){
        require(phangorn)
        tr <- upgma(dist, method = "complete") %>% root(outgroup = "diploid") %>% drop.tip(excl)
    } else { stop(paste("Unrecognized tree initialization function:",tree_init)) }
    return(tr)
}

plot_ggtree <- function(tr, show_labels=F){
    require(ggtree)
    # ggplot tree vertically
    p <- ggtree(tr, ladderize=F, size=1) + 
        geom_tiplab(size=0,align=TRUE) +
        hexpand(.01) +
        theme(
            panel.spacing.x = unit(0, 'mm'),
            panel.spacing.y = unit(1, 'mm'),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0, 'cm')
        )
    if(show_labels) p <- p + geom_tiplab(size=3,align=TRUE)
    return(p)
}

plot_cn_grid <- function(mtx, excl_clone="0", order=NULL, chr_pos=NULL, labels=T){
  brks <- c(-Inf,0,1,2,3,4,20,50,Inf)
  cols <- c("darkblue","blue","white","red","darkred","brown","green","yellow")
  cols_vec2 <- setNames(as.character(cut(0:200, breaks = brks, labels = cols)),nm=0:200)
  
  clone.idx=T
  if(!is.null(excl_clone)) clone.idx <- which(colnames(mtx)!=excl_clone) # Drop selected clones
  df <- as.data.frame(mtx[,clone.idx]) %>% 
    # mutate(idx=paste0(disc.seg$chr,":",disc.seg$start_bp,"-",disc.seg$end_bp)) %>% 
    mutate(idx=seq_along(.[[1]])) %>% 
    gather(clone, cn, -idx)
  if(!is.null(order)) df$clone <- factor(df$clone, levels=rev(order))
  p <- ggplot(df, aes(idx, clone)) + 
    geom_tile(aes(fill=as.factor(cn)), show.legend=F) +
    # scale_fill_manual(values=str_replace_all(cols_vec2, "grey", "white")) +
    scale_fill_manual(values=cols_vec2) +
    labs(x=NULL,y=NULL) +
    theme_classic() +
    theme(legend.position = "none",
          panel.border = element_rect(size = 1, color = 'black', fill = NA),
          panel.spacing.x = unit(0, 'mm'),
          panel.spacing.y = unit(1, 'mm'),
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(t = 0, r = 1, b = 0, l = 0, 'pt')
    ) +
    scale_x_continuous(expand=c(0,0))
  if(!is.null(chr_pos)) p <- p + geom_vline(xintercept=chr_pos$line_bin, color="grey")
  if(labels) p <- p + theme(axis.text.x=element_text(angle=90,size=8,vjust=0.5,hjust=1))
  return(p)
}

# Segs: All "common breakpoint segments" to use for scoring: start, end required (in bin positions)
# Clones (df): clone, dna_library_id
# Counts (matrix): bin filtered raw counts x bins
# excl_clone: clone/clones to exclude from QC
calc_clone_cor <- function(segs, counts, clones, excl_clone=c("0","S","G2M")){
  require(Rfast)
  if(!any(class(clones) %in% "data.frame") & !all(c("dna_library_id","clone") %in% colnames(clones))) stop("clones needs to be a data.frame with columns dna_library_id and clone")
  clones.list <- split(clones$dna_library_id, clones$clone)
  if(!is.null(excl_clone)) {
    clones.excl <- grep("_G2M$|_S$|_0$",names(clones.list),value = T)
    clones.list <- clones.list[!grepl("_G2M$|_S$|_0$",names(clones.list))]
  }
  # Make sure counts is in same order as clones
  counts <- counts[,clones$dna_library_id]
  # Aggregate each clone per segment
  clone.counts <- lapply(clones.list, function(x){
    lapply(1:nrow(segs), function(i){
      sum(counts[segs$start[[i]]:segs$end[[i]],x]) # sum of counts for that segment and those cells/clone  
    }) %>% unlist
  }) %>% bind_cols
  # Aggregate each cell per segment
  seg.counts <- apply(counts, 2, function(x){
    sapply(1:nrow(segs), function(i){
      sum(x[segs$start[[i]]:segs$end[[i]]])  
    })
  })
  
  # Normalize clones and cells
  clone.counts.norm <- t(t(clone.counts) / colSums(clone.counts)) # OBS! Have to transpose since R does -first division
  # Remove unwanted clones
  clone.counts.filt <- clone.counts.norm[,!colnames(clone.counts.norm) %in% excl_clone]
  
  seg.counts.norm <- t(t(seg.counts) / colSums(seg.counts))
  segs.len <- segs$end - segs$start
  
  # THIS is the best one, I think. counts*log2(length) / length. So it's average read density * log(length). Makes intuitive sense that long fragments are more important but in a log function rather than linear.
  a <- dista(xnew=t(seg.counts.norm*log2(segs.len)/segs.len), 
             x=t(clone.counts.norm*log2(segs.len)/segs.len), k=1, index=TRUE)
  a.filt <- dista(xnew=t(seg.counts.norm*log2(segs.len)/segs.len), 
                  x=t(clone.counts.filt*log2(segs.len)/segs.len), k=1, index=TRUE)
  
  # Store results
  clones$new_clone  <- colnames(clone.counts.norm)[a]
  clones$new_clone_filt <- colnames(clone.counts.filt)[a.filt]
  
  # Cross tables
  clones.tab <- table("original"=clones$clone,
                      "new"=clones$new_clone)
  clones.tab.filt <- table("original"=clones$clone,
                           "new"=clones$new_clone_filt)
  # o.hc <- hclust(dist(clones.tab))
  # o.ord <- o.hc$labels[o.hc$order]
  # n.hc <- hclust(dist(t(clones.tab)))
  # n.ord <- n.hc$labels[n.hc$order]
  clones.tab.df <- as.data.frame(clones.tab) %>% group_by(original) %>% 
    mutate(n=sum(Freq), original=as.character(original), new=as.character(new))
           # original=factor(as.character(original),levels=c(setdiff(n.ord,clones.excl),clones.excl)),
           # new=factor(as.character(new),levels=c(setdiff(o.ord,clones.excl),clones.excl)))
  clones.tab.df.filt <- as.data.frame(clones.tab.filt) %>% group_by(original) %>% mutate(n=sum(Freq), original=as.character(original), new=as.character(new))
  
  p.tab <- ggplot(clones.tab.df, aes(original, new, fill=Freq/n)) + geom_tile() + scale_fill_viridis_c() +
    geom_text(aes(label = ifelse(original==new, round(Freq/n, 2), NA))) + labs(subtitle="All clones") +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  p.tab.filt <- ggplot(clones.tab.df.filt, aes(original, new, fill=Freq/n)) + geom_tile() + scale_fill_viridis_c() +
    geom_text(aes(label = ifelse(original==new, round(Freq/n, 2), NA))) + labs(subtitle="Filtered clones") +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  
  # Return as list
  return(list("clones"=clones, "plot_all"=p.tab, "plot_filtered"=p.tab.filt))
}

smooth_bins_mean <- function(mat, w, bins=NULL){
  # Faster version for mean across bin windows 
  newbins <- c(sapply(1:ceiling(nrow(mat)/w),function(x){rep(x, w)}))[1:nrow(mat)]
  mat2 <- rowsum(mat, newbins) / tabulate(newbins) # Sum / nr of bins = mean
  if(!is.null(bins)){
    aggr.bin.start <- bins %>% dplyr::slice(which(row_number() %% w == 1)) %>% pull(bin)
    rownames(mat2) <- aggr.bin.start
    bins.m <- aggregate(as.integer(bins$chr_int), by=list(newbins), FUN=min)
    bins.m <- bins.m[order(as.integer(bins.m$Group.1)),] %>%
      mutate(chr_short=case_when(x == 23 ~ "X", x == 24 ~ "Y", TRUE ~ as.character(x)))
    return(list(mat=mat2, bins=bins.m))
  } else {
    rownames(mat2) <- 1:nrow(mat2)
    return(list(mat=mat2, bins=NULL))
  }
}

smooth_bins_subset <- function(m, w, bins=NULL){
  sub <- m[seq(1, nrow(m), w), ]
  # Debug: Expand back out: rep(sub,each=10)
  # Faster version for mean across bin windows 
  if(!is.null(bins)){
    bins.sub <- bins[seq(1, nrow(bins), w),]
    rownames(sub) <- bins.sub$bin
    return(list(mat=sub, bins=bins.sub))
  } else {
    rownames(sub) <- 1:nrow(sub)
    return(list(mat=sub, bins=NULL))
  }
}

# Tree/dist fx
make_tree <- function(cn, dist_method="manhattan", tree_method="nj", rootid="root", ladderize=T, return_dist=F, cores=1, is_dist=F){
  require(ape)
  if(is_dist){
    # Skip dist calculation
    dist <- cn
  } else {
    cat("Calc dist matrix...\n")
    cn2 <- cbind(cn, 2) # Add diploid root
    colnames(cn2)[ncol(cn2)] <- rootid
    dist <- amap::Dist(t(cn2), method=dist_method, nbproc = cores)
    if(return_dist){
      dist.mtx <- as.matrix(dist)
      diag(dist.mtx) <- NA
      return(dist.mtx)
    }
  }
  cat("Calc tree...\n")
  if(tree_method=="nj"){
    cat("Tree with neighbour joining...\n")
    tree <- nj(dist)
  } else if(tree_method=="fast"){
    cat("Tree with minimum evolution algo...\n")
    tree <- fastme.bal(dist)
  } else {
    stop("tree method not recognized")
  }
  tree <- root.phylo(tree, outgroup = which(tree$tip.label == rootid), resolve.root = TRUE)
  tree <- drop.tip(tree, tip = rootid)
  if(ladderize) tree <- ladderize(tree)
  return(tree)
}

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

# TODO: Common DNA functions
# - plot legend for CN colors
# - annotate segments with any genes (bed format)

plot_clone <- function(df.list=NULL, show_detail=F, cn.trunc=4, labels=NULL, df.chr=NULL, region=NULL){
  brks <- c(-Inf,0,1,2,3,4,20,50,Inf)
  cols <- c("darkblue","blue","grey","red","darkred","brown","green","yellow")
  cols_vec2 <- setNames(as.character(cut(0:200, breaks = brks, labels = cols)),nm=0:200)
  if(is.null(names(df.list))) warning("List without names, naming clones sequentially")
  df <- bind_rows(df.list, .id="clone")
  df$bin <- as.integer(df$bin)
  if(!is.null(labels)) labels <- labels[unique(df$clone)]
  chr_hjust=0.5
  if(!is.null(region)){
    if(all(grepl("chr",region))){
      df <- df %>% filter(chr %in% region)
      df.chr <- df.chr %>% filter(chr %in% unique(df$chr))
    } else if(all(is.numeric(region))){
      df <- df %>% filter(bin >= min(region) & bin <= max(region))
      bp <- paste0(unique(df$chr),":",round(df$bp_chr_start[1]/1e6,1),"-",round(df$bp_chr_end[nrow(df)]/1e6,1),"Mb")
      region.anno.df <- tibble(clone=names(labels),chr=df$chr[nrow(df)],
                               label=paste0("Bins: ",min(df$bin),"-",max(df$bin), 
                                            "(",bp,")"))
      # new_label_pos <- df$bp[which(df$bin==round(min(region)+((max(region)-min(region))/2)))][1]
      df.chr <- df.chr %>% filter(chr %in% unique(df$chr)) %>% 
        mutate(line=NULL, label_pos=df$bp[which(df$bin==min(region))][1])
      chr_hjust=-0.5
    } else {
      stop("Unclear region: expecting either chromosome (eg chr1) or numeric vector of two c(100,120)")
    }
  }
  if(!is.null(cn.trunc)){
    ymax <- max(min(cn.trunc+0.5, max(df$cn.pruned)+1), 4) # 1) If cn.trunc > max in set, stop at max 2) Never go lower than 4
  } else {
    ymax <- max(df$cn.pruned)+1
  }
  
  p <- ggplot(df,aes(x=bp)) + 
    geom_hline(yintercept=2,color="grey",linetype="solid",size=0.1)
  if(show_detail) p <- p + geom_point(aes(y=raw),pch=19,size=0.1,alpha=0.1,color="grey")
  p <- p +  geom_point(aes(y=gc,color=as.factor(cn.pruned)),pch=19,size=0.1,alpha=0.5,show.legend=F)
  if(show_detail){
    p <- p + geom_point(aes(y=seg),pch=19,size=0.15,col="orange",alpha=0.5) +
      geom_point(aes(y=cn),color="black",pch=19,size=0.3,alpha=0.5,show.legend=F)
  }
  p <- p + geom_point(aes(y=cn.pruned, color=as.factor(cn.pruned)),pch=19,size=0.5,show.legend=F) +
    scale_color_manual(values=cols_vec2) +
    # scale_x_continuous(expand=c(0.01,0.01)) +
    scale_x_continuous(expand=c(0,0)) +
    lims(y=c(0,ymax)) + # Max ymax, min 4, otherwise dependent on actual copy number
    theme_dntr(axis_ticks = T) +
    theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), plot.margin=margin(5,5,0,5)) +
    labs(x=NULL, y="copy number") +
    geom_text(data=df.chr, aes(x=label_pos,y=0,label=chr_short), check_overlap = T, hjust=chr_hjust, vjust=-0.2)
  if(exists("region.anno.df")) p <- p + geom_text(data=region.anno.df, aes(x=Inf,y=Inf,label=label), inherit.aes = F, hjust=1, vjust=1.5)
  
  if(any(df$cn.pruned>cn.trunc)) p <- p + geom_point(data=filter(df,cn.pruned>cn.trunc), aes(y=cn.trunc+0.5, fill=as.factor(cn.pruned)), color="black", pch=24,size=2.5,show.legend=F)
  if(length(unique(df$chr))>1 && length(unique(df$chr))<24){
    clone.anno.df <- tibble(clone=names(labels), label=labels, chr=unique(df$chr)[1])
    p <- p + 
      geom_text(data=clone.anno.df, aes(x=-Inf,y=Inf,label=label), inherit.aes = F, hjust=-0.01, vjust=1.5) +
      facet_grid(factor(clone,levels=str_sort(unique(clone),numeric=T)) ~ chr, space="free", scales="free") + 
      theme(strip.text.x=element_blank(), panel.spacing.x=unit(0, "lines"))
  } else {
    p <- p + geom_vline(xintercept=df.chr$line, color="grey") + 
      facet_wrap(vars(clone), ncol=1, labeller=labeller(clone=labels))
  }
  p
}
