suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
library(caTools)
library(Rcpp)
library(ggplot2)
library(patchwork)
source(snakemake@params[["src_general"]])
sourceCpp(snakemake@params[["src_cpp"]])

options(scipen=999)

counts.file <- snakemake@input[["cell"]]
map <- as.numeric(readLines(snakemake@input[["map"]]))
gc <- as.numeric(readLines(snakemake@input[["gc"]]))
bins <- read.table(snakemake@input[["bins"]],col.names=c("chr","start_coord","end_coord"),stringsAsFactors=F) # Use this for splitting
if(!is.null(snakemake@input[["badbins"]])) badbins <- read.table(snakemake@input[["badbins"]],header=T,stringsAsFactors=F) else badbins <- NULL
gc_min <- snakemake@params[["gc_min"]]
map_min <- snakemake@params[["map_min"]]
max.w <- snakemake@params[["max_w"]]
min.w <- snakemake@params[["min_w"]]
min.scale <- snakemake@params[["min_scale"]]
max.scale <- snakemake@params[["max_scale"]]
reads.per.window <- snakemake@params[["reads_per_window"]]
outfile.qc <- snakemake@output[["pdf_qc"]]
outfile.cnv <- snakemake@output[["cnv_png"]]
outfile.cnv2 <- snakemake@output[["cnv2_png"]]
outfile.seg <- snakemake@output[["seg"]]
outfile.seg_cn <- snakemake@output[["seg_cn"]]
outfile.cn <- snakemake@output[["cn"]]
outfile.scalefactors <- snakemake@output[["scalefactors"]]
outfile.normcounts <- snakemake@output[["normcounts"]]
cn.trunc <- snakemake@params[["cn_trunc"]]

bins$chr <- factor(bins$chr, levels=unique(bins$chr)) # Maintain natural sort order for chromosomes
bins$bin_unfilt <- seq_along(bins$chr)

### Start
cell_id <- read_lines(counts.file, n_max=1)
counts <- as.numeric(read_lines(counts.file, skip=1)) # Skip header line
good.bins <- gc > gc_min & map > map_min & ! bins$bin_unfilt %in% badbins$bin_unfilt

# Filtered
bins.f <- bins[good.bins,]
bins.f$bin <- seq_along(bins.f$chr)
#counts.f <- counts[good.bins]
counts.f <- counts
gc.f <- gc[good.bins]
map.f <- map[good.bins]

# Split by chromosome
counts.fs <- split(counts.f, bins.f[[1]])
gc.fs <- split(gc.f, bins.f[[1]])
map.fs <- split(map.f, bins.f[[1]])
bins.fs <- split(bins.f, bins.f[[1]])
bins.f$bin_local <- unlist(lapply(bins.fs, function(x){1:nrow(x)})) # Bin IDs for chromosome-level

# Calculate quadratic fit to gc. Do once per cell!
gc.coefs <- summary(lm(counts.f ~ gc.f + I(gc.f^2)))$coefficients[,1]
# Calculate window size. Do once per cell!
w <- as.numeric(max(min.w, min(round(reads.per.window/mean((counts.f))), max.w)))

# Functions
gc.map.correct <- function(bin.counts, gc, map, gc.coefs) {
    predq <- function(qest, gc, scale=1) {
        (qest[1]+gc*qest[2]+gc^2*qest[3])/scale
    }
    gc.scale <- predq(gc.coefs, 0.45, 1)
    weight <- map*predq(gc.coefs, gc, gc.scale)
    bin.counts/weight
}
call.breakpointsCPP <- function(bin.counts, gc, map, gc.coefs, w, peak.size.cutoff=10, plotQC=F, chr=NULL) {
    empty_df <- data.frame(start=numeric(0),end=numeric(0),peak.size=numeric(0),local.max=numeric(0),origin=character(0),llr=numeric(0))
    if(is.null(chr)) chr <- ""
    predq <- function(qest, gc, scale=1) {
        (qest[1]+gc*qest[2]+gc^2*qest[3])/scale
    }
    gc.scale <- predq(gc.coefs, 0.45, 1)
    
    if(w>length(bin.counts)/2-1) w <- round(length(bin.counts)/2-1) # Window can't be larger than half the chromosome boundaries (in practice only applies to chrY) #MARTIN: window should be smaller than half.
    map.norm <- map/mean(map)
    
    # cat("bins:",length(bin.counts),"\n")
    # cat("gc:", length(gc),"; map:",length(map),"\n")
    # cat("w:", w, "\n")
    
    llr <- sapply(w:(length(bin.counts)-w), function(i) { # TODO fix all the edge cases (w:... instead of 1:..)
        calcLLR(bin.counts, i-1, w, map.norm, gc, gcscale=gc.scale, gc.coefs) #MARTIN -1 is to translate between 1-base (R) and 0-based (C++) indices
    })
    # Call peaks by cutoff
    above.cut <- llr > 5
    #cat("Initial peaks with LLR >5:",sum(above.cut),"\n")
    if(sum(above.cut,na.rm=T)==0) { 
        #cat("^^No initial peaks detected\n")
        if(plotQC) {
            require(caTools)
            par(mfrow=c(3,1))
            plot(runmean(bin.counts, k=10), type='p', pch=19, cex=0.2, main=paste0(chr, " Raw counts; run mean, k=10"))
            plot(w:(length(llr)+w-1), llr, type='l', col='blue', lwd=1.5, main=paste0("LLR with window=", w, " pruned by size (green) and inter-peak interval (red)"))
            abline(h=5)
        }
        return(empty_df) 
    }
    
    crosses <- which(above.cut[-1] != above.cut[-length(above.cut)])+1+w # Add w to fix initial shift. I have a feeling I'll introduce a diff by 1 somewhere here. 
    if(above.cut[1]) {
        crosses <- c(1,crosses)
    }
    if(above.cut[length(above.cut)]) { # This will probably ruin sth lower down if it ever triggers
        #crosses <- c(crosses,length(bin.counts)) # MARTIN change length(above cuts)+w*2+1, etc, for length(bin.count)
        # VZ: This doesn't seem right. Will go out of index if triggered. Upper limit should be length-w
        crosses <- c(crosses,length(bin.counts)-w)
    }
    
    peaks <- as.data.frame(matrix(crosses, ncol=2, byrow=T))
    #cat("Initial peaks table:",nrow(peaks),"\n")
    
    colnames(peaks) <- c("start", "end")
    peaks$start <- peaks$start+1
    peaks$peak.size <- peaks$end-peaks$start+1
    peaks$small.peaks <- peaks$peak.size < peak.size.cutoff
    peaks$local.max <- apply(peaks, 1, function(x) {
        st <- x[1]-w  # MARTIN: remove -1-w
        en <- min(x[2]-w, length(llr)) # MARTIN: remove -1-w and add +w to llr length
        x[1]+which.max(llr[st:en])
    }) # Subtract w to fix shift
    peaks$origin <- "Initial"
    longpeaks <- peaks[!peaks$small.peaks,] # Equivalent to B0 breakpoint set (10X)
    #cat("Initial peaks, small peaks removed:",nrow(longpeaks),"\n")
    if(nrow(longpeaks)==0) { 
        #cat("^^No initial peaks left after small peak removal\n")
        if(plotQC) {
            par(mfrow=c(3,1))
            require(caTools)
            plot(runmean(bin.counts, k=10), type='p', pch=19, cex=0.2, main=paste0(chr, " Raw counts; run mean, k=10"))
            plot(w:(length(llr)+w-1), llr, type='l', col='blue', lwd=1.5, main=paste0("LLR with window=",w,"pruned by size (green) and inter-peak interval (red)"))
            abline(h=5)
        }
        return(empty_df)
    }
    
    #cat("Pruning breakpoints by LLR\n")
    maxs.l <- c(1,longpeaks$local.max-1, length(bin.counts)-1) # Maybe first peak should be at 1 instead of 0? Probably not though.
    longpeaks$llr <- sapply(2:(length(maxs.l)-1), function(i) {
        lr <- c(maxs.l[i]-maxs.l[i-1],maxs.l[i+1]-maxs.l[i])
        llr <- calcLLR(bin.counts, i=maxs.l[i], w=lr, mapnorm=map.norm, gc=gc, gcscale=gc.scale, qest=gc.coefs) # MARTIN: remove -w-i for maxs.l[i]
        #        cat("i, lr for pruning:", maxs.l[i], lr, ", llr:", llr, "\n") # MARTIN: remove -w-1 for maxs.l[i]
        llr
    })
    
    good.peaks <- longpeaks[longpeaks$llr > 5,] # Equivalent to B1 breakpoint set (10X)
    #cat("Nr of good peaks: ",nrow(good.peaks),"\n")
    final.peaks <- empty_df # MARTIN: declare final peaks before loop.
    
    if(nrow(good.peaks)>=1){
        #cat("Interval to scan:",w+1,"-",length(bin.counts)-w-1,"\n") 
        # Sliding left + right window going by bin between two breaks, but minimum w as window size
        perbase <- sapply((w+1):(length(bin.counts)-w-1), function(i) { # MARTIN: Pretty sure this shoudl be w+1/-w-1 rather than just w
            l <- max(c(w, good.peaks$local.max[good.peaks$local.max <= i])) # MARTIN: And this should be w, not w+1
            r <- min(c(good.peaks$local.max[good.peaks$local.max >= i]), length(bin.counts))
            #cat("w interval: ", i-l, r-i, "\n")
            if(i-l < w | r-i < w) {
                return(0)
            }
            return(calcLLR(bin.counts, i=i-1, w=c(i-l,r-i), mapnorm=map.norm, gc=gc, gcscale=gc.scale, qest=gc.coefs))
        })
        #perbin_df <- data.frame(idx=seq_along(perbase), bin=c((w+1):(length(perbase)+w)), llr=perbase)
        #cat("Perbase\n")
        #print(length(perbase))
        #print(range(perbase))
        # Call peaks by cutoff again
        above.cut <- perbase > 5
        #cat("putative interpeak breaks: ",sum(above.cut),"\n")
        if(sum(above.cut>0)){
            crosses <- which(above.cut[-1] != above.cut[-length(above.cut)])+w+1 # MARTIN get back to correct coords again with w+1. Should be w+2? VZ: No
            #crosses <- perbin_df[which(above.cut[-1] != above.cut[-length(above.cut)]),"bin"] # MARTIN get back to correct coords again with w+1. Should be w+2?
            if(above.cut[1]) {
                crosses <- c(0,crosses)
            }
            if(above.cut[length(above.cut)]) {
                #crosses <- c(crosses,length(bin.counts)) # MARTIN not sure about this
                crosses <- c(crosses,length(bin.counts)-w) # VZ has to be max bins - window I think. Maybe -1 too? 
            }
            peaks2 <- as.data.frame(matrix(crosses, ncol=2, byrow=T))
            colnames(peaks2) <- c("start", "end")
            peaks2$start <- peaks2$start+1 # VZ: Why? Was already added above
            peaks2$peak.size <- peaks2$end-peaks2$start+1 # VZ: Why? Was already added above
            peaks2$small.peaks <- peaks2$peak.size < peak.size.cutoff
            peaks2$local.max <- apply(peaks2, 1, function(x) {x[1]+which.max(perbase[(x[1]:x[2])-w-1])}) # MARTIN: Fix indices
            peaks2$llr <- apply(peaks2, 1, function(x) {
                max(perbase[(x[[1]]-w):(x[[2]]-w-1)])
                #max(perbase[(x[1]:x[2])-w-1])
            }) # MARTIN: -w to -w-1 #VZ: Throws error for chr4:6566-8549
            
            peaks2$origin <- "Interpeak"
            longpeaks2 <- peaks2[!peaks2$small.peaks,]
            #cat("Added interpeak breaks:",nrow(longpeaks2),"\n") #DEB
            final.peaks <- rbind(good.peaks, longpeaks2)
            final.peaks <- final.peaks[order(final.peaks$local.max),-4] # Sort, remove "small peak" entry
        } else {
            final.peaks <- good.peaks[,-4] # All putative interpeak breaks where filterd/removed. Return empty
        }
    }
    
    if(plotQC) {
        require(caTools)
        par(mfrow=c(4,1))
        plot(runmean(bin.counts, k=10), type='p', pch=19, cex=0.2, main=paste0(chr, " Raw counts; run mean, k=10. Final peaks in red"))
        abline(v=final.peaks$local.max, col="red", lwd=2, lty=1)
        plot(w:(length(llr)+w-1), llr, type='l', col='blue', lwd=1.5, main="LLR with window ~ 200 reads, pruned by size (green) and inter-peak interval (red)")
        abline(h=5)
        abline(v=longpeaks$local.max, col='green')
        abline(v=longpeaks$local.max[longpeaks$llr > 5], col='red')
        plot(longpeaks$local.max, longpeaks$llr, xlim=c(0, length(bin.counts)))
        if(exists("perbase")){
            plot(w:(length(perbase)+w-1), perbase, type='l', ylim=c(-2,40), main="LLR with variable windows, between adjacent peaks, pruned by size (red)")
            abline(h=5)
            if(exists("longpeaks2")) abline(v=longpeaks2$local.max, col="red") # Dont plot if no added peaks  
        } else {
            plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.5, paste("No added interpeak breakpoints (no peaks left after pruning)"))
        }
    }
    #cat("n peaks:",nrow(final.peaks),"\n")
    return(final.peaks)
    #return(list(gp=good.peaks, peaks=peaks2, pb=perbase, fp=final.peaks))
}

euclid.dist <- function(factor, segs) {
    sum(((segs*factor)-round(segs*factor))^2)
}
manhattan.dist <- function(factor, segs) {
    sum(abs((segs*factor)-round(segs*factor)))
}

pdf(outfile.qc)
breakpoints.byChr <- lapply(names(counts.fs), function(x){ 
    # cat(">>>", x,"\n");
    call.breakpointsCPP(counts.fs[[x]], gc.fs[[x]], map.fs[[x]],gc.coefs, w=w, peak.size.cutoff=10, plotQC=T, chr=x) 
})
dev.off()

call.segs <- function(breakpoints, counts, w){
    # 1. Prune breakpoints on inter-breakpoint distance (as nr of bins)
    # In no-break chrs, this will create a vector of inter-chr distance 
    # ie chrY has 650 bins at 20k, bp.dists=c(649)
    if(nrow(breakpoints)==0){
        #cat("No breakpoints found. Skip pruning...\n")
        seg=data.frame(start=1, end=length(counts), start.llr=Inf, end.llr=Inf, size=length(counts))
        return(seg)
    }
    bp.dists <- c(breakpoints$local.max, length(counts))-c(1, breakpoints$local.max)
    o <- which(bp.dists < w)
    # bp.dists[i] distance between breakpoints[i-1] to breakpoints[i]
    bp <- breakpoints
    #print(bp)
    # Go through breakpoints near eachother (< w x 1) and figure out which to merge (by higher LLR)
    for(i in o) {
        # Compare LLR to breakpoint before
        topllr <- which.max(bp$llr[(i-1):i])
        # Skip if only two breaks detected, but too close to eachother (< w)
        if(i-2+topllr>0){ 
            bp[i-2+(c(2,1)[topllr]),] <- bp[i-2+topllr,]
        }
    }
    bp <- unique(bp)
    #cat("Before:",nrow(breakpoints),"After:",nrow(bp),"\n")
    #cat("Pruned ",nrow(breakpoints)-nrow(bp)," breaks (distance < ",w," bins)\n",sep="")
    
    # 2. Make segments
    seg <- data.frame(start=c(1, bp$local.max), end=c(bp$local.max, length(counts)), start.llr=c(Inf, bp$llr), end.llr=c(bp$llr, Inf))
    seg$size <- seg$end-seg$start
    seg$size[seg$start==1] <- seg$size[seg$start==1] + 1 # +1 if segments starts at chr boundary
    return(seg)
}

# Call segments
segs.byChr <- lapply(1:length(breakpoints.byChr), function(x){
    #cat("^^Chr",x,"\n")
    call.segs(breakpoints=breakpoints.byChr[[x]], counts=counts.fs[[x]], w=w)
})

# Combine to genome level
names(segs.byChr) <- names(counts.fs)
segs <- dplyr::bind_rows(segs.byChr, .id="chr")
segs <- mutate_at(segs, c("start","end","size"), as.integer)

# Figure out global bin indices
segs2 <- segs %>% 
    left_join(select(bins.f,chr,start_global=bin,bin_local,start_coord), by=c("chr", "start"="bin_local")) %>% # Add global start bin and coordinates
    left_join(select(bins.f,chr,end_global=bin,bin_local,end_coord), by=c("chr", "end"="bin_local")) %>% # Add global end bin
    mutate(size_bp=end_coord-start_coord)

# Scale at genome-level and return to bin-level matrix
normcounts <- gc.map.correct(counts.f, gc.f, map.f, gc.coefs)
segs_bin_coords <- select(segs2,start_global,end_global)
segs2$meanCov <- apply(segs_bin_coords, 1, function(x) {mean(normcounts[x[1]:x[2]])}) # Have to convert to global bin indices
segs2$meanCov.norm <- segs2$meanCov/mean(normcounts) # normalized to read depth

# 3. Return to bin-level values
seg.counts <- rep(segs2$meanCov, segs2$size)/mean(normcounts) # Cannot use median here, bins are too small. Maybe median seg$meanCov?

# Whole genome level: scaling factor and CNV-profile
s <- seq(min.scale, max.scale, length.out = 1000) # Reasonable ploidies/scaling factors
s.tests <- sapply(s, function(x) {manhattan.dist(x, seg.counts)})
if(all(is.na(s.tests))){ # If no scaling factor can be found
    best.s <- NA 
} else {
    best.s <- s[which.min(s.tests)]
}
CNV.profile <- round(seg.counts*best.s)

write_lines(seg.counts*best.s, outfile.seg_cn)
write_lines(normcounts, outfile.normcounts)
write_lines(s.tests, outfile.scalefactors)

# Add scaling factor and copy number integers to segments for saving
segs2$scaling_factor <- best.s
segs2$copy_number <- round(segs2$meanCov.norm*best.s)

# Add quality scores
calcQ <- function(counts.f, gc.f, map.f, l, r, S, C) {
    if(is.na(S)) return(NA) # If no scaling factor found, return NA
    S.10x <- mean(counts.f)/S
    norm.coef <- gc.map.correct(rep(1, length(l:r)), gc.f[l:r], map.f[l:r], gc.coefs)
    mu <- mean(counts.f[l:r]*norm.coef)
    stdev <- sd(counts.f[l:r]*norm.coef)
    val1 <- round(mean((C-0.5)*S.10x*norm.coef))
    val2 <- round(mean((C+0.5)*S.10x*norm.coef))
    #    return(-log10(1-(pnorm(q=val2, mean=mu, sd=stdev, lower.tail=TRUE) - pnorm(q=val1, mean=mu, sd=stdev, lower.tail=TRUE))))
    if(stdev > mu)
        return(-log10(1-sum(dnbinom(x=c(val1:val2), mu=mu, size=stdev))))
    else
        return(-log10(1-sum(dpois(x=c(val1:val2), mu))))
}

a <- apply(segs2[,c(7,9,15)], 1, function(x) { calcQ(counts.f, gc.f, map.f, x[1], x[2], best.s, x[3]) })
a[is.infinite(a) | is.na(a)] <- 0 # Added 250130(VZ): 0-count segments will produce infinite or NA scores and fail plotting
segs2$score <- a

seg.scores <- rep(segs2$score, segs2$size)

# Global coordinates for bins
bins.global <- bins %>% mutate(end_cum=cumsum(as.numeric(end_coord)-as.numeric(start_coord)))
norm.df <- cbind(bins.f,"counts.scaled.runmean"=runmean(normcounts,k=10)/mean(normcounts)*best.s) %>% 
    left_join(select(bins.global, bin_unfilt, end_cum), by="bin_unfilt")

# Plot at genome level
png(outfile.cnv,width=2000,height=1400,res=150)
mat <- matrix(c(1,2,2,3,3,4))
layout(mat)
par(mar=c(4,4,1.5,1.5))
#par(mfrow=c(3,1))
plottitle <- paste0(cell_id,"; avg reads/bin=",round(mean(counts.f),1)," (",basename(snakemake@input[["bins"]]),"); w=",w,"; scaling factor=",round(best.s,2),"; segs=",nrow(segs2))
plot(runmean(counts.f, k=10), type='p', pch=19, cex=0.2, main=plottitle, ylab="Raw counts", xlab="Bins")
if(!is.na(best.s)){
    # Plot 2
    plot(runmean(normcounts, k=10)/mean(normcounts)*best.s, type='p', cex=.1, col=alpha("black", 0.3), xlab="Bins", ylab="Copy number")
    points(CNV.profile, col='red', pch=19, cex=0.3)
    points(seg.counts*best.s, pch=19, cex=0.1, col='green')
    abline(v=unlist(lapply(split(bins.f$bin, bins.f$chr), max)), col="grey")
    text(x=unlist(lapply(split(bins.f$bin, bins.f$chr), function(x){round(mean(x))})),
         y=0, labels=sub("chr","",levels(bins.f$chr)), col="black")
    # Plot 3
    plot(norm.df$end_cum, norm.df$counts.scaled.runmean, pch=19, cex=.1, col=alpha("black", 0.3), xlab="Coordinate", ylab="Copy number")
    points(norm.df$end_cum, CNV.profile, col='red', pch=19, cex=0.3)
    points(norm.df$end_cum, seg.counts*best.s, pch=19, cex=0.1, col='green')
    abline(v=unlist(lapply(split(norm.df$end_cum, norm.df$chr), max)), col="grey")
    text(x=unlist(lapply(split(norm.df$end_cum, norm.df$chr), function(x){round(mean(x))})),
         y=0, labels=sub("chr","",levels(bins.f$chr)), col="black")
    # Plot 4
    plot(norm.df$end_cum, seg.scores, pch=19, cex=.1, col="orange", xlab="Coordinate", ylab="-log10(score)")
} else {
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("Scaling factor is NA. Scaled plot not possible"))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("Scaling factor is NA. Copy integer plot not possible"))
}
dev.off()

# Cleaner plot2 with ggplot
if(!is.na(best.s)){
    df <- tibble(bin=norm.df$bin, bp=norm.df$end_cum, bp_local=norm.df$end_coord, chr=norm.df$chr,
                 raw=norm.df$counts.scaled.runmean,
                 seg=seg.counts*best.s,
                 cn=CNV.profile)
    df.chr <- norm.df %>% 
        group_by(chr) %>% 
        summarize(line=max(end_cum),label_pos=round(mean(end_cum))) %>% 
        mutate(chr=sub("chr","",chr))
    ymax <- max(min(cn.trunc+0.5, max(df$cn)+1), 4) # Upper limit will be at least 4, maximum 8.5 and if always add 1 CN space above integer number to not cut off underlying raw data points
    
    brks <- c(-Inf,0,1,2,3,4,20,50,Inf)
    cols <- c("darkblue","blue","grey","red","darkred","brown","green","yellow")
    cols_vec2 <- setNames(as.character(cut(0:100, breaks = brks, labels = cols)),nm=0:100)
    
    p <- ggplot(df,aes(x=bp)) + 
        geom_vline(xintercept=df.chr$line, color="grey") +
        geom_hline(yintercept=2,color="grey",linetype="solid",size=0.1) +
        geom_point(aes(y=raw),pch=19,size=0.1,alpha=0.1) +
        geom_point(aes(y=seg),pch=19,size=0.15,col="orange") +
        geom_point(aes(y=cn,color=as.factor(cn)),pch=19,size=0.3,show.legend=F) +
        scale_color_manual(values=cols_vec2) +
        scale_fill_manual(values=cols_vec2) +
        scale_x_continuous(expand=c(0.01,0.01)) +
        lims(y=c(0,ymax)) + # Max ymax, min 4, otherwise dependent on actual copy number
        theme_dntr(axis_ticks = T) +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), plot.margin=margin(5,5,0,5)) +
        labs(x=NULL,y="copy number",
             subtitle=paste0(cell_id," (",paste0(round(sum(counts/1e3),0),"k read-pairs; ", as.integer(sub(".*-(.*).bed","\\1",basename(snakemake@input[["bins"]])))/1e3,"kb bins)")))
    if(any(df$cn>cn.trunc)) p <- p + geom_point(data=filter(df,cn>cn.trunc), aes(y=cn.trunc+0.5, fill=as.factor(cn)), color="black", pch=24,size=2.5,show.legend=F)
    pchr <- ggplot(df.chr) + geom_text(aes(x=label_pos,y=0,label=chr), check_overlap = T, hjust=0.5, vjust=0.2) +
        scale_x_continuous(expand=c(0.01,0.01), limits=c(min(df$bp),max(df$bp))) + theme_void() + theme(plot.margin=margin(0,0,0,0))
    
    pcomb <- p + pchr + plot_layout(ncol=1, heights = c(15,1)) + plot_annotation(theme = theme(plot.margin = margin(0,0,0,0)))
    
    # Outliers
    oi <- which(CNV.profile>cn.trunc)
    if(length(oi)>5){
        oi.breaks1 <- which(oi[2:length(oi)]!=oi[-length(oi)]+1)
        oi.breaks2 <- which(norm.df$chr[oi[2:length(oi)]]!=norm.df$chr[oi[-length(oi)]+1])
        oi.breaks <- union(oi.breaks1,oi.breaks2)
        oi.df <- data.frame(start=c(min(oi),oi[oi.breaks+1]),
                            end=c(oi[oi.breaks],max(oi)))
        oi.df$len <- (oi.df$end-oi.df$start)+1
        oi.df$avg_cn <- apply(oi.df, 1, function(x){ mean(CNV.profile[x[["start"]]:x[["end"]]]) })
        oi.df$start_pad <- oi.df$start-round(oi.df$len*0.1)
        oi.df$end_pad <- oi.df$end+round(oi.df$len*0.1)
        # TODO: Can't use padding unless implement a check for chr boundaries
        oi.s <- apply(oi.df, 1, function(x){ df[x["start"]:x["end"],] })
        
        # If massive number of outliers, will overwhelm plotting. Sample up to 20 of them instead
        cap <- NULL
        if(length(oi.s)>20){
            cap <- paste0(length(oi.s), " outliers total")
            oi.s <- oi.s[sample(1:length(oi.s),20)]
        }
        
        # Plot each outlier in list
        plist <- lapply(oi.s, function(df2){
            df2_lab <- paste0(df2$chr[1],":",df2$bp_local[1],"-",df2$bp_local[nrow(df2)])
            ggplot(df2,aes(x=bp_local)) + 
                geom_point(aes(y=raw),pch=19,size=0.5,alpha=0.5) +
                geom_point(aes(y=seg),size=0.5,col="orange") +
                geom_point(aes(y=cn,color=as.factor(cn)),size=1,show.legend=F) +
                scale_color_manual(values=cols_vec2) +
                scale_x_continuous(expand=c(0.01,0.01),breaks=pretty_breaks(2), labels=comma_format(scale=1e-6, suffix="Mb", big.mark="")) +
                scale_y_continuous(limits=c(max(cn.trunc,min(df2$raw)),max(df2$raw)),expand=c(0.05,0.05),
                                   breaks=pretty_breaks(5)) +
                theme_dntr(axis_ticks = T,font_size = 10) +
                labs(x=NULL,y=NULL,subtitle=df2_lab)
        })
        pout <- wrap_plots(plist,ncol=3)
        pfinal <- cowplot::plot_grid(pcomb + labs(caption=cap), pout,align="v",axis="lr",nrow=2, rel_heights = c(1.5,ceiling(length(oi.s)/3)))
        ggsave(snakemake@output$cnv2_png, pfinal, width=1200,height=300+(ceiling(length(oi.s)/3)*200), unit="px",dpi=100)
    } else {
        # No outliers
        ggsave(snakemake@output$cnv2_png, pcomb, width=1200,height=300,units = "px", dpi=100)
    }
} else {
    # No data to plot
    png(outfile.cnv2,width=1200,height=300,res=100)
    par(mar=c(2,4,0.5,1.5))
    plottitle <- paste0(cell_id,"; ",round(sum(counts/1e3),0),"k reads","; ",as.integer(sub(".*-(.*).bed","\\1",basename(snakemake@input[["bins"]])))/1e3,"kb bins; scaling factor=",round(best.s,2),"; segs=",nrow(segs2))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("Scaling factor is NA. Scaled plot not possible"))
    dev.off()
}

# Save profiles
cat("Save seg and CN profile...\n")
head(CNV.profile)
write_lines(x=c(cell_id,CNV.profile), path=outfile.cn)
head(segs2)
write_tsv(x=segs2, path=outfile.seg)
cat("Done.\n")
