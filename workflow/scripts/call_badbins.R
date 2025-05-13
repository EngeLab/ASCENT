# Find bad bins
library(readr)
library(dplyr)
library(caTools)
library(data.table)
options(scipen=999)

counts.file <- snakemake@input[["counts"]]
map <- as.numeric(readLines(snakemake@input[["map"]]))
gc <- as.numeric(readLines(snakemake@input[["gc"]]))
bins <- read.table(snakemake@input[["bins"]],col.names=c("chr","start_coord","end_coord"),stringsAsFactors=F) # Use this for splitting
gc_min <- snakemake@params[["gc_min"]]
map_min <- snakemake@params[["map_min"]]
badbins.bed <- snakemake@output[["badbins"]]
badbins.pdf <- snakemake@output[["pdf"]]
cen <- NULL
gap <- NULL

bins$chr <- factor(bins$chr, levels=unique(bins$chr)) # Maintain natural sort order for chromosomes
bins$bin_unfilt <- seq_along(bins$chr)
bins.global <- bins %>% mutate(end_cum=cumsum(as.numeric(end_coord)-as.numeric(start_coord)))

# Drop bins with low mappability and extreme GC values
# Or use all of them (badbins and low map/high GC will overlap then)
good.bins <- gc > gc_min & map > map_min
bins.f <- bins[good.bins,]
map.f <- map[good.bins]
gc.f <- gc[good.bins]

# Load countsfile with normal cells
counts <- as.matrix(data.table::fread(counts.file,data.table=F))
counts.f <- counts[good.bins,]
dim(counts.f)

# Calculate global bin coordinates for bins/gaps etc
bnds.g <- unlist(lapply(split(bins$bin_unfilt, bins$chr), max))
bnds.coord <- unlist(lapply(split(bins.global$end_cum,bins.global$chr),max))
bnds_col <- c(rep(c("#FFFFFF","#DDDDDD"),12),"#FFFFFF")

cum_addon <- c(0,bnds.coord[-length(bnds.coord)])
names(cum_addon) <- paste0("chr",c(1:22,"X","Y"))

# Centromeres and gaps
if(!is.null(snakemake@params[["centromeres"]])){
  cen <- read_tsv(snakemake@params[["centromeres"]], col_names=c("chr","start","end"))  
  cen$start_cum <- cen$start + c(0,bnds.coord[-length(bnds.coord)])
  cen$end_cum <- cen$end + c(0,bnds.coord[-length(bnds.coord)])
}
if(!is.null(snakemake@params[["genome_gaps"]])){
  gap <- read_tsv(snakemake@params[["genome_gaps"]], col_names=c("chr","start","end"))
  gap$start_cum <- gap$start + cum_addon[gap$chr]
  gap$end_cum <- gap$end + cum_addon[gap$chr]
}

# Poisson method for bad bins
ppois.per.bin <- function(bcounts, map, gc, qest, gc.scale, lower.tail=TRUE) {
  weight <- map*predq(qest, gc, gc.scale)
  expected <- mean(bcounts)/weight
  return(ppois(q=bcounts, lambda=expected, lower.tail=lower.tail, log.p=TRUE))
}
predq <- function(qest, gc, scale=1) {
  (qest[1]+gc*qest[2]+gc^2*qest[3])/scale
}
compute_badbins <- function(counts, map, gc){
  map.norm <- map/mean(map)
  cat("Calc upper tail...")
  b <- apply(counts, 2, function(x) {
    gc.coefs <- summary(lm(x ~ gc + I(gc^2)))$coefficients[,1]
    gc.scale <- predq(gc.coefs, 0.45, 1)
    ppois.per.bin(x, map.norm, gc, gc.coefs, gc.scale, lower.tail=F)
  })
  cat("done.\n")
  cat("Calc lower tail...")
  bb <- apply(counts, 2, function(x) {
    gc.coefs <- summary(lm(x ~ gc + I(gc^2)))$coefficients[,1]
    gc.scale <- predq(gc.coefs, 0.45, 1)
    ppois.per.bin(x, map.norm, gc, gc.coefs, gc.scale, lower.tail=T)
  })
  cat("done.\n")
  bsub <- rowSums(bb-b,na.rm=T)
  
  # Stats
  b.mean <- mean(-rowSums(b,na.rm=T))
  b.sd <- sd(-rowSums(b,na.rm=T))
  bsub.mean <- mean(bsub)
  bsub.sd <- sd(bsub)
  
  df <- data.frame("upper"=-rowSums(b,na.rm=T) > b.mean+b.sd*2 | -rowSums(b,na.rm=T) < b.mean-b.sd*2,
                   "twotailed"=bsub > bsub.mean+bsub.sd*2 | bsub < bsub.mean-bsub.sd*2,
                   "upper25sd"=-rowSums(b,na.rm=T) > b.mean+b.sd*2.5 | -rowSums(b,na.rm=T) < b.mean-b.sd*2.5,
                   "twotailed25sd"=bsub > bsub.mean+bsub.sd*2.5 | bsub < bsub.mean-bsub.sd*2.5)
  
  return(list("b"=b,"bb"=bb,"bsub"=bsub,
              "b.mean"=b.mean,"b.sd"=b.sd,"bsub.mean"=bsub.mean,"bsub.sd"=bsub.sd,
              "stats"=df))
}

# Split into auto, X and Y chrs
auto.idx <- bins.f$chr %in% paste0("chr",1:22)
x.idx <- bins.f$chr=="chrX"
idx.list <- list("autosomal"=auto.idx, "X"=x.idx)

# Autosomal bins (both male and female)
bb.auto <- compute_badbins(counts.f[auto.idx,], map.f[auto.idx], gc.f[auto.idx])
# chrX
bb.x <- compute_badbins(counts.f[x.idx,], map.f[x.idx], gc.f[x.idx])

pdf(badbins.pdf,width=14,height=8)
par(mfrow=c(3,2),mar=c(2,5,2,2))
bl <- list("autosomal"=bb.auto,"X"=bb.x)
for(i in names(bl)){
  b.mean <- bl[[i]][["b.mean"]]
  b.sd <- bl[[i]][["b.sd"]]
  bsub.mean <- bl[[i]][["bsub.mean"]]
  bsub.sd <- bl[[i]][["bsub.sd"]]
  # Upper tail
  drop <- paste0(sum(bl[[i]][["stats"]][,"upper25sd"])," of ",nrow(bl[[i]][["stats"]]))
  hist(-rowSums(bl[[i]][["b"]]), breaks=200, main=paste0(i," -log(p) lower.tail=F (",drop,")"))
  abline(v=b.mean,col="grey")
  abline(v=c(b.mean-b.sd*2,b.mean+b.sd*2),col="black",lty="dashed")
  abline(v=c(b.mean-b.sd*2.5,b.mean+b.sd*2.5),col="red",lty="dashed")
  abline(v=c(b.mean-b.sd*3,b.mean+b.sd*3),col="black",lty="dashed")
  # Two-tailed
  drop2 <- paste0(sum(bl[[i]][["stats"]][,"twotailed25sd"])," of ",nrow(bl[[i]][["stats"]]))
  hist(bl[[i]][["bsub"]], breaks=200,main=paste0(i," -log(p) bb-b (",drop2,")"))
  abline(v=bsub.mean,col="grey")
  abline(v=c(bsub.mean-bsub.sd*2,bsub.mean+bsub.sd*2),col="black",lty="dashed")
  abline(v=c(bsub.mean-bsub.sd*2.5,bsub.mean+bsub.sd*2.5),col="red",lty="dashed")
  abline(v=c(bsub.mean-bsub.sd*3,bsub.mean+bsub.sd*3),col="black",lty="dashed")
}

# Plot on genome coords with centromeres and gaps
par(mfrow=c(3,1),mar=c(2,5,2,2))
for(i in names(bl)){
  binset <- bins.f$bin_unfilt[idx.list[[i]]]
  badbins <- bl[[i]][["stats"]][,"twotailed25sd"]
  bsub.mean <- bl[[i]][["bsub.mean"]]
  bsub.sd <- bl[[i]][["bsub.sd"]]
  maxy=max(bl[[i]][["bsub"]],na.rm=T)*1.05
  miny=min(bl[[i]][["bsub"]],na.rm=T)*1.05
  x_breaks <- pretty(n=10, c(1,max(bnds.coord)))
  x_labels <- scales::number(x_breaks, big.mark = "", accuracy = 0.1, suffix="Gb", scale=1e-9, trim=T)
  title <- paste0(i, " exp vs obs counts/bin (bb-b, low.tail=T; n=",ncol(bl[[i]][["b"]])," diploid cells; n=",length(bl[[i]][["bsub"]])," bins)")
  plot(NULL, xlim=c(1,max(bnds.coord)), ylim=c(miny,maxy), ylab="-log(p)", xaxt="n", xlab=NA, main=title)
  axis(side=1,at=x_breaks,labels=x_labels)
  rect(xleft=c(1,bnds.coord[-length(bnds.coord)]),xright=bnds.coord,ybottom=miny,ytop=maxy,col=rep(c("#FFFFFF","#EEEEEE"),12), lwd=0)
  if(!is.null(gap)) rect(xleft=gap$start_cum, xright=gap$end_cum, ybottom=miny, ytop=maxy,col="lightblue",lwd=0) # gaps
  if(!is.null(cen)) rect(xleft=cen$start_cum, xright=cen$end_cum, ybottom=miny, ytop=maxy,col="orange",lwd=0) # centromeres
  points(x=bins.global[binset,"end_cum"], y=bl[[i]][["bsub"]], cex=0.1)
  # Show filtered bins
  points(x=bins.global[binset[badbins],"end_cum"], y=bl[[i]][["bsub"]][badbins],cex=0.1,col="red")
  abline(h=c(bsub.mean-bsub.sd*2.5,bsub.mean+bsub.sd*2.5),col="red",lty="dashed")
  abline(h=c(bsub.mean-bsub.sd*3,bsub.mean+bsub.sd*3),col="purple",lty="dashed")
}
dev.off()

# Merge to one set of bad bins auto, X, and Y 
bad.auto <- bins.f %>% filter(chr %in% paste0("chr",1:22)) %>% filter(bb.auto$stats$twotailed25sd)
bad.x <- bins.f %>% filter(chr=="chrX") %>% filter(bb.x$stats$twotailed25sd)
bad.all <- bind_rows(bad.auto,bad.x) %>% distinct() %>% arrange(bin_unfilt)

write_tsv(bad.all, badbins.bed)

# Only filtered/good bins: 
good.bins.final <- gc > gc_min & map > map_min & !bins$bin_unfilt %in% bad.all$bin_unfilt
bins.final <- bins[good.bins.final, ] %>% mutate(bin=1:nrow(.))
write_tsv(bins.final, snakemake@output[["goodbins"]])
