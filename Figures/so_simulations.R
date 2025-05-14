library(tidyverse)
library(rtracklayer)

frac<-0.005

files <- system("ls /wrk/data/human_genomic/aligned/VZA02704*/*bed.gz", intern=T) #Single cell fragment file 
files <- files[!grepl(".scp.", files)]
numfrags <- sapply(files, function(file) {
    length(import(file))
})
good.files <- files[numfrags > 70000]

set.seed(123)
rand.files <- sample(good.files, 100, replace=F)

rand.files.frags <- lapply(rand.files, function(file) {
    lengths(import(file))
})

makeDopedCovStats <- function(frags, seg.ploidy, duprate=0.0, niter=20, seqlen=10E6) {
    my.gr <- makeGranges(frags=frags, ploidy=seg.ploidy, duprate=duprate, seqlen=seqlen)
    cov2 <- sapply(1:niter, function(i) {
        my.cov <- coverage(sample(my.gr, round(seqlen/mean(frags)*frac), replace=F))[[1]]
        c(sum(width(my.cov[my.cov==1])), sum(width(my.cov[my.cov==2])), sum(width(my.cov[my.cov==3])))
    })
    rowMeans(cov2)/seqlen
}

makeGranges <- function(frags, ploidy, duprate=0, seqlen=10E6) {
    my.gr <- GRanges()
    for(i in 1:ploidy) {
        fragsites <- cumsum(sample(frags, replace=T, seqlen/mean(frags)))
        my.gr <- c(my.gr, GRanges(seqnames='1', ranges=IRanges(start=c(1, fragsites[-length(fragsites)]), end=fragsites-1)))
        num.dups <- round(duprate * length(my.gr))
#        cat(length(my.gr), "\t")
        my.gr <- c(my.gr, sample(my.gr, num.dups, replace=T))
#        cat(length(my.gr), "\n")
    }
    my.gr
}

runploidy <- function(fragl.list, seqlen, dr, ploidies=1:4) {
    res <- lapply(1:4, function(seg.ploidy) {
        cat(paste("\nploidy", seg.ploidy, "\n"))
        unlist(t(as.data.frame(lapply(fragl.list, function(frags) {
            cat(".")
            res <- makeDopedCovStats(frags=frags, seg.ploidy=seg.ploidy, duprate=dr, seqlen=seqlen)
        })))[,2])
    })
    cat("\n")
    res
}

gsize.1e9.dr0.1 <- runploidy(rand.files.frags, seqlen=1e9, dr=0.1, ploidies=1:4)
gsize.1e8.dr0.1 <- runploidy(rand.files.frags, seqlen=1e8, dr=0.1, ploidies=1:4)
gsize.1e7.dr0.1 <- runploidy(rand.files.frags, seqlen=1e7, dr=0.1, ploidies=1:4)
gsize.1e6.dr0.1 <- runploidy(rand.files.frags, seqlen=1e6, dr=0.1, ploidies=1:4)


gsize.1e9.dr0.15 <- runploidy(rand.files.frags, seqlen=1e9, dr=0.15, ploidies=1:4)
gsize.1e8.dr0.15 <- runploidy(rand.files.frags, seqlen=1e8, dr=0.15, ploidies=1:4)
gsize.1e7.dr0.15 <- runploidy(rand.files.frags, seqlen=1e7, dr=0.15, ploidies=1:4)
gsize.1e6.dr0.15 <- runploidy(rand.files.frags, seqlen=1e6, dr=0.15, ploidies=1:4)

gsize.1e9.dr0.05 <- runploidy(rand.files.frags, seqlen=1e9, dr=0.05, ploidies=1:4)
gsize.1e8.dr0.05 <- runploidy(rand.files.frags, seqlen=1e8, dr=0.05, ploidies=1:4)
gsize.1e7.dr0.05 <- runploidy(rand.files.frags, seqlen=1e7, dr=0.05, ploidies=1:4)
gsize.1e6.dr0.05 <- runploidy(rand.files.frags, seqlen=1e6, dr=0.05, ploidies=1:4)

save(list=c("rand.files.frags", ls()[grep("gsize*", ls())]), file="results_01.rda")

beeswarm(c(gsize.1e6.dr0.05, gsize.1e7.dr0.05, gsize.1e8.dr0.05, gsize.1e9.dr0.05), pch=19, cex=0.2, ylim=c(0,1.5e-5))

boxplot(c(gsize.1e6.dr0.05, gsize.1e7.dr0.05, gsize.1e8.dr0.05, gsize.1e9.dr0.05), pch=19, cex=0.2, ylim=c(0,1.5e-5))

