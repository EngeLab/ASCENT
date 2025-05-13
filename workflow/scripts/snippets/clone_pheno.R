library(tidyverse)
library(RColorBrewer)

seur <- readRDS("./rna_seurat.Rds")
annot <- read.table("enge13-clones-20000.txt", header=T)
normal.dna <- c("enge13_1")
blast.dna <- c("enge13_10", "enge13_11", "enge13_12", "enge13_13",
               "enge13_14", "enge13_6", "enge13_7", "enge13_8", "enge13_9")

seur <- readRDS("./enge20/rna_seurat.Rds")
annot <- read.table("enge20/enge20-clones-20000.txt", header=T)
normal.dna <-c("enge20_1", "enge20_2")
tumor.dna <-c("enge20_3", "enge20_4", "enge20_5", "enge20_6", "enge20_7", "enge20_8",
                                   "enge20_9")

seur <- readRDS("./src20/rna_seurat.Rds")
annot <- read.table("src20/src20-clones-20000.txt", header=T)
normal.dna<-c("src20_1")
blast.dna<-c("src20_2", "src20_3", "src20_4", "src20_5", "src20_6", "src20_7")

dta <- as.matrix(seur@assays$RNA@counts)
dta.log <- as.matrix(seur@assays$RNA@data)

clones <- annot$dna_class_merge
names(clones) <- gsub("D_", ".", annot$dna_library_id)
dta.log <- dta.log[,colnames(dta.log) %in% names(clones)]
clones <- clones[colnames(dta.log)]

clones <- clones[clones %in% c(normal.dna, blast.dna)]
#clones <- clones[clones %in% blast.dna]
dta.log <- dta.log[,names(clones)]

clone.aggrs <- apply(2^dta.log, 1, function(x) {
    sapply(split(x, clones), mean)
})

mxs <- apply(clone.aggrs, 2, max)
mean(mxs > 10)
o <- order(mxs, decreasing=T)[1:5000]
o.n <- names(sort(mxs, decreasing=T)[1:2000])

cl.dst <- dist(clone.aggrs[,o])
cl.dst <- cor(t(clone.aggrs[,o]), method="spearman")
cl.dst <- cor(t(log2(clone.aggrs[,o]+1)), method="pearson")
cl.dst <- cor(t(clone.aggrs[,seur@assays$RNA@var.features]), method="spearman")

# Variable features + highly expressed
cl.dst <- cor(t(clone.aggrs[,unique(o.n, seur@assays$RNA@var.features)]), method="spearman")

#aa <- t(clone.aggrs[,all.seur@assays$RNA@var.features])

png(width=1000, height=1000)
#image(cl.dst)
heatmap(cl.dst, scale='none', col=colorRampPalette(c("white", "red"))(100))
dev.off()

pca <- prcomp(1-cl.dst)

cl.csc <- cmdscale(as.dist(1-cl.dst), 2) # For correlation
#cl.csc <- cmdscale(cl.dst, 2)
library(RColorBrewer)
my.cols <- brewer.pal(12, "Set3")[as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]}))]
plot(cl.csc[,2:1], col=my.cols, pch=19, cex=1.5, main="Pearson cor log2(vals)")
#identify(cl.dst, labels=rownames(cl.dst))
legend("topright", fill=brewer.pal(12, "Set3")[unique(as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]})))], legend=unique(as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]}))))


dta.filt <- 2^dta.log

sampsize  <- 5
numsamp <- 200
# Test 1. Randomly select genes within loop.
samp.clone.aggrs <- apply(dta.filt, 1, function(x) {
    sapply(split(x, clones), function(x) {
        sapply(1:numsamp, function(y) {
            mean(sample(x, sampsize))
        })
    })
})
rownames(samp.clone.aggrs) <- paste(rep(rownames(clone.aggrs), each=numsamp), 1:numsamp, sep='.')
tmp <- clone.aggrs
rownames(tmp) <- paste0(rownames(clone.aggrs), ".0")
samp.clone.aggrs <- rbind(tmp, samp.clone.aggrs)
rm(tmp)

# Remove normal cells?
#samp.clone.aggrs <- samp.clone.aggrs[!gsub("\\..+", "", rownames(samp.clone.aggrs)) %in% normal.dna,]

#samp.clone.aggrs.e13 <- samp.clone.aggrs
#samp.clone.aggrs.final <- samp.clone.aggrs[!gsub("\\..+", "", rownames(samp.clone.aggrs)) %in% normal.dna,]
o <- match(colnames(samp.clone.aggrs), colnames(samp.clone.aggrs.e13))
samp.clone.aggrs.final <- rbind(samp.clone.aggrs[,!is.na(o)], samp.clone.aggrs.e13[,o[!is.na(o)]])
samp.clone.aggrs.final <- samp.clone.aggrs.final[!gsub("\\..+", "", rownames(samp.clone.aggrs.final)) %in% c("enge13_1", "src20_1"),]

mxs <- apply(samp.clone.aggrs.final, 2, max)
#mxs <- apply(cl.tp.aggrs, 2, max)
mean(mxs > 10)
o <- order(mxs, decreasing=T)[1:5000]
o.n <- names(sort(mxs, decreasing=T)[1:2000])

cl.dst <- dist(samp.tp.aggrs[,o])
cl.dst <- cor(t(samp.tp.aggrs[,o]), method="pearson")
cl.dst <- cor(t(samp.clone.aggrs.final[,o]), method="spearman")
cl.dst <- cor(t(samp.tp.aggrs[,seur@assays$RNA@var.features]), method="spearman")
cl.dst <- cor(t(samp.clone.aggrs.final[,unique(o.n, seur@assays$RNA@var.features)]), method="spearman")

cl.csc <- cmdscale(as.dist(1-cl.dst), 2) # For correlation
#cl.csc <- cmdscale(cl.dst, 2)
#my.cols <- paste0(brewer.pal(12, "Set3")[as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]}))], "11")
transp <- 'FF'
transp <- '22'
categs <- as.factor(gsub("\\..+", "", rownames(cl.csc)))
my.cols <- paste0(brewer.pal(12, "Set3")[categs], transp) # Make transparent colors
#my.cols <- paste0(brewer.pal(12, "Set3"))
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
if(length(levels(categs)) > 12) {
    my.cols <- paste0(jet.colors(length(levels(categs)))[as.factor(gsub("\\..+", "", rownames(cl.csc)))], transp) # Make transparent colors
}

#labels <- as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]}))  # Patient
labels <- as.factor(gsub(".+_", "", sapply(strsplit(rownames(cl.csc), "\\."), function(x) {x[1]}))) # Clone
labels <- as.factor(sapply(strsplit(rownames(cl.csc), "\\."), function(x) {x[1]})) # Clone+patient

#my.cols <- paste0(sample(jet.colors(length(levels(labels))))[labels], "11")
#my.pch <- as.factor(sapply(strsplit(rownames(cl.csc), "\\."), function(x) {x[2]}))
my.pch <- 19

title <- gsub("_.+", "", normal.dna[1])
pdf(paste0("clouds_", title, ".pdf"), width=7, height=7.5)
plot(cl.csc, col=my.cols, pch=as.integer(my.pch), cex=0.5, lwd=0, main=title)

#identify(cl.csc, labels=rownames(cl.dst))
#legend("topright", fill=brewer.pal(12, "Set3")[unique(as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]})))], legend=unique(as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]}))))

colnames(cl.csc) <- c('x', 'y')
mids <- cl.csc %>% as.data.frame %>% mutate(labels = labels) %>% group_by(labels) %>% summarize(mx=median(x), my=median(y))
text(x=mids$mx, y=mids$my, labels=mids$labels, cex=1)
#mids <- cl.csc[1:length(levels(as.factor(clones))),] %>% as.data.frame %>% mutate(labels=labels[1:length(levels(as.factor(clones)))])
#points(x=mids$x, y=mids$y, pch=19, cex=0.3)
#text(x=mids$x, y=mids$y, labels=mids$labels, cex=0.5)

dev.off()

# Make ridge clouds

library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggrepel)

pdf("Clouds2_enge13.pdf",width=5, height=5)

colnames(cl.csc) <- c('x', 'y')
cl.2 <- cl.csc %>% as.data.frame %>% mutate(labels = labels)

scatterplot <- cl.2 %>% filter(!(labels %in% normal.dna)) %>% ggplot(aes(x, y, color = labels)) +
#  geom_point() +
    geom_density_2d(binwidth=NULL) +
    theme_classic() +
theme(axis.text = element_blank(), axis.title = element_blank())

scatterplot <- scatterplot + theme(legend.position = "none")
mids <- cl.2 %>% as.data.frame %>% filter(!(labels %in% normal.dna)) %>% mutate(labels = labels) %>% group_by(labels) %>% summarize(mx=median(x), my=median(y))
scatterplot <- scatterplot + geom_point(data=mids, aes(x=mx, y=my, col='black')) + geom_label_repel(data=mids, aes(x=mx, y=my, label=labels), box.padding=0.5)


#scatterplot <- scatterplot + geom_point(x=mids$mx, y=mids$my, labels=mids$labels, cex=0.5)
scatterplot

dev.off()

cl.csc3 <- cmdscale(as.dist(1-cl.dst), 3) # For correlation

library("scatterplot3d")
scatterplot3d(cl.csc3, color=my.cols, angle=55)
legend('topleft', legend=levels(labels), fill=unique(my.cols))
