library(tidyverse)
library(RColorBrewer)

seur<-readRDS("/wrk/resources/dntr2/results/DNA_150_SKE_newRNA/rna_seurat.Rds")
clones<-readRDS("/wrk/resources/dntr2/results/DNA_150_SKE_newRNA/clones/enge13final_clone_object.Rds")
clones <- clones$cells %>% 
  filter(!is.na(clone_final)) %>% 
  select(cell=dna_library_id, clone=clone_final)
clones
annot<-as.data.frame(clones)
annot
rownames(annot) <- (annot$cell)
rownames(annot)<-gsub("D_", "\\.", rownames(annot))
normal.dna<-c("5_1")
blast.dna<-setdiff(annot$clone, "5_1")
#blast.dna<-c("ALL66_2", "ALL66_3", "ALL66_4")

# 
# annot<-read.table("/wrk/resources/dntr/results/ALL_DNA_37_240904/clones/ALL71-clones-20000.txt", header=TRUE)
# normal.dna<-c("ALL71_1")
# blast.dna<-c("ALL71_2", "ALL71_3", "ALL71_4", "ALL71_5", "ALL71_6")

# seur<-readRDS("/wrk/resources/dntr/results/DNA_150_SKF_240903/rna_seurat.Rds")
# annot<-read.table("/wrk/resources/dntr/results/DNA_150_SKF_240903/clones/mel308_OP-clones-20000.txt", header=TRUE)
# normal.dna<-c("mel308_OP_1")
# blast.dna<-c("mel308_OP_2", "mel308_OP_3", "mel308_OP_4", "mel308_OP_5", "mel308_OP_6", "mel308_OP_7")




dta <- as.matrix(seur@assays$RNA@counts)
dta.log <- as.matrix(seur@assays$RNA@data)

# clones <- annot$dna_class
# names(clones) <- gsub("D_", ".", annot$dna_library_id)
dta.log <- dta.log[,colnames(dta.log) %in% rownames(annot)]
annot <- annot[colnames(dta.log),]
dim(annot)
dim(dta.log)

#clones <- clones[clones %in% c(normal.dna, blast.dna)]
clones <- clones[clones %in% blast.dna]
dta.log <- dta.log[,rownames(annot)]

clones<-annot$clone
names(clones)<-rownames(annot)

clone.aggrs <- apply(2^dta.log, 1, function(x) {
  sapply(split(x, clones), mean)
})

mxs <- apply(clone.aggrs, 2, max)
mean(mxs > 10)
o <- order(mxs, decreasing=T)[1:5000]
o.n <- names(sort(mxs, decreasing=T)[1:2000])

# cl.dst <- dist(clone.aggrs[,o])
# cl.dst <- cor(t(clone.aggrs[,o]), method="spearman")
# cl.dst <- cor(t(log2(clone.aggrs[,o]+1)), method="pearson")
# cl.dst <- cor(t(clone.aggrs[,seur@assays$RNA@var.features]), method="spearman")
# 
# # Variable features + highly expressed
# cl.dst <- cor(t(clone.aggrs[,unique(o.n, seur@assays$RNA@var.features)]), method="spearman")
# 
# #aa <- t(clone.aggrs[,all.seur@assays$RNA@var.features])
# 
# #png(width=1000, height=1000)
# #image(cl.dst)
# heatmap(cl.dst, scale='none', col=colorRampPalette(c("white", "red"))(100))
# #dev.off()
# 
# pca <- prcomp(1-cl.dst)
# 
# cl.csc <- cmdscale(as.dist(1-cl.dst), 2) # For correlation
# #cl.csc <- cmdscale(cl.dst, 2)
# library(RColorBrewer)
# my.cols <- brewer.pal(12, "Set3")[as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]}))]
# plot(cl.csc[,2:1], col=my.cols, pch=19, cex=1.5, main="Pearson cor log2(vals)")
# #identify(cl.dst, labels=rownames(cl.dst))
# legend("topright", fill=brewer.pal(12, "Set3")[unique(as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]})))], legend=unique(as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]}))))

#Un logging the data..
dta.filt <- 2^dta.log

sampsize  <- 5
numsamp <- 100

#Here including the normal cells..

# Test 1. Randomly select genes within loop - 
#What it's doing is per gene it splits by clone - then it samples 200 times the mean of that gene from 5 cells per clone
#This takes timmeeee 
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
dim(samp.clone.aggrs) #But this is 2010 - - not 2200 - because we skip a clone?? or what? Yeah we skip enge13_1 
#samp.clone.aggrs
# Remove normal cells?
#samp.clone.aggrs <- samp.clone.aggrs[!gsub("\\..+", "", rownames(samp.clone.aggrs)) %in% normal.dna,] #Martin had commented this out 


# samp.clone.aggrs.e13 <- samp.clone.aggrs #Martin had commented this out 
# samp.clone.aggrs.final <- samp.clone.aggrs[!gsub("\\..+", "", rownames(samp.clone.aggrs)) %in% normal.dna,] #Martin had commented this out 
# o <- match(colnames(samp.clone.aggrs), colnames(samp.clone.aggrs.e13))
# samp.clone.aggrs.final <- rbind(samp.clone.aggrs[,!is.na(o)], samp.clone.aggrs.e13[,o[!is.na(o)]])
# samp.clone.aggrs.final <- samp.clone.aggrs.final[!gsub("\\..+", "", rownames(samp.clone.aggrs.final)) %in% c("enge13_1", "src20_1"),]

#Should we remove normal cells? 

samp.clone.aggrs.final <- samp.clone.aggrs.final[!gsub("\\..+", "", rownames(samp.clone.aggrs.final)) %in% c("5_1"),]

#And here we are to two times the original - why is that? 
#Lets skip that and just put it differently 
#samp.clone.aggrs.final<-samp.clone.aggrs

dim(samp.clone.aggrs.final) #[1]  4020 26577
dim(samp.clone.aggrs)
mxs <- apply(samp.clone.aggrs.final, 2, max)
head(mxs)
#What is the highest value per gene from these subsampled means?
#mxs <- apply(cl.tp.aggrs, 2, max)
mean(mxs > 10)
o <- order(mxs, decreasing=T)[1:5000] #Take the top 5000 expressed ish 
#o.n <- names(sort(mxs, decreasing=T)[1:2000])

#cl.dst <- dist(samp.tp.aggrs[,o])
##cl.dst <- cor(t(samp.tp.aggrs[,o]), method="pearson")
#cl.dst <- cor(t(samp.clone.aggrs.final[,o]), method="spearman") #Then if we use this method here 
#cl.dst <- cor(t(samp.tp.aggrs[,seur@assays$RNA@var.features]), method="spearman") 
cl.dst <- cor(t(samp.clone.aggrs.final[,unique(o, seur@assays$RNA@var.features)]), method="spearman") 
cl.dst
class(cl.dst)
heatmap(cl.dst[1:11, 1:11])

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

my.cols <- paste0(sample(jet.colors(length(levels(labels))))[labels], "11")
#my.pch <- as.factor(sapply(strsplit(rownames(cl.csc), "\\."), function(x) {x[2]}))
my.pch <- 19

title <- gsub("_.+", "", normal.dna[1])
#pdf(paste0("clouds_", title, ".pdf"), width=7, height=7.5)
plot(cl.csc, col=my.cols, pch=as.integer(my.pch), cex=0.5, lwd=0, main=title)

#But then like how do we get to a value of like what 
#Because the clouds are very .. small 
#Like can I get them to . like.. 
#eh

cl.csc[1:2,1:2]

#I think what we were saying is if I find the pairwise distance between two cells in two separate clones based on their expression - except i check like 200 times..
#But then on different gene sets? 

# dta.filt[o,]
# samp.clone.aggrs.final
# #So thats each clone....
# dta.log
# Idents(seur)<-seur$cell_id
# ids<-seur$cell_id[seur$cell_id%in%names(clones)]
# ids%in%names(clones)
# Idents(object = seur)
# levels(x = pbmc)
# 
# seur<-subset(seur, idents=unlist(ids))
# 
# DimPlot(seur, group.by="plate_id")
# #But yeah ok so if we do this for random subset 
# downsampled.obj <- seur[, sample(colnames(seur), size =20, replace=F)]
# downsampled.obj
# downsampled.obj <- FindVariableFeatures(downsampled.obj)
# downsampled.obj <- RunPCA(downsampled.obj)
# downsampled.obj <- FindNeighbors(downsampled.obj, dims=1:30)
# downsampled.obj <- FindClusters(downsampled.obj, resolution = 0.5)
# downsampled.obj <- RunUMAP(downsampled.obj, dims = 1:30, n.neighbors = 5)
# DimPlot(downsampled.obj)
# #
# #OK no - wait 
# dta.log <- as.matrix(downsampled.obj@assays$RNA@data)
# #Then I need like a dataframe with distances - so if i do just pairwise distances?? 
# cl.dst <- cor((dta.log[o,]), method="spearman")#Or then use one of martins correlation things?? 
# cl.dst
# 
# kNNdist(cl.dst, 10)
# 
# nn_indices <- get.knnx(cl.dst, cl.dst, k = 10)$nn.index
# nn_indices
# 
# #Ok so - what if I would then take always the same amount of cells per clone 
# # Split the data by the `enge13_*` values
# split_data <- split(names(clones), clones)
# 
# # Function to randomly sample 6 IDs for each group
# random_sample <- function(ids) {
#   if (length(ids) > 6) {
#     return(sample(ids, 6))
#   } else {
#     return(ids)
#   }
# }
# 
# # Apply the function to each group
# subsampled_data <- lapply(split_data, random_sample)
# random<-unlist(subsampled_data, use.names=FALSE)
# dta.log<-as.matrix(seur@assays$RNA@data)
# cl.dst <- cor((dta.log[o,random]), method="spearman")#
# nn_indices <- get.knnx(cl.dst, cl.dst, k = 6)$nn.index
# head(nn_indices)
# head(cl.dst)
# 
# col_names <- colnames(cl.dst)
# # Use the indices from nn_indices to get the corresponding column names
# nn_colnames <- apply(nn_indices, 2, function(idx) col_names[idx])
# # Display the mapped column names
# head(nn_colnames)
# 
# #ok and then count how many are the correct clone?
# head(nn_colnames) 
# head(clones)
# 
# nn_clones <- apply(nn_colnames, 2, function(x) clones[x])
# print(nn_clones) # 
# #And then count the proportion of the first one 
# 
# count_matches <- function(matrix_data) {
#   # Apply a function to each row of the matrix
#   apply(matrix_data, 1, function(row) {
#     # Compare columns 2:6 with the first column and sum up the matches
#     sum(row[2:6] == row[1])/5*100
#   })
# }
# 
# count_matches(nn_clones)
# 
# grouped_proportion_matches <- function(matrix_data) {
#   # Count matches for each row
#   match_counts <- apply(matrix_data, 1, function(row) {
#     sum(row[2:6] == row[1])
#   })
#   
#   # Create a data frame with the first column and the match counts
#   df <- data.frame(first_col = matrix_data[, 1], match_counts = match_counts)
#   
#   # Calculate the total possible matches for each group
#   df$total_comparisons <- 5
#   
#   # Calculate the sum of match counts and total comparisons, grouped by the first column
#   grouped_sums <- aggregate(cbind(match_counts, total_comparisons) ~ first_col, data = df, sum)
#   
#   # Calculate the proportion of matches
#   grouped_sums$proportion <- grouped_sums$match_counts / grouped_sums$total_comparisons
#   
#   # Return the grouped proportions
#   return(grouped_sums[, c("first_col", "proportion")])
# }
# grouped_proportion_matches(nn_clones)
# 
# 
# #Yeah so maybe if I do this 100 times 
# 
# #So if I do a umap of 100 cells - based on variable features? and then check for each cell the k nearest neighbors 
# 
# 
# #identify(cl.csc, labels=rownames(cl.dst))
# #legend("topright", fill=brewer.pal(12, "Set3")[unique(as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]})))], legend=unique(as.factor(sapply(strsplit(rownames(cl.csc), "_"), function(x) {x[1]}))))
# 
# colnames(cl.csc) <- c('x', 'y')
# mids <- cl.csc %>% as.data.frame %>% mutate(labels = labels) %>% group_by(labels) %>% summarize(mx=median(x), my=median(y))
# text(x=mids$mx, y=mids$my, labels=mids$labels, cex=1)
# #mids <- cl.csc[1:length(levels(as.factor(clones))),] %>% as.data.frame %>% mutate(labels=labels[1:length(levels(as.factor(clones)))])
# #points(x=mids$x, y=mids$y, pch=19, cex=0.3)
# #text(x=mids$x, y=mids$y, labels=mids$labels, cex=0.5)
# 
# dev.off()

# Make ridge clouds

library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggrepel)

#pdf("Clouds2_enge13.pdf",width=5, height=5)
#dev.off()
colnames(cl.csc) <- c('x', 'y')
cl.2 <- cl.csc %>% as.data.frame %>% mutate(labels = labels)

label_distances <- cl.2 %>%
  group_by(labels) %>%
  summarise(
    center_x = median(x),
    center_y = median(y),
    max_distance = max(sqrt((x - median(x))^2 + (y - median(y))^2))  # Euclidean distance
  )

# Print the result
print(label_distances[order(label_distances$max_distance, decreasing=TRUE),])



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
#Lets figure out what the ridge values mean - like can we get some sort of distance measure that would be similar to entropy?

#dev.off()

#cl.csc3 <- cmdscale(as.dist(1-cl.dst), 3) # For correlation

#library("scatterplot3d")
#scatterplot3d(cl.csc3, color=my.cols, angle=55)
#legend('topleft', legend=levels(labels), fill=unique(my.cols))

