library(scater)
library(tidyverse)
library(dplyr)
library(readr)

target<-as.numeric(snakemake@params[["downsample"]])
cell_bins<-read.table(snakemake@input[["cell"]], header=TRUE)
good_bins<-read.table(snakemake@input[["good_bins"]], header=TRUE)

cell_bins_f<-cell_bins[good_bins$bin_unfilt,]
sums<-sum(cell_bins_f)
props <- pmin(1, target / sums)  # Ensure we don't exceed available counts

downsampled_matrix <- downsampleMatrix(as.matrix(cell_bins_f), prop = props, bycol = TRUE)
counts.f<-as.matrix(downsampled_matrix)
dim(counts.f)
sum(counts.f[,1])
colnames(counts.f)<-colnames(cell_bins)
write_tsv(as.data.frame(counts.f), file=snakemake@output[["cell_out"]])



  
