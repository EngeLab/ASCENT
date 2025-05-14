################################
#Downsampling - find difference from reference and plot
#Performed here for llr
library(ggbeeswarm)
library(dplyr)
library(tidyverse)

bins_500<-read.table("resources/goodbins-500000.bed", header=TRUE)
bins_200<-read.table("resources/goodbins-200000.bed", header=TRUE)
bins_80<-read.table("resources/goodbins-80000.bed", header=TRUE)
bins_40<-read.table("resources/goodbins-40000.bed", header=TRUE)
bins_10<-read.table("resources/goodbins-10000.bed", header=TRUE)

#Read in the bins that are found across all resolutions - as only these will be used to make the comparison fair
bins_500_final<-read.table("../Figures/Bins500_common.tsv", header=TRUE)
bins_200_final<-read.table("../Figures/Bins200_common.tsv", header=TRUE)
bins_80_final<-read.table("../Figures/Bins80_common.tsv", header=TRUE)
bins_40_final<-read.table("../Figures/Bins40_common.tsv", header=TRUE)
bins_10_final<-read.table("../Figures/Bins10_common.tsv", header=TRUE)


k500 <- list.files("results/Downsample_llr/ALL40/",
                   pattern = "ALL40-copynumber_llr-500000",
                   full.names = TRUE)
k500_expanded <- list()
for (file in k500) {
  a <- read.table(file, header = TRUE)
  bin_counts <- as.data.frame(table(bins_500_final$bin))  
  colnames(bin_counts) <- c("bin", "count")
  all_bins <- data.frame(bin = unique(bins_500$bin)) 
  bin_counts_full <- merge(all_bins, bin_counts, by = "bin", all.x = TRUE)
  bin_counts_full$count[is.na(bin_counts_full$count)] <- 0
  k500_cn_expanded <- a[rep(1:nrow(a), times = bin_counts_full$count), ]
  k500_expanded[[file]]<-k500_cn_expanded
}

k200 <- list.files("results/Downsample_llr/ALL40/",
                   pattern = "ALL40-copynumber_llr-200000",
                   full.names = TRUE)
k200_expanded <- list()
for (file in k200) {
  a <- read.table(file, header = TRUE)
  bin_counts <- as.data.frame(table(bins_200_final$bin))  
  colnames(bin_counts) <- c("bin", "count")
  all_bins <- data.frame(bin = unique(bins_200$bin)) 
  bin_counts_full <- merge(all_bins, bin_counts, by = "bin", all.x = TRUE)
  bin_counts_full$count[is.na(bin_counts_full$count)] <- 0
  k200_expanded[[file]] <- a[rep(1:nrow(a), times = bin_counts_full$count), ]
}

k80 <- list.files("results/Downsample_llr/ALL40/",
                   pattern = "ALL40-copynumber_llr-80000",
                   full.names = TRUE)


k80_expanded <- list()
for (file in k80) {
  a <- read.table(file, header = TRUE)
  bin_counts <- as.data.frame(table(bins_80_final$bin))  
  colnames(bin_counts) <- c("bin", "count")
  all_bins <- data.frame(bin = unique(bins_80$bin)) 
  bin_counts_full <- merge(all_bins, bin_counts, by = "bin", all.x = TRUE)
  bin_counts_full$count[is.na(bin_counts_full$count)] <- 0
  k80_expanded[[file]] <- a[rep(1:nrow(a), times = bin_counts_full$count), ]
}

k40 <- list.files("results/Downsample_llr/ALL40/",
                   pattern = "ALL40-copynumber_llr-40000",
                   full.names = TRUE)

k40_expanded <- list()
for (file in k40) {
  a <- read.table(file, header = TRUE)
  bin_counts <- as.data.frame(table(bins_40_final$bin))  
  colnames(bin_counts) <- c("bin", "count")
  all_bins <- data.frame(bin = unique(bins_40$bin)) 
  bin_counts_full <- merge(all_bins, bin_counts, by = "bin", all.x = TRUE)
  bin_counts_full$count[is.na(bin_counts_full$count)] <- 0
  k40_expanded[[file]] <- a[rep(1:nrow(a), times = bin_counts_full$count), ]
  
}

k10 <- list.files("results/Downsample_llr/ALL40/",
                   pattern = "ALL40-copynumber_llr-10000",
                   full.names = TRUE)

k10_expanded <- list()
for (file in k10) {
  a <- read.table(file, header = TRUE)
  bin_counts <- as.data.frame(table(bins_10_final$bin))  
  colnames(bin_counts) <- c("bin", "count")
  all_bins <- data.frame(bin = unique(bins_10$bin)) 
  bin_counts_full <- merge(all_bins, bin_counts, by = "bin", all.x = TRUE)
  bin_counts_full$count[is.na(bin_counts_full$count)] <- 0
  k10_cn_expanded <- a[rep(1:nrow(a), times = bin_counts_full$count), ]
  k10_expanded[[file]]<-k10_cn_expanded
}


combined<-c(k500_expanded, k80_expanded, k40_expanded, k10_expanded, k200_expanded)
names(combined)<-gsub(".*_llr-|\\.tsv\\.gz$", "", names(combined))
combined <- combined[!grepl("-900000", names(combined))]
gold<-read_tsv("../Figures/ALL40_10k_reference.tsv")

#And then to split 
d<-readRDS("results/ALL/ALL40/clones/ALL40-final_clone_object-g10-10000.Rds")
clones_final <- tibble(dna_library_id=d$cells$dna_library_id, 
                       clone=d$cells[["clone_final"]], 
                       dna_reads=d$cells$bam_read_pairs, 
                       rna_phase=d$cells$cell_phase,
                       timepoint=d$cells$timepoint)
cells_used<-as.data.frame(colnames(combined$`500000-100000`))
colnames(cells_used)<-"dna_library_id"
info<-left_join(cells_used, clones_final)
table(info$clone)

#So we need to split the combined object into combined_tumor and combined_normal
normal<-info%>%filter(clone=="2_1")%>%select(dna_library_id)
tumor<-info%>%filter(!clone=="2_1")%>%select(dna_library_id)
subset_df <- function(df, ids) {
  df[, colnames(df) %in% ids, drop = FALSE]
}

# Create combined_normal and combined_tumor
combined_normal <- lapply(combined, subset_df, ids = normal$dna_library_id)
combined_tumor  <- lapply(combined, subset_df, ids = tumor$dna_library_id)


gold_use<-gold[,colnames(combined$`500000-100000`)]
gold_use_normal<-gold[,colnames(combined_normal$`500000-100000`)]
gold_use_tumor<-gold[,colnames(combined_tumor$`500000-100000`)]


#Do for tumor and normal separately (change code here)
for (name in names(combined_tumor)) {
  df <- combined_tumor[[name]]
  diff_counts[[name]] <- colSums(df != gold_use_tumor, na.rm = TRUE)
}
diff_counts_df <- as.data.frame(diff_counts)
rownames(diff_counts_df) <- colnames(gold_use_tumor)
print(diff_counts_df)


#Prepare for plotting 
long_df <- diff_counts_df %>%
  pivot_longer(cols = starts_with("X"), names_to = "condition", values_to = "count")

# Merge with the summary statistics
summary_stats <- diff_counts_df %>%
  pivot_longer(cols = starts_with("X"), names_to = "condition", values_to = "count")%>%
  group_by(condition) %>%
  summarise(
    median_count = median(count, na.rm = TRUE),
    sd_count = sd(count, na.rm = TRUE)
  )
merged_df <- left_join(long_df, summary_stats, by = "condition")
merged_df$resolution<-sub("\\..*", "", merged_df$condition)
merged_df$resolution<-gsub("X", "", merged_df$resolution)
merged_df$downsample <- sub(".*\\.", "", merged_df$condition)

table(merged_df$downsample)
table(merged_df$resolution)
library(ggplot2)
merged_df$downsample<-as.numeric(merged_df$downsample)/1000
merged_df$downsample <- factor(merged_df$downsample, levels = c(20, 30,40,50,60,70,80,90,100, 200, 300, 400))
merged_df$resolution<-as.numeric(merged_df$resolution)/1000
merged_df$resolution<-factor(merged_df$resolution, levels=c(10, 40,80, 200,500))
table(merged_df$resolution)
summary_df <- merged_df %>%
  group_by(downsample, resolution) %>%
  mutate(log_count=log(count+1))%>%
  summarise(
    median_count = median(count),
    sd_count = sd(count),
    sd_log=sd(log_count),
    .groups = 'drop'
  )

summary_df

ggplot(merged_df, aes(x = downsample, y = log(count+1), color = resolution)) + 
  geom_boxplot(aes(group = interaction(downsample, resolution)), outliers = FALSE)+
  geom_quasirandom(aes(group = interaction(downsample, resolution)), width = 0.1, dodge.width = 0.75, alpha=0.25, size=0.4) +
  theme_bw() +
  labs(x = "Downsampling", y="log(Basepairs different from gold standard)",  color = "resolution")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_color_brewer(palette="Dark2")+
  ylim(0,12)+
  ggtitle("llr segmentation vs 10k Reference")

