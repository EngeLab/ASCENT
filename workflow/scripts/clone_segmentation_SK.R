library(tidyverse)
library(GenomicRanges)
library(copynumber)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ape)
options(scipen = 999)
try(RhpcBLASctl::blas_set_num_threads(1))
try(RhpcBLASctl::omp_set_num_threads(1))
source("workflow/scripts/clone_functions.R")
# library(mgsub)

# Parameters
run <- "results/ALL-latest/"
patient_id <- "ALL64"
read_length <- 37

sc_gamma <- 1
sc_method <- "umap"
binsize <- 20000
threads <- 4

diploid_clone <- "1_1"

if (read_length == 150) {
  normal_counts_file <- "resources/normals_scaling-normal_150_novogene-20000.tsv.gz"
} else {
  normal_counts_file <- "resources/normals_scaling-normals_37bp_231206-20000.tsv.gz"
}

# Paths
clone_file <- file.path(paste0(run,"clones/", patient_id,"-clones-",sc_method,"-g",sc_gamma,"-",binsize,".txt"))
counts_file <- file.path(run,paste0(patient_id,"-bincounts-",binsize,".tsv.gz"))
sfs <- read.table(paste0(run, "clones/", patient_id, "-scalefactors-g", sc_gamma, "-20000.txt"), col.names = c("dna_library_id", "scale_factor")) #Always g1 because thats where scp where calculated from 
scp <- read.table(paste0(run, "clones/", patient_id, "-scps-", binsize, ".txt"), header=TRUE)

# Load data
bins_all <- load_bins(bins_file = paste0("resources/fixed-",binsize,".bed"),
                      map_file = paste0("resources/fixed-",binsize,".map.txt"),
                      gc_file = paste0("resources/fixed-",binsize,".gc.txt"),
                      cytoband_file = "/wrk/resources/genomes/hg38-iGenome-ERCC/cytoBand.txt.gz") #Do we need a new cytoband file? na?? 
bins_good <- read_tsv(paste0("resources/goodbins-",binsize,".bed"), skip=1,
                      col_names=c("chr","start","end","id"), col_select=1:4, show_col_types = F) # NOTE: Fixed at the 37bp good bin indices
counts <- as.matrix(data.table::fread(counts_file))
clones <- read_tsv(clone_file) %>% select(dna_library_id, clone=clone_final)

# Add more cell-level metadata, optional
metadata_files <- list("meta"=file.path(run, "metadata_long.tsv"), 
                       "rna"=file.path(run, "rna_phases.txt"))
cell_data <- load_and_join(metadata_files, join_by="cell_id") %>% 
  filter(dna_library_id %in% colnames(counts)) %>% 
  left_join(clones, by="dna_library_id") %>% 
  select(dna_library_id, clone, everything())

### SCP and scale factor processing
# #Here - read in scp data and decide on acceptible scale factors per clone 
# a<-left_join(scp, clones)
# df <- left_join(a, scalefactors) %>%filter(ploidy < 5) %>% mutate(copy = factor(ploidy))
# 
# summary_stats <- summarize_stats(df, diploid_clone)
# 
# missing_clones_df <- tibble(
#   clone = unique(df$clone[!df$clone%in%summary_stats$clone]),
#   mean_coverage_at_2x = NA,
#   mean_scale_factor = NA,
#   sd_cov = NA,
#   diploid_coverage = NA,
#   sd_cov_dip = NA,
#   scalefactor = "good",
#   min_scale_factor = 1.5,
#   max_scale_factor = 4.5
# )
# 
# summary_stats <- bind_rows(summary_stats, missing_clones_df)
# summary_stats
# ggplot(df, aes(x = copy, y = 100000*(depth2), fill = as.factor(clone))) +
#   geom_boxplot() +
#   labs(x = "Copy", y = "Fraction Covered * 10e6", fill = "Clone") +
#   ggtitle("Boxplots of Depth 2 by Clone") +
#   ylim(0,2)+
#   theme_minimal()
# 
# 
# # Add more cell-level metadata, optional
# metadata_files <- list("meta"=file.path(run, "metadata_long.tsv"), 
#                        "rna"=file.path(run, "rna_phases.txt"))
# 
# #If I update the cell_data object??? 
# 
# scp_2<-scp%>%filter(ploidy==2)%>%select(depth2, dna_library_id)

# Load cnv data
d <- create_pseudobulk_analysis(counts_matrix = counts, 
                                bins_info = bins_all, 
                                good_bins = bins_good$id, 
                                cell_metadata = cell_data,
                                normal_counts = normal_counts_file, 
                                params = list(patient_id=patient_id, gamma=sc_gamma, min_cells=8, dim_reduction=sc_method))

# Processing
d <- normalize_counts(d, methods=c("ft_lowess", "gcmap"))
d <- call_segments(d, gamma=0.25, norm_segments = "ft_lowess_normal", norm_ratio="gcmap_normal", verbose=T)
d <- merge_small_segments(d, current="initial", revision="merged", min_bins_filter=10, boundary_filter=50, update_clones=T)
# d <- calc_cn_integers_SK(d, diploid_clone=diploid_clone)
d <- calc_cn_integers(d,scale_range=c(1.5, 4.5))
d <- mask_high_residuals(d, max_residual = 0.3, clone_filter_fraction=0.3, update_clones=T, diploid_clone=diploid_clone) #Need to be able to split the clones that are flagged here... 
d <- split_mixed_clones(d, residual_threshold = 0.3, improvement_threshold = 0.8, update_clones = T, verbose=F, plot=T, diploid_clone=diploid_clone) 
d <- remove_bad_clones(d, diploid_clone=diploid_clone)
d <- refine_segments_from_cn(d, diploid_clone=diploid_clone)
d<- merge_duplicate_clones(d, diploid_clone=diploid_clone)
d <- merge_almost_duplicate_clones(d, diploid_clone=diploid_clone)
d<-remove_small_clones(d, min_size_clone = 4)
plot_clone_heatmap(d, show_chr = T)

#  14   q  21860000  22460000       23      36  97266  97288    FALSE          NA            NA                   NA
#So - which genes are in there? 
#32  14   q  21860000  22460000       23      36  97266  97288    FALSE          NA            NA                   NA
#need to figure out which bins unfiltered these are 110669 110698

# 
# plot(109576:114928, gc[109576:114928], col = "blue", xlab = "Index", ylab = "Values")
# 
# # Overlay points (use the same x-values)
# points(110669:110698, gc[110669:110698], col = "red", pch = 19)

#gamma 1 would be fine - but lets try 0.75 

pdf(paste0(run, "clones/", patient_id, "-chromosome.pdf"))
plot_clone_detail(d, region="chr1")
plot_clone_detail(d, region="chr2") 
plot_clone_detail(d, region="chr3") 
plot_clone_detail(d, region="chr4") 
plot_clone_detail(d, region="chr5") 
plot_clone_detail(d, region="chr6") #Correct calls - maybe a breakpoint missing 
plot_clone_detail(d, region="chr7") 
plot_clone_detail(d, region="chr8") 
plot_clone_detail(d, region="chr9") 
plot_clone_detail(d, region="chr10")
plot_clone_detail(d, region="chr11")
plot_clone_detail(d, region="chr12") 
plot_clone_detail(d, region="chr13") 
plot_clone_detail(d, region="chr14") 
plot_clone_detail(d, region="chr15") 
plot_clone_detail(d, region="chr16") #missing breakpoint 
plot_clone_detail(d, region="chr17") 
plot_clone_detail(d, region="chr18") 
plot_clone_detail(d, region="chr19") #missing breakpoints for sure 
plot_clone_detail(d, region="chr20") 
plot_clone_detail(d, region="chr21") 
plot_clone_detail(d, region="chr22") 
dev.off()



#But then here - i am re-calculating the ones in zero - but maybe thats just good - and totally on purpose..

d <- calc_cell_cn(d,diploid_clone=diploid_clone) 

#So lets make sure that the cn data is there now!! 

plot_cell_heatmap(d) #6_1, 6_3 - 7_1ma and 9_1 #91 is different 

plot_clone_heatmap(d, show_chr=T)
d2 <- recall_cells(d, mad_cutoff = , plot=F)

#So we need to remove that recall_prob dependency here. . . ? if I don't want that to be what is being plotted in the plot_clone_detail stuff... 
plot_cell_heatmap(d2, smooth_bins = 50)
plot_clone_heatmap(d2, show_chr = T, clone_types=T)

#So then I would save this d2 object . . 
saveRDS(d2, file=paste0(run, "clones/", patient_id, "final_clone_object.Rds"))
png(paste0(run, "clones/", patient_id, "-final-clone-heatmap.png"),  width=2400,height=1400,units="px",res=150)
plot_clone_heatmap(d2, show_chr=T)
dev.off()
png(paste0(run, "clones/", patient_id, "-final-cell-heatmap.png"),  width=2400,height=1400,units="px",res=150)
plot_cell_heatmap(d2, smooth_bins = 50)
dev.off()

#d<-readRDS(paste0(run, "clones/", patient_id, "final_clone_object.Rds"))


#I mean it's only that specific clone - lets look at a different tumor... 
#plot(d[["cells"]][["cn"]][[7]]) one of the cells that get the wrong copy numbers.... why...... 

#Should I make like a separate function for this??? 
#Or what should we do regarding this fuckin recall -like 

#So then this is good but then need the traditional one as well - - next stage figure that out 
# 
# 
# png(paste0("analysis/",patient_id,"-cell_heatmap.png"), width=2000, height=1100)
# plot_cell_heatmap(d, smooth_bins = 5, annotate="timepoint", annotate_colors=timepoint_cols)
# dev.off()
#     

#Beaut 


# 241202: Make matrix and annot df to Martin
# Run cn profiles of replicating cells 
# repl_idx <- grepl("S|G2M",d$cells$cell_phase)
# d <- calc_cell_cn(d, cell_idx=repl_idx, scale_range=c(1.5, 4.5))
# all_idx <- !is.na(d$cells$clone_initial) | repl_idx
# cell_bins.l <- lapply(d$cells$cn[all_idx], function(x) expand_gaps_vec(x, length.out=nrow(d$bins$all), idx=d$bins$good$id))
# cell_mtx <- do.call(cbind, cell_bins.l)
# colnames(cell_mtx) <- d$cells$dna_library_id[all_idx]
# annot <- tibble(cell=d$cells$dna_library_id[all_idx], 
#                 rna_phase=d$cells$cell_phase[all_idx],
#                 clone_match=d$cells$clone_recall[all_idx],
#                 scale_factor=unlist(d$cells$scale_factor[all_idx])
# )
# clone_bins.l <- lapply(d$clones$cn, function(x) expand_gaps_vec(x, length.out=nrow(d$bins$all), idx=d$bins$good$id))
# clone_mtx <- do.call(cbind, clone_bins.l)
# colnames(clone_mtx) <- d$clone$clone_id
# save(list=c("cell_mtx","clone_mtx", "annot"), file="analysis/ALL40_wSMphase.Rda")


# TODO list
# - Filter clones with bad residuals all over -> to clone 0
# - plot_clones()/plot_clone_detail() # Before splitting/merger
# - Add LOH/SCP logic to scale correctly
