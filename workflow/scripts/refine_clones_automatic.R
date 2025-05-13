#Script for automatic clone refinement 

library(tidyverse)
library(GenomicRanges)
library(copynumber)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ape)
library(ggbeeswarm)
library(patchwork)
library(mgsub)
options(scipen = 999)
try(RhpcBLASctl::blas_set_num_threads(1))
try(RhpcBLASctl::omp_set_num_threads(1))
source("workflow/scripts/clone_functions_forPaper.R")
threads=6
# Setup

clone_min_bins = 10
clone_boundary_filter = 30


# Paths
clone_file <- file.path(snakemake@input[["clones"]])
counts_file <- file.path(snakemake@input[["counts"]])
normal_counts_file <- snakemake@input[["normal_cells"]] 
sf<-read.table(snakemake@input[["sf"]], col.names = c("dna_library_id", "scale_factor"))
logodds<-read_tsv(snakemake@input[["logodds"]])

bins_all <- load_bins(bins_file = snakemake@input[["bins"]],
                      map_file = snakemake@input[["map"]],
                      gc_file = snakemake@input[["gc"]],
                      cytoband_file = snakemake@input[["cytoband"]])
bins_good <- read_tsv(snakemake@input[["good_bins"]], skip=1,
                      col_names=c("chr","start","end","id"), col_select=1:4, show_col_types = F) # NOTE: Fixed at the 37bp good bin indices
counts <- as.matrix(data.table::fread(counts_file))
clones <- read_tsv(clone_file) %>% select(dna_library_id, clone=clone_final)

#Create single cell level scale factor 
logodds<-left_join(sf, logodds)
logodds$correct_scalefactor<-logodds$scale_factor*logodds$multiplication
logodds <- logodds %>%
  mutate(correct_scalefactor = coalesce(correct_scalefactor, scale_factor))

clones<-left_join(clones, logodds)

meta_seed <- read_tsv(snakemake@input[["meta"]])
meta_qc_dna <- read_tsv(snakemake@input[["qc_dna"]]) %>% 
  mutate(cell_id=dna_to_cellid(dna_library_id)) %>% select(-dna_library_id, -starts_with("dna_frac"), -starts_with("fq"))

cell_data <- meta_seed %>% 
  #left_join(meta_phase, by="cell_id") %>% 
  #left_join(meta_qc_rna, by="cell_id") %>% 
  left_join(meta_qc_dna, by="cell_id") %>% 
  left_join(clones, by="dna_library_id") %>% 
  left_join(logodds)%>%
  select(dna_library_id, clone, everything()) %>% 
  filter(!is.na(dna_library_id)) # Added 250128: Hit-picking in some plates will mean different number of RNA vs DNA libraries in a plate 

# Load cnv data

d <- create_pseudobulk_analysis(counts_matrix = counts, 
                                bins_info = bins_all, 
                                good_bins = bins_good$id, 
                                cell_metadata = cell_data,
                                normal_counts = normal_counts_file)

# Processing
d <- normalize_counts(d, methods=c("ft_lowess", "gcmap"))
d <- call_segments(d, gamma=0.5, norm_segments = "ft_lowess_normal", norm_ratio="gcmap_normal", verbose=T)
d <- merge_small_segments(d, current="initial", revision="merged", min_bins_filter=10, boundary_filter=40, update_clones=T)
d <- calc_cn_integers(d) 
plot_clone_heatmap(d)
d <- split_mixed_clones(d, residual_threshold = 0.3, improvement_threshold = 0.8, update_clones = T, verbose=F, plot=T)
d <- mask_high_residuals(d, max_residual = 0.3, clone_filter_fraction=0.3, update_clones=T)
d <- remove_bad_clones(d)
d <- refine_segments_from_cn(d)
d <- merge_duplicate_clones(d)
d <- remove_small_clones(d, min_size_clone = 5)
plot_clone_heatmap(d)
### Re-calculate single cell copy numbers based on the refined segments 
d <- calc_cell_cn(d)

#Remove cells that don't fit well 
d <- remove_bad_cells(d, max_diff_bins = 1000, update_clones=T)

plot_clone_detail(d, region="chr12")

#Single cell copynumber heatmap
png(file=snakemake@output[["sc_heatmap"]], width=2400,height=1400,units="px",res=150)
plot_cell_heatmap(d, filtered=TRUE, smooth_bins = 10)
dev.off()

#Clone heatmap and per clone chromosome profiles 
#This is key to analyse if your clones are high quality and if the breakpoints detected are accurate 
pdf(file=snakemake@output[["chr_heatmap"]], width=12, height=max(8, nrow(d$clones)*1.3))
plot_clone_heatmap(d, show_chr=T, clone_types=T, only_divergent = F)
for(chr in levels(d$bins$all$chr)){
  plot_clone_detail(d, region = chr)
}
dev.off()



# Final clones
clone_slot <- names(d$cells)[max(grep("clone_", names(d$cells)))]
clones_final <- tibble(dna_library_id=d$cells$dna_library_id,
                       clone=d$cells[[clone_slot]],
                       dna_reads=d$cells$bam_read_pairs,
                       #rna_counts=d$cells$rna_counts,
                       #rna_phase=d$cells$cell_phase,
                       timepoint=d$cells$timepoint)
# clones_final <- unnest(d$clones, cols="cells") %>% select(dna_library_id=cells, clone_id, revision)
write_tsv(clones_final, snakemake@output[["final_clones"]])

# Save Rds
write_rds(d, snakemake@output[["final_clone_object"]])
