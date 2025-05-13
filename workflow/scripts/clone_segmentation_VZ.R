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
options(scipen = 999)
try(RhpcBLASctl::blas_set_num_threads(1))
try(RhpcBLASctl::omp_set_num_threads(1))
source("workflow/scripts/clone_functions_VZ.R")

# Setup
patient_id <- "ALL40"
threads = 12

# Defaults and load patient-specific parameters
default_params <- list(
  binsize = 40000,
  run = "/wrk/resources/dntr2_chr_arm/results/Downsampling_ALL40/",
  sc_gamma = 0.5,
  sc_method = "umap",
  clone_gamma = 0.1,
  clone_min_bins = 10,
  clone_boundary_filter = 30
)

#p <- load_parameters(patient_id, yaml_file="config/ALL_patient_params.yaml", defaults = default_params)
p<-default_params
print(unlist(p))
default_params$run

# Paths
clone_file <- file.path(p$run,patient_id,"clones", paste0(patient_id,"-clones-",p$sc_method,"-g",p$sc_gamma,"-",p$binsize,".txt"))
counts_file <- file.path(p$run,patient_id,paste0(patient_id,"-bincounts-",p$binsize,".tsv.gz"))
normal_counts_file <- "resources/normals_scaling-normals_37bp_231206-20000.tsv.gz"

# Load data
bins_all <- load_bins(bins_file = paste0("resources/fixed-",p$binsize,".bed"),
                      map_file = paste0("resources/fixed-",p$binsize,".map.txt"),
                      gc_file = paste0("resources/fixed-",p$binsize,".gc.txt"),
                      cytoband_file = "/wrk/resources/genomes/hg38-iGenome-ERCC/cytoBand.txt.gz")
bins_good <- read_tsv(paste0("resources/goodbins-",p$binsize,".bed"), skip=1,
                      col_names=c("chr","start","end","id"), col_select=1:4, show_col_types = F) # NOTE: Fixed at the 37bp good bin indices
counts <- as.matrix(data.table::fread(counts_file))
clones <- read_tsv(clone_file) %>% select(dna_library_id, clone=clone_final)

# Add more cell-level metadata, optional
meta_seed <- read_tsv(file.path(p$run, patient_id, paste0(patient_id,"-metadata_long.tsv")))
if(file.exists(file.path(p$run, patient_id, paste0(patient_id,"-rna_phases.txt")))){
  meta_phase <- read_tsv(file.path(p$run, patient_id, paste0(patient_id,"-rna_phases.txt")))
} else {
  meta_phase <- data.frame(cell_id=character(0), cell_phase=character(0))
}
if(file.exists(file.path(p$run, patient_id, "qc", paste0(patient_id, "-qc_rna.tsv")))) {
  meta_qc_rna <- read_tsv(file.path(p$run, patient_id, "qc", paste0(patient_id, "-qc_rna.tsv"))) %>% 
    mutate(cell_id=rna_to_cellid(rna_library_id)) %>% select(-rna_library_id, -starts_with("fq"))
} else {
  meta_qc_rna <- data.frame(cell_id=character(0))
}
meta_qc_dna <- read_tsv(file.path(p$run, patient_id, "qc", paste0(patient_id, "-qc_dna.tsv"))) %>% 
  mutate(cell_id=dna_to_cellid(dna_library_id)) %>% select(-dna_library_id, -starts_with("dna_frac"), -starts_with("fq"))

cell_data <- meta_seed %>% 
  left_join(meta_phase, by="cell_id") %>% 
  left_join(meta_qc_rna, by="cell_id") %>% 
  left_join(meta_qc_dna, by="cell_id") %>% 
  left_join(clones, by="dna_library_id") %>% 
  select(dna_library_id, clone, everything()) %>% 
  filter(!is.na(dna_library_id)) # Added 250128: Hit-picking in some plates will mean different number of RNA vs DNA libraries in a plate 

# Load cnv data
#Fix some shit here - i have done that before so should be able to find that... 

d <- create_pseudobulk_analysis(counts_matrix = counts, 
                                bins_info = bins_all, 
                                good_bins = bins_good$id, 
                                cell_metadata = cell_data,
                                normal_counts = normal_counts_file, 
                                params = p)

# Processing
d <- normalize_counts(d, methods=c("ft_lowess", "gcmap"))
d <- call_segments(d, gamma=p$clone_gamma, norm_segments = "ft_lowess_normal", norm_ratio="gcmap_normal", verbose=T)
d <- merge_small_segments(d, current="initial", revision="merged", min_bins_filter=p$clone_min_bins, boundary_filter=p$clone_boundary_filter, update_clones=T)
d <- calc_cn_integers(d, scale_range=c(1.5, 4.5))
d <- split_mixed_clones(d, residual_threshold = 0.3, improvement_threshold = 0.8, update_clones = T, verbose=F, plot=T)
d <- mask_high_residuals(d, max_residual = 0.3, clone_filter_fraction=0.3, update_clones=T)
d <- remove_bad_clones(d)
d <- refine_segments_from_cn(d)
d <- merge_duplicate_clones(d)

### Initial plots
d <- calc_cell_cn(d, scale_range=c(1.5, 4.5))
png(file=file.path(p$run, patient_id, paste0(patient_id, "-clones_cells_initial.png")), width=2400,height=1400,units="px",res=150)
plot_cell_heatmap(d, filtered=TRUE, smooth_bins = 10, annotate="timepoint", annotate_colors=timepoint_cols)
dev.off()

pdf(file=file.path(p$run, patient_id, paste0(patient_id, "-clones_by_chr.pdf")), width=12, height=max(8, nrow(d$clones)*1.3))
plot_clone_heatmap(d, show_chr=T, clone_types=T, only_divergent = F)
for(chr in levels(d$bins$all$chr)){
  plot_clone_detail(d, region = chr)
}
dev.off()

### Optional: Targeted split based on review of initial heatmap
# ALL42 split
# d <- split_mixed_clones(d, clones="1_1", residual_threshold = 0.3, improvement_threshold = 0.8, total_improvement_threshold = 3, update_clones = T, verbose=T, plot=T)
# ALL4 split
# d <- split_mixed_clones(d, clones=c("1_1","5_1"), residual_threshold = 0.2, improvement_threshold = 0.8, total_improvement_threshold = 1, update_clones = T, verbose=T, plot=T)

### Optional mask specific segment
# d <- mask_segment(d, mask=NULL)
# ALL3 masking:
# plot_clone_detail(d, region = "chr14:21000000-24000000")
d <- mask_segment(d, mask=67) # ALL3

### Optional: fuzzy merge
d <- fuzzy_merge_clones(d, max_diff_frac = 0.1, max_residual = 0.4)
# d <- fuzzy_merge_clones(d, max_diff_frac = 0.1, max_residual = 0.2) # ALL35

# Rerun masking and cell-level cn?
d <- mask_high_residuals(d, max_residual = 0.3, clone_filter_fraction=0.3, update_clones=T)
d <- calc_cell_cn(d, scale_range=c(1.5, 4.5))

# Recall from 0 and/or previously excluded
pdf(file=file.path(p$run, patient_id, paste0(patient_id, "-clones_recall.pdf")), width=10)
# d <- recall_cells(d, from="0", mad_cutoff=7, plot=T, update_clones=T)
# d <- recall_cells(d, from=c(d$clones$clone_id,"0",NA), mad_cutoff=7, plot=T, update_clones=T)
d <- recall_cells(d, from="0", mad_cutoff=7, plot=T, update_clones=T)
dev.off()

# Remove obvious outliers
# d <- remove_bad_cells(d, max_diff_bins = 1000, clone_slot = "clone_recall", update_clones=T)
d <- remove_bad_cells(d, max_diff_bins = 1000, update_clones=T)

# Drop small/orphaned clones. TODO: Check that this works properly?
# d <- remove_small_clones(d, min_size_clone = 5)

# Re-plot detailed clones (after recall + removal of bad cells)
pdf(file=file.path(p$run, patient_id, paste0(patient_id, "-clones_by_chr-final.pdf")), width=12, height=max(8, nrow(d$clones)*1.3))
plot_clone_heatmap(d, show_chr=T, clone_types=T, only_divergent = F)
for(chr in levels(d$bins$all$chr)){
  plot_clone_detail(d, region = chr)
}
dev.off()

# Cell-level heatmaps
png(file=file.path(p$run, patient_id, paste0(patient_id, "-clones_cells_final.png")), width=2400,height=1400,units="px",res=150)
plot_cell_heatmap(d, filtered=TRUE, smooth_bins = 10, annotate="timepoint", annotate_colors=timepoint_cols)
dev.off()

png(file=file.path(p$run, patient_id, paste0(patient_id, "-clones_cells_all.png")), width=2400,height=1400,units="px",res=150)
plot_cell_heatmap(d, filtered=FALSE, smooth_bins = 10, annotate="timepoint", annotate_colors=timepoint_cols)
dev.off()

png(file=file.path(p$run, patient_id, paste0(patient_id, "-clones_cells_all_qc.png")), width=2400,height=1400,units="px",res=150)
plot_cell_heatmap(d, filtered=FALSE, smooth_bins = 10, annotate="timepoint", annotate_colors=timepoint_cols, group_by="pass_dna_qc")
dev.off()

# Segments + cn profiles
segs_slot <- length(d$segments)
cn <- as.matrix(sapply(d$clones$cn, function(x) x[d$segments[[segs_slot]]$start]))
colnames(cn) <- d$clones$clone_id
cn.segs <- cbind(d$segments[[segs_slot]][,1:4], cn)
write_tsv(cn.segs, file.path(p$run, patient_id, paste0(patient_id, "-cn_segments.txt")))

# Final clones
clone_slot <- names(d$cells)[max(grep("clone_", names(d$cells)))]
clones_final <- tibble(dna_library_id=d$cells$dna_library_id, 
                       clone=d$cells[[clone_slot]], 
                       dna_reads=d$cells$bam_read_pairs, 
                       rna_counts=d$cells$rna_counts, 
                       rna_phase=d$cells$cell_phase,
                       timepoint=d$cells$timepoint)
# clones_final <- unnest(d$clones, cols="cells") %>% select(dna_library_id=cells, clone_id, revision)
write_tsv(clones_final, file.path(p$run, patient_id, paste0(patient_id, "-clones_final.txt")))

# Save Rds
write_rds(d, file.path(p$run, patient_id, paste0(patient_id, "-clones.Rds")))