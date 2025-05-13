#Clone refinement 
#Clone and breakpoint refinement 

setwd("/wrk/resources/dntr2_extrascp/")
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
source("workflow/scripts/clone_functions_forPaper.R")

# Setup
patient_id <- "ALL40" #Change this  
threads = 12

# Defaults and load patient-specific parameters
#These are decided based on visual observation of different gammas from the pipeline output
#It is possible that different sc_gammas need to be run until desired outcome is achieved
default_params <- list(
  binsize = 10000,
  run = "/wrk/resources/dntr2_extrascp/results/DNA_37_SK_250428/", 
  sc_gamma = 5,
  sc_method = "umap",
  clone_gamma = 0.5,
  clone_min_bins = 10,
  clone_boundary_filter = 30
)

p<-default_params
print(unlist(p))
default_params$run

# Paths
clone_file <- file.path(p$run,patient_id,"clones", paste0(patient_id,"-clones-",p$sc_method,"-g",p$sc_gamma,"-",p$binsize,".txt"))
counts_file <- file.path(p$run,patient_id,paste0(patient_id,"-bincounts-",p$binsize,".tsv.gz"))
normal_counts_file <- paste0("resources/normals_scaling-normals_37bp_231206-", p$binsize, ".tsv.gz") 
sf<-read.table(paste0(p$run,"/", patient_id, "/clones/", patient_id, "-scalefactors-g", p$sc_gamma, "-", p$binsize, ".txt"), col.names = c("dna_library_id", "scale_factor"))
logodds<-read_tsv(paste0(p$run, "/", patient_id, "/clones/", patient_id, "-log_odds_df-scCN-", p$binsize, "-g", p$sc_gamma, ".tsv"))
scps <-read_tsv(paste0(p$run, "/", patient_id, "/clones/", patient_id, "-scps-scCN-", p$binsize, "-g", p$sc_gamma, ".txt"))

bins_all <- load_bins(bins_file = paste0("resources/fixed-",p$binsize,".bed"),
                      map_file = paste0("resources/fixed-",p$binsize,".map.txt"),
                      gc_file = paste0("resources/fixed-",p$binsize,".gc.txt"),
                      cytoband_file = "/wrk/resources/genomes/hg38-iGenome-ERCC/cytoBand.txt.gz")
bins_good <- read_tsv(paste0("resources/goodbins-",p$binsize,".bed"), skip=1,
                      col_names=c("chr","start","end","id"), col_select=1:4, show_col_types = F) # NOTE: Fixed at the 37bp good bin indices
counts <- as.matrix(data.table::fread(counts_file))
clones <- read_tsv(clone_file) %>% select(dna_library_id, clone=clone_final)

#Create single cell level scale factor 
logodds<-left_join(sf, logodds)
logodds$correct_scalefactor<-logodds$scale_factor*logodds$multiplication
logodds <- logodds %>%
  mutate(correct_scalefactor = coalesce(correct_scalefactor, scale_factor))

clones<-left_join(clones, logodds)

#Visual inspection of outliers 
boxplot(clones$correct_scalefactor~clones$clone)
#Find outliers
outliers <- clones %>%
  group_by(clone) %>%
  mutate(mean_sf = mean(correct_scalefactor, na.rm = TRUE)) %>%
  filter(abs(correct_scalefactor - mean_sf) > 0.5) %>%
  select(dna_library_id, clone, correct_scalefactor, mean_sf)

existing_clones <- unique(clones$clone)
numeric_part <- as.numeric(gsub("_1", "", existing_clones[grepl("^[0-9]+_1$", existing_clones)]))
max_clone <- max(numeric_part, na.rm = TRUE)

#If outliers are five or more we add them to their own clone, if they are fewer we move them to clone zero
if (length(outliers$dna_library_id) > 4) {
  new_clone <- paste0(max_clone + 1, "_1")
} else {
  new_clone <- "0"
}

clones$clone[clones$dna_library_id %in% outliers$dna_library_id] <- new_clone


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
  left_join(logodds)%>%
  select(dna_library_id, clone, everything()) %>% 
  filter(!is.na(dna_library_id)) # Added 250128: Hit-picking in some plates will mean different number of RNA vs DNA libraries in a plate 

# Load cnv data

d <- create_pseudobulk_analysis(counts_matrix = counts, 
                                bins_info = bins_all, 
                                good_bins = bins_good$id, 
                                cell_metadata = cell_data,
                                normal_counts = normal_counts_file, 
                                params = p)

# Processing
d <- normalize_counts(d, methods=c("ft_lowess", "gcmap"))
d <- call_segments(d, gamma=0.5, norm_segments = "ft_lowess_normal", norm_ratio="gcmap_normal", verbose=T)
d <- merge_small_segments(d, current="initial", revision="merged", min_bins_filter=10, boundary_filter=40, update_clones=T)
d <- calc_cn_integers(d) 
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

#Single cell copynumber heatmap
png(file=paste0(p$run,"/", patient_id, "/clones/", patient_id, "-final_clone_heatmap.png"), width=2400,height=1400,units="px",res=150)
plot_cell_heatmap(d, filtered=TRUE, smooth_bins = 10, annotate="timepoint", annotate_colors=timepoint_cols)
dev.off()

#Clone heatmap and per clone chromosome profiles 
#This is key to analyse if your clones are high quality and if the breakpoints detected are accurate 
pdf(file=file.path(p$run, patient_id, paste0(patient_id, "-clones_by_chr.pdf")), width=12, height=max(8, nrow(d$clones)*1.3))
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
                       rna_counts=d$cells$rna_counts,
                       rna_phase=d$cells$cell_phase,
                       timepoint=d$cells$timepoint)
# clones_final <- unnest(d$clones, cols="cells") %>% select(dna_library_id=cells, clone_id, revision)
write_tsv(clones_final, file.path(p$run, patient_id, paste0(patient_id, "-clones_final.txt")))

# Save Rds
write_rds(d, file.path(paste0(p$run,"/", patient_id, "/clones/", patient_id, "-final_clone_object.Rds")))
