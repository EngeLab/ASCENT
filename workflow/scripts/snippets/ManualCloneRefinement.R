#For complicated samples or with new methods automatic clone refinement within the pipelines sometimes gives unsatisfying results
#In those cases it is good to run clone_refinement interactively and finetune some parameters to the case in question 
#This is a script to make that easier 

setwd("/wrk/resources/ASCENT/")
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
library(mgsub)

# Setup
patient_id <- "ALL40" #Change this  
threads = 12

# Defaults and load patient-specific parameters
#These are decided based on visual observation of different gammas from the pipeline output
#It is possible that different sc_gammas need to be run until desired outcome is achieved
default_params <- list(
  binsize = 10000,
  binsize_refine=10000,
  run = "results/ALL/", 
  sc_gamma = 10,
  sc_method = "umap",
  clone_min_bins = 10,
  clone_boundary_filter = 30
)

p<-default_params
print(unlist(p))
default_params$run

# Paths
clone_file <- file.path(p$run,patient_id,"clones", paste0(patient_id,"-clones-",p$sc_method,"-g",p$sc_gamma,"-",p$binsize,".txt"))
counts_file <- file.path(p$run,patient_id,paste0(patient_id,"-bincounts-",p$binsize_refine,".tsv.gz"))
normal_counts_file <- paste0("resources/normals_scaling-normals_37bp_231206-", p$binsize_refine, ".tsv.gz") 
sf<-read.table(paste0(p$run,"/", patient_id, "/clones/", patient_id, "-scalefactors-g", p$sc_gamma, "-", p$binsize, ".txt"), col.names = c("dna_library_id", "scale_factor"))
logodds<-read_tsv(paste0(p$run, "/", patient_id, "/clones/", patient_id, "-log_odds_df-scCN-", p$binsize, "-g", p$sc_gamma, ".tsv"))
scps <-read_tsv(paste0(p$run, "/", patient_id, "/clones/", patient_id, "-scps-scCN-", p$binsize, "-g", p$sc_gamma, ".txt"))

bins_all <- load_bins(bins_file = paste0("resources/fixed-",p$binsize_refine,".bed"),
                      map_file = paste0("resources/fixed-",p$binsize_refine,".map.txt"),
                      gc_file = paste0("resources/fixed-",p$binsize_refine,".gc.txt"),
                      cytoband_file = "/wrk/resources/genomes/hg38-iGenome-ERCC/cytoBand.txt.gz")

bins_good <- read_tsv(paste0("resources/goodbins-",p$binsize_refine,".bed"), skip=1,
                      col_names=c("chr","start","end","id"), col_select=1:4, show_col_types = F) 

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

# Normal panel scaling
norm<-"gcmap_norm"

use_normal_scaling <- !is.null(normal_counts_file) && file.exists(normal_counts_file)
if (!use_normal_scaling) {
  cat("Normal cell panel not provided or not found. Running without normal scaling.\n")
  normal_counts_file <- NULL
  norm<-"gcmap"
}



# Load cnv data
d <- create_pseudobulk_analysis(counts_matrix = counts, 
                                bins_info = bins_all, 
                                good_bins = bins_good$id, 
                                cell_metadata = cell_data,
                                normal_counts = normal_counts_file)

#The first tunable parameter, and the most important, is the gamma parameter of mpcf, 
#if too many breakpoints were found, increase it, if breakpoints are missing - decrease it 

# Processing
d <- normalize_counts(d, methods=c("ft_lowess", "gcmap"))
if(is.null(normal_counts_file)){
  d <- call_segments(d, gamma=2.5, norm_segments = "ft_lowess", norm_ratio="gcmap", verbose=T)
} else {
  d <- call_segments(d, gamma=0.5, norm_segments = "ft_lowess_normal", norm_ratio="gcmap_normal", verbose=T)
}
d <- call_segments(d, gamma=0.5, norm_segments = "ft_lowess_normal", norm_ratio="gcmap_normal", verbose=T)
d <- merge_small_segments(d, current="initial", revision="merged", min_bins_filter=10, boundary_filter=40, update_clones=T)
d <- calc_cn_integers(d) 
#Visualize the segments found: 
plot_clone_heatmap(d, show_chr = T) 
#It is good to look at individual chromosomes to understand if the gamma needs to be increased or decreased 
plot_clone_detail(d, region="chr1", norm=norm)
d <- split_mixed_clones(d, residual_threshold = 0.3, improvement_threshold = 0.8, update_clones = T, verbose=F, plot=T)
plot_clone_heatmap(d, show_chr = T)
d <- remove_bad_clones(d)
d <- refine_segments_from_cn(d) 
plot_clone_heatmap(d, show_chr = T)
d <- mask_high_residuals(d, max_residual = 0.3, clone_filter_fraction=0.3, update_clones=T) #Mask high residuals again? 
plot_clone_heatmap(d, show_chr = T)

#Sometimes it is desired to ignore specific segments (for example TCR regions on chr14 and chr7 sometimes also show up in diploid immune cells, and often we would like to mask those regions so they don't affect clones)
#Look for segments on chr that are problematic like this: 

d<-mask_segment(d, mask=42)
d<- merge_duplicate_clones(d)
d <- remove_small_clones(d, min_size_clone = 2)
#Fuzzy merge clones merges clones that are identical except for in regions where residuals are high - 
#This could indicate regions where copy numbers are not called accurately and we often would like to merge these clones 
d<-fuzzy_merge_clones(d)


### Re-calculate single cell copy numbers based on the refined segments 
d <- calc_cell_cn(d)

#Then here for the recall - 

#We have now additionally included a recall function in which cells from clone zero are tested against all other clones to see if we can assign them to any clone
#This is optional and not run automatically in the pipeline 

#d <- recall_cells(d, mad=5)


#Remove cells that don't fit well - the max_diff_bins parameter here can be changed depending on sample quality + complexity 
d <- remove_bad_cells(d, max_diff_bins = 1000, update_clones=T)
d <- remove_small_clones(d, min_size_clone = 5)



png(file=paste0(p$run,"/", patient_id, "/clones/", patient_id, "-final_clone_heatmap_manual.png"), width=2400,height=1400,units="px",res=150)
plot_cell_heatmap(d, filtered=TRUE, smooth_bins = 10, annotate="timepoint", annotate_colors=timepoint_cols)
dev.off()

#Clone heatmap and per clone chromosome profiles 
#This is key to analyse if your clones are high quality and if the breakpoints detected are accurate 
pdf(file=file.path(p$run, patient_id, "/clones/", paste0(patient_id, "-clones_by_chr_manual.pdf")), width=12, height=max(8, nrow(d$clones)*1.3))
plot_clone_heatmap(d, show_chr=T, clone_types=T, only_divergent = F)
for(chr in levels(d$bins$all$chr)){
  plot_clone_detail(d, region = chr, norm=norm)
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
write_tsv(clones_final, file.path(p$run, patient_id, "/clones", paste0(patient_id, "-clones_final_manual.txt")))

# Save Rds
write_rds(d, file.path(paste0(p$run,"/", patient_id, "/clones/", patient_id, "-final_clone_object_manual.Rds")))

#This rds object can then be used as the input for haplotyping + medicc2 if desired ! 

