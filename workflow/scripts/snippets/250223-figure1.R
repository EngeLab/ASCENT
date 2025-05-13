# 250223: Big plot ALL3
library(tidyverse)
library(Seurat)
source("workflow/scripts/general.R")
make_cn_colorscale <- function(max_cn) {
  if (max_cn <= 5) {
    return(circlize::colorRamp2(c(0, 2, 4, 5), c("blue","white","red", "brown")))
  } else if(max_cn <= 20) {
    return(circlize::colorRamp2(c(0, 2, 4, 5, 20), c("blue","white","red", "brown", "green")))
  } else {
    return(circlize::colorRamp2(c(0, 2, 4, 5, 20, 50), c("blue","white","red", "brown", "green", "yellow")))
  }
}
# source("workflow/scripts/clone_functions_VZ.R")
run <- "results/ALL-latest"
patient_id <- "ALL3"
normal <- c("1_1","9_1")
exclude <- "0"

# d <- read_rds(file.path(run, patient_id, paste0(patient_id, "-clones.Rds")))

clones <- read_tsv(file.path(run, patient_id, paste0(patient_id, "-clones_final.txt"))) %>% 
  filter(!is.na(clone) & !clone %in% c(normal,exclude))

# Load non-diploid clones and aberrant segments only
segs <- read_tsv(file.path(run, patient_id, paste0(patient_id, "-cn_segments.txt"))) %>% 
  select(-all_of(normal)) %>% 
  filter(rowSums(across(5:ncol(.), ~. != 2)) > 0 & ! chr %in% c("chrY"))

idx.subclonal <- apply(segs[,5:ncol(segs)], 1, function(x) length(unique(x)) > 1)
idx.clonal <- apply(segs[,5:ncol(segs)], 1, function(x) length(unique(x))==1)

str(clones)
str(segs)

### Test 1: Complex Heatmap
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Create timepoint-clone combinations with relative frequencies
clone_freqs <- clones %>%
  group_by(clone, timepoint) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 0) %>%
  # Calculate relative frequency within each timepoint
  group_by(timepoint) %>%
  mutate(rel_freq = count / sum(count)) %>%
  # Sort within each timepoint by frequency
  arrange(timepoint, rel_freq) %>%
  ungroup()

# Reshape segments data to long format with better segment names
segs_long <- segs %>%
  pivot_longer(cols = c("2_1", "3_1", "4_1", "5_1", "6_1", "7_1"),
               names_to = "clone",
               values_to = "copy_number") %>%
  mutate(segment = sprintf("%s%s:%.1f-%.1fMb", 
                           gsub("chr", "", chr),
                           arm,
                           start.pos/1e6,
                           end.pos/1e6))

# Create the matrix maintaining timepoint grouping
matrix_data <- clone_freqs %>%
  # Join with segments data to get copy numbers
  left_join(segs_long, by = "clone") %>%
  # Create the matrix
  select(segment, timepoint, clone, copy_number) %>%
  pivot_wider(id_cols = segment,
              names_from = c(timepoint, clone),
              values_from = copy_number) %>%
  column_to_rownames("segment") %>%
  as.matrix()

# Get column info for proper ordering and splits
col_info <- tibble(
  full_name = colnames(matrix_data),
  timepoint = sub("_.*$", "", full_name),
  clone = sub("^[^_]*_", "", full_name)
) %>%
  left_join(clone_freqs, by = c("timepoint", "clone")) %>%
  arrange(timepoint, rel_freq)

# Reorder matrix columns
matrix_data <- matrix_data[, col_info$full_name]
colnames(matrix_data) <- col_info$clone

# Clone frequencies 
freq_anno <- col_info$rel_freq
names(freq_anno) <- colnames(matrix_data)

# Colors
col_fun <- make_cn_colorscale(max(matrix_data, na.rm=TRUE))
freq_col_fun <- colorRamp2(c(0, 1), c("white", "black"))

# Create column splits for timepoints
column_splits <- col_info$timepoint %>%
  factor(levels = unique(.))

# Annotation
ha <- HeatmapAnnotation(
  "Relative Frequency" = freq_anno,
  col = list("Relative Frequency" = freq_col_fun),
  show_legend = FALSE
)

# Heatmap
ht <- Heatmap(matrix_data,
              name = "Copy Number",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              column_names_rot = 45,
              column_split = column_splits,
              column_title_gp = gpar(fontsize = 10),
              column_title_side = "top",
              border = TRUE,
              column_gap = unit(5, "mm"),
              row_title = "Segments",
              top_annotation = ha,
              show_heatmap_legend = FALSE)

draw(ht)

### Layer 2: RNA data
dat <- read_rds(file.path(run, patient_id, paste0(patient_id, "-rna_seurat.Rds")))
dat$clone <- clones$clone[match(dat$dna_library_id, clones$dna_library_id)]
clone_cols <- get_clones_col(dat$clone)
timepoints_short <- c("d0", "d2", "d3", "d5", "d15", "d29", "rel", "rel2","other")
timepoint_cols_short <- c("d0"="#cfcfcf","d2"="#fee08b","d3"="#d9ef8b","d5"="#91cf60","d15"="#1a9850","d29"="#136b39","rel"="#d73027", "rel2"="#002984", "other"="#002984")
dat@meta.data <- dat@meta.data %>% 
  mutate(
    timepoint_short=factor(case_when(timepoint=="diagnosis" ~ "d0",
                                     timepoint=="day 2" ~ "d2",
                                     timepoint=="day 3" ~ "d3",
                                     timepoint=="day 5" ~ "d5",
                                     timepoint=="day 15" ~ "d15",
                                     timepoint=="day 29" ~ "d29",
                                     timepoint %in% c("relapse","relapse 1") ~ "rel",
                                     timepoint=="relapse 2" ~ "rel2",
                                     TRUE ~ "other"),levels=timepoints_short),
    clone_time=case_when(
      !is.na(clone) ~ paste0(clone, "_", timepoint_short),
      TRUE ~ NA
    )
  )

DimPlot(dat, group.by="clone", cols=clone_cols)
DimPlot(dat, group.by="timepoint_short", cols=timepoint_cols_short)
DimPlot(dat, group.by="clone_time", order=T)

genes <- read_tsv("/wrk/resources/genomes/grch38-gdc/star_index/geneInfo.tab",skip=1,col_names=c("gene_id","gene_name","gene_type"))
g <- setNames(genes$gene_id, nm=genes$gene_name)

VlnPlot(dat, c(g["DNTT"],g["HDAC9"]), group.by="clone_time")
VlnPlot(dat, g[c("SDK2","STAT1","DDIT4")], group.by="clone_time")
VlnPlot(dat, "HALLMARK_MYC_TARGETS_V2_UCell", group.by="clone_time")

# Signatures with UCell
# https://bioconductor.org/packages/release/bioc/vignettes/UCell/inst/doc/UCell_Seurat.html
# remotes::install_github("carmonalab/SignatuR")
library(UCell)
library(SignatuR)
library(msigdbr)

ig <- lapply(GetSignature(SignatuR$Hs$Compartments$Immunoglobulins), function(x) g[x])
dat <- AddModuleScore_UCell(dat, features=ig, name=NULL)

hs <- as.data.frame(msigdbr(species = "Homo sapiens", category = "H"))
names <- unique(hs$gs_name)
h.markers <- list()
for (h in names) {
  sub <- hs[hs$gs_name==h,"ensembl_gene"]
  h.markers[[h]] <- unique(sub)
}
dat <- AddModuleScore_UCell(dat, features = h.markers, ncores=14)
signature.names <- paste0(names(h.markers),"_UCell")

# dat <- SetIdent(dat, value=dat$clone_time)
# m <- FindMarkers(dat, ident.1 = "5_1_rel", ident.2 = "5_1_rel2", group.by="clone_time")

# Setup gene/pathway names
sets_to_plot <- c(
  "HALLMARK_APOPTOSIS_UCell",
  "HALLMARK_MYC_TARGETS_V2_UCell",
  "nCount_RNA"
)
genes_to_plot <- c(
"ENSG00000069188",
"ENSG00000115415",
"ENSG00000168209"
)

# Extract expression data and calculate means per clone+timepoint
expression_data <- data.frame(
  clone = dat$clone,
  timepoint = dat$timepoint
)

# Add expression data for each feature
for(set in sets_to_plot) {
  expression_data[[set]] <- dat@meta.data[[set]]
}
for(gene in genes_to_plot) {
  expression_data[[names(g)[g==gene]]] <- dat@assays$RNA@scale.data[gene,]
}

# Calculate means per clone+timepoint
all_features <- colnames(expression_data)[3:ncol(expression_data)]
expression_means <- expression_data %>%
  filter(!is.na(clone)) %>%
  group_by(clone, timepoint) %>%
  summarise(across(all_of(all_features), mean, na.rm = TRUE), .groups = "drop") %>% 
  # Z-scale each feature across all groups
  mutate(across(all_of(all_features), scale))

