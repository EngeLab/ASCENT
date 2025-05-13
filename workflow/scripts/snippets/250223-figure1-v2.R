library(tidyverse)
library(ComplexHeatmap)
library(circlize)

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

# Original copy number processing code
clone_freqs <- clones %>%
  group_by(clone, timepoint) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1) %>%
  group_by(timepoint) %>%
  mutate(rel_freq = count / sum(count)) %>%
  arrange(timepoint, rel_freq) %>%
  ungroup()

segs_long <- segs %>%
  pivot_longer(cols = c("2_1", "3_1", "4_1", "5_1", "6_1", "7_1"),
               names_to = "clone",
               values_to = "copy_number") %>%
  mutate(segment = sprintf("%s%s:%.1f-%.1fMb", 
                           gsub("chr", "", chr),
                           arm,
                           start.pos/1e6,
                           end.pos/1e6))

# Create the copy number matrix
cn_matrix <- clone_freqs %>%
  left_join(segs_long, by = "clone") %>%
  select(segment, timepoint, clone, copy_number) %>%
  pivot_wider(id_cols = segment,
              names_from = c(timepoint, clone),
              values_from = copy_number) %>%
  column_to_rownames("segment") %>%
  as.matrix()

# Get column info for proper ordering and splits
col_info <- tibble(
  full_name = colnames(cn_matrix),
  timepoint = sub("_.*$", "", full_name),
  clone = sub("^[^_]*_", "", full_name)
) %>%
  left_join(clone_freqs, by = c("timepoint", "clone")) %>%
  arrange(timepoint, rel_freq)

# Reorder matrix columns
cn_matrix <- cn_matrix[, col_info$full_name]
colnames(cn_matrix) <- col_info$clone

# Create expression matrix matching the column order
expr_matrix <- matrix(NA, nrow = length(all_features), 
                      ncol = ncol(cn_matrix),
                      dimnames = list(all_features, colnames(cn_matrix)))

for(i in seq_len(ncol(cn_matrix))) {
  current_clone <- colnames(cn_matrix)[i]
  current_timepoint <- col_info$timepoint[i]
  
  # Find corresponding expression data
  expr_data <- expression_means %>%
    filter(clone == current_clone, timepoint == current_timepoint)
  
  if(nrow(expr_data) > 0) {
    expr_matrix[, i] <- as.numeric(expr_data[, all_features])
  }
}

# Create frequency annotation data
freq_anno <- col_info$rel_freq
names(freq_anno) <- colnames(cn_matrix)

# Create color functions
cn_col_fun <- make_cn_colorscale(max(cn_matrix, na.rm=TRUE))
freq_col_fun <- colorRamp2(c(0, 1), c("white", "black"))
expr_col_fun <- colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027"))

# Create top annotation
ha <- HeatmapAnnotation(
  "Relative Frequency" = freq_anno,
  col = list("Relative Frequency" = freq_col_fun),
  show_legend = FALSE
)

# Adjust text sizes and dimensions
ht_list <- Heatmap(cn_matrix,
                   name = "Copy Number",
                   col = cn_col_fun,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_row_names = TRUE,
                   row_names_gp = gpar(fontsize = 6),  # smaller row names
                   show_column_names = TRUE,
                   column_names_gp = gpar(fontsize = 6),  # smaller column names
                   column_names_rot = 45,
                   column_split = col_info$timepoint,
                   # column_title = NULL,
                   column_title_gp = gpar(fontsize = 7),  # smaller title
                   border = TRUE,
                   column_gap = unit(2, "mm"),  # reduced gap
                   row_title = "Segments",
                   row_title_gp = gpar(fontsize = 7),  # smaller row title
                   top_annotation = ha,
                   show_heatmap_legend = FALSE
) %v%
  Heatmap(expr_matrix,
          name = "Expression",
          col = expr_col_fun,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_gp = gpar(fontsize = 6),  # smaller row names
          show_column_names = TRUE,
          column_names_gp = gpar(fontsize = 6),  # smaller column names
          column_names_rot = 45,
          column_split = col_info$timepoint,
          column_title = col_info$timepoint %>% unique(),
          column_title_gp = gpar(fontsize = 7),  # smaller title
          border = TRUE,
          column_gap = unit(2, "mm"),  # reduced gap
          show_heatmap_legend = FALSE
  )

# Save with specific dimensions
pdf("results/ALL-latest/ALL3_fig1test.pdf",width = 3, height=3)
draw(ht_list)
dev.off()
