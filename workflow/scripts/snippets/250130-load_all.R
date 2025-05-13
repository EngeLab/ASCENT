# Load clones and annotate RNA data
library(tidyverse)
library(scales)
library(Seurat)
library(patchwork)

# Basics
source("workflow/scripts/general.R")
source("workflow/scripts/snippets/250131-clone_clouds.R")

# Cohort
mrd <- c("ALL5","ALL8","ALL47","ALL64","ALL66","ALL67","ALL68","ALL71")
rel <- c("ALL1","ALL3","ALL4","ALL6", "ALL35","ALL40","ALL42")
excl <- c("ALL2","ALL7","ALL55")

# Load DNA data
diploids <- read_tsv("results/ALL-latest/analysis/all_diploid_clones-250201.txt")
clone_files <- list.files("results/ALL-latest/", pattern = "*-clones_final.txt", recursive = T, full.names = T)

# Modify timepoint data
timepoint_cols["diagnosis"] <- "#cfcfcf"
timepoints_short <- c("d0", "d2", "d3", "d5", "d15", "d29", "rel", "rel2","other")
timepoint_cols_short <- c("d0"="#cfcfcf","d2"="#fee08b","d3"="#d9ef8b","d5"="#91cf60","d15"="#1a9850","d29"="#136b39","rel"="#d73027", "rel2"="#002984", "other"="#002984")

cl <- tibble(
  file = clone_files,
  patient_id = str_extract(basename(file), "ALL\\d+"),
  data = map(file, ~read_tsv(.x))
) %>% select(-file) %>% 
  unnest(data)

cl.annot <- cl %>% left_join(mutate(diploids, dna_type="diploid"), by=c("patient_id", "clone"="diploid_clone")) %>% 
  mutate(dna_type=case_when(
    is.na(clone) | clone=="0" ~ NA,
    is.na(dna_type) ~ "aneuploid",
    .default = dna_type
  ))

# Load RNA data
load_pids <- unique(cl.annot$patient_id)
dl <- sapply(load_pids, function(x) {
  cat(sprintf("Loading seurat object for: %s", x),"\n")
  x <- read_rds(paste0("results/ALL-latest/",x,"/",x,"-rna_seurat.Rds"))
  x$clone <- cl.annot$clone[match(x$dna_library_id, cl.annot$dna_library_id)]
  x$dna_type <- cl.annot$dna_type[match(x$dna_library_id, cl.annot$dna_library_id)]
  x$dna_reads <- cl.annot$dna_reads[match(x$dna_library_id, cl.annot$dna_library_id)]
  # Shorten timepoint labels
  x@meta.data <- x@meta.data %>% 
    mutate(
      timepoint_short=factor(case_when(timepoint=="diagnosis" ~ "d0",
                                       timepoint=="day 2" ~ "d2",
                                       timepoint=="day 3" ~ "d3",
                                       timepoint=="day 5" ~ "d5",
                                       timepoint=="day 15" ~ "d15",
                                       timepoint=="day 29" ~ "d29",
                                       timepoint %in% c("relapse","relapse 1") ~ "rel",
                                       timepoint=="relapse 2" ~ "rel2",
                                       TRUE ~ "other"),levels=timepoints_short)
    )
  return(x)
})
str(dl, max.level=3)

# Plots
# pid <- "ALL3"
pid <- "ALL6"

# pdf("results/ALL-latest/clone_rna_overview-250130.pdf", width=12, height=6)
pdf("results/ALL-latest/ALL4_overview-250201.pdf", width=12, height=6)
for(pid in names(dl)){
# All cells
da <- dl[[pid]]
# Blasts only
db <- dl[[pid]][,dl[[pid]]$dna_type %in% "aneuploid"]
db$clone_time <- paste(db$clone, db$timepoint_short, sep="_")
clone_cols <- get_clones_col(da$clone)
clone_time_cols <- setNames(timepoint_cols_short[sub(".*_(d.*)$","\\1", db$clone_time)],nm=db$clone_time)
clone_time_mid_cols <- setNames(clone_cols[sub("(.*)_d.*$","\\1", db$clone_time)],nm=db$clone_time)

p1 <- DimPlot(da, group.by="dna_type", order=T, pt.size=.6) + theme_dntr(axis_labels = F, font_size = 10)
p2 <- DimPlot(da, group.by="clone",order=T, pt.size=.6) + theme_dntr(axis_labels = F, font_size = 10) + scale_color_manual(values=clone_cols)
p3 <- DimPlot(da, group.by="timepoint", order=T, pt.size=.6) + theme_dntr(axis_labels = F, font_size = 10) + scale_color_manual(values=timepoint_cols)
p4 <- DimPlot(da, group.by="cell_phase", order=T, pt.size=.6) + theme_dntr(axis_labels = F, font_size = 10) + scale_color_manual(values=phase_cols)
p5 <- FeaturePlot(da, "ENSG00000107447", order=T, pt.size=.6, min.cutoff = "q5", max.cutoff = "q95") + theme_dntr(axis_labels = F, font_size = 10)
# p6 <- FeaturePlot(da, "ENSG00000156738", order=T, pt.size=.6, min.cutoff = "q5", max.cutoff = "q95") + theme_dntr(axis_labels = F, font_size = 10)
p6a <- plot_cloud(db, clone_slot="clone", threads=20) + scale_color_manual(values=clone_cols)
pall <- p1 + p2 + p3 + p4 + p5 + p6a + plot_layout(nrow=1)

p1 <- DimPlot(db, group.by="dna_type", order=T, pt.size=.6) + theme_dntr(axis_labels = F, font_size = 10)
p2 <- DimPlot(db, group.by="clone",order=T, pt.size=.6) + theme_dntr(axis_labels = F, font_size = 10) + scale_color_manual(values=clone_cols)
p3 <- DimPlot(db, group.by="timepoint", order=T, pt.size=.6) + theme_dntr(axis_labels = F, font_size = 10) + scale_color_manual(values=timepoint_cols)
p4 <- DimPlot(db, group.by="cell_phase", order=T, pt.size=.6) + theme_dntr(axis_labels = F, font_size = 10) + scale_color_manual(values=phase_cols)
p5 <- FeaturePlot(db, "ENSG00000107447", order=T, pt.size=.6, min.cutoff = "q5", max.cutoff = "q95") + theme_dntr(axis_labels = F, font_size = 10)
pbc <- plot_cloud(db, clone_slot="clone_time", threads=20)
p6b <-  pbc + scale_color_manual(values=clone_time_cols) + scale_fill_manual(values=clone_time_mid_cols)
pblasts <- p1 + p2 + p3 + p4 + p5 + p6b + plot_layout(nrow=1)

pcomb <- pall / pblasts + plot_annotation(title=paste0(pid, " (",ncol(da)," cells, ", ncol(db)," blasts)"))
plot(pcomb)
}
dev.off()

# Merge objects
# d <- dl[[pid]][,dl[[pid]]$dna_type %in% "aneuploid"]


# Clouds for specific patient
pid = ALL66
