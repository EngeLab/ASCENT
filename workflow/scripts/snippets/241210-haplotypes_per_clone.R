# 241210: clone and segment-level major/minor haplotype
library(tidyverse)
library(parallel)
library(GenomicRanges)
library(caTools)
library(patchwork)
options(scipen=999)
source("workflow/scripts/clone_functions_VZ.R")
threads <- 12
segs_slot <- "refined"
# clone_slot <- "recall"

# Get object from clone_segmentation_VZ.R
patient_id <- "ALL40"
run <- "results/DNA_37_SK_newRNA/"
cn_obj <- file.path(run, patient_id, paste0(patient_id,"-clones.Rds"))
cn_obj<-"/wrk/data/solrun/AnalysisPipelinePaper/ALL40/ALL40-bins-10-clones-goldstandard-useclonecopynumber_removeBadCells.Rds"
d <- read_rds(cn_obj)

clones <- d$cells %>% 
  filter(!is.na(clone_final)) %>% 
  select(cell=dna_library_id, clone=clone_final)


snps.all <- read_tsv("/wrk/resources/dntr2/results/DNA_37_SK_newRNA/phasing/ALL40-baf.txt.gz", col_names=c("chr","pos","cell","hapA_count","hapB_count")) %>% 
  left_join(clones, by="cell") %>% 
  mutate(chr=factor(chr, levels=paste0("chr",1:22)))

# summed_snps <- snps.all %>%
#   group_by(chr, pos) %>%
#   summarize(
#     # # Check if position is heterozygous in diploid clone
#     # diploid_A = sum(hapA_count[clone == diploid]),
#     # diploid_B = sum(hapB_count[clone == diploid]),
#     # Total counts across all clones
#     total_A = sum(hapA_count),
#     total_B = sum(hapB_count)
#   )
# filtered_snps <- summed_snps %>%
#   # Keep position if:
#   # - It has both alleles in diploid clone (A≥1 AND B≥1)
#   # - OR it has stronger evidence across all clones (A≥2 AND B≥2)
#   filter(
#     # (diploid_A >= 1 | diploid_B >= 1)
#     (total_A >= 5 & total_B >= 5)
#   )
# snps <- snps.all %>%
#   semi_join(filtered_snps, by = c("chr", "pos"))

snps <- snps.all
snps.list <- split(snps, snps$chr)
block_size <- 50000

### Sum per position. Parallelized per chr. Group into phased blocks 50-100kb
pos <- mclapply(snps.list, function(x){
  x %>% 
    group_by(chr, pos) %>%
    summarise(
      hapA = sum(hapA_count),
      hapB = sum(hapB_count),
      total = hapA + hapB,
      BAF = hapB / total,
      mBAF = pmin(BAF, 1-BAF),
      n_cells = n(),
      mean_coverage = (sum(hapA_count) + sum(hapB_count)) / n(),
      .groups = 'drop'
    ) %>% arrange(pos) %>% 
    mutate(block_id = ceiling(pos / block_size)) %>% 
    # Stricter filtering
    filter(n_cells >= 2)
}, mc.cores = threads)

### Block summary
blocks <- mclapply(pos, function(x){
  x %>%  
    group_by(block_id) %>% 
    summarise(
      pos=mean(pos),
      hapA = sum(hapA),
      hapB = sum(hapB),
      total = hapA + hapB,
      BAF = hapB / total,
      mBAF = pmin(BAF, 1-BAF),
      n_snps = n(),
      snp_dens_kb=1000*n_snps/block_size,
      mean_coverage = (sum(hapA) + sum(hapB)) / n(),
      major = ifelse( hapB / total > 0.5, hapB, hapA),
      minor = ifelse( hapB / total > 0.5, hapA, hapB),
      switch = major != hapA,
      .groups = 'drop'
    )
  # filter(n_snps >= min_snps_per_block, total >= min_reads_per_block)
}, mc.cores=threads)


### Add back to cell-level data (snps.list)
snps.phased <- mclapply(seq_along(snps.list), function(i){
  snps.list[[i]] %>% 
    left_join(select(pos[[i]], chr, pos, block_id), by=c("chr","pos")) %>% 
    left_join(select(blocks[[i]], block_id, switch), by="block_id") %>% 
    mutate(cA=ifelse(switch, hapB_count, hapA_count), 
           cB=ifelse(switch, hapA_count, hapB_count)) %>% 
    filter(!is.na(switch))
}, mc.cores=threads)
names(snps.phased) <- names(snps.list)

### Group by clone and block
res <- mclapply(snps.phased, function(x){
  xl <- split(x, x$clone)
  lapply(xl, function(xclone){
    xclone %>% 
        group_by(block_id) %>% 
      summarize(
        chr=unique(chr),
        start=min(pos),
        mid=mean(pos),
        end=max(pos),
        A=sum(cA),
        B=sum(cB),
        total=A+B,
        BAF=B/total,
        clone=unique(clone),
        n_snps=n(),
        snp_dens_kb=1000*n_snps/block_size
      )
  })
}, mc.cores=threads)

### 250223: mpcf segmentation


### Get CN and segment data
# source("scripts/snippets/241205-all3data.R")
segs <- d$segments[[segs_slot]]
segs <- d$segments[[segs_slot]] %>% 
  dplyr::rename(id=seg_idx) %>% 
  mutate(local_index=1:nrow(segs))

### Convert to clone-level segments
# 1. Restrict to arm-level boundaries
chr_arm <- factor(paste0(segs$chr, segs$arm), 
                  levels=unique(paste0(segs$chr, segs$arm)))
chr_arm.l <- split(1:length(chr_arm), chr_arm)

# 2. Find segments per clone
cn_clones <- do.call(cbind, d$clones$cn)[d$segments[[segs_slot]]$start,]
colnames(cn_clones) <- d$clones$clone_id
dim(cn_clones)

clone_breaks <- apply(cn_clones, 2, function(clone_cn){
  lapply(chr_arm.l, function(id) {
    s <- clone_cn[id]
    c(TRUE, diff(s) != 0)
  }) %>% unlist
})
head(clone_breaks)

clone_segs <- lapply(colnames(clone_breaks), function(i){
  x <- clone_breaks[,i]
  change_idx <- which(x)
  starts <- segs$start.pos[change_idx]
  ends <- segs$end.pos[c(change_idx[-1] - 1, nrow(segs))]
  bin_starts <- segs$start[change_idx]
  bin_ends <- segs$end[c(change_idx[-1] - 1, nrow(segs))]
  r <- tibble(chr=segs$chr[change_idx], 
              arm=segs$arm[change_idx],
              start.pos=starts,
              end.pos=ends,
              start=as.integer(bin_starts),
              end=as.integer(bin_ends),
              cn=cn_clones[change_idx,i],
              # id=segs$id[change_idx])
              id=segs$local_index[change_idx]) # 250222(VZ)
  r$gcmap <- calc_segment_means(d$clones$gcmap_normal[[which(d$clones$clone_id==i)]], r)
  return(r)
})
names(clone_segs) <- colnames(clone_breaks)
sapply(clone_segs, nrow)

clone_segs.gr <- lapply(clone_segs, function(x){
  GRanges(
    seqnames = x$chr,
    ranges = IRanges(start = x$start.pos, end = x$end.pos),
    arm = x$arm,
    id = x$id,
    cn = x$cn,
    gcmap = x$gcmap,
  )
})

# Process clone blocks by finding overlaps with segments
# Outer loop by clone (df of blocks). Inner loop by chr
comb <- lapply(names(clone_segs.gr), function(clone) {
  segs_gr <- clone_segs.gr[[clone]]
  # cat("\n\n",clone)
  inner <- lapply(names(res), function(chr){
    # cat(chr,"..")
    segs_gr.f <- segs_gr[seqnames(segs_gr) == chr]
    
    blocks <- res[[chr]][[clone]]
    blocks_gr <- GRanges(
      seqnames = blocks$chr,
      ranges = IRanges(start = blocks$start, end = blocks$end)
    )
    
    hits <- findOverlaps(segs_gr.f, blocks_gr)
    
    # Convert to indices
    seg_idx <- as.integer(queryHits(hits))
    block_idx <- as.integer(subjectHits(hits))
    
    # Create result tibble
    result <- tibble(
      chr = as.character(seqnames(segs_gr.f)),
      start = start(segs_gr.f),
      end = end(segs_gr.f),
      arm = segs_gr.f$arm,
      bins = segs_gr.f$size,
      id = segs_gr.f$id,
      cn = segs_gr.f$cn,
      gcmap = segs_gr.f$gcmap
    ) %>%
      group_by(chr, start, end, arm) %>%
      mutate(
        alleles = list(blocks[block_idx[seg_idx == cur_group_id()],]), # ?
        A = map_dbl(alleles, ~sum(.$A)),
        B = map_dbl(alleles, ~sum(.$B)),
        total = map_dbl(alleles, ~sum(.$total)),
        BAF = B/total,
        mBAF = pmin(BAF, 1-BAF),
        avgBAF = map_dbl(alleles, ~mean(.x$BAF,na.rm=T)),
        avg_mBAF = pmin(avgBAF, 1-avgBAF)
      ) %>%
      # select(-alleles) %>%
      ungroup()
    result
  })
  
  inner %>% bind_rows
})
names(comb) <- names(clone_segs.gr)

### Plot individual chr
chr="chr12"
snp_dens_limit = 0.2
plist <- list()
for(cl in names(comb)){
  # seg_idx <- segs$chr == chr
  idx <- comb[[cl]]$chr==chr
  seg_bnd <- unique(c(clone_segs[[cl]]$start.pos[idx], clone_segs[[cl]]$end.pos[idx]))
  # raw_baf <- bind_rows(comb[[cl]][idx,"alleles"])
  raw_baf <- do.call(rbind, comb[[cl]]$alleles[idx])
  plist[[cl]] <- comb[[cl]][idx,] %>% 
    ggplot() +
    geom_point(data=filter(raw_baf, snp_dens_kb > snp_dens_limit), aes(mid, runmean(BAF,k=5)), alpha=0.1, size=0.1) +
    geom_line(data=raw_baf, aes(mid, runmean(snp_dens_kb,k=20)), col="orange", alpha=0.2) +
    geom_segment(aes(x = start, xend = end, y = mBAF, yend = mBAF)) +
    geom_segment(aes(x = start, xend = end, y = avg_mBAF, yend = avg_mBAF),color="grey") +
    theme_dntr(axis_ticks = T) +
    labs(x = "Position", y = "mBAF") +
    scale_y_continuous(limits=c(0,0.5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(5),
                       labels = function(x) paste0(x/1e6,"Mb")) + 
    geom_hline(yintercept=c(0.2, 0.25, 0.33), linetype="dotted") +
    geom_hline(yintercept=snp_dens_limit, color="orange", linetype="solid", alpha=0.2) +
    geom_vline(xintercept=seg_bnd, color="red") +
    labs(subtitle=paste(chr,", clone ",cl))
}
wrap_plots(plist, ncol=2)

### Plot BAF x CN/gcmap
pl2 <- list()
id_colors <- sample(scales::hue_pal()(max(sapply(comb, function(x) max(x$id)))))
# id_colors <- sample(scales::hue_pal()(max(sapply(comb, function(x) max(x$id, na.rm=T)))))
names(id_colors) <- 1:max(sapply(comb, function(x) max(x$id)))
ymax <- max(sapply(comb, function(x) max(x$gcmap)))*1.1
for(cl in names(comb)){
  raw_baf <- pmap_dfr(
    list(
      alleles = comb[[cl]]$alleles,
      parent_id = comb[[cl]]$id,
      parent_cn = comb[[cl]]$cn,
      parent_gcmap = comb[[cl]]$gcmap
    ),
    function(alleles, parent_id, parent_cn, parent_gcmap) {
      alleles %>%
        mutate(
          id = parent_id,
          cn = parent_cn,
          gcmap = parent_gcmap
        )
    }
  ) %>%
    filter(snp_dens_kb > 0.08,
           total > 15, 
           BAF > 0, 
           cn > 0)
  
  pl2[[cl]] <- comb[[cl]] %>% 
    filter(cn > 0) %>%
    ggplot() +
    geom_hline(yintercept=seq(0,ymax,0.5), linetype="dotted", color="grey") +
    geom_vline(xintercept=seq(0.1, 0.5, 0.1), linetype="dotted", color="grey") +
    geom_point(data=raw_baf, aes(x=BAF, y=gcmap, 
                                 fill=ifelse(gcmap > 0.9 & gcmap < 1.1, "CN", id),
                                 size=snp_dens_kb), 
               # pch=21, alpha=0.1, position=position_jitter(height=0.1, width=0), 
               pch=21, alpha=0.1, stroke=0, show.legend = F) +
    geom_point(aes(x=mBAF, y=gcmap, 
                   fill=ifelse(gcmap > 0.9 & gcmap < 1.1, "CN", id)),
                   size=3, stroke=1, color="black", pch=21, show.legend = F) +
    scale_fill_manual(values=c("CN"="grey", id_colors)) +
    xlim(c(0,0.5)) +
    theme_dntr(axis_ticks = T) +
    scale_y_continuous(limits=c(0,ymax), breaks = scales::pretty_breaks(3)) +
    labs(subtitle=paste0(cl))
}


pl2[[4]]
comb

#It's f

### Get major/minor cn
baf.l <- lapply(comb, function(x){
  x[,c(1:7,11,12)]
})

map_to_original_segs <- function(reduced_segs, original_segs, value_cols) {
  # Lookup table for original vs reduced/clone-level segments
  # segment_mapping <- sapply(original_segs$id, function(orig_id) {
    # orig_seg <- original_segs[original_segs$id == orig_id, ]
  segment_mapping <- sapply(original_segs$local_index, function(orig_id) {
    orig_seg <- original_segs[original_segs$local_index == orig_id, ]
    containing_seg <- which(
      reduced_segs$chr == orig_seg$chr &
        reduced_segs$start <= orig_seg$start.pos &
        reduced_segs$end >= orig_seg$end.pos
    )
    
    if(length(containing_seg) == 0) return(NA)
    if(length(containing_seg) > 1) warning(sprintf("Original segment %d matches multiple reduced segments", orig_id))
    return(containing_seg[1])
  })
  
  # Add specific values back to original seg file
  # result <- original_segs %>% select(chr,arm,start.pos,end.pos,n.probes,id,start,end,filtered)
  result <- original_segs
  for(col in value_cols) {
    result[[col]] <- reduced_segs[[col]][segment_mapping]
  }
  
  return(result)
}

bafs <- lapply(baf.l, function(x){
  r <- map_to_original_segs(x, segs, c("BAF","cn","gcmap","total"))
  r %>% 
    mutate(
      BAF=ifelse(cn==0, 0, BAF), # fix to BAF 0 with 0 copies, for medicc
      cn_a=round((1-BAF)*cn),
      cn_b=round(BAF*cn)
    )
})

# Add to cnv object
d$clones$baf <- lapply(bafs, function(x) x$BAF[rep(1:nrow(x), x$n.probes)])
d$clones$cn_a <- lapply(bafs, function(x) x$cn_a[rep(1:nrow(x), x$n.probes)])
d$clones$cn_b <- lapply(bafs, function(x) x$cn_b[rep(1:nrow(x), x$n.probes)])

# Medicc2 input
# sample_id, chrom, start, end, cn_a, cn_b
medicc <- bafs %>% bind_rows(.id="sample_id") %>% 
  select(sample_id, chrom=chr, start=start.pos, end=end.pos, cn_a, cn_b) %>% 
  mutate(chrom=factor(chrom, levels=paste0("chr",c(1:22,"X","Y")))) %>% 
  arrange(chrom, start, end) %>% 
  filter(chrom %in% paste0("chr",1:22) &
         (!is.na(cn_a) & !is.na(cn_b)))

# Optional: remove 
medicc <- medicc %>% filter()
medicc_file <- file.path(run,patient_id,paste0(patient_id,"-medicc.txt"))
write_tsv(medicc, file=medicc_file)

# Run medicc2 (in bash)
# conda activate medicc
# medicc2 -n 1_1 results/ALL-latest/ALL3/ALL3-medicc.txt results/ALL-latest/ALL3/medicc

# Plot BAFs. Sort from diploid 
library(ape)
library(ComplexHeatmap)
segs.m <- do.call(cbind, lapply(d$clones$cn, function(x) x[d$segments$refined$start]))
colnames(segs.m) <- d$clones$clone_id
o <- sort_from_diploid_root(segs.m, na.rm = T)

pdf(file.path(run,patient_id,paste0(patient_id,"-clone_heatmaps.pdf")), width=6, height=3)
plot_clone_heatmap(d, order=o, clone_types=get_clone_types(d, diploid="1_1"))
plot_clone_heatmap(d, cn_slot="baf", show_chr=T, clone_types=NULL, order=o)
plot_clone_heatmap(d, cn_slot="cn_a", show_chr=T, clone_types=NULL, highlight_dups = F, order=o)
plot_clone_heatmap(d, cn_slot="cn_b", show_chr=T, clone_types=NULL, highlight_dups = F, order=o)
dev.off()

write_rds(d, file=cn_obj)