library(tidyverse)
library(parallel)
library(GenomicRanges)
library(caTools)
library(patchwork)
options(scipen=999)
segs_slot <- "refined"
source("/wrk/resources/dntr2_extrascp/workflow/scripts/clone_functions_forPaper.R")

patient_id <- "ALL40"
run <- "/wrk/resources/dntr2_extrascp/results/DNA_37_SK_250405"

cn_obj<-"final_clone_object.Rds"
d <- read_rds(cn_obj)

clones <- d$cells %>% 
  filter(!is.na(clone_final)) %>% 
  select(cell=dna_library_id, clone=clone_final)


snps.all <- read_tsv(paste0(run, "/phasing/", patient_id, "-baf.txt.gz"), col_names=c("chr","pos","cell","hapA_count","hapB_count")) %>% 
  left_join(clones, by="cell") %>% 
  mutate(chr=factor(chr, levels=paste0("chr",1:22)))

snps <- snps.all
snps.list <- split(snps, snps$chr)
block_size <- 50000

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
}, mc.cores = 16)

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
}, mc.cores=16)


### Add back to cell-level data (snps.list)
snps.phased <- mclapply(seq_along(snps.list), function(i){
  snps.list[[i]] %>% 
    left_join(select(pos[[i]], chr, pos, block_id), by=c("chr","pos")) %>% 
    left_join(select(blocks[[i]], block_id, switch), by="block_id") %>% 
    mutate(cA=ifelse(switch, hapB_count, hapA_count), 
           cB=ifelse(switch, hapA_count, hapB_count)) %>% 
    filter(!is.na(switch))
}, mc.cores=16)
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
}, mc.cores=16)

### Get CN and segment data
segs <- d$segments$refined #%>% dplyr::rename(id=seg_idx)
segs <- segs %>% 
  mutate(
    local_index=1:nrow(segs)
  )

### Convert to clone-level segments
# 1. Restrict to arm-level boundaries
chr_arm <- factor(paste0(segs$chr, segs$arm), 
                  levels=unique(paste0(segs$chr, segs$arm)))
chr_arm.l <- split(1:length(chr_arm), chr_arm)

# 2. Find segments per clone
cn_clones <- do.call(cbind, d$clones$cn)[d$segments$refined$start,]
colnames(cn_clones) <- d$clones$clone_id

clone_breaks <- apply(cn_clones, 2, function(clone_cn){
  lapply(chr_arm.l, function(id) {
    s <- clone_cn[id]
    s[is.na(s)] <- (-1)
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
  print(i)
  r <- tibble(chr=segs$chr[change_idx], 
              arm=segs$arm[change_idx],
              start.pos=starts,
              end.pos=ends,
              start=as.integer(bin_starts),
              end=as.integer(bin_ends),
              cn=cn_clones[change_idx,i],
              #id=segs$id[change_idx])
              id=segs$local_index[change_idx])
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
chr="chr3"
snp_dens_limit = 0.05
plist <- list()
for(cl in names(comb)){
  # seg_idx <- segs$chr == chr
  idx <- comb[[cl]]$chr==chr
  seg_bnd <- unique(c(clone_segs[[cl]]$start.pos[idx], clone_segs[[cl]]$end.pos[idx]))
  # raw_baf <- bind_rows(comb[[cl]][idx,"alleles"])
  raw_baf <- do.call(rbind, comb[[cl]]$alleles[idx])
  plist[[cl]] <- comb[[cl]][idx,] %>% 
    ggplot() +
    geom_point(data=filter(raw_baf, snp_dens_kb > snp_dens_limit), aes(mid, runmean(BAF,k=5)), alpha=0.5, size=0.1) +
    geom_point(data=(raw_baf), aes(mid, runmean(BAF,k=5)), alpha=0.5, size=0.5) +
    #geom_line(data=raw_baf, aes(mid, runmean(snp_dens_kb,k=5)), col="orange", alpha=1) +
    geom_segment(aes(x = start, xend = end, y = mBAF, yend = mBAF), color="grey", linewidth=2) +
    geom_segment(aes(x = start, xend = end, y = avg_mBAF, yend = avg_mBAF),color="black", linewidth=2) +
    theme_dntr(axis_ticks = T) +
    labs(x = "Position", y = "mBAF") +
    scale_x_continuous(breaks = scales::pretty_breaks(5),
                       labels = function(x) paste0(x/1e6,"Mb")) + 
    #geom_hline(yintercept=c(0.2, 0.25, 0.33), linetype="dotted") +
    #geom_hline(yintercept=snp_dens_limit, color="orange", linetype="solid", alpha=0.2) +
    geom_vline(xintercept=seg_bnd, color="red") +
    labs(subtitle=paste(chr,", clone ",cl))
}

wrap_plots(plist, ncol=2)

#
baf.l <- lapply(comb, function(x){
  x[,c(1:7,11,12)]
})

map_to_original_segs <- function(reduced_segs, original_segs, value_cols) {
  # Lookup table for original vs reduced/clone-level segments
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

#If we have less then 5 snps then we set a as major 

bafs_f <- lapply(baf.l, function(x) {
  r <- map_to_original_segs(x, segs, c("BAF", "cn", "gcmap", "total"))
  r <- r %>%
    mutate(
      BAF = ifelse(cn == 0, 0, BAF), 
      cn_a = ifelse(total < 5, 
                    case_when(
                      cn == 1 ~ 1,
                      cn == 2 ~ 1,
                      cn == 3 ~ 2,
                      cn == 4 ~ 2,
                    ),
                    round((1 - BAF) * cn)),
      cn_b = ifelse(total < 5, 
                    case_when(
                      cn == 1 ~ 0,
                      cn == 2 ~ 1,
                      cn == 3 ~ 1,
                      cn == 4 ~ 2,
                    ),
                    round(BAF * cn)))
  return(r)
})


bafs_forPlot <- lapply(baf.l, function(x) {
  r <- map_to_original_segs(x, segs, c("BAF", "cn", "gcmap", "total"))
  r <- r %>%
    mutate(
      BAF = ifelse(cn == 0, 0, BAF), 
      cn_a = ifelse(total < 5, 
                    case_when(
                      cn == 0 ~ 0.1,
                      cn == 1 ~ 1.1,
                      cn == 2 ~ 1.1,
                      cn == 3 ~ 2.1,
                      cn == 4 ~ 2.1,
                    ),
                    round((1 - BAF) * cn)),
      cn_b = ifelse(total < 5, 
                    case_when(
                      cn == 0 ~ 0.1,
                      cn == 1 ~ 0.1,
                      cn == 2 ~ 1.1,
                      cn == 3 ~ 1.1,
                      cn == 4 ~ 2.1,
                    ),
                    round(BAF * cn)))
  return(r)
})

# Assign a unique number to each combination, but set to NA if it includes NA/NaN

assign_baf_cn <- function(df, cn_combination_mapping) {
  df <- df %>%
    mutate(baf_cn = ifelse(
      is.na(cn_a) | is.na(cn_b) | is.nan(cn_a) | is.nan(cn_b), 
      NA,  # If either cn_a or cn_b is NA or NaN, set baf_cn to NA
      cn_combination_mapping[paste(cn_a, cn_b)]  # Use the mapping to assign baf_cn
    ))
  return(df)
}

combined_cns <- do.call(rbind, lapply(bafs_forPlot, function(df) df[, c("cn_a", "cn_b")]))
unique_combinations <- unique(paste(combined_cns$cn_a, combined_cns$cn_b))

cn_combination_mapping <- setNames(seq_along(unique_combinations), unique_combinations)
library(scales)
gray_colors <- grey_pal()(seq(0, 1, length.out = 6))


bafs_forPlot <- lapply(bafs_forPlot, assign_baf_cn, cn_combination_mapping = cn_combination_mapping)


# Apply the function to all elements in bafs
d$clones$allelespecific <- lapply(bafs_forPlot, function(x) x$baf_cn[rep(1:nrow(x), x$n.probes)])
library(ape)
#OK figure out plotting and so 

pdf(paste0(run, "/", patient_id, "/clones/", patient_id, "-allelespecificCloneHeatmap.pdf" width=20, height=10)
plot_clone_heatmap_2(d, cn_slot="allelespecific", show_chr=T)
dev.off()

medicc<-bafs_f %>% 
  bind_rows(.id = "sample_id") %>% 
  select(sample_id, chrom = chr, start = start.pos, end = end.pos, cn_a, cn_b) %>% 
  mutate(chrom = factor(chrom, levels = paste0("chr", c(1:22, "X", "Y")))) %>% 
  arrange(chrom, start, end) %>% 
  filter(chrom %in% paste0("chr", 1:22)) %>% 
  group_by(start) %>% 
  filter(!any(is.na(cn_a) | is.na(cn_b))) %>%  # Remove groups with any NA in cn_a or cn_b
  ungroup()

#
write_tsv(medicc, file=paste0(run, "/", patient_id, "/", patient_id, "-mediccinput.tsv")

# Run medicc2 (in bash)
# conda activate medicc
# medicc2 -n diploid_clone mediccinput.tsv results/medicc

