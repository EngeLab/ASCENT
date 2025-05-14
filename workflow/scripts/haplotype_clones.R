#Get haplotype info within pipeline
library(tidyverse)
library(parallel)
library(GenomicRanges)
library(caTools)
library(patchwork)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(grid)
library(ape)
library(ggbeeswarm)


options(scipen=999)
segs_slot <- "refined"
source("workflow/scripts/clone_functions_forPaper.R")


cn_obj<-snakemake@input[["cn_obj"]]
d <- read_rds(cn_obj)

clones <- d$cells %>% 
  filter(!is.na(clone_final)) %>% 
  select(cell=dna_library_id, clone=clone_final)

snps.all <- read_tsv(snakemake@input[["phasing"]], col_names=c("chr","pos","cell","hapA_count","hapB_count")) %>% 
  left_join(clones, by="cell") %>% 
  mutate(chr=factor(chr, levels=paste0("chr",1:22)))

snps <- snps.all
snps.list <- split(snps, snps$chr)
block_size <- 50000

# Functions
block_snps <- function(snps, block_size = 50000, min_snps_global=0, min_snps_clone=0) {
  # Create block_id based on genomic position
  snps_with_blocks <- snps %>%
    mutate(block_id = ceiling(pos / block_size)) %>% 
    group_by(block_id) %>% 
    mutate(n_snps_all=n()) %>% 
    ungroup()
  
  # Summarize allele counts by chromosome, block, and clone
  blocks_by_clone <- snps_with_blocks %>%
    group_by(chr, block_id, clone) %>%
    summarize(
      n_snps_all = unique(n_snps_all),
      start_pos = min(pos),
      end_pos = max(pos),
      mid_pos = mean(pos),
      A = sum(hapA_count),
      B = sum(hapB_count),
      total = A + B,
      BAF = if_else(total > 0, B / total, NA_real_),
      n_snps = n(),
      .groups = 'drop'
    )
  
  blocks_filtered <- blocks_by_clone %>% filter(n_snps_all >= min_snps_global, n_snps >= min_snps_clone)
  return(blocks_filtered)
}

# 1. Group snps into blocks
snps.f <- snps.all %>% filter(! clone %in% c(0, NA)) # Filter out clone 0 and cells without clone assignment
snps.l <- split(snps.f, snps.f$chr)
blocks <- mclapply(snps.l, function(x){
  block_snps(x, block_size=50000, min_snps_global=0, min_snps_clone=0)
})

# 2. Get segments + cn integers. Make local index for segments for matching with blocks
# segs <- d$segments[[segs_slot]] %>% select(1:4) %>% mutate(seg_idx=1:nrow(.))
segs <- d$segments[[segs_slot]] %>% mutate(seg_idx=1:nrow(.))
cn_clones <- do.call(cbind, d$clones$cn)[d$segments[[segs_slot]]$start,]
colnames(cn_clones) <- d$clones$clone_id
dim(cn_clones)
nrow(segs) == nrow(cn_clones)
head(cn_clones)
head(segs)
segs.l <- split(segs, segs$chr)[names(blocks)] # subset to same chr with allelic info

# segments with no SNPs are left out

# 3. Integrate segments with block-level data
blocks.segs <- lapply(names(blocks), function(x){
  chr_blocks <- blocks[[x]]
  chr_segs <- segs.l[[x]]
  
  segs_gr <- GRanges(
    seqnames = chr_segs$chr,
    ranges = IRanges(start = chr_segs$start.pos, end = chr_segs$end.pos),
    seg_id = chr_segs$seg_idx,
  )
  
  blocks_gr <- GRanges(
    seqnames = chr_blocks$chr,
    ranges = IRanges(start = chr_blocks$start_pos, end = chr_blocks$end_pos),
    block_id = chr_blocks$block_id,
  )
  
  hits <- findOverlaps(blocks_gr, segs_gr)
  block_idx <- queryHits(hits)
  seg_idx <- subjectHits(hits)
  
  chr_blocks$seg_id <- NA
  chr_blocks$seg_id[block_idx] <- segs_gr$seg_id[seg_idx]
  
  idx <- cbind(chr_blocks$seg_id, match(chr_blocks$clone, colnames(cn_clones)))
  chr_blocks$cn <- cn_clones[idx]
  
  # Drop SNPs outside segments (not informative in this analysis)
  chr_blocks %>% filter(!is.na(seg_id))
})
names(blocks.segs) <- names(blocks)

# 4.  Rule based assignment of major/minor allele per segment
# Rules:
# 1. CN=1 -- no SNP filter applied. Remaining allele is major allele (=fipped to consistent minor/B allele=0)
# 2. CN>1 and significant skew in abs(BAF-0.5) --> likely CNLOH. Filter for n snps/block first
# 3. CN>2 -- BAF > 0.5 after filtering sets major allele. Filter for n snps/block first
# 4. None of the above, left as is from local phasing (eg random switches ~BAF 0.5)

# At clone level: 
# If more than one clone with same CN -- pick the one with highest coverage to fix major/minor

get_seg_priority <- function(x, min_snp_block=0, min_snp_segment=0, min_blocks_segment=0){
  x %>%
    filter(!is.na(cn)) %>% 
    mutate(snp_filter = n_snps >= min_snp_block) %>% # Allow filter by priority level/context later
    group_by(seg_id) %>% 
    mutate(
      seg_reads_all=sum(total),
      n_blocks_all=length(unique(block_id))
    ) %>% 
    group_by(seg_id, clone) %>%
    summarize(
      chr = as.character(chr[1]),
      cn = cn[1],
      n_blocks_all = n_blocks_all[1],
      seg_reads_all = seg_reads_all[1],
      n_blocks = n(),
      seg_snps = sum(n_snps),
      seg_reads = sum(total),
      seg_snps_filtered = sum(n_snps[snp_filter]),
      seg_reads_filtered = sum(total[snp_filter]),
      reads_frac = seg_reads/seg_reads_all,
      avg_BAF = mean(BAF, na.rm = TRUE),
      avg_mBAF = mean(pmin(BAF, 1-BAF), na.rm = TRUE),
      avg_BAF_deviation = mean(abs(BAF-0.5)),
      # Filtered versions
      avg_BAF_filtered = mean(BAF[snp_filter], na.rm = TRUE),
      avg_mBAF_filtered = mean(pmin(BAF[snp_filter], 1-BAF[snp_filter]), na.rm = TRUE),
      avg_BAF_deviation_filtered = mean(abs(BAF[snp_filter]-0.5)),
      .groups = 'drop'
    ) %>%
    # priority order
    mutate(
      filter_seg = seg_snps < min_snp_segment | n_blocks < min_blocks_segment,
      priority = case_when(
        filter_seg ~ 4, # If too low SNP coverage in segment, do not use
        cn == 1 ~ 1,
        cn > 1 & avg_mBAF_filtered < 0.1 ~ 2,  # CN-LOH candidate TODO: Does this need to be coverage normalized? 
        cn > 2 ~ 3, # Any gain
        TRUE ~ 4  # Normal CN=2
      )
    ) %>% 
    group_by(seg_id, priority)
}

# Add optional filter by priority level/context
min_snp_block <- 2
min_snp_segment <- 20 
min_blocks_segment <- 0
# Make filtered + unfiltered version for plotting overview
prio <- lapply(blocks.segs, get_seg_priority, min_snp_block=min_snp_block, min_snp_segment=min_snp_segment, min_blocks_segment=min_blocks_segment)
prio.u <- lapply(blocks.segs, get_seg_priority)

# Combine for plots
pa <- bind_rows(prio) %>% select(-avg_mBAF) %>% dplyr::rename(avg_mBAF=avg_mBAF_filtered)
pa.u <- bind_rows(prio.u)
pa.u$priority <- pa$priority
pa$dataset <- "filtered"
pa.u$dataset <- "unfiltered"
combined_data <- rbind(pa, pa.u)
# Precalculate vertical jittering for plots
combined_data <- combined_data %>%
  group_by(clone, seg_id) %>%
  mutate(
    jit_x = avg_mBAF,
    jit_y = cn + runif(n(), -0.15, 0.15)
  )

# Sum and order clones
clone_counts <- clones %>%
  count(clone) %>%
  arrange(desc(n)) %>%
  mutate(clone_label = paste0(clone, " (n=", n, ")"))
clone_order <- clone_counts$clone
combined_data$clone_ordered <- factor(combined_data$clone, levels = clone_order)
clone_labels <- setNames(clone_counts$clone_label, clone_counts$clone)

pc <- combined_data %>% 
  # filter(chr=="chr10") %>%
  ggplot(aes(jit_x, jit_y, col=factor(priority))) + 
  geom_line(aes(group=interaction(clone, seg_id)), alpha=0.2) +
  geom_point(aes(shape=dataset), alpha=0.7, size=3) +
  geom_vline(xintercept = 0.1, linetype="dotted") +
  scale_shape_manual(values=c(16, 21)) +
  scale_color_manual(values=c("1"="red","2"="orange","3"="darkgreen", "4"="grey")) +
  labs(x = "avg_mBAF", y = "cn") +
  theme_dntr(axis_ticks = T) + 
  facet_wrap(~ clone_ordered, labeller = labeller(clone_ordered = clone_labels))
pc

# Convert to lists by segment 
# a) segment level, prioritized by hierarchy above (1-4)
pa <- bind_rows(prio)
pl <- split(pa, pa$seg_id)
# b) block level BAF
ba <- bind_rows(blocks.segs) %>% mutate(chr_block_id=paste(chr, block_id, sep="_"))
bl <- split(ba, ba$seg_id)

# Only run informative segs (priority < 4)
informative_segs <- sapply(pl, function(x) any(x$priority < 4))
pl.f <- pl[informative_segs]
switch_list <- lapply(pl.f, function(x){
  s <- x$seg_id[1]
  prio <- min(x$priority)
  i <- x$priority == prio # Create index with highest prio
  # If more than one clone in segment with priority (length(i) > 1) -> pick the highest coverage one
  cl <- x$clone[i][which.max(x$reads_frac[i])] 
  
  # Get priority segment + clone --> determine major allele and create flip T/F vector
  x2 <- bl[[as.character(s)]] %>% filter(clone==cl)
  xfix <- x2 %>% mutate(
    major = ifelse(BAF > 0.5, B, A), #####
    minor = total - major,
    switch = major != A,
  )
  # Return flip vector (block id + switch T/F)
  xfix[,c("chr","block_id","switch")]
})
switches <- bind_rows(switch_list)
nrow(switches)
switches %>% head

# Join block level data with switch vector and priority info (from segment level)
ba.switched <- ba %>% 
  left_join(switches, by=c("chr","block_id")) %>% 
  left_join(select(pa, clone, seg_id, priority), by=c("clone", "seg_id")) %>% 
  # Corrected (c) A/B counts, BAF
  # For non-informative segments, cB is NA and BAF is used
  mutate(
    cA = case_when(switch ~ B,
                   TRUE ~ A), # Default (switch=F or NA) to A 
    cB = case_when(switch ~ A,
                   TRUE ~ B),
    cBAF = cB / total
  )
ba.sl <- split(ba.switched, ba.switched$chr)

switched.segs <- ba.switched %>% 
  # mutate(snp_filter = priority == 1 | n_snps >= min_snp_block ) %>% # Allow filter by priority level/context
  group_by(seg_id) %>% 
  mutate(
    seg_reads_all=sum(total),
    n_blocks_all=length(unique(block_id))
  ) %>%
  group_by(seg_id, clone) %>%
  summarize(
    # snp_filter = priority[1] == 1 | sum(total) >= min_snp_segment, # Filter: total reads / seg
    snp_filter = priority[1] == 1 | sum(n_snps) >= min_snp_segment, # Filter: total SNPs / seg
    chr = as.character(chr[1]),
    start_pos = min(start_pos),
    end_pos = max(end_pos),
    cn = cn[1],
    n_blocks_all = n_blocks_all[1],
    seg_reads_all = seg_reads_all[1],
    n_blocks = n(),
    seg_snps = sum(n_snps),
    seg_reads = sum(total),
    seg_snps_filtered = sum(n_snps[snp_filter]),
    seg_reads_filtered = sum(total[snp_filter]),
    reads_frac = seg_reads/seg_reads_all,
    priority=priority[1],
    
    # 250508: Recalculate from actual counts instead of means
    # A: Filtered versions first
    f_BAF_avg = mean(BAF[snp_filter], na.rm = TRUE),
    f_mBAF_avg = mean(pmin(BAF[snp_filter], 1-BAF[snp_filter]), na.rm = TRUE),
    f_A = sum(A[snp_filter]), 
    f_B = sum(B[snp_filter]),
    f_BAF = f_B / (f_A+f_B),
    f_mBAF = pmin(f_BAF, 1-f_BAF),
    
    f_cBAF_avg = mean(cBAF[snp_filter], na.rm=TRUE), 
    f_cmBAF_avg = mean(pmin(cBAF[snp_filter], 1-cBAF[snp_filter]), na.rm=TRUE),
    f_cA = sum(cA[snp_filter]), 
    f_cB = sum(cB[snp_filter]),
    f_cBAF = f_cB / (f_cA + f_cB),
    f_cmBAF = pmin(f_cBAF, 1-f_cBAF),    
    
    # B: Unfiltered versions
    BAF_avg = mean(BAF, na.rm = TRUE),
    mBAF_avg = mean(pmin(BAF, 1-BAF), na.rm = TRUE),
    A = sum(A), 
    B = sum(B), 
    BAF = B / (A+B),
    mBAF = pmin(BAF, 1-BAF),
    
    cBAF_avg = mean(cBAF, na.rm=TRUE), 
    cmBAF_avg = mean(pmin(cBAF, 1-cBAF), na.rm=TRUE),
    cA = sum(cA), 
    cB = sum(cB),
    cBAF = cB / (cA + cB),
    cmBAF = pmin(cBAF, 1-cBAF),
    
    .groups = 'drop'
  )

plot_baf_comp <- function(x, var1, var2){
  x %>% 
    ggplot(aes(!!sym(var1), !!sym(var2), color=factor(cn))) + 
    geom_abline(slope=1, linetype="dotted", color="grey30") +
    geom_point() + facet_wrap(~clone, nrow=1) +
    theme_dntr(axis_ticks = T, legend = "right") + 
    scale_x_continuous(breaks=c(0,0.5,1), limits = function(x) c(0, ifelse(max(x) <= 0.5, 0.5, 1))) +
    scale_y_continuous(breaks=c(0,0.5,1), limits = function(y) c(0, ifelse(max(y) <= 0.5, 0.5, 1)))
}

plot_baf_comp(switched.segs, "BAF", "BAF_avg") / 
  plot_baf_comp(switched.segs, "cBAF", "cBAF_avg") / 
  plot_baf_comp(switched.segs, "mBAF", "mBAF_avg") / 
  plot_baf_comp(switched.segs, "cmBAF", "cmBAF_avg")

plot_baf_comp(switched.segs, "f_BAF", "f_BAF_avg") / 
  plot_baf_comp(switched.segs, "f_cBAF", "f_cBAF_avg") / 
  plot_baf_comp(switched.segs, "f_mBAF", "f_mBAF_avg") / 
  plot_baf_comp(switched.segs, "f_cmBAF", "f_cmBAF_avg")



### Output plots
# 1. Plot chromosome level (or zoom to region)
plot_chr <- "chr3"
seg_breaks <- blocks.segs[[plot_chr]] %>% 
  group_by(seg_id) %>% 
  summarize(start_block=min(block_id), 
            end_block=max(block_id),
            start_pos=min(start_pos),
            end_pos=max(end_pos))
seg_cn <- switched.segs %>% 
  filter(chr == plot_chr) %>% 
  mutate(clone_ordered=factor(clone, clone_order)) %>% 
  left_join(select(seg_breaks, seg_id, start_block, end_block), by="seg_id")


s1 <- ba.sl[[plot_chr]] %>%
  # filter(block_id %in% f) %>%
  mutate(clone_ordered = factor(clone, levels=clone_order)) %>% 
  # ggplot(aes(block_id, cBAF, size=n_snps, color=switch)) +
  # ggplot(aes(block_id, cBAF, color=switch)) +
  ggplot(aes(start_pos, cBAF)) +
  geom_vline(xintercept=union(seg_breaks$start_pos-1, seg_breaks$end_pos), 
             color="red") +
  geom_hline(yintercept=c(1/3,0.5,2/3), linetype="solid", color="grey", alpha=0.5) +
  geom_point(alpha=.5, size=.5) + 
  geom_point(data=. %>% filter(switch), pch=21, alpha=.5, size=.7, color="orange") + 
  geom_segment(data=seg_cn, aes(y=cBAF, yend=cBAF, x=start_pos, xend=end_pos, color=factor(cn)), 
               linewidth=1.5, alpha=0.7, inherit.aes = F) +
  # geom_line(aes(y=BAF_runmean), color="orange", size=1) +
  # geom_line(aes(y=deviation_runmean), color="orange", linewidth=1) +
  # scale_color_manual(values=c("TRUE"="red","FALSE"="black", "NA"="grey")) +
  scale_color_manual(values=c("1"="blue","2"="black","3"="red","4"="darkred", "5"="darkred"), "") +
  facet_wrap(~ clone_ordered, labeller = labeller(clone_ordered = clone_labels), ncol=1) +
  theme_dntr(axis_ticks = T, legend = "none") + 
  scale_y_continuous(expand=c(0.05, 0.05), limits=c(0,1)) + labs(title=paste0(plot_chr, ", haplotyped/fixed, no snp filter"), breaks=pretty_breaks())
s1
# Show next to cn plot of same chr
plot_clone_detail(d, region="chr12p", cn.trunc = 8) + s1

# 2. Plot summary cn x correctedBAF, per clone. Each point is a segment
switched.segs$clone_ordered <- factor(switched.segs$clone, levels = clone_order)
s2 <- switched.segs %>% 
  # ggplot(aes(abs(cBAF-0.5), cn, col=factor(priority), shape=seg_snps >= min_snp_segment)) +
  ggplot(aes(cBAF, cn, col=factor(priority), shape=seg_snps >= min_snp_segment)) +
  # geom_vline(xintercept=seq(0, 0.5, by = 0.1), linetype="dotted", color="grey") +
  geom_vline(xintercept=c(1/3,1/4), linetype="dotted", color="grey") +
  geom_vline(xintercept=c(1/2), linetype="solid", color="grey") +
  geom_jitter(alpha=0.7, height=0.2, width=0, size=3) +
  scale_shape_manual(values=c("TRUE"=19, "FALSE"=13)) +
  scale_color_manual(values=c("1"="red","2"="orange","3"="darkgreen", "4"="grey", "black"="black")) +
  theme_dntr(axis_ticks = T) + 
  scale_y_continuous(expand=c(0.1, 0.1), breaks=unique(switched.segs$cn)) +
  facet_wrap(~ clone_ordered, labeller = labeller(clone_ordered = clone_labels), ncol=1) +
  theme(legend.position = "right", legend.direction = "vertical") + 
  # labs(x="phase-corrected, mirrored BAF [min(BAF, 1-BAF)]")
  labs(x="phase-corrected BAF")
# labs(x="phase-corrected |BAF-0.5|")
s2

s3 <- switched.segs %>% 
  # ggplot(aes(abs(f_cBAF-0.5), cn, col=factor(priority), shape=seg_snps >= min_snp_segment)) +
  ggplot(aes(f_cmBAF, cn, col=factor(priority), shape=seg_snps >= min_snp_segment)) +
  # geom_vline(xintercept=seq(0, 0.5, by = 0.1), linetype="dotted", color="grey") +
  geom_vline(xintercept=c(1/2), linetype="solid", color="grey") +
  geom_vline(xintercept=c(1/3,1/4), linetype="dotted", color="grey") +
  geom_jitter(alpha=0.7, height=0.2, width=0, size=3) +
  scale_color_manual(values=c("1"="red","2"="orange","3"="darkgreen", "4"="grey", "black"="black")) +
  scale_shape_manual(values=c("TRUE"=19, "FALSE"=13)) +
  theme_dntr(axis_ticks = T) + 
  scale_y_continuous(expand=c(0.1, 0.1), breaks=unique(switched.segs$cn)) +
  facet_wrap(~ clone_ordered, labeller = labeller(clone_ordered = clone_labels), ncol=5) +
  theme(legend.position = "right", legend.direction = "vertical") + 
  labs(x="phase-corrected BAF", caption=paste0("filtered for >=",min_snp_block," SNPs/block"))
# labs(x="phase-corrected |BAF-0.5|", caption=paste0("filtered for >=",min_snp_block," SNPs/block"))

pdf(snakemake@output[["cmBAF"]], width=15, height=3)
s3
dev.off()

#Probability based genotype calls
source("workflow/scripts/model_gt_wplots.R")
#Tunable filters 
error_rate <- 0.02
error_multiplier <- 5
min_gt_prob <- 0.97 # or something (max 1)


switched.segs$chr <- factor(switched.segs$chr, levels=levels(segs$chr))
ss.l <- split(switched.segs, switched.segs$clone)
baf_segs.l <- lapply(ss.l, function(x) {
  cl <- x$clone[1]
  x %>% 
    full_join(segs, by = c("seg_id" = "seg_idx", "chr")) %>% 
    arrange(seg_id) %>% 
    select(chr, arm, start.pos, end.pos, n.probes, everything(), -clone_ordered) %>% 
    mutate(
      clone = ifelse(is.na(clone), cl, clone), # Segments without SNP blocks will have clone NA
      cn = ifelse(is.na(cn), cn_clones[seg_id, cl], cn),
      filter_seg = seg_snps < min_snp_segment | is.na(seg_snps),
      BAF_final = case_when(
        cn == 0 | chr == "chrY" ~ 0, 
        filter_seg & !is.na(cn) & is.na(BAF) ~ floor(cn/2)/cn, 
        # chr == "chrX" & !is.na(cn) & is.na(BAF) ~ floor(cn/2)/cn, 
        # cn > 0 & filter_seg ~ floor(cn/2)/cn,
        TRUE ~ cBAF
      ),
      # Default allele-specific copy numbers
      cn_a_orig = round((1 - BAF_final) * cn),
      cn_b_orig = cn - cn_a_orig,
      pass = !is.na(cn) & cn > 0 & !filter_seg & !is.na(cA) & !is.na(cB)
    ) %>% 
    mutate(
      gt_test = pmap(list(cA, cB, cn, pass), function(A, B, cn, pass) {
        if(pass) calc_gt_prob(A, B, cn, err = error_rate, dynamic_error = TRUE, max_multiplier = error_multiplier)
      }),
      gt_res = map(gt_test, ~ if(!is.null(.x)) get_gt(.x) else NULL),
      gt = map_chr(gt_res, ~ if(!is.null(.x)) .x$best_genotype else NA),
      gt_prob = map_dbl(gt_res, ~ if(is.null(.x)) NA_real_ else .x$probability),
      cn_a = map2_int(gt, gt_prob, function(gt, prob) if(is.na(gt) | prob < min_gt_prob) NA else stringr::str_count(gt, "A")),
      cn_b = map2_int(gt, gt_prob, function(gt, prob) if(is.na(gt) | prob < min_gt_prob) NA else stringr::str_count(gt, "B"))
    )
})

bafs <- bind_rows(baf_segs.l)

### Plotting cn_a/cn_b
plot_data <- bafs %>%
  mutate(
    cn_a_plot = cn_a + 0.1,
    cn_a_orig_plot = cn_a_orig + 0.1,
    clone_ordered = factor(clone, levels=clone_order)
  ) %>%
  group_by(clone, chr) %>%
  mutate(
    order_index = order(start.pos),
    seg_number = seq_along(start.pos)
  ) %>%
  ungroup()

# Create the plot
plot_data %>% 
  # filter(chr=="chr11") %>% 
  filter(chr %in% c("chr1", "chr3", "chr4", "chr10", "chr11","chr12", "chr14", "chr15", "chr16", "chr18", "chr20", "chr22")) %>% 
  ggplot() +
  geom_segment(aes(x = start.pos, xend = end.pos, 
                   y = cn_a_plot, yend = cn_a_plot),
               color = "gray40", size = 2) +
  geom_segment(aes(x = start.pos, xend = end.pos, 
                   y = cn_b, yend = cn_b),
               color = "blue", size = 2) +
  geom_segment(aes(x = start.pos, xend = end.pos, 
                   y = cn_a_orig_plot, yend = cn_a_orig_plot),
               color = "lightgray", size = 1.5, alpha = 0.8) +
  geom_segment(aes(x = start.pos, xend = end.pos, 
                   y = cn_b_orig, yend = cn_b_orig),
               color = "lightblue", size = 1.5, alpha = 0.8) +
  facet_grid(clone_ordered ~ chr, scales = "free_x", space = "free_x") +
  theme_dntr(axis_ticks=T) + 
  # theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) + 
  labs(
    x = "position (Mb)",
    y = "copy number",
  ) +
  scale_y_continuous(limits = c(-0.1, NA), breaks=function(y) c(0:max(y))) +
  scale_x_continuous(breaks=pretty_breaks(4), labels=function(x) x/1e6, name = "position (Mb)")


#Fill in information for segments with haplotyped copynumbers, and segments without haplotyped copynumbers
baf_segs.l <- map(baf_segs.l, function(df) {
  df$cn_ab <- map2_chr(df$cn_a, df$cn_b, function(a, b) {
    if (is.na(a) || is.na(b)) {
      NA_character_
    } else {
      paste0(strrep("A", a), strrep("B", b))
    }
  })
  df
})

baf_segs.l <- map(baf_segs.l, function(df) {
  df$cn_ab_noNA <- map2_dbl(df$cn_a_orig, df$cn_b_orig, function(a, b) {
    if (is.na(a) || is.na(b)) {
      NA_real_  
    } else {
      a + b  
    }
  })
  df$final_ab <- ifelse(
    is.na(df$cn_ab), 
    df$cn_ab_noNA,
    df$cn_ab  
  )
  return(df) 
})

#Add to clone object 

d$clones$final_ab <- lapply(baf_segs.l[d$clones$clone_id], function(x) x$final_ab[rep(1:nrow(x), x$n.probes)])
d$clones$baf <- lapply(baf_segs.l[d$clones$clone_id], function(x) x$BAF_final[rep(1:nrow(x), x$n.probes)])
d$clones$cn_a <- lapply(baf_segs.l[d$clones$clone_id], function(x) x$cn_a[rep(1:nrow(x), x$n.probes)])
d$clones$cn_b <- lapply(baf_segs.l[d$clones$clone_id], function(x) x$cn_b[rep(1:nrow(x), x$n.probes)])


pdf(snakemake@output[["heatmap"]], width=20, height=10)
o <- intersect(clone_order, d$clones$clone_id)
plot_clone_heatmap(d, show_chr=T, order=o)
plot_clone_heatmap(d, cn_slot="baf", show_chr=T, clone_types=NULL, order=o)
plot_clone_heatmap(d, cn_slot="cn_a", show_chr=T, clone_types=NULL, highlight_dups = F, order=o)
plot_clone_heatmap(d, cn_slot="cn_b", show_chr=T, clone_types=NULL, highlight_dups = F, order=o)
plot_clone_heatmap(d, cn_slot="final_ab", show_chr=T)
dev.off()

write_rds(d, snakemake@output[["cn_obj_out"]])


#Prepare medicc object 

bafs<-do.call(rbind,baf_segs.l)
bafs<-bafs%>%select(chr, clone, start.pos, end.pos, cn_a, cn_b, final_ab)

bafs <- bafs %>%
  mutate(
    to_fix = is.na(cn_a) & is.na(cn_b) & !is.na(final_ab),
    cn_a = ifelse(to_fix,
                  case_when(
                    final_ab == "1" ~ 1L,
                    final_ab == "0" ~ 0L,
                    final_ab == "2" ~ 1L,
                    final_ab == "3" ~ 2L,
                    TRUE ~ NA_integer_
                  ),
                  cn_a),
    cn_b = ifelse(to_fix,
                  case_when(
                    final_ab == "1" ~ 0L,
                    final_ab == "0" ~ 0L,
                    final_ab == "2" ~ 1L,
                    final_ab == "3" ~ 1L,
                    TRUE ~ NA_integer_
                  ),
                  cn_b)
  ) %>%
  select(-to_fix)

bafs

medicc <- bafs %>% 
  select(clone, chrom = chr, start = start.pos, end = end.pos, cn_a, cn_b) %>% 
  mutate(chrom = factor(chrom, levels = paste0("chr", c(1:22, "X", "Y")))) %>% 
  arrange(chrom, start, end) %>%
  group_by(start) %>% 
  filter(!any(is.na(cn_a) | is.na(cn_b))) %>%  # Remove groups with any NA in cn_a or cn_b
  ungroup()

colnames(medicc)<-gsub("clone", "sample_id", colnames(medicc))

write_tsv(medicc, snakemake@output[["medicc"]])


