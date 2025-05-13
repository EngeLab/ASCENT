library(tidyverse)
options(scipen=999)

all <- lapply(res, function(chr){
  bind_rows(chr, .id = "clone") %>%
    # select(block_id, chr, start, end, BAF, clone)
    select(block_id, chr, start, end, A, B, total, BAF, clone)
}) %>% bind_rows

bin_width <- 5e6  
baf_mtx <- all %>%
  mutate(
    start = floor(start / bin_width) * bin_width,
    end = start + bin_width
  ) %>%
  group_by(chr, start, end, clone) %>%
  # Opt 1. Take mean BAF per bin per clone
  # Opt 2. Re-calculate BAF from sum A + B of larger bin
  summarize(
    # BAF = mean(BAF, na.rm = TRUE),
    A = sum(A,na.rm=T),
    B = sum(B,na.rm=T),
    total = A + B,
    BAF = B/total,
    .groups = 'drop'
  ) %>%
  mutate(
    mBAF = pmin(BAF, 1-BAF)
  ) %>% 
  pivot_wider(
    id_cols = c(chr, start, end),
    names_from = clone,
    values_from = mBAF
  ) %>%
  arrange(chr, start)

# Check result
head(baf_mtx)

# Get arm annotations to same bin size
arms <- d$bins$all %>%
  mutate(
    bin_start = floor(start / bin_width) * bin_width,
    bin_end = bin_start + bin_width
  ) %>%
  # Count how many 20kb segments of each arm are in each 5Mb bin
  group_by(chr, bin_start, bin_end, arm) %>%
  summarize(
    count = n(),
    .groups = 'drop'
  ) %>%
  # For bins that span arms, keep the one with more 20kb segments
  group_by(chr, bin_start, bin_end) %>%
  slice_max(count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(chr, start=bin_start, end=bin_end, arm)

# Combine
baf_mpcf <- baf_mtx %>%
  left_join(arms, by = c("chr", "start", "end")) %>% 
  arrange(chr, start)

baf_mpcf.input <- data.frame(
  chr = as.integer(sub("chr", "", baf_mpcf$chr)),
  pos = baf_mpcf$start,
  baf_mpcf[,grep("[0-9]_.*",colnames(baf_mpcf))],
  check.names = F
)
baf_mpcf.input[is.na(baf_mpcf.input)] <- 0

# Drop single-bin arms
arm_counts <- table(paste(baf_mpcf$chr, baf_mpcf$arm))
single_bins <- names(which(arm_counts == 1))
rows_to_remove <- which(paste(baf_mpcf$chr, baf_mpcf$arm) %in% single_bins)

res_mpcf <- copynumber::multipcf(
  data = baf_mpcf.input[-rows_to_remove,],
  # Y = mpcf.inputY,
  pos.unit = "bp",
  arms = baf_mpcf$arm[-rows_to_remove],
  gamma = 0.1,
  normalize = FALSE,
  fast = FALSE,
  verbose = TRUE,
  return.est = TRUE,
  # w = weights
)

copynumber::plotHeatmap(res_mpcf,upper.lim = 0.5, lower.lim = -0.01)
