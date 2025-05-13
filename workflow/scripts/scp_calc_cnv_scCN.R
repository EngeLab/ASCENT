suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
library(parallel)
options(scipen=999)
#Inputs: 
mtx <- as.matrix(read_tsv(snakemake@input[["cn"]]))
# ss <- colnames(mtx)


bins <- read_tsv(snakemake@input[["bins"]], skip=1, col_names = c("chr","start","end","idx","bin")) %>%
  mutate(
    chr_short=sub("chr","",chr),
    chr_int=as.integer(case_when(chr_short == "X" ~ "23", chr_short == "Y" ~ "24", TRUE ~ chr_short)),
    chr=factor(chr, levels=unique(chr))
  )

bins.all <- read_tsv(snakemake@input[["all_bins"]],col_names=c("chr","start","end")) %>% 
  mutate(idx=seq_along(chr),
         bin=seq_along(chr),
         chr_short=sub("chr","",chr),
         chr_int=as.integer(case_when(chr_short == "X" ~ "23", chr_short == "Y" ~ "24", TRUE ~ chr_short)),
         chr=factor(chr, levels=unique(chr)))

bins.bed<- import(snakemake@input[["all_bins"]])
segments <- read_tsv(snakemake@input[["sc_segments"]])
frag_files <- setNames(snakemake@input[["frag_files"]], nm=snakemake@params[["cells"]])


#Calculation 
#1: Expand copynumbers to all bins 
mtx.e <- matrix(
  NA,
  nrow = max(bins.all$idx),
  ncol = ncol(mtx),
  dimnames = list(
    1:max(bins.all$idx),
    colnames(mtx)
  )
)

#Need to mask out short segments! 
mtx_s<-mtx[segments$start,]
mtx_s[segments$size < snakemake@params[["min_binsize"]], ] <- NA
mtx_s[1:5,1:5]
# split(clones$dna_library_id, clones$clone_final)

# Calculate median cn per segment per clone

cn_clone.segs<-mtx_s

# Return to bin-level and expand to include all bins
cn_clone.bins <- cn_clone.segs[rep(1:nrow(cn_clone.segs), segments$size), ]
mtx.e[bins$idx, ] <- cn_clone.bins

# Ensure same order
frag_files.o <- frag_files[colnames(mtx.e)]

single_cell_SCP <- function(cn=NULL, frag_file=NULL, iter=20, frac=0.005, bins=bins.bed) {
  cn_states <- table(cn)
  cat("^^^Running cell: ", names(frag_file), "\n")
  # gr <- granges(import(frag_file)) # granges call is to remove metadata
  gr <- import(frag_file) # granges call is to remove metadata
  #resize(ga_b, width=width(ga_b) - 10, fix="center")    
  # ga_bf <- ga_b[width(ga_b) > 50] #Maybe change this as well... 
  q <- quantile(width(gr), c(0.05, 0.95))
  gr.f <- gr[width(gr) >=  q[1] & width(gr) <= q[2]]
  gr.r <- resize(gr.f, width=width(gr.f) - 10, fix="center")
  
  #Here it removes when there are only few bins which exist at that ploidy 
  cn_states.f <- cn_states[cn_states > 1000]
  
  #Here we check if total cnstates bins are more than 60k bins 
  
  res <- lapply(as.numeric(names(cn_states.f)), function(state) {
    print(paste(names(frag_file), "cn state:", state))
    #bins.pl <- bins[!is.na(cn) & cn == state & !bins@seqnames %in% c("chrX","chrY")]
    bins.pl <- bins[!is.na(cn) & cn == state & !bins@seqnames %in% c("chrY")]
    gr.pl <- subsetByOverlaps(gr.r, bins.pl, type='within') # Get subset in pl bins
    totBinSize <- sum(width(bins.pl))
    meanFragsize <- mean(width(gr.pl))
    numReads <- length(gr.pl)
    if(numReads==0){
      return(NA)
    }
    f <- frac*totBinSize/meanFragsize # Means to aim for 0.5% coverage
    if(f > numReads) {
      print(paste0("Ploidy: ", state))
      print(paste0("meanFragsize: ", round(meanFragsize,1), ", totBinSize: ", round(totBinSize/1e6), " Mb"))
      print(paste0("Samplesize: ", round(f), ", Readnum: ", numReads))
      print("")
      return(NA)
    }
    gr.pl.cov <- sapply(1:iter, function(i) {
      cov <- coverage(gr.pl[sample(1:length(gr.pl), f, replace=FALSE)])
      sl  <- IRanges::slice(cov, lower=2)
      lns <- sum(unlist(lapply(sl, length)), na.rm=T)
      formn <- unlist(lapply(sl, width))
      mns <- 0
      if(length(formn) > 0) {
        mns <- mean(formn, na.rm=T)
      }        
      #        setNames(sapply(1:4, function(k){ sum(sum(cov==k))/totBinSize }), nm=1:4)
      setNames(c(sapply(1:4, function(k){ sum(sum(cov==k))/totBinSize }), lns, mns), nm=c(1:4, 'numOverlaps', 'meanOverlapSize'))
    })
    # gr.pl.cov
    setNames(c(rowMeans(gr.pl.cov, na.rm=T), meanFragsize, numReads, f), nm=c(rownames(gr.pl.cov), "meanFragSize", "numReads", "sampleSize")) # Keep only mean
  })
  names(res) <- names(cn_states.f)
  res
}

cell.scps <- mclapply(colnames(mtx.e), function(i){
  single_cell_SCP(cn=mtx.e[,i], frag_file=frag_files.o[i], bins=bins.bed, frac=0.005, iter=20)
},mc.cores = snakemake@threads)
names(cell.scps) <- colnames(mtx.e)


result_df <- enframe(unlist(cell.scps)) %>%
  separate(name, into = c("dna_library_id", "ploidy", "depth_num"), sep = "\\.") %>%
  filter(!is.na(value)) %>%
  group_by(dna_library_id, ploidy) %>%
  mutate(depth_num = paste0("depth", row_number()),
         ploidy=as.integer(ploidy)) %>%
  pivot_wider(
    names_from = depth_num,
    values_from = value
  ) %>%
  ungroup()

allCell.ploidies <- mclapply(colnames(mtx.e), function(i) {
  table(mtx.e[,i])
},mc.cores = snakemake@threads)
names(allCell.ploidies) <- colnames(mtx.e)

df <- map_dfr(names(allCell.ploidies), ~ tibble(
  dna_library_id = .x,
  ploidy = as.integer(names(allCell.ploidies[[.x]])),
  nr_bins = as.numeric(allCell.ploidies[[.x]])
))

final_df <- left_join(result_df, df)
colnames(final_df) <- c(colnames(final_df)[1:6], "numOverlaps", "meanOverlapSize",  "meanFragSize", "numReads", "sampleSize", "nr_bins")

write_tsv(final_df, snakemake@output[["scps"]])
