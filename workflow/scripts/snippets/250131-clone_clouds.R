library(RColorBrewer)
library(ggrepel)
library(parallel)

# pid <- "ALL4"
# dat <- dl[[pid]][,dl[[pid]]$dna_type %in% "aneuploid"]
# dat$clone_time <- paste(dat$clone, dat$timepoint, sep="_")
# DimPlot(dat, group.by = c("clone","timepoint"))
# DimPlot(dat, group.by = c("clone_time"))
# 
# plot_cloud(dat)
# plot_cloud(dat, clone_slot="clone_time")

plot_cloud <- function(dat, clone_slot="clone", sampsize=5, numsamp=100, threads=20){
  dta <- as.matrix(dat@assays$RNA@counts)
  dta.log <- as.matrix(dat@assays$RNA@data)

  dim(dta.log)
  
  clones <- dat@meta.data[[clone_slot]]
  clone_sizes <- table(clones)
  print(clone_sizes)
  valid_clones <- names(clone_sizes)[clone_sizes >= sampsize]
  if(length(setdiff(unique(clones), valid_clones))>0) {
    cat("Excluded (low cell count):", setdiff(unique(clones), valid_clones),"\n")
    clones <- clones[clones %in% valid_clones]
    dta.log <- dta.log[,which(clones %in% valid_clones)]
  }
  length(clones) == ncol(dta.log)
  print(dim(dta.log))
  
  clone.aggrs <- apply(2^dta.log, 1, function(x) {
    sapply(split(x, clones), mean)
  })
  
  if(length(valid_clones)>1){
    mxs <- apply(clone.aggrs, 2, max)  
  } else {
    # Single clone
    mxs <- max(clone.aggrs)
  }
  
  mean(mxs > 10)
  o <- order(mxs, decreasing=T)[1:5000]
  o.n <- names(sort(mxs, decreasing=T)[1:2000])
  
  dta.filt <- 2^dta.log

    # Run 1
  cat(sprintf("Sampling expression for %s:\n %d clones, sample size %d, %d iterations", unique(dat$patient_id), length(unique(clones)), sampsize, numsamp),"\n")
  samp.clone.aggrs.list <- mclapply(1:nrow(dta.filt), function(i) {
    x <- dta.filt[i,]
    result <- sapply(split(x, clones), function(x) {
      sapply(1:numsamp, function(y) {
        mean(sample(x, sampsize))
      })
    })
    # Preserve the exact same structure as the apply version
    colnames(result) <- names(split(x, clones))
    rownames(result) <- 1:numsamp
    result
  }, mc.cores = threads)
  names(samp.clone.aggrs.list) <- rownames(dta.filt)
  
  # Convert to matrix preserving names
  samp.clone.aggrs <- do.call(rbind, lapply(1:length(valid_clones), function(i) {
    do.call(cbind, lapply(samp.clone.aggrs.list, function(x) x[,i]))
  }))

  if(length(valid_clones)>1){
    rownames(samp.clone.aggrs) <- paste(rep(rownames(clone.aggrs), each=numsamp), 1:numsamp, sep='.')
    tmp <- clone.aggrs
    rownames(tmp) <- paste0(rownames(clone.aggrs), ".0")
  } else {
    rownames(samp.clone.aggrs) <- paste(rep(valid_clones, each=numsamp), 1:numsamp, sep='.')
    tmp <- t(as.matrix(clone.aggrs))
    rownames(tmp) <- paste0(valid_clones, ".0")
  }
  samp.clone.aggrs <- rbind(tmp, samp.clone.aggrs)
  rm(tmp)
  
  dim(samp.clone.aggrs)
  # 303 x 27173. (n clones x samples) + (n clones x 1) x n genes
  
  mxs <- apply(samp.clone.aggrs, 2, max)
  head(mxs)
  
  mean(mxs > 10)
  o <- order(mxs, decreasing=T)[1:5000] #Take the top 5000 expressed ish 
  
  cl.dst <- cor(t(samp.clone.aggrs[,unique(o, dat@assays$RNA@var.features)]), method="spearman") 
  class(cl.dst)
  
  cl.csc <- cmdscale(as.dist(1-cl.dst), 2) # For correlation
  labels <- as.factor(gsub(".+_", "", sapply(strsplit(rownames(cl.csc), "\\."), function(x) {x[1]}))) # Clone
  labels <- as.factor(sapply(strsplit(rownames(cl.csc), "\\."), function(x) {x[1]})) # Clone+patient
  
  # Ridge clouds
  colnames(cl.csc) <- c('x', 'y')
  cl.2 <- cl.csc %>% as.data.frame %>% mutate(labels = labels)
  
  label_distances <- cl.2 %>%
    group_by(labels) %>%
    summarise(
      center_x = median(x),
      center_y = median(y),
      max_distance = max(sqrt((x - median(x))^2 + (y - median(y))^2))  # Euclidean distance
    )
  
  # Print the result
  print(label_distances[order(label_distances$max_distance, decreasing=TRUE),])
  mids <- cl.2 %>% as.data.frame %>% 
    mutate(labels = labels) %>% 
    group_by(labels) %>% 
    summarize(mx=median(x), my=median(y))
  
  scatterplot <- cl.2 %>% 
    ggplot(aes(x, y, color = labels)) +
    # geom_density_2d(binwidth=NULL) +
    # Try density plot, if fails just continue without it
    {
      tryCatch(
        geom_density_2d(binwidth=NULL),
        error = function(e) {
          message("Could not compute density for some clones, showing only centers")
          geom_blank()
        }
      )
    } +
    theme_classic() +
    theme(axis.text = element_blank(), axis.title = element_blank()) + 
    geom_point(data=mids, aes(x=mx, y=my), fill='black', pch=21) + 
    geom_label_repel(data=mids, aes(x=mx, y=my, label=labels), max.overlaps = 20, box.padding=0.5) +
    theme_dntr(axis_labels = F, legend="none", font_size = 10) + 
    labs(title=unique(dat$patient_id), caption=paste0("slot: ", clone_slot, " (sampled ",numsamp," x ",sampsize," cells per clone)"))
  scatterplot
  
  return(scatterplot) 
}
