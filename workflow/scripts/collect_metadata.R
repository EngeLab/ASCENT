library(tidyverse)
options(scipen=999)

dna_to_cellid <- function(x){ sub("([A-Z]{3}[0-9]{5})(D?)_([A-Z][0-9]{2}).*","\\1.\\3", x) }
rna_to_cellid <- function(x){ sub("([A-Z]{3}[0-9]{5})([DR]?).([A-Z][0-9]{2}).*","\\1.\\3", x) }
cellid_to_plateid <- function(x){ sub("([A-Z]{3}[0-9]{5}).*","\\1", x) }
cellid_to_well <- function(x){ sub("([A-Z]{3}[0-9]{5}).([A-Z][0-9]{2})","\\2", x) }
cellid_to_row <- function(x){ sub("([A-Z]{3}[0-9]{5}).([A-Z])[0-9]{2}","\\2", x) }
cellid_to_col <- function(x){ sub(".*[A-Z]{1}([0-9]{2})$","\\1", x) }

d.input <- snakemake@input[["dnafile"]]
r.input <- snakemake@input[["rnafile"]]
outfile <- snakemake@output[["meta_per_cell"]]
outfile_sort <- snakemake@output[["meta_wsort"]]
outfile_facs <- snakemake@output[["meta_wfacs"]]

if(!is.null(d.input)){
  d <- read_tsv(d.input, col_names = c("dna_library_id","dna_fq1","dna_fq2"), col_types="ccc", skip=1) %>% 
    filter(dna_library_id %in% snakemake@params[["cells_dna"]]) %>%
    mutate(cell_id=dna_to_cellid(dna_library_id),
           dna_bam=paste0(snakemake@config[["path"]][["dnaout"]],"/",dna_library_id,"/",
                          dna_library_id,".dedup.bam")) %>% 
    group_by(cell_id) %>% 
    summarize(across(everything(), ~paste(unique(.x),collapse=";"))) %>% # Combine multiple sequencing runs to one row per cell_id (dna_fq1 ; delim)
    select("cell_id",everything())
} else { d <- tibble(cell_id = character(), dna_library_id = character()) }
if(!is.null(r.input)){
  r <- read_tsv(r.input, col_names = c("rna_library_id","rna_fq1","rna_fq2"), col_types="ccc", skip=1) %>% 
    filter(rna_library_id %in% snakemake@params[["cells_rna"]]) %>%
    mutate(cell_id=rna_to_cellid(rna_library_id),
           rna_bam=paste0(snakemake@config[["path"]][["rnaout"]],"/",rna_library_id,"/",
                          rna_library_id,".sorted.bam"),
           rna_library_id=sub("-",".",rna_library_id)) %>% 
    group_by(cell_id) %>% 
    summarize(across(everything(), ~paste(unique(.x),collapse=";"))) %>% # Combine multiple sequencing runs to one row per cell_id (dna_fq1 ; delim)
    select("cell_id",everything())
} else { r <- tibble(cell_id = character(),rna_library_id = character()) }

# Load metadata
expand_cols <- function(x,sep="-"){
  if(is.na(x)) return(1:24)
  xs <- unlist(strsplit(x,sep))
  if(length(xs)==1) return(as.integer(xs))
  return(seq(xs[1],xs[2]))
}

if(!is.null(snakemake@input[["metafile"]])){
  m <- read_tsv(snakemake@input[["metafile"]],show_col_types=F)
  m$col <- lapply(m$subset_cols, expand_cols)
  ml <- m %>% select(-subset_cols) %>% unnest(col) %>% mutate(col=sprintf("%02d",col)) 
} else {
  cat("No metadata file provided, using plate prefixes as patient_id\n")
  ml <- tibble(
    plate_id=character(),
    patient_id=character(),
    col=character()
  )
}

j <- full_join(r, d, by="cell_id") %>% 
  mutate(plate_id=cellid_to_plateid(cell_id),
         col=cellid_to_col(cell_id),
         row=cellid_to_row(cell_id)) %>% 
  left_join(ml, by=c("plate_id","col")) %>% 
  select("cell_id","plate_id", "patient_id",everything())

if(sum(is.na(j$patient_id))>0){
  warning(paste0("Several cells without patient_id (n=",sum(is.na(j$patient_id)),")\nFirst/last five:\n"),
          paste(head(j$cell_id[is.na(j$patient_id)]),collapse=" "),"\n",
          paste(tail(j$cell_id[is.na(j$patient_id)]),collapse=" "))
  j <- j %>% mutate(patient_id = ifelse(is.na(patient_id), str_sub(cell_id, 1, 6), patient_id))
}

# FACS data
sr.l <- list.files("facs_data/sort_result", paste0(unique(j$plate_id),collapse="|"), full.names=T)
if(length(sr.l)>=1){
  names(sr.l) <- substr(basename(sr.l),1,8)
  sr <- map_dfr(sr.l, ~read_delim(.x, delim=",", show_col_types=F),.id="plate_id") %>% 
    mutate(index=paste0(str_sub(`Well Number`,1,1), str_pad(str_sub(`Well Number`, 2,-1),2,"left",0)), # Pad well numbers (H1 --> H01)
           cell_id = paste0(plate_id, ".", index)) %>% 
    select(cell_id, everything(), -plate_id, -`Well Number`, -index) %>% 
    filter(`Sorted Count`>0) %>% 
    mutate(sort_gate_short=case_when(grepl("Non-blasts|Non-leukemia|non blasts|non-blasts|b-cells|t-cells", `Sort Gate`) ~ "normal",
                                     grepl("true blasts|blast gate|blasts|Blasts|Mostly blasts|Generic blast gate", `Sort Gate`) ~ "blasts",
                                     grepl("Singlet|True singlets", `Sort Gate`) ~ "singlets",
                                     TRUE ~ as.character(`Sort Gate`)))
}

is.l <- list.files("facs_data/index_sort", paste0(unique(j$plate_id),collapse="|"), full.names=T)
if(length(is.l)>=1){
  names(is.l) <- substr(basename(is.l),1,8)
  is.ll <- lapply(is.l, function(x){ 
    read.csv(x, stringsAsFactors=F) %>% filter(!is.na(Index)) %>% 
      rename_with(~paste0(substr(.x,1,4),"FACS"), matches("CD[0-9]+"))
  })
  is <- bind_rows(is.ll, .id="plate_id") %>% 
    filter(!is.na(Index)) %>% 
    mutate(index=paste0(str_sub(Index,1,1), str_pad(str_sub(Index, 2,-1),2,"left",0)), # Pad well numbers (H1 --> H01)
           cell_id = paste0(plate_id, ".", index)) %>% 
    select(cell_id, everything(), -Index, -index, -plate_id) %>% 
    filter(! cell_id %in% cell_id[duplicated(cell_id)]) # Drop 5-cell wells b/c multiple FACS lines per cell
}

if(exists("sr")){
  js <- j %>% 
    left_join(sr, by="cell_id")
  jf <- j %>% 
    left_join(sr, by="cell_id") %>% 
    left_join(is, by="cell_id")
} else { js=j; jf=j }

write_tsv(j, outfile)
write_tsv(js, outfile_sort)
write_tsv(jf, outfile_facs)