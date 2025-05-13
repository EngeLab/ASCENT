options(scipen=999)
suppressPackageStartupMessages(library(tidyverse))
library(cowplot)

# Functions and plot theme
cellid_to_plateid <- function(x){ sub("([A-Z]{3}[0-9]{5}).*","\\1", x) }
rna_to_cellid <- function(x){ sub("([A-Z]{3}[0-9]{5})([DR]?).([A-Z][0-9]{2}).*","\\1.\\3", x) }
cellid_to_well <- function(x){ sub("([A-Z]{3}[0-9]{5}).([A-Z][0-9]{2})","\\2", x) }
cellid_to_row <- function(x){ sub("([A-Z]{3}[0-9]{5}).([A-Z])[0-9]{2}","\\2", x) }
cellid_to_col <- function(x){ sub(".*[A-Z]{1}([0-9]{2})$","\\1", x) }
theme_dntr <- function (font_size = 14, font_family = "Helvetica", line_size = 0.5) {
  half_line <- font_size/2
  small_rel <- 0.857
  small_size <- small_rel * font_size
  theme_grey(base_size = font_size, base_family = font_family) %+replace% 
    theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black", size=1), 
          legend.justification = "top", legend.background = element_blank(), legend.key = element_blank(),
          legend.key.size = unit(1, "lines"), strip.background = element_blank(), strip.text=element_text(hjust=0,size=12),
          rect = element_rect(fill = "transparent", color = "black", size = 1, linetype = "solid"), axis.line = element_blank(),
          text = element_text(family = font_family, face = "plain", color = "black", size = font_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9, margin = margin(), debug = FALSE),
          axis.text.x = element_text(color="black", margin = margin(t = small_size/4), vjust = 1), 
          axis.text.y = element_text(color="black", margin = margin(r = small_size/4), hjust = 1), 
          axis.title.x = element_text(margin = margin(t = small_size/4, b = small_size/4)), 
          axis.title.y = element_text(angle = 90, margin = margin(r = small_size/4, l = small_size/4)), 
          axis.ticks = element_line(colour = "black", size = line_size))
}

d <- data.table::fread(snakemake@input[["counts"]], sep="\t", stringsAsFactors=F, header=T, data.table=F)
dim(d)

mtgenes <- readLines(snakemake@params[["mt_genes"]])

# Write out basic RNA QC too
ercc.idx <- grepl("^ERCC-",d$gene_id)
htseq.idx <- grepl("^__",d$gene_id)
mt.idx <- d$gene_id %in% mtgenes # MT-genes
ercc.frac <- colSums(d[ercc.idx,-1])/colSums(d[!htseq.idx,-1])
mt.frac <- colSums(d[mt.idx,-1])/colSums(d[!htseq.idx & !ercc.idx,-1])
filtered <- as.matrix(d[!htseq.idx & !ercc.idx,-1])

dups <- read_tsv(snakemake@input[["dupstats"]]) %>% 
  mutate(rna_library_id=sub("(.*)-picard-mark_duplicates.txt","\\1",LIBRARY)) %>% 
  select(rna_library_id, 
         rna_read_pairs=READ_PAIRS_EXAMINED,
         rna_unmapped=UNMAPPED_READS,
         rna_dup_frac=PERCENT_DUPLICATION,
         rna_unpaired=UNPAIRED_READS_EXAMINED
  )

trim.df <- tibble(
  cell_id=snakemake@params[["cells"]],
  trim1=snakemake@input[["trim1"]],
  trim2=snakemake@input[["trim2"]]
)
d.trim <- map_dfr(seq_along(trim.df$cell_id), function(x){
  r <- readLines(trim.df$trim1[[x]])
  tot <- parse_number(r[grep("Total reads processed",r)])
  adap <- parse_number(r[grep("Reads with adapters",r)])
  filt <- parse_number(r[grep("Reads written",r)])
  totbp <- parse_number(r[grep("Total basepairs processed",r)])
  r2 <- readLines(trim.df$trim2[[x]])
  tot2 <- parse_number(r2[grep("Total reads processed",r2)])
  adap2 <- parse_number(r2[grep("Reads with adapters",r2)])
  filt2 <- parse_number(r2[grep("Reads written",r2)])
  totbp2 <- parse_number(r2[grep("Total basepairs processed",r2)])
  data.frame(rna_library_id=trim.df$cell_id[[x]],
             fq1_trimfile=trim.df$trim1[[x]],
             fq2_trimfile=trim.df$trim2[[x]],
             fq1_reads=tot,fq1_reads_adap=adap,fq1_trimfilter=filt,fq1_bp=totbp,
             fq2_reads=tot2,fq2_reads_adap=adap2,fq2_trimfilter=filt2,fq2_bp=totbp2)
})

rna_qc <- tibble(
  rna_library_id=colnames(filtered),
  rna_counts=colSums(filtered),
  rna_features=colSums(filtered > 0), # Nr of features with at least 1 read 
  ercc_frac=ercc.frac,
  mtDNA_frac=mt.frac,
) %>% 
  left_join(dups, by="rna_library_id") %>% 
  left_join(d.trim, by="rna_library_id") %>%
  mutate(across(ends_with("_frac"), ~round(.x,3)))

# Plots
df <- rna_qc %>% 
  mutate(cell_id=rna_to_cellid(rna_library_id),
         col=cellid_to_col(cell_id),
         row=cellid_to_row(cell_id),
         plate_id=cellid_to_plateid(cell_id))
if(length(unique(df$plate_id))>40){
  # too many plates, or something wrong -- force into common group
  warning(paste0("n=",length(unique(df$plate_id)), " plates. Too many to plot individually, summarizing to 3-letter prefix"))
  df$plate_id <- str_sub(df$plate_id,1,3)
}

n_plates <- length(unique(df$plate_id))
sets <- max(1, n_plates/8)
if(any(sapply(unique(df$col), nchar) > 2) || any(sapply(unique(df$row), nchar) > 1)) r_scales = "free" else r_scales="fixed"

p_rna_counts <- df %>% 
  ggplot(aes(col,factor(row,levels=rev(LETTERS[1:16])),fill=log10(rna_counts))) + geom_tile() + scale_fill_viridis_c() +
  theme_dntr() + labs(x=NULL,y=NULL, caption="RNA counts") + facet_wrap(~plate_id, scales = r_scales, ncol=2) + 
  guides(fill=guide_colorbar(title.position="top")) + theme(legend.position="bottom")
ggsave(snakemake@output[["counts_plate"]], p_rna_counts, height=6.5*sets, width=8.5, limitsize=F)

p_rna_genes <- df %>% 
  ggplot(aes(col,factor(row,levels=rev(LETTERS[1:16])),fill=rna_features)) + geom_tile() + scale_fill_viridis_c() +
  theme_dntr() + labs(x=NULL,y=NULL, caption="RNA genes") + facet_wrap(~plate_id, scales = r_scales, ncol=2) + 
  guides(fill=guide_colorbar(title.position="top")) + theme(legend.position="bottom")
ggsave(snakemake@output[["genes_plate"]], p_rna_genes, height=6.5*sets, width=8.5, limitsize=F)

p_rna_ercc <- df %>% 
  ggplot(aes(col,factor(row,levels=rev(LETTERS[1:16])),fill=ercc_frac)) + geom_tile() + scale_fill_viridis_c() +
  theme_dntr() + labs(x=NULL,y=NULL, caption="RNA ERCC fraction") + facet_wrap(~plate_id, scales = r_scales, ncol=2) + 
  guides(fill=guide_colorbar(title.position="top")) + theme(legend.position="bottom")
ggsave(snakemake@output[["ercc_plate"]], p_rna_ercc, height=6.5*sets, width=8.5, limitsize=F)

p_rna_dups <- df %>% 
  ggplot(aes(log10(rna_counts),rna_dup_frac)) + geom_point(alpha=0.3) +
  theme_dntr() + labs(x="log10(read counts)",y="RNA duplicate rate") +
  facet_wrap(~plate_id,scales="free",ncol=3)
ggsave(snakemake@output[["reads_vs_dups"]], p_rna_dups, height=4.5*sets, width=8.5, limitsize=F)

# Boxplots
p_rna_counts_box <- df %>% 
  ggplot(aes(plate_id,log10(rna_counts))) + geom_boxplot(outlier.size=0.5) +
  theme_dntr() + labs(x=NULL,y="log10(counts)") +
  theme(axis.text.x=element_blank())
p_rna_genes_box <- df %>% 
  ggplot(aes(plate_id,rna_features)) + geom_boxplot(outlier.size=0.5) +
  theme_dntr() + labs(x=NULL,y="Genes") +
  theme(axis.text.x=element_blank())
p_rna_ercc_box <- df %>% 
  ggplot(aes(plate_id,ercc_frac)) + geom_boxplot(outlier.size=0.5) +
  theme_dntr() + labs(x=NULL,y="ERCC frac") +
  theme(axis.text.x=element_blank())
p_rna_dup_box <- df %>% 
  ggplot(aes(plate_id,rna_dup_frac)) + geom_boxplot(outlier.size=0.5) +
  theme_dntr() + labs(x=NULL,y="Duplicates") +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=10))

p_boxes <- plot_grid(p_rna_counts_box, p_rna_genes_box, p_rna_ercc_box, p_rna_dup_box, align="v", ncol=1, rel_heights = c(1,1,1,1.5))
ggsave(snakemake@output[["boxplots"]], p_boxes, height=8)


# If detailed qc == TRUE
if(!is.null(snakemake@input[["distr"]])){
  # rseqc read distributions
  distr_files <- snakemake@input[["distr"]]
  names(distr_files) <- snakemake@params[["cells"]]
  
  get_distr_totals <- function(x){
    x[1:3] %>%
      str_match("(.*?)\\s+(\\d+)") %>%
      as.data.frame() %>%
      select(Metric = V2, Value = V3) %>%
      mutate(Value = as.numeric(Value)) %>%
      pivot_wider(names_from = Metric, values_from = Value)
  }
  get_distr_data <- function(x){
    x[6:15] %>%
      str_match("(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)") %>%
      as.data.frame() %>%
      select(Group = V2, Total_bases = V3, Tag_count = V4, Tags_per_Kb = V5) %>%
      mutate(across(c(Total_bases, Tag_count, Tags_per_Kb), as.numeric))
  }
  
  distr.raw <- lapply(distr_files, read_lines)
  distr.tot <- lapply(distr.raw, get_distr_totals)
  distr.dat <- lapply(distr.raw, get_distr_data)
  distr.tot.df <- bind_rows(distr.tot, .id="rna_library_id")
  intergenic <- distr.tot.df$`Total Tags`-distr.tot.df$`Total Assigned Tags`
  intergenic.df <- tibble(
    rna_library_id=distr.tot.df$rna_library_id,
    Group="Intergenic",
    Tag_count=intergenic
  )
  distr.dat.df <- bind_rows(distr.dat, .id="rna_library_id") %>% 
    bind_rows(intergenic.df) %>% 
    arrange(rna_library_id, Group)
  distr.df <- distr.tot.df %>% 
    left_join(distr.dat.df, by="rna_library_id") %>% 
    mutate(Tag_frac=Tag_count/`Total Tags`) %>% 
    arrange(rna_library_id)
  # distr.df.wide.count <- distr.df %>% 
  #   select(-Total_bases, -Tags_per_Kb, -Tag_frac) %>% 
  #   spread(key="Group",value=c("Tag_count"))
  distr.df.wide.frac <- distr.df %>% 
    select(-Total_bases, -Tags_per_Kb, -Tag_count) %>% 
    spread(key="Group",value=c("Tag_frac"))
  
  distr.plot.df <- distr.df.wide.frac %>% 
    left_join(enframe(mt.frac,name="rna_library_id",value="mtDNA"), by="rna_library_id") %>% 
    mutate(across(starts_with("Total "), ~log10(.x))) %>% 
    pivot_longer(!rna_library_id, names_to="Group", values_to="val") %>% 
    # filter(Group %in% c("CDS_Exons","5'UTR_Exons","3'UTR_Exons","Introns","Intergenic", "mtDNA")) %>% 
    mutate(stat_group=factor(
      case_when(Group %in% c("CDS_Exons","5'UTR_Exons","3'UTR_Exons") ~ "On-target", 
                Group %in% c("Introns","Intergenic","mtDNA") ~ "Off-target",
                TRUE ~ "Totals (log10)"),
      levels=c("On-target","Off-target","Totals (log10)")),
      Group=case_when(Group=="CDS_Exons" ~ "Exons",
                      Group=="5'UTR_Exons" ~ "5'UTR",
                      Group=="3'UTR_Exons" ~ "3'UTR",
                      TRUE ~ Group),
      Group=factor(Group, levels=c("Exons","5'UTR","3'UTR","Introns","mtDNA","Intergenic","Total Reads","Total Tags","Total Assigned Tags")))

  rna_qc <- rna_qc %>% left_join(distr.df.wide.frac,by="rna_library_id")

  pdistr <- distr.plot.df %>% 
    left_join(select(df, rna_library_id, plate_id),by="rna_library_id") %>% 
    filter(!is.na(Group)) %>% 
    # filter(Group %in% c("CDS_Exons","5'UTR_Exons","3'UTR_Exons","Introns","Intergenic", "mtDNA")) %>% 
    ggplot(aes(Group, val)) + 
    geom_boxplot() +
    scale_color_brewer(palette="Set1") +
    labs(x="", y="RNA fraction") +
    facet_wrap(plate_id~stat_group, nrow=1, scales="free") +
    theme(axis.title.x=element_blank(), strip.text=element_text(hjust=0,size=10),
          panel.border = element_rect(fill = NA, colour = "black", size=1),
          rect = element_blank(), axis.line=element_blank(),
          legend.position="bottom", legend.justification = "right", legend.title = element_blank(),
          legend.direction = "vertical", legend.margin=margin(t=-.5,unit="lines"))
  
  ggsave(snakemake@output[["rseqc_distr"]], pdistr, height=8)
}

write_tsv(rna_qc, snakemake@output[["qc"]])