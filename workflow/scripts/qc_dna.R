suppressPackageStartupMessages(library(tidyverse))
# Functions
cellid_to_plateid <- function(x){ sub("([A-Z]{3}[0-9]{5}).*","\\1", x) }
dna_to_cellid <- function(x){ sub("([A-Z]{3}[0-9]{5})(D?)_([A-Z][0-9]{2}).*","\\1.\\3", x) }
cellid_to_well <- function(x){ sub("([A-Z]{3}[0-9]{5}).([A-Z][0-9]{2})","\\2", x) }
cellid_to_row <- function(x){ sub("([A-Z]{3}[0-9]{5}).([A-Z])[0-9]{2}","\\2", x) }
cellid_to_col <- function(x){ sub("([A-Z]{3}[0-9]{5}).[A-Z]([0-9]{2})","\\2", x) }
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

# Get Picard dup metrics
dupfiles <- snakemake@input[["dupstats"]]
d.dup <- read_tsv(dupfiles, col_types="cnnnnnnnnnn") %>% 
  select(dna_library_id=LIBRARY, 
         bam_read_pairs=READ_PAIRS_EXAMINED, 
         bam_unmapped=UNMAPPED_READS, 
         bam_dup_frac=PERCENT_DUPLICATION, 
         bam_unpaired=UNPAIRED_READS_EXAMINED) %>% 
  distinct()
cat("Dup files:",nrow(d.dup),"\n")

# Get FASTQ counts 
# TODO: Add back and account for multiple fastq per cell (when pick-up)
# fqfiles <- snakemake@input[["fqfiles"]]
# cellids <- sub("(.*)_S.*","\\1",basename(fqfiles))
# trimfiles <- paste0(snakemake@params[["dnaout"]],"/",cellids,"/trimmed/",basename(fqfiles),"_trimming_report.txt")
# names(trimfiles) <- cellids

# fqfiles2 <- sub("_R1","_R2",fqfiles)
# trimfiles2 <- paste0(snakemake@params[["dnaout"]],"/",cellids,"/trimmed/",basename(fqfiles2),"_trimming_report.txt")
# names(trimfiles2) <- cellids

# Get Picard insert size metrics
insertfiles <- snakemake@input[["insfiles"]]
d.ins <- read_tsv(insertfiles, col_types="cnn",col_names=c("filename","bp","count")) %>% 
  filter(!is.na(count)) %>%
  mutate(dna_library_id=sub("(.*)\\-picard.*","\\1",basename(filename))) %>% 
  select(dna_library_id,everything(),-filename) %>% 
  distinct()
cat("Insert files:",length(unique(d.ins$dna_library_id)),"\n")

d.ins.means <- d.ins %>% 
  group_by(dna_library_id) %>% 
  summarize(dna_insert_mean=sum(bp*count)/sum(count),
            dna_insert_sd=sqrt(sum((bp-dna_insert_mean)**2*count)/(sum(count)-1)))

d.ins.plate_means <- d.ins %>% 
  mutate(plate_id=cellid_to_plateid(dna_library_id)) %>% 
  group_by(plate_id) %>% 
  summarize(mean=sum(bp*count)/sum(count),
            sd=sqrt(sum((bp-mean)**2*count)/(sum(count)-1)))

trim.df <- tibble(
  cell_id=snakemake@params[["cells"]],
  trim1=snakemake@input[["trim1"]],
  trim2=snakemake@input[["trim2"]]
)
# cellids <- seedfile$cell
# trimfiles <- snakemake@input[["trim1"]]
# trimfiles2 <- snakemake@input[["trim2"]]
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
  data.frame(dna_library_id=trim.df$cell_id[[x]],
             fq1_trimfile=trim.df$trim1[[x]],
             fq2_trimfile=trim.df$trim2[[x]],
             fq1_reads=tot,fq1_reads_adap=adap,fq1_trimfilter=filt,fq1_bp=totbp,
             fq2_reads=tot2,fq2_reads_adap=adap2,fq2_trimfilter=filt2,fq2_bp=totbp2)
})

# Get chromosome counts per cell
chrfiles <- snakemake@input[["chrfiles"]]
d.cc <- read_tsv(chrfiles,col_names=c("filename","chr","n"),col_types="ccn") %>% 
  mutate(dna_library_id=sub("(.*)\\-chr_counts.*","\\1",filename)) %>% select(-filename) %>% 
  distinct()

cat("Chr read-count files:",length(unique(d.cc$dna_library_id)),"\n")
d.cc.frac <- d.cc %>% 
  group_by(dna_library_id) %>% 
  mutate(total_n=sum(n)) %>% 
  ungroup() %>% 
  mutate(frac=n/total_n,
         chr=factor(chr,levels=paste0("chr",c(1:22,"X","Y","M"))))
d.cc.w <- d.cc.frac %>% 
  select(dna_library_id,chr,frac) %>% 
  rename(val=frac,dna_frac=chr) %>% 
  spread(dna_frac,val,sep="_")
d.cc.w.plate <- d.cc.w %>% 
  mutate(plate_id=cellid_to_plateid(dna_library_id)) %>% 
  group_by(plate_id) %>% 
  summarise(across(starts_with("dna_frac_chr"),
                   mean, na.rm=T, .names = "mean_{.col}"))

# d.scp <- read_tsv(snakemake@input[["scp"]], 
#                   col_names=c("dna_library_id","scp_cov_bin","scp_depth","scp"),
#                   col_types = "cinn", skip = 1) %>% 
#   filter(scp_cov_bin==as.numeric(snakemake@params[["scp_cov_bin"]]))

d.out <- full_join(d.dup, d.trim, by="dna_library_id") %>% 
  left_join(d.ins.means, by="dna_library_id") %>% 
  left_join(d.cc.w, by="dna_library_id")
  # left_join(d.scp, by="dna_library_id")
head(d.out)

df <- d.out %>% 
  mutate(bam_mapped_frac=1-(bam_unmapped/((bam_read_pairs*2)+bam_unpaired+bam_unmapped)),
         cell_id=dna_to_cellid(dna_library_id),
         col=cellid_to_col(cell_id),
         row=cellid_to_row(cell_id),
         plate_id=cellid_to_plateid(cell_id))
n_plates <- length(unique(df$plate_id))
sets <- max(1, n_plates/4)

p_dna_fq <- df %>% 
  ggplot(aes(col,factor(row,levels=rev(LETTERS[1:16])),fill=log10(fq1_reads))) + geom_tile() + scale_fill_viridis_c() +
  theme_dntr() + labs(x=NULL,y=NULL, caption="DNA FASTQ read-pairs") + facet_wrap(~plate_id, ncol=2) + 
  guides(fill=guide_colorbar(title.position="top")) + theme(legend.position="bottom")
ggsave(snakemake@output[["fq_plate"]], p_dna_fq, height=6.5*sets, width=8.5,limitsize=F)

p_dna_fq_box <- df %>% 
  ggplot(aes(plate_id,log10(fq1_reads))) + geom_boxplot(outlier.size=0.5) +
  theme_dntr() + labs(x=NULL,y="DNA FASTQ reads, log10") +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=10))
ggsave(snakemake@output[["fq_boxplot"]], p_dna_fq_box, height=4)

p_dna_mapped <- df %>% 
  ggplot(aes(col,factor(row,levels=rev(LETTERS[1:16])),fill=bam_mapped_frac)) + geom_tile() + scale_fill_viridis_c() +
  theme_dntr() + labs(x=NULL,y=NULL, caption="DNA mapped fraction") + facet_wrap(~plate_id, ncol=2) + 
  guides(fill=guide_colorbar(title.position="top")) + theme(legend.position="bottom")
ggsave(snakemake@output[["map_plate"]], p_dna_mapped, height=6.5*sets, width=8.5,limitsize=F)

p_dna_dup <- df %>% 
  ggplot(aes(col,factor(row,levels=rev(LETTERS[1:16])),fill=bam_dup_frac)) + geom_tile() + scale_fill_viridis_c() +
  theme_dntr() + labs(x=NULL,y=NULL, caption="DNA duplicate fraction") + facet_wrap(~plate_id, ncol=2) + 
  guides(fill=guide_colorbar(title.position="top")) + theme(legend.position="bottom")
ggsave(snakemake@output[["dup_plate"]], p_dna_dup, height=6.5*sets, width=8.5,limitsize=F)

p_dna_dup_box <- df %>% 
  ggplot(aes(plate_id,bam_dup_frac)) + geom_boxplot(outlier.size=0.5) +
  theme_dntr() + labs(x=NULL,y="DNA duplicate rate") +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=10))
ggsave(snakemake@output[["dup_boxplot"]], p_dna_dup_box, height=4)

p_dna_dup_v_reads <- df %>%
  ggplot(aes(log10(bam_read_pairs),bam_dup_frac)) + geom_point() +
  theme_dntr() + labs(x="log10(read pairs)",y="DNA duplicate rate") +
  geom_vline(xintercept=log10(snakemake@config$dna$min_count),linetype="dotted")+
  geom_hline(yintercept=snakemake@config$dna$max_dup_frac,linetype="dotted")+
  facet_wrap(~plate_id,scales="free",ncol=3)
ggsave(snakemake@output[["reads_vs_dups"]], p_dna_dup_v_reads, height=4.5*sets, width=8.5, limitsize=F)

p_dna_ins <- df %>% 
  ggplot(aes(col,factor(row,levels=rev(LETTERS[1:16])),fill=dna_insert_mean)) + geom_tile() + scale_fill_viridis_c() +
  theme_dntr() + labs(x=NULL,y=NULL, caption="DNA insert size means (bp)") + facet_wrap(~plate_id, ncol=2) + 
  guides(fill=guide_colorbar(title.position="top")) + theme(legend.position="bottom")
ggsave(snakemake@output[["ins_plate"]], p_dna_ins, height=6.5*sets, width=8.5, limitsize=F)

p_dna_ins_hist <- d.ins %>% 
  mutate(plate_id=cellid_to_plateid(dna_library_id)) %>% 
  ggplot(aes(x=bp, y=count,group=dna_library_id)) + geom_line(alpha=0.2) + theme_dntr() +
  xlim(c(0,500)) +
  geom_vline(data=d.ins.plate_means, aes(xintercept=mean),linetype="dotted",color="red") +
  geom_vline(data=d.ins.plate_means, aes(xintercept=mean-sd),linetype="dotted",color="black") +
  geom_vline(data=d.ins.plate_means, aes(xintercept=mean+sd),linetype="dotted",color="black") +
  geom_text(data=d.ins.plate_means, aes(x=Inf,y=Inf,label=round(mean)),color="red",hjust=1.2,vjust=1.2,inherit.aes=F) +
  facet_wrap(~plate_id,scales = "free_y",ncol=2) +
  theme(axis.text.y = element_blank())
ggsave(snakemake@output[["ins_hist"]], p_dna_ins_hist, height=6.5*sets, width=8.5,limitsize=F)

p_dna_ins_box <- df %>% 
  ggplot(aes(plate_id,dna_insert_mean)) + geom_boxplot(outlier.size=0.5) +
  theme_dntr() + labs(x=NULL,y="DNA insert size") +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=10))
ggsave(snakemake@output[["ins_boxplot"]], p_dna_ins_box, height=4)

p_dna_chrfrac <- d.cc.frac %>% 
  mutate(plate_id=cellid_to_plateid(dna_library_id)) %>% 
  ggplot(aes(chr,log(frac))) + 
  geom_boxplot(outlier.size = 0.5) +
  theme_dntr() +
  facet_wrap(~plate_id,scales="free_y",ncol=2) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=8)) +
  labs(x=NULL,y="log fraction")
ggsave(snakemake@output[["counts_boxplot"]], p_dna_chrfrac, height=6.5*sets, width=8.5,limitsize=F)

mfrac <- d.cc.w %>% 
  mutate(cell_id=dna_to_cellid(dna_library_id),
         col=cellid_to_col(cell_id),
         row=cellid_to_row(cell_id),
         plate_id=cellid_to_plateid(cell_id))

p_dna_mfrac <- mfrac %>% 
  ggplot(aes(plate_id,log(dna_frac_chrM))) + geom_boxplot(outlier.size=0.5) +
  theme_dntr() + labs(x=NULL,y="DNA chrM fraction, log scale") +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=10))
ggsave(snakemake@output[["counts_boxplot_chrM"]], p_dna_mfrac, height=4)

p_dna_mfrac_plate <- mfrac %>% 
  ggplot(aes(col,factor(row,levels=rev(LETTERS[1:16])),fill=log(dna_frac_chrM))) + geom_tile() + scale_fill_viridis_c() +
  theme_dntr() + labs(x=NULL,y=NULL, caption="DNA chrM fraction, log scale") + facet_wrap(~plate_id, ncol=2) + 
  guides(fill=guide_colorbar(title.position="top")) + theme(legend.position="bottom",axis.text=element_text(size=8))
ggsave(snakemake@output[["counts_plate_chrM"]], p_dna_mfrac_plate, height=6.5*sets, width=8.5,limitsize=F)

# p_dna_scp_box <- df %>% 
#   # ggplot(aes(plate_id,scp)) + geom_boxplot(outlier.size=0.5) +
#   ggplot(aes(plate_id,scp)) + ggbeeswarm::geom_beeswarm() +
#   theme_dntr() + labs(x=NULL,subtitle=paste0("SCP (fraction ",unique(d.scp$scp_cov_bin),"x cov at model ",unique(d.scp$scp_depth),"x mapped)")) +
#   theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=10))
# ggsave(snakemake@output[["scp_boxplot"]], p_dna_scp_box, height=4)

# Print out joint qc table per cell
write_tsv(d.out, snakemake@output[["qc_tsv"]])