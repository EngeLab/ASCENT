rna_to_cellid <- function(x){
  sub("([A-Z]{3}[0-9]{5})(R?).([A-Z][0-9]{2}).*","\\1.\\3", x)
}
# From load_all.R
d <- dl$ALL6
d2 <- dl$ALL5

# Chr counts metadata
# qc <- read_tsv("results/ALL-latest/ALL6/qc/ALL6-qc_dna.tsv") %>% 
qc <- read_tsv("results/ALL-latest/ALL5/qc/ALL5-qc_dna.tsv") %>% 
  select(dna_library_id, starts_with("dna_frac")) %>% as.data.frame
rownames(qc) <- dna_to_cellid(qc$dna_library_id)
d <- AddMetaData(d, qc)

# Fusions
fu <- read_tsv("/wrk/resources/dntr/results/ALL-latest/rna_fusions/ALL6/star-fusion.fusion_predictions.tsv.samples_deconvolved.abridged.tsv")
etv6_fused <- fu %>% filter(grepl("ETV6", `#FusionName`)) %>% pull(Cell) %>% rna_to_cellid()
msh6_fused <- fu %>% filter(grepl("MSH6", `#FusionName`)) %>% pull(Cell) %>% rna_to_cellid()

d$etv6_fusion <- d$cell_id %in% etv6_fused
d$msh6_fusion <- d$cell_id %in% msh6_fused

ptime <- DimPlot(d, group.by="timepoint_short", cols=timepoint_cols_short, order=T, pt.size=0.5) + theme_dntr(axis_labels=F)
pdntt <- FeaturePlot(d, "ENSG00000107447", order=T) + theme_dntr(axis_labels=F)
pxist <- FeaturePlot(d, "ENSG00000229807", order=T) + theme_dntr(axis_labels=F)
p1 <- DimPlot(d, group.by = "etv6_fusion", order=T) + scale_color_manual(values=c("FALSE"="grey","TRUE"="red")) + theme_dntr(axis_labels = F)
p2 <- DimPlot(d, group.by = "msh6_fusion", order=T) + scale_color_manual(values=c("FALSE"="grey","TRUE"="red")) + theme_dntr(axis_labels = F)
ptime + pdntt + p1 + p2 + plot_layout(nrow=1)

# clones <- cl.annot %>% filter(patient_id=="ALL6")

# DNTT: ENSG00000107447
# XIST: ENSG00000229807

# Add X/Y chr set. Add module score
gtf <- read_tsv("/wrk/resources/genomes/grch38-gdc/gencode.v36.filtered.gtf",
                comment = "#",
                col_names = c("chr", "source", "type", "start", "end", "score", 
                              "strand", "phase", "attributes"))
gtf_genes <- gtf %>%
  filter(type == "gene") %>%
  mutate(
    gene_id = str_extract(attributes, 'gene_id "([^"]+)"') %>%
      str_replace('gene_id "', '') %>%
      str_replace('"', '')
  )

y_genes <- gtf_genes %>% filter(chr=="chrY") %>% pull(gene_id)
x_genes <- gtf_genes %>% filter(chr=="chrX") %>% pull(gene_id)
gene_lists <- list("y_genes"=y_genes, "x_genes"=x_genes)

dj <- merge(d, d2)

dj <- AddModuleScore(dj, gene_lists, name=names(gene_lists))
# d <- AddModuleScore(d, gene_lists, name=names(gene_lists))
# d2 <- AddModuleScore(d2, gene_lists, name=names(gene_lists))

par(mfrow=c(1,3))
xr <- quantile(dj$y_genes1, c(.01, .99))
yr <- quantile(dj$dna_frac_chrY, c(.01, .99),na.rm = T)
plot(dj$y_genes1[dj$patient_id=="ALL5"], 
     dj$dna_frac_chrY[dj$patient_id=="ALL5"], 
     main="ALL5 (XY) vs ALL6 (XX)", ylab="DNA read frac at Y",xlab="Y-chr gene set expr", 
     sub="red = relapse1",
     col="black", xlim=xr, ylim=yr, cex=.8)
points(dj$y_genes1[dj$patient_id=="ALL6" & dj$timepoint!="relapse 1"], 
     dj$dna_frac_chrY[dj$patient_id=="ALL6" & dj$timepoint!="relapse 1"],
     col="blue", cex=.8)
points(dj$y_genes1[dj$patient_id=="ALL6" & dj$timepoint=="relapse 1"], 
       dj$dna_frac_chrY[dj$patient_id=="ALL6" & dj$timepoint=="relapse 1"],
       col="red", cex=.8)

boxplot(split(dj@assays$RNA@data["ENSG00000229807",dj$patient_id=="ALL5"], dj$plate_id[dj$patient_id=="ALL5"]), las=2, main="XIST RNA, ALL5")
boxplot(split(dj@assays$RNA@data["ENSG00000229807",dj$patient_id=="ALL6"], dj$plate_id[dj$patient_id=="ALL6"]), las=2, main="XIST RNA, ALL6")


yr <- c(0, quantile(dj@assays$RNA@data["ENSG00000229807",], c(.99),na.rm = T))
plot(dj$dna_frac_chrY[dj$patient_id=="ALL5"], 
     dj@assays$RNA@data["ENSG00000229807", dj$patient_id=="ALL5"],
     main="ALL5 (male, black) vs ALL6 (red, female)", 
     ylab="XIST RNA",xlab="Y-chr gene set expr", 
     col="black", xlim=xr, ylim=yr, cex=.8, pch=19)
points(dj$dna_frac_chrY[dj$patient_id=="ALL6" & dj$timepoint!="relapse 1"],
       dj@assays$RNA@data["ENSG00000229807", dj$patient_id=="ALL6" & dj$timepoint!="relapse 1"], 
       col="blue", cex=.8)
points(dj$dna_frac_chrY[dj$patient_id=="ALL6" & dj$timepoint=="relapse 1"],
       dj@assays$RNA@data["ENSG00000229807", dj$patient_id=="ALL6" & dj$timepoint=="relapse 1"], 
       col="red", cex=.8)


VlnPlot(d, features = "ENSG00000229807", group.by = "clone")
