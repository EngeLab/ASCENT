#Script to calculate log-odds ratios of current copy number
library(tidyverse)
library(uwot)
library(dbscan)
library(ComplexHeatmap)
library(circlize)
library(furrr)
library(patchwork)
library(gridExtra)
library(ggplot2)
library(ggbeeswarm)


scps_sc<-read_tsv(snakemake@input[["scps_sc"]])
bins<-read_tsv(snakemake@input[["bins"]])
binsize <- as.integer(snakemake@wildcards[["binsize"]]) 

scps_sc_f<-scps_sc%>%
  group_by(dna_library_id)%>%
  mutate(sum_bins=sum(nr_bins))%>%
  filter(sum_bins>dim(bins)[1]/2)

ploidy.stat <- function(ploidy, frac, dup.rate) {
  # prob. of picking the same number twice * number of possible pairs (n over k)
  # n over k: n!/(k!*(n-k)!). k is 2 (=pairs)orial(2)*factorial(ploidy+ploidy*dup.rate-2))
  (frac/(ploidy+ploidy*dup.rate/2))^2 * factorial(ploidy+ploidy*dup.rate)/(factorial(2)*factorial(ploidy+ploidy*dup.rate-2))
}

ploidy.rates <- setNames(sapply(1:20, function(ploidy) {
  ploidy.stat(ploidy, frac=0.005, dup.rate=0.1)
}), nm=1:20)

pval.ploidy <- function(xmat, ploidy.rates) {
  overlapSize = sum(xmat$meanFragSize*xmat$numReads)/(sum(xmat$numReads)*2)
  cur.p <- dpois(x=round(sum(xmat$numOverlaps)), lambda=sum(ploidy.rates[xmat$ploidy] * xmat$nr_bins*binsize/overlapSize), log=T)
  p.double <- dpois(x=round(sum(xmat$numOverlaps)), lambda=sum(ploidy.rates[xmat$ploidy*2] * xmat$nr_bins*binsize/overlapSize), log=T)
  p.half <- dpois(x=round(sum(xmat$numOverlaps)), lambda=sum(ploidy.rates[round(xmat$ploidy*0.5001)] * xmat$nr_bins*binsize/overlapSize), log=T)
  
  #lods.up <- cur.p-p.double
  lods.up <- p.double-cur.p
  #lods.down <- cur.p-p.half
  lods.down <- p.half-cur.p
  setNames(c(cur.p, lods.up, lods.down), nm=c("logp.dens", "lods.up", "lods.down"))
}

ret.mat <- sapply(unique(scps_sc_f$dna_library_id), function(x) {
  xmat <- scps_sc_f %>% filter(dna_library_id == x)
  xmat<-xmat%>%filter(!ploidy==0) #Skip ploidy == 0 (happens very rarely -1 cell so far but messes with calculations)
  pval.ploidy(xmat, ploidy.rates)
})

df_wide <- as.data.frame(t(ret.mat)) %>%
  rownames_to_column(var = "dna_library_id")

#Output plots

p1<-ggplot(df_wide, aes(x=1, y=lods.up))+
  geom_quasirandom(width=0.2, alpha=0.5)+
  theme_bw()+
  geom_hline(aes(yintercept=4.6))

p2<-ggplot(df_wide, aes(x=1, y=lods.down))+
  geom_quasirandom(width=0.2, alpha=0.5)+
  theme_bw()+
  geom_hline(aes(yintercept=4.6))


df_wide<-df_wide%>%mutate(multiplication=case_when(
  lods.down>(4.6)~0.5, #This means 1% chance that it is by chance 
  lods.up>(4.6)~2, # 
  TRUE ~ 1
))


scps_plot<-left_join(scps_sc_f, df_wide)

p3<-ggplot(scps_plot, aes(x=ploidy, y=depth2, color=as.factor(multiplication)))+
  geom_quasirandom(width=0.2)+
  theme_bw()+
  scale_color_brewer(palette="Dark2", na.value = "grey50")

png(snakemake@output[["log_odds_sc"]], width=2400, height=3000, units="px", res=150)
grid.arrange(p1, p2, p3)
dev.off()

write_tsv(df_wide, snakemake@output[["log_odds_df_sc"]])








