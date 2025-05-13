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
scalefactors<-read.table("results/DNA_37_SK_250502_Downsample/ALL40/clones/ALL40-scalefactor_minima-g10-10000.txt")
scalefactors_2<-scalefactors[,-1]
rownames(scalefactors_2)<-scalefactors$V1
scps_sc<-read.table("results/DNA_37_SK_250502_Downsample/ALL40/clones/ALL40-scps-scCN-10000-g10.txt", header=TRUE)
#scalefactors<-t(as.data.frame(local_min_s))
binsize<-10000
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

# pval.ploidy <- function(xmat, ploidy.rates) {
#   overlapSize = sum(xmat$meanFragSize*xmat$numReads)/(sum(xmat$numReads)*2)
#   cur.p <- dpois(x=round(sum(xmat$numOverlaps)), lambda=sum(ploidy.rates[xmat$ploidy] * xmat$nr_bins*binsize/overlapSize), log=T)
#   p.double <- dpois(x=round(sum(xmat$numOverlaps)), lambda=sum(ploidy.rates[xmat$ploidy*2] * xmat$nr_bins*binsize/overlapSize), log=T)
#   p.half <- dpois(x=round(sum(xmat$numOverlaps)), lambda=sum(ploidy.rates[round(xmat$ploidy*0.5001)] * xmat$nr_bins*binsize/overlapSize), log=T)
# 
#   lods.up <- p.double-cur.p
#   lods.down <- p.half-cur.p
#   setNames(c(cur.p, lods.up, lods.down), nm=c("logp.dens", "lods.up", "lods.down"))
# }

pois.ploidy<-function(xmat, ploidy.rates, test){
  overlapSize = sum(xmat$meanFragSize*xmat$numReads)/(sum(xmat$numReads)*2)
  cur.p <- dpois(x=round(sum(xmat$numOverlaps)), lambda=sum(ploidy.rates[round(xmat$ploidy*test)] * xmat$nr_bins*binsize/overlapSize), log=T)
}

# #So we need to test 
# head(scalefactors)
# 
# ret.mat <- sapply(unique(scps_sc_f$dna_library_id), function(x) {
#   xmat <- scps_sc_f %>% filter(dna_library_id == x)
#   xmat<-xmat%>%filter(!ploidy==0) #Skip ploidy == 0 (happens very rarely -1 cell so far but messes with calculations)
#   pois.ploidy(xmat, ploidy.rates, 1) #Then here should sapply through the different ones 
# })
# 
# 
# scps_sc_f2<-scps_sc_f%>%filter(dna_library_id%in%rownames(scalefactors))
ret.mat <- sapply(unique(scps_sc$dna_library_id), function(x) {
  xmat <- scps_sc %>% filter(dna_library_id == x) %>% filter(ploidy != 0)

  # Apply the pois.ploidy function for each column in scalefactors (6 columns)
  scale_vals <- sapply(1:ncol(scalefactors_2), function(i) {
    # Compute scale_val and ensure it's at least 0.5001
    #scale_val <- pmax(scalefactors[x, i] / 2, 0.5001)
    scale_val <- (scalefactors_2[x, i] / 2)
    pois.ploidy(xmat, ploidy.rates, scale_val)
  })

  # Add a final column with scale_val = 1
  scale_vals_with_one <- c(scale_vals, pois.ploidy(xmat, ploidy.rates, 1))

  return(scale_vals_with_one)
})


ret.mat <- sapply(unique(scps_sc$dna_library_id), function(x) {
  xmat <- scps_sc %>% filter(dna_library_id == x) %>% filter(ploidy != 0)
  scale_vals <- sapply(setdiff(1:6, 2), function(i) {
    scale_val <- (i + 0.0001) / 2  # e.g., (1.0001)/2 = 0.50005, (3.0001)/2 = 1.50005, etc.
    pois.ploidy(xmat, ploidy.rates, scale_val)
  })
  scale_vals_with_one <- c(scale_vals, pois.ploidy(xmat, ploidy.rates, 1))
  
  return(scale_vals_with_one)
})

df_wide <- as.data.frame(t(ret.mat)) %>%
  rownames_to_column(var = "dna_library_id")

colnames(df_wide)<-c("dna_library_id", "n0.5", "n1.5", "n2", "n2.5", "n3", "Initial")
head(df_wide)
df_wide$lods.0.5<-df_wide$n0.5-df_wide$Initial
df_wide$lods.1.5<-df_wide$n1.5-df_wide$Initial
df_wide$lods.2<-df_wide$n2-df_wide$Initial
df_wide$lods.2.5<-df_wide$n2.5-df_wide$Initial
df_wide$lods.3<-df_wide$n3-df_wide$Initial

df_wide
#Output plots

p1<-ggplot(df_wide, aes(x=1, y=n0.5))+
  geom_quasirandom(width=0.2, alpha=0.5)+
  theme_bw()+
  geom_hline(aes(yintercept=4.6))

p2<-ggplot(df_wide, aes(x=1, y=n1.5))+
  geom_quasirandom(width=0.2, alpha=0.5)+
  theme_bw()+
  geom_hline(aes(yintercept=4.6))

p3<-ggplot(df_wide, aes(x=1, y=n2))+
  geom_quasirandom(width=0.2, alpha=0.5)+
  theme_bw()+
  geom_hline(aes(yintercept=4.6))

p4<-ggplot(df_wide, aes(x=1, y=n2.5))+
  geom_quasirandom(width=0.2, alpha=0.5)+
  theme_bw()+
  geom_hline(aes(yintercept=4.6))

p5<-ggplot(df_wide, aes(x=1, y=n3))+
  geom_quasirandom(width=0.2, alpha=0.5)+
  theme_bw()+
  geom_hline(aes(yintercept=4.6))

p6<-ggplot(df_wide, aes(x=1, y=Initial))+
  geom_quasirandom(width=0.2, alpha=0.5)+
  theme_bw()+
  geom_hline(aes(yintercept=4.6))


grid.arrange(p1, p2, p3, p4, p5, p6)

df_wide%>%filter(lods.1.5>4.6)

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








