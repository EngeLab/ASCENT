library(tidyverse)
library(rtracklayer)
###
#Code snippet to calculate noise level for absolute copy number calling
#This code is for data-adaptive estimation of noise level to be used for absolute copy number calling
#The prerequisite to do this is to have cells that are known to be correctly scaled (so that copy state=2 is indeed 2, but not 4 for example). 
#The scps-scCN-{binsize}-g{gamma}.txt file can be used as input and filtered for cells that have a known ploidy

###Functions: 
ploidy.stat <- function(ploidy, frac, dup.rate) {
  # prob. of picking the same number twice * number of possible pairs (n over k)
  # n over k: n!/(k!*(n-k)!). k is 2 (=pairs)orial(2)*factorial(ploidy+ploidy*dup.rate-2))
  (frac/(ploidy+ploidy*dup.rate/2))^2 * factorial(ploidy+ploidy*dup.rate)/(factorial(2)*factorial(ploidy+ploidy*dup.rate-2))
}

calc.noise.from.known.diploid <- function(duprates.test, sample.SCP) {
  res.mat2.closed <- sapply(duprates.test, function(x) {
    ploidy.stat(2, frac, x)
  })
  
  my.lm.calc <- lm(res.mat2.closed ~ I(duprates.test))
  
  ret <- (median(sample.SCP$depth2[sample.SCP$ploidy == 2])-coefficients(my.lm.calc)[1])/(coefficients(my.lm.calc)[2])
  names(ret) <- NULL
  return(ret)
}

#Calculations: 

frac <- 0.005
duprates.test <- seq(0, 0.4, by=0.1)

sample <- read.table("results/ALL/ALL40/clones/ALL40-scps-scCN-10000-g10.txt", header=T)
sample <- sample %>% filter(sample$ploidy %in% c(2,3,4))

calculated_noise<-calc.noise.from.known.diploid(duprates.test, sample)
calculated_noise

#This calculated noise can replace the dup.rate in workflow/scripts/scp_log_odds.R 
