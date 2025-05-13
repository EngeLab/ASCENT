library(tidyverse)
library(dplyr)
library(readr)

arms <- read.table(snakemake@input[["arms"]], header = FALSE, col.names = c("chr_arm", "start_arm", "end_arm"))
bed <- read.table(snakemake@input[["bed"]], header = FALSE, col.names = c("chr_arm", "start", "end"))

arms$chr <- gsub("_p|_q", "", arms$chr_arm)

p_arm_ends <- arms %>%
  filter(grepl("_p$", chr_arm)) %>%
  select(chr, end_arm)
colnames(p_arm_ends)<-c("chr","p_end")

bed_transformed <- bed %>%
  left_join(arms, by = c("chr_arm")) %>%  # Join by chromosome name
  left_join(p_arm_ends, by = "chr") %>%  # Add p-arm end for q-arm adjustment
  mutate(
    new_start = ifelse(grepl("_q$", chr_arm), start + p_end, start),  # Adjust start if in q-arm
    new_end = ifelse(grepl("_q$", chr_arm), end + p_end, end)  # Adjust end if in q-arm
  ) 
print<-bed_transformed%>%select(chr, new_start, new_end)


# Write output

write.table(print, snakemake@output[["bed"]], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

