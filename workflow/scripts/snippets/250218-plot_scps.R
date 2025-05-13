# SCPs ALL
cellid_to_plateid <- function(x){
  sub("([A-Z]{3}[0-9]{5}).*","\\1", x)
}
library(tidyverse)
library(ggbeeswarm)
options(scipen=999)

# Metadata
m_files <- list.files("results/ALL-latest","metadata_long\\.tsv$",recursive=T,full.names=T)
m <- map_dfr(m_files, ~read_tsv(.x, show_col_types = FALSE))

# QC 
q_files <- list.files("results/ALL-latest","qc_dna\\.tsv$",recursive=T,full.names=T)
q <- map_dfr(q_files, ~read_tsv(.x, show_col_types = FALSE))

# Clones
cl_files <- list.files("results/ALL-latest","clones_final\\.txt$",recursive=T,full.names=T)
cl <- map_dfr(cl_files, ~read_tsv(.x, show_col_types = FALSE))

do <- m %>% 
  left_join(q) %>% 
  left_join(select(cl, -timepoint))

# do6 <- do %>% filter(patient_id=="ALL6")
# mpcf <- read_tsv("results/ALL-latest/ALL6/clones/ALL6-mpcf-g1-20000.txt.gz")

scp_files <- list.files(
  "results/ALL-latest",
  pattern = "scps-20000\\.txt$",  # Escaped dot for literal match
  recursive = TRUE,
  full.names = TRUE
)
d <- map_dfr(scp_files, ~read_tsv(.x, show_col_types = FALSE)) %>% 
  left_join(m, by="dna_library_id") %>% 
  left_join(select(cl, -timepoint), by="dna_library_id")


d %>% 
  # filter(patient_id=="ALL68") %>% 
  filter(ploidy==2 & !is.na(clone) & clone!=0) %>% 
  ggplot(aes(x=plate_id, y=depth2*1e5, color=patient_id)) +
  geom_quasirandom(alpha=0.5, show.legend = F) +
  theme(axis.text.x = element_text(angle=90)) + 
  facet_grid(.~patient_id, scales = "free_x", space = "free")


# ALL6 different scp run comparisons
d1 <- d %>% filter(patient_id=="ALL6")
d2 <- read_tsv("~/scps/ALL6-scps-dupremoval_total.txt") %>% 
  mutate(plate_id=cellid_to_plateid(dna_library_id))

d1 %>% 
  filter(ploidy==2 & !is.na(clone) & clone!=0) %>% 
  ggplot(aes(x=plate_id, y=depth2*1e5)) +
  geom_quasirandom(alpha=0.5, show.legend = F) +
  geom_quasirandom(data=d2, col="red", alpha=0.5) +
  theme(axis.text.x = element_text(angle=90))
