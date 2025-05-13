library(tidyverse)
library(scales)

mrd <- c("ALL5","ALL8","ALL47","ALL64","ALL66","ALL67","ALL68","ALL71")
rel <- c("ALL1","ALL3","ALL4","ALL6", "ALL35","ALL40","ALL42")
excl <- c("ALL2","ALL7","ALL55")

diploids <- read_tsv("results/ALL-latest/analysis/all_diploid_clones-250129.txt")
clone_files <- list.files("results/ALL-latest/", pattern = "*-clones_final.txt", recursive = T, full.names = T)

timepoints_short <- c("d0", "d2", "d3", "d5", "d15", "d29", "rel", "rel2", "rel2+")
timepoint_cols_short <- c("d0"="#fc8d59","d2"="#fee08b","d3"="#d9ef8b","d5"="#91cf60","d15"="#1a9850","d29"="#136b39","rel"="#d73027", "rel2"="#002984", "rel2+"="#002050")

cl <- tibble(
  file = clone_files,
  patient_id = str_extract(basename(file), "ALL\\d+"),
  data = map(file, ~read_tsv(.x))
) %>% select(-file) %>% 
  unnest(data)

cl.annot <- cl %>% left_join(
  mutate(diploids, dna_type="diploid"), by=c("patient_id", "clone"="diploid_clone")
) %>% 
  mutate(dna_type=ifelse(is.na(dna_type), "aneuploid", dna_type))

cl.f <- cl %>% 
  # Filter out diploid clones
  anti_join(
    diploids,
    by = c("patient_id", "clone" = "diploid_clone")
  ) %>%
  # Filter out small clones
  group_by(patient_id, clone) %>%
  filter(n() >= 5, !is.na(clone), ! clone %in% "0") %>%
  ungroup() %>% 
  group_by(patient_id, timepoint) %>% 
  filter(n() >= 5) %>% 
  ungroup

cl2 <- cl.f %>% 
  group_by(patient_id) %>%
  mutate(
    clone_short = dense_rank(clone),
  ) %>% ungroup() %>% 
  mutate(timepoint_short=factor(case_when(timepoint=="diagnosis" ~ "d0",
                                          timepoint=="day 2" ~ "d2",
                                          timepoint=="day 3" ~ "d3",
                                          timepoint=="day 5" ~ "d5",
                                          timepoint=="day 15" ~ "d15",
                                          timepoint=="day 29" ~ "d29",
                                          timepoint %in% c("relapse","relapse 1") ~ "rel",
                                          timepoint=="relapse 2" ~ "rel2",
                                          grepl("relapse 2",timepoint) ~ "rel2+",
                                          TRUE ~ "other"),levels=timepoints_short)) %>% 
  filter(!is.na(timepoint_short)) # drops 1d29

clone_cols <- get_clones_col(cl2$clone_short)

patient_counts <- cl2 %>%
  group_by(patient_id) %>%
  summarize(
    n_cells = n(),
    n_clones = n_distinct(clone_short)
  ) %>%
  mutate(label = paste0(patient_id, "\n(n=", n_cells, "; ", n_clones, " clones)"))

p.rel <- cl2 %>% 
  filter(patient_id %in% rel) %>% 
  # filter(both_pass, rna_blast, cell_phase=="G1", clone_short!=0) %>% 
  ggplot(aes(timepoint_short, fill=as.factor(clone_short))) + geom_bar(position="fill", show.legend = T) + 
  scale_fill_manual(values=clone_cols) +
  scale_y_continuous(breaks=pretty_breaks(1), expand=c(0,0)) +
  labs(y="Frac",x=NULL) +
  theme_dntr(axis_ticks = T) +
  theme(legend.title=element_blank(), strip.text=element_text(size=8)) +
  theme(legend.title=element_blank()) + 
  facet_grid(.~patient_id, scales="free", space="free") +
  facet_grid(.~patient_id, scales="free", space="free",
             labeller = labeller(patient_id = setNames(patient_counts$label[patient_counts$patient_id %in% rel],
                                                       patient_counts$patient_id[patient_counts$patient_id %in% rel])))

p.mrd <- cl2 %>% 
  filter(patient_id %in% mrd) %>% 
  # filter(both_pass, rna_blast, cell_phase=="G1", clone_short!=0) %>% 
  ggplot(aes(timepoint_short, fill=as.factor(clone_short))) + geom_bar(position="fill", show.legend = T) + 
  scale_fill_manual(values=clone_cols) +
  scale_y_continuous(breaks=pretty_breaks(1), expand=c(0,0)) +
  labs(y="Frac",x=NULL) +
  theme_dntr(axis_ticks = T) +
  theme(legend.title=element_blank(), strip.text=element_text(size=8)) +
  theme(legend.title=element_blank()) +
  facet_grid(.~patient_id, scales="free", space="free") +
  facet_grid(.~patient_id, scales="free", space="free",
             labeller = labeller(patient_id = setNames(patient_counts$label[patient_counts$patient_id %in% mrd],
                                                       patient_counts$patient_id[patient_counts$patient_id %in% mrd])))

library(patchwork)
pdf("results/ALL-latest/clone-fractions-250131.pdf", width=9, height=6)
p.mrd / p.rel
dev.off()
