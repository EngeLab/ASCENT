# Overlaps
cellid_to_plateid <- function(x){
  sub("([A-Z]{3}[0-9]{5}).*","\\1", x)
}

files <- list.files("~/scps/fq_stats/overlaps",full.names = T)
names(files) <- sub(".*/(.*)\\.overlaps.txt","\\1",files)

d <- map_dfr(files, ~read_tsv(.x, show_col_types=F, col_names = c("frag_len","overlap")), .id="cell")

d2 <- d %>% mutate(plate_id=cellid_to_plateid(cell))

r <- d2 %>% group_by(cell, plate_id) %>% 
  summarize(avg_frag_len=mean(frag_len), avg_overlap=mean(overlap[overlap>0]), 
            n=n(), n_overlap=sum(overlap>0), overlap_frac=n_overlap/n)

r <- d2 %>% group_by(cell, plate_id) %>% 
  summarize(avg_frag_len=mean(frag_len), avg_overlap=mean(overlap), 
            n=n(), n_overlap=sum(overlap>0), overlap_frac=n_overlap/n)

r %>% 
  filter(n>50000) %>% 
  ggplot(aes(avg_frag)) +
  geom_histogram(bins=100) +
  facet_wrap(~plate_id, ncol=1)

# Single plate
d2 %>% 
  filter(plate_id=="VZA02601") %>%
  slice_sample(n=1000) %>%
  ggplot(aes(frag_len, overlap)) +
  geom_point(alpha=0.2, size=.8)

range <- c(0,500)
p1 <- d2 %>% 
  filter(plate_id=="VZA02601") %>%
  slice_sample(n=10000) %>%
  ggplot(aes(frag_len, overlap)) +
  geom_density2d_filled(show.legend = F) +
  labs(caption="sample 10k reads per plate", title="VZA02601") + 
  lims(x=range) +
  theme_minimal()

p2 <- d2 %>% 
  filter(plate_id=="VZA02603") %>%
  slice_sample(n=10000) %>%
  ggplot(aes(frag_len, overlap)) +
  geom_density2d_filled(show.legend = F) +
  labs(caption="sample 10k reads per plate", title="VZA02602") + 
  lims(x=range) +
  theme_minimal()


p1 <- d2 %>% 
  ggplot(aes(frag_len)) + geom_histogram(bins=100) + facet_wrap(~plate_id, ncol=1, scales="free_y") +
  lims(x=c(0,300))

p2 <- d2 %>% 
  filter(overlap>0) %>% 
  ggplot(aes(overlap)) + geom_histogram(bins=30) + facet_wrap(~plate_id, ncol=1, scales="free_y")

p1 + p2

p1 <- r %>% 
  filter(n>50000) %>% 
  # ggplot(aes(avg_frag, n_overlap_frac, color=plate_id)) +
  ggplot(aes(avg_frag_len, avg_overlap, color=plate_id)) +
  geom_point() + 
  guides(color=guide_legend(override.aes=list(size=4)))

p2 <- r %>% 
  filter(n>50000) %>% 
  ggplot(aes(avg_frag_len, overlap_frac, color=plate_id)) +
  # ggplot(aes(avg_frag_len, avg_overlap, color=plate_id)) +
  geom_point() +
  guides(color=guide_legend(override.aes=list(size=4)))
