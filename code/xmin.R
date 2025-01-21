library(poweRlaw)
library(tidyverse)
# this method follows Clauset et al. 20096
by estimating the minimum
# size for which the data follow a power law using by minimizing the K-S statistic.
# i.e., xmin is model-based

# load data
dat = readRDS("data/sizes_full_2.rds") %>% 
  group_by(site, month) %>% 
  sample_n(2000, weight = counts, replace = T)

dat_list = dat %>% group_by(site, month) %>% group_split()

xmin_list = list()

for(i in 1:length(dat_list)){
  powerlaw = conpl$new(dat_list[[i]]$pg_dm)
  xmin_list[[i]] = tibble(xmin_clauset = estimate_xmin(powerlaw)$xmin,
                          month = unique(dat_list[[i]]$month),
                          site = unique(dat_list[[i]]$site))
}

xmins_clauset = bind_rows(xmin_list)

dat_clauset_xmins = readRDS("data/sizes_full_2.rds") %>% 
  left_join(xmins_clauset) %>% 
  group_by(site, month) %>% 
  filter(pg_dm >= xmin_clauset) %>%
  mutate(xmin = xmin_clauset)

saveRDS(dat_clauset_xmins, file = "data/dat_clauset_xmins.rds")

dat_clauset_xmins <- read_rds("data/dat_clauset_xmins.rds")

readRDS("data/sizes_full_2.rds") %>% 
  left_join(xmins_clauset) %>% 
  group_by(site, month) %>% 
  filter(pg_dm >= xmin_clauset) %>% 
  distinct(xmin, xmin_clauset, site, month) %>% 
  ggplot(aes(x = xmin, y = xmin_clauset)) +
  geom_point()


readRDS("data/sizes_full_2.rds") %>% 
  left_join(xmins_clauset) %>% 
  group_by(site, month) %>% 
  mutate(id = cur_group_id()) %>% 
  filter(id == 1) %>% 
  ggplot(aes(x = pg_dm)) + 
  geom_density() +
  scale_x_log10() +
  geom_vline(aes(xintercept = xmin_clauset)) +
  geom_vline(aes(xintercept = xmin))
