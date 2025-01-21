library(rstan)
library(brms)
library(tidyverse)
library(tidybayes)
library(dplyr)
library(isdbayes)
library(ggplot2)

#sizes
sizes_full_2 = readRDS(file = "data/sizes_full.rds") %>%
  mutate(xmax = max(pg_dm),
         counts = 1/samvol_ul) %>% 
  mutate(chl_s = (chl - mean(chl))/sd(chl)) %>% 
  group_by(site, month) %>% 
  mutate(xmin = min(pg_dm))

saveRDS(sizes_full_2, file = "data/sizes_full_2.rds")

#general model
brm_river= brm(pg_dm | vreal(counts, xmin, xmax) ~ 1 + site + month + site:month,
                    data = sizes_full,
                    stanvars = stanvars,  # keep don't change
                    family = paretocounts(),  # keep don't change
                    prior = c(prior(normal(-1.1, 0.5), class = "Intercept"),
                              prior(normal(0, 0.5), class = "b")),
                              #prior(exponential(6), class = "sd")), # regularizing priors
                    chains = 1, iter = 10,  # keep don't change
                    cores = 4, # keep don't change
                    file = 'models/brm_river.rds',  # keep don't change
                    file_refit = "on_change")       # keep don't change

brm_river_update = update(brm_river,  chains = 4, iter = 2000)
saveRDS(brm_river_update, file = "models/brm_river_update.rds")

conditional_effects(brm_river_update)

cond_river = conditional_effects((brm_river_update))


site_river = sizes_full %>% ungroup %>% distinct(site, river)

cond_data = cond_river$`site:month` %>% left_join(site_river) 

#new x_min model (with sites and rivers)
brm_river_xmin= brm(pg_dm | vreal(counts, xmin, xmax) ~ 1 + site + month + site:month + (1|river),
               data = dat_clauset_xmins,
               stanvars = stanvars,  # keep don't change
               family = paretocounts(),  # keep don't change
               prior = c(prior(normal(-1.1, 0.5), class = "Intercept"),
                         prior(normal(0, 0.5), class = "b")),
               #prior(exponential(6), class = "sd")), # regularizing priors
               chains = 1, iter = 10,  # keep don't change
               cores = 4, # keep don't change
               file = 'models/brm_river_xmin.rds',  # keep don't change
               file_refit = "on_change")       # keep don't change

brm_river_xmin_update = update(brm_river_xmin,  chains = 4, iter = 2000)
saveRDS(brm_river_xmin_update, file = "models/brm_river_xmin_update.rds")

cond_xmin = conditional_effects(brm_river_xmin_update)

plot(conditional_effects(brm_river_xmin_update, effects = "site:month")) + 
                   select(site == c('3', '36'))

#new x_min with rivers instead of sites
brm_river_xmin_2= brm(pg_dm | vreal(counts, xmin, xmax) ~ 1 + river + month + river:month + (1|site),
                    data = dat_clauset_xmins,
                    stanvars = stanvars,  # keep don't change
                    family = paretocounts(),  # keep don't change
                    prior = c(prior(normal(-1.1, 0.5), class = "Intercept"),
                              prior(normal(0, 0.5), class = "b")),
                    #prior(exponential(6), class = "sd")), # regularizing priors
                    chains = 1, iter = 10,  # keep don't change
                    cores = 4, # keep don't change
                    file = 'models/brm_river_xmin_2.rds',  # keep don't change
                    file_refit = "on_change")       # keep don't change

brm_river_xmin_2_update = update(brm_river_xmin_2,  chains = 4, iter = 2000)
saveRDS(brm_river_xmin_2_update, file = "models/brm_river_xmin_2_update.rds")

conditional_effects(brm_river_xmin_2_update)


#looking at posteriors :)
posts = dat_clauset_xmins %>% 
  distinct(site, month, river, xmin, xmax) %>% 
  mutate(counts = 1) %>% 
  add_epred_draws(brm_river_xmin_update)


posts %>% 
  ggplot(aes(x = river, y = .epred, color = month, group = interaction(month, site))) + 
  stat_pointinterval() 


posts %>% 
  group_by(site) %>% 
  mutate(median = median(.epred)) %>% 
  ggplot(aes(y = reorder(site, median), x = .epred, color = month, group = interaction(month, site),
             fill = month)) + 
  stat_slab()


#checking posteriors
fit2 <- readRDS("models/brm_river_xmin_update.rds")
pp_check(fit2) + scale_x_log10()
pp_check(brm_river_update) + scale_x_log10()

# chl model
brm_river_chl = brm(pg_dm | vreal(counts, xmin, xmax) ~ 1 + chl_s + 
                 (1|site) + (1|river) + (1|month),
               data = sizes_full,
               stanvars = stanvars,  # keep don't change
               family = paretocounts(),  # keep don't change
               prior = c(prior(normal(-1.1, 0.5), class = "Intercept"),
                         prior(normal(0, 0.5), class = "b"),
               prior(exponential(6), class = "sd")), # regularizing priors
               chains = 1, iter = 10,  # keep don't change
               cores = 4, # keep don't change
               file = 'models/brm_river_chl.rds',  # keep don't change
               file_refit = "on_change")       # keep don't change

brm_river_chl_update = update(brm_river_chl,  chains = 4, iter = 1000)
saveRDS(brm_river_chl_update, file = "models/brm_river_chl_update.rds")
  
cond_chl = plot(conditional_effects(brm_river_chl_update))


#updated xmin for chlorophyll
brm_river_xmin_chl = brm(pg_dm | vreal(counts, xmin, xmax) ~ 1 + chl_s + 
                      (1|site) + (1|river) + (1|month),
                    data = dat_clauset_xmins,
                    stanvars = stanvars,  # keep don't change
                    family = paretocounts(),  # keep don't change
                    prior = c(prior(normal(-1.1, 0.5), class = "Intercept"),
                              prior(normal(0, 0.5), class = "b"),
                              prior(exponential(6), class = "sd")), # regularizing priors
                    chains = 1, iter = 10,  # keep don't change
                    cores = 4, # keep don't change
                    file = 'models/brm_river_xmin_chl.rds',  # keep don't change
                    file_refit = "on_change")       # keep don't change

brm_river_xmin_chl_update = update(brm_river_xmin_chl,  chains = 4, iter = 2000)
saveRDS(brm_river_xmin_chl_update, file = "models/brm_river_xmin_chl_update.rds")

cond_chl_xmin = plot(conditional_effects(brm_river_xmin_chl_update))


# set up grid
post_dots = sizes_full %>% 
  distinct(chl_s, site, river, month, xmin, xmax) %>% 
  mutate(counts = 1) %>% 
  add_epred_draws(brm_river_update, re_formula = NULL)

cond_chl$chl_s$data %>% 
  ggplot(aes(x = chl_s)) +
  geom_lineribbon(aes(y = estimate__, ymin = lower__, ymax = upper__)) +
  stat_pointinterval(data = post_dots, aes(y = .epred))

post_diffs = post_dots %>% 
  # filter(.draw == 54) %>% 
  # filter(site == "33") %>%
  ungroup %>% 
  select(site, .draw, .epred, month) %>% 
  arrange(.draw) %>% 
  pivot_wider(names_from = month, values_from = .epred) %>% 
  mutate(diff = june - july)

post_diffs %>% 
  group_by(site) %>%
  median_qi(diff, na.rm = T) %>% 
  View()


# isd plot
isd_dots = brm_river_update$data %>% 
  group_by(site, month) %>% 
  arrange(-pg_dm) %>% 
  mutate(order = row_number(),
         max = max(order),
         order = order/max)

isd_dots %>% 
  ggplot(aes(x = pg_dm, y =order,color =month)) +
  geom_point(aes(size = pg_dm)) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~site)
