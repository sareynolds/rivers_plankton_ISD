library(tidyverse)
library(rstan)
library(brms)
library(tidybayes)

sizes_full_2 <- sizes_full_2 %>% 
  ungroup() %>% 
  mutate(mean_pgdm = mean(pg_dm, na.rm = T)) %>% 
  mutate(pg_dm_s = pg_dm / mean_pgdm) #divided by global mean (mean centered)

sizes_full_2 %>% 
  summarize(mean = mean(pg_dm_s))


#model (standardized)
brm_average_size <- brm(pg_dm_s ~ 1 + river + month + river:month +
                        (1 + river + month + river:month|site),
                        data = sizes_full_2,
                        family = Gamma (link = "log"),
                        prior = c(prior(normal(0, 1), class = "Intercept"),
                                  prior(normal(0, 1), class = "b"),
                                  prior(exponential(1), class = "sd")),
                        chain = 1, iter = 10
                        )

brm_average_size_update <- update(brm_average_size, chains = 4, iter = 2000,
                                  newdata = sizes_full_2, 
                                  data2 = list(mean_pgdm = unique(sizes_full_2$mean_pgdm)))

saveRDS(brm_average_size_update, "models/brm_average_size_update.rds")

pp_check(brm_average_size_update, type = "boxplot") + scale_y_log10()
pp_check(brm_average_size_update) + scale_x_log10()

cond_aver_size <- plot(conditional_effects(brm_average_size_update, effect = "river:month"))

raw_means2b = sizes_full_2 %>% 
  group_by(site, month) %>% 
  reframe(mean_size = mean(pg_dm, na.rm = T))

#posteriors
post_average_size <- brm_average_size_update$data %>% 
  select(-pg_dm_s) %>% 
  distinct() %>% 
  #group_by(site) %>% 
  add_epred_draws(brm_average_size_update, re_formula = NULL) %>% 
  filter(.draw < 1000) %>% 
  ungroup() %>% 
  # mutate(.epred_uns = .epred * mean(sizes_full_2$pg_dm)) %>% 
  left_join(raw_means2b) %>% 
  mutate(.epred_uns = .epred * mean_size)

post_as_qi <- post_average_size %>% 
  group_by(site, month) %>% 
  median_qi(.epred_uns)

max(post_as_qi$.epred_uns)
min(post_as_qi$.epred_uns)

#plot (sites :))
plot_average_size <- post_average_size %>% 
  ggplot(aes(x = site, y = .epred_uns, color = month)) +
  #geom_point(data = sizes_full_2, aes(x = site, y = pg_dm, color = month),
             #position = position_jitter(width = 0.1),
             #position = position_dodge(width = 0.5)) +
  stat_pointinterval(position = position_dodge(width = 0.5)) +
  # scale_y_log10() +
  facet_wrap(~river,
             scales = "free",
             labeller = labeller(river = 
                                   c("big" = "Big Sioux",
                                     "james" = "James",
                                     "mis" = "Missouri",
                                     "verm" = "Vermillion"))) +
  # geom_boxplot(data = sizes_full_2,aes(x = site, y = pg_dm, color = month)) +
  #ylim(20000, 100000)
  theme_classic() +
  scale_color_discrete(name = "Month", labels = c("June", "July")) +
  labs(y = "Dry Mass (pg)",
       x = "River")+
  theme(axis.text = element_text(color = "black", family = "serif", size = 10),
        text = element_text(family = "serif", size = 12, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        legend.title = element_text(color = "black", family = "serif", size = 12),
        legend.text = element_text(color = "black", family = "serif", size = 10),
        strip.background = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 5)))

plot_average_size

ggsave(plot_average_size, file = "plots/plot_average_size_sites.jpg", width = 5.5, height = 4.5, dpi =500)

#plot for rivers
plot_average_size_river <- post_average_size %>% 
  ggplot(aes(x = river, y = .epred_uns, color = month)) +
  #geom_point(data = sizes_full_2, aes(x = river, y = pg_dm, color = month),
             #position = position_dodge(width = 0.5),
             #position = position_jitter(width = 0.1))
  #stat_boxplot(position = position_dodge(width = 1))
  stat_pointinterval(position = position_dodge(width = 0.5)) +
  labs(y = "Dry Mass (pg)",
       x = "River") +
  theme_classic() +
  scale_x_discrete(labels = c('Big Sioux', 'James', 'Missouri', 'Vermillion')) +
  scale_color_discrete(name = "Month", labels = c("June", "July")) +
  theme(axis.text = element_text(color = "black", family = "serif", size = 10),
        text = element_text(family = "serif", size = 12, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        legend.title = element_text(color = "black", family = "serif", size = 12),
        legend.text = element_text(color = "black", family = "serif", size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 5)))

plot_average_size_river

ggsave(plot_average_size_river, file = "plots/plot_average_size_river.jpg", width = 5.5, height = 4.5, dpi =500)

plot_average_size_river_poster <- post_average_size %>% 
  ggplot(aes(x = river, y = .epred_uns, color = month)) +
  stat_pointinterval(position = position_dodge(width = 0.5),
                     point_size = 10) +
  labs(y = "Dry Mass (pg)",
       x = "River") +
  theme_classic() +
  scale_x_discrete(labels = c('Big Sioux', 'James', 'Missouri', 'Vermillion')) +
  scale_color_discrete(name = "Month", labels = c("June", "July")) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 30))
  
plot_average_size_river_poster

ggsave(plot_average_size_river_poster, file = "plots/plot_average_size_river_poster.jpg", width = 9, height = 7.3, dpi =500)

  


post_average_size %>% 
  group_by(river, month) %>% 
  summarize(mean = mean(.epred_uns),
            sd = sd(.epred_uns))

#summarization of models

river_average_size_summary = post_average_size %>% 
  group_by(river, month)  %>% 
  group_by(river, month, .draw) %>% 
  reframe(.epred_uns = mean(.epred_uns)) %>% 
  group_by(river, month) %>%
  summarize(mean = mean(.epred_uns),
            sd = sd(.epred_uns))
  #group_by(river, month) %>% 
  #mean_qi(.epred_uns)
river_average_size_summary

site_average_size_summary = post_average_size %>% 
  #group_by(site, month) %>% 
  group_by(site, month, .draw) %>% 
  reframe(.epred_uns = mean(.epred_uns)) %>% 
  #group_by(site, month) %>% 
  #summarize(mean = mean(.epred_uns),
            #sd = sd(.epred_uns)) 
  group_by(site, month) %>% 
  mean_qi(.epred_uns)

site_average_size_summary

##probability of difference :)
posts_as = brm_average_size_update$data %>% 
  distinct(river, site, month) %>% 
  add_epred_draws(brm_average_size_update, re_formula = NA)

#between rivers
post_average_size %>% 
  group_by(river, .draw) %>% 
  reframe(.epred_uns = mean(.epred_uns)) %>% 
  ungroup %>% 
  pivot_wider(names_from = river, values_from = .epred_uns) %>% 
  #select( .draw, `big`, `james`, `mis`, `verm`) %>% 
  mutate(diffbig_james = `big` - `james`,
         diffbig_mis = `big` - `mis`,
         diffbig_verm = `big` - `verm`,
         diffjames_mis = `james` - `mis`,
         diffjames_verm = `james` - `verm`,
         diffmis_verm = `mis` - `verm`) %>% 
  reframe(probbig_james = sum(diffbig_james > 0) / max(.draw),
          probbig_mis = sum(diffbig_mis > 0) / max(.draw),
          probbig_verm = sum(diffbig_verm > 0) / max(.draw),
          probjames_mis = sum(diffjames_mis > 0) / max(.draw),
          probjames_verm = sum(diffjames_verm > 0) / max(.draw),
          probmis_verm = sum(diffmis_verm > 0) / max(.draw))

#between months (rivers)
post_average_size %>% 
  group_by(river, month, .draw) %>% 
  reframe(.epred_uns = mean(.epred_uns)) %>% 
  ungroup %>% 
  pivot_wider(names_from = month, values_from = .epred_uns) %>% 
  select(river, .draw, `june`, `july`) %>% 
  mutate(diffjuly_june = `july` - `june`) %>% 
  group_by(river) %>% 
  reframe(probjuly_june = sum(diffjuly_june > 0) / max(.draw))


#between sites

#big sioux (sites 11, 12, 13)
post_average_size %>% 
  group_by(site) %>% 
  ungroup() %>% 
  select(-.row, -.chain, -.iteration, -river, -.epred, -mean_size) %>% 
  pivot_wider(names_from = site, values_from = .epred_uns) %>% 
  select(month, .draw, `11`, `12`, `13`) %>% 
  mutate(diff11_12 = `11` - `12`,
         diff11_13 = `11` - `13`,
         diff12_13 = `12` - `13`) %>% 
  group_by(month) %>% 
  reframe(prob11_12 = sum(diff11_12 > 0) / max(.draw),
          prob11_13 = sum(diff11_13 > 0) / max(.draw),
          prob12_13 = sum(diff12_13 > 0) / max(.draw))

#probability difference between vermillion (9, 8, 6) sites
post_average_size %>% 
  group_by(site) %>% 
  ungroup %>% 
  select(-.row, -.chain, -.iteration, -river, -.epred, -mean_size) %>% 
  pivot_wider(names_from = site, values_from = .epred_uns) %>% 
  select(month, .draw, `9`, `8`, `6`) %>% 
  mutate(diff9_8 = `9` - `8`,
         diff9_6 = `9` - `6`,
         diff8_6 = `8` - `6`) %>% 
  group_by(month) %>% 
  reframe(prob9_8 = sum(diff9_8 > 0) / max(.draw),
          prob9_6 = sum(diff9_6 > 0) / max(.draw),
          prob8_6 = sum(diff8_6 > 0) / max(.draw))


#probability difference between james (22, 18, 17, 16, 15) sites 
post_average_size %>% 
  group_by(site) %>% 
  ungroup %>% 
  select(-.row, -.chain, -.iteration, -river, -.epred, -mean_size) %>% 
  pivot_wider(names_from = site, values_from = .epred_uns) %>% 
  select(month, .draw, `22`, `18`, `17`, `16`, `15`) %>% 
  mutate(diff22_18 = `22` - `18`,
         diff22_17 = `22` - `17`,
         diff22_16 = `22` - `16`,
         diff22_15 = `22` - `15`,
         diff18_17 = `18` - `17`,
         diff18_16 = `18` - `16`,
         diff18_15 = `18` - `15`,
         diff17_16 = `17` - `16`,
         diff17_15 = `17` - `15`,
         diff16_15 = `16` - `15`)%>% 
  group_by(month) %>% 
  reframe(#prob22_18 = sum(diff22_18 > 0) / max(.draw),
          #prob22_17 = sum(diff22_17 > 0) / max(.draw),
          #prob22_16 = sum(diff22_16 > 0) / max(.draw),
          #prob22_15 = sum(diff22_15 > 0) / max(.draw),
          #prob18_17 = sum(diff18_17 > 0) / max(.draw),
          #prob18_16 = sum(diff18_16 > 0) / max(.draw),
          #prob18_15 = sum(diff18_15 > 0) / max(.draw),
          #prob17_16 = sum(diff17_16 > 0) / max(.draw))
          prob17_15 = sum(diff17_15 > 0) / max(.draw),
          prob16_15 = sum(diff16_15 > 0) / max(.draw))

#probability difference between missouri (36, 34, 33, 14, 3, 26) sites
post_average_size %>% 
  group_by(site) %>% 
  ungroup %>% 
  select(-.row, -.chain, -.iteration, -river, -.epred, -mean_size) %>% 
  pivot_wider(names_from = site, values_from = .epred_uns) %>% 
  select(month, .draw, `36`, `34`, `33`, `14`, `3`, `26`) %>% 
  mutate(diff36_34 = `36` - `34`,
         diff36_33 = `36` - `33`,
         diff36_14 = `36` - `14`,
         diff36_3 = `36` - `3`,
         diff36_26 = `36` - `26`,
         diff34_33 = `34` - `33`,
         diff34_14 = `34` - `14`,
         diff34_3 = `34` - `3`,
         diff34_26 = `34` - `26`,
         diff33_14 = `33` - `14`,
         diff33_3 = `33` - `3`,
         diff33_26 = `33` - `26`,
         diff14_3 = `14` - `3`,
         diff14_26 = `14` - `26`,
         diff3_26 = `3` - `26`)%>% 
  group_by(month) %>% 
  reframe(#prob36_34 = sum(diff36_34 > 0) / max(.draw),
          #prob36_33 = sum(diff36_33 > 0) / max(.draw),
          #prob36_14 = sum(diff36_14 > 0) / max(.draw),
          #prob36_3 = sum(diff36_3 > 0) / max(.draw),
          #prob36_26 = sum(diff36_26 > 0) / max(.draw),
          #prob34_33 = sum(diff34_33 > 0) / max(.draw),
          #prob34_14 = sum(diff34_14 > 0) / max(.draw))
          prob34_3 = sum(diff34_3 > 0) / max(.draw),
          prob34_26 = sum(diff34_26 > 0) / max(.draw),
          prob33_14 = sum(diff33_14 > 0) / max(.draw),
          prob33_3 = sum(diff33_3 > 0) / max(.draw),
          prob33_26 = sum(diff33_26 > 0) / max(.draw),
          prob14_3 = sum(diff14_3 > 0) / max(.draw),
          prob14_26 = sum(diff14_26 > 0) / max(.draw),
          prob3_26 = sum(diff3_26 > 0) / max(.draw))

#between months (sites)
post_average_size %>% 
  group_by(site, month) %>% 
  ungroup %>% 
  select(-.row, -.chain, -.iteration, -river, -.epred, -mean_size) %>% 
  pivot_wider(names_from = month, values_from = .epred_uns) %>% 
  select(site, .draw, `june`, `july`) %>% 
  mutate(diffjuly_june = `july` - `june`) %>% 
  group_by(site) %>% 
  reframe(probjuly_june = sum(diffjuly_june > 0) / max(.draw))
