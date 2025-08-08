library(rstan)
library(brms)
library(tidyverse)
library(tidybayes)
library(isdbayes)


#isd model
brm_river_xmin_update <- readRDS("models/brm_river_xmin_update.rds")

##getting posteriors
# set up grid
posterior = brm_river_xmin_update$data %>% 
  distinct(river, site, month, xmin, xmax) %>% 
  mutate(counts = 1) %>% 
  add_epred_draws(brm_river_xmin_update, re_formula = NULL)

#Average lambda for each river
posterior_rivers = posterior %>% 
  group_by(river, month) %>% 
  mutate(xmin = min(xmin),
         xmax = max(xmax)) %>% 
  group_by(river, month, xmin, xmax, .draw) %>% 
  reframe(.epred = mean(.epred))
  
post_weighted <- posterior %>% 
  group_by(river, .draw, month) %>% 
  reframe(.epred = mean(.epred)) %>% 
  ggplot(aes(x = river, y = .epred, fill = month)) +
  stat_halfeye() +
  theme_ggdist() +
  stat_pointinterval(data = posterior, aes(group = interaction(month, site), color = month))

ggsave(post_weighted, file = "plots/post_weighted.jpg", width = 5.5, height =3.5, dpi =500)

#custom order of river sites in data frames for plots
posterior$site = factor(posterior$site, levels = c('11', '12', '13', '22', '18', '17', '16', '15', '9', '8', '6', '36', '34', '33', '14', '3', '26'))

posterior$month <- factor(posterior$month, levels = c('june', 'july'))

posterior %>% 
  arrange(site)

sizes_full_2$site = factor(sizes_full_2$site, levels = c('11', '12', '13', '22', '18', '17', '16', '15', '9', '8', '6', '36', '34', '33', '14', '3', '26'))

sizes_full_2 %>% 
  arrange(site)

sizes_full_2$month <- factor(sizes_full_2$month, levels = c('june', 'july'))

#rivers plot using posterior
river_plot <- posterior %>% 
  ggplot(aes(x = river, y = .epred, color = month)) +
  stat_pointinterval(position = "dodge") +
  theme_classic() +
  labs(y = "Lambda (λ)",
       x = "River") +
  scale_x_discrete(labels = c('Big Sioux', 'James', 'Missouri', 'Vermillion')) +
  scale_color_discrete(name = "Month", labels = c("June", "July")) +
  #geom_hline(yintercept = -2, linetype = "dashed", color = "darkgray") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.text = element_text(color = "black", family = "serif", size = 10),
        text = element_text(family = "serif", size = 12, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        legend.title = element_tex(color = "black", family = "serif", size = 12),
        legend.text = element_text(color = "black", family = "serif", size = 10))

river_plot

ggsave(river_plot, file = "plots/river_plot.jpg", width = 5.5, height = 4.5, dpi =500)

#summary 
mean_lambda <- posterior %>%
  group_by(counts) %>% 
  summarise(mean = mean(.epred, na.rm = TRUE),
            sd = sd(.epred, na.rm = TRUE))

mean_lambda

#site plot using posteriors
site_plot <- posterior %>% 
  ggplot(aes(x = site, y = .epred, color = month)) +
  facet_wrap(~river,
             scales = "free",
             labeller = labeller(river = 
                                   c("big" = "Big Sioux",
                                     "james" = "James",
                                     "mis" = "Missouri",
                                     "verm" = "Vermillion"))) +
  stat_pointinterval(position = "dodge") +
  theme_classic() +
  labs(y = "Lambda (λ)",
       x = "Site") +
  ylim(-3, -1) +
  scale_color_discrete(name = "Month", labels = c("June", "July")) +
  theme(axis.text = element_text(color = "black", family = "serif", size = 10),
        text = element_text(family = "serif", size = 12, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        legend.title = element_text(color = "black", family = "serif", size = 12),
        legend.text = element_text(color = "black", family = "serif", size = 10),
        strip.background = element_blank()) +
  #geom_hline(yintercept = -2, linetype = "dashed", color = "darkgray")+
  guides(color = guide_legend(override.aes = list(size = 5)))

site_plot

ggsave(site_plot, file = "plots/site_plot.jpg", width = 5.5, height = 4.5, dpi =500)


# lambda table (sites)
cond_site = plot(conditional_effects(brm_river_xmin_update,
                                       effects = "site:month"))

cond_site_plot = cond_site$`site:month` +
  theme_ggdist() +
  labs(y = "Lambda (λ)",
       x = "Site")+
  scale_color_discrete(name = "Month", labels = c("July", "June"))+
  guides(fill = "none")+
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 10,
                                 color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.line=element_line(linewidth = 2))

cond_site_plot$data %>% 
  select("site", "month", "estimate__", "se__", "lower__", "upper__")


chl_s = (chl - mean(chl))/sd(chl))

#chlorophyll plot (and unstandardizing :()

chl_unstandarizing = sizes_full_2 %>% 
  select(river, chl, pg_dm) %>% 
  mutate(chl_sd = sd(chl),
         chl_mean = mean(chl)) %>% 
  select(-chl) %>% 
  distinct()

post_chl = brm_river_xmin_chl_update$data %>% 
  select(-pg_dm) %>% 
  distinct() %>% 
  add_epred_draws(brm_river_xmin_chl_update, re_formula = NULL) %>% 
  left_join(chl_unstandarizing) %>% 
  mutate(chl = (chl_s * chl_sd) + chl_mean)

cond_chl_xmin$chl_s$data %>% 
  left_join(chl_unstandarizing)
  

post_chl %>% 
  group_by(chl) %>% 
  mean_qi(.epred)

plot_chl = post_chl %>% 
  ggplot(aes(x = chl_s, y = .epred)) +
  geom_lineribbon(data = cond_chl_xmin$chl$data,
                  aes(y = estimate__,
                      ymin = lower__,
                      ymax = upper__),
                  color = "black") +
  theme_classic() +
  labs(y = "Lambda (λ)",
       x = "Chlorophyll (μg/L)") +
  theme(axis.text = element_text(color = "black", family = "serif", size = 10),
        text = element_text(family = "serif", size = 12, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"))

plot_chl

ggsave(plot_chl, file = "plots/plot_chl.jpg", width = 5.5, height = 4.5, dpi = 500)

#summary of chl at each river
chl %>% 
  group_by(river, month) %>% 
  summarize(mean = mean(chl),
            sd = sd(chl))

#lambda table (for rivers :))
cond_river = plot(conditional_effects(brm_river_xmin_update,
                                      effects = "site:month"))

cond_river_plot = cond_river$`site:month` +
  theme_ggdist() +
  labs(y = "Lambda (λ)",
       x = "Site")+
  scale_color_discrete(name = "Month", labels = c("July", "June"))+
  guides(fill = "none")+
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 10,
                                 color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.line=element_line(linewidth = 2))

cond_river_plot$data %>% 
  select("river", "month", "estimate__", "se__", "lower__", "upper__")

# simulate body sizes -----------------------------------------------------

post_preds = brm_river_xmin_update$data %>% 
  distinct(river, site, month, xmax) %>% 
  mutate(xmin = min(brm_river_xmin$data$xmin)) %>% 
  mutate(counts = 1) %>% 
  add_predicted_draws(brm_river_xmin, re_formula = NULL)

post_preds %>%
  group_by(site) %>% 
  mutate(median = median(.prediction)) %>% 
  ggplot(aes(x = reorder(site,median), y = .prediction, color = month)) +
  geom_jitter(width = 0.1, size = 0.6) +
  facet_wrap(~month)

#mass verse site
mass_site <- sizes_full_2 %>% 
  group_by(site, month) %>% 
  mutate(median = median(pg_dm)) %>% 
  ggplot(aes(x= reorder(site, median), y = pg_dm/100000, color = month)) + 
  labs(x = "Site",
       y = "Dry Mass (pg) x 100,000") +
  theme_ggdist() +
  geom_jitter(shape = 1, size = 2, width = 0.2) +
  facet_wrap(month ~ river, scales = "free_x") +
  theme(axis.title = element_text(size = 10),
        strip.text = element_blank()) +
  scale_x_discrete(limits=c("36","34","33", "14", "3", "26"))
ggsave(mass_site, file = "plots/mass_site.jpg", width = 5.5, height = 4.5, dpi =500)

mass_site

#mass verse site and month
sizes_full_2 %>% 
  group_by(site) %>% 
  mutate(median = median(pg_dm)) %>% 
  ggplot(aes(x= reorder(site, median), y = pg_dm, color = month)) + 
  geom_jitter(shape = 1, size = 0.2, width = 0.2) +
  facet_wrap(river~month) 

#simulated data (comparison of lambdas)
sim_data = tibble(site_17 = rparetocounts(n = 1000, lambda = -2, xmin = min(sizes_full_2$pg_dm),
              xmax = max(sizes_full_2$pg_dm)),
              site_18 = rparetocounts(n = 1000, lambda = -1.5, xmin = min(sizes_full_2$pg_dm),
                                      xmax = max(sizes_full_2$pg_dm))) 

sim_data %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(x = name, y = value)) + 
  geom_jitter() +
  scale_y_log10()

#simulated slope -2
sim_slope2 <- sim_data %>% 
  arrange(-site_17) %>% 
  mutate(order = 1:nrow(.)) %>% 
  ggplot(aes(x = site_17, y = order)) + 
  geom_point(aes(size = site_17)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth() +
  theme_ggdist() +
  labs(y = "Abundance",
       x = "Body Sizes (pg)") +
  theme(axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"))
sim_slope2

ggsave(sim_slope2, file = "plots/sim_slope2.jpg", width = 5.5, height =3.5, dpi =500)


#simulated slope -1.5
sim_slope1.5 <- sim_data %>% 
  arrange(-site_18) %>% 
  mutate(order = 1:nrow(.)) %>% 
  ggplot(aes(x = site_18, y = order)) + 
  geom_point(aes(size = site_18)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth() +
  theme_ggdist() +
  labs(y = "Abundance",
       x = "Body Sizes (pg)") +
  theme(axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"))

ggsave(sim_slope1.5, file = "plots/sim_slope1.5.jpg", width = 3, height = 3, dpi =500)

#chlorophyl isd

brm_river_xmin_chl_update <- readRDS("models/brm_river_xmin_chl_update.rds")

cond_chl_xmin = plot(conditional_effects(brm_river_xmin_chl_update))

cond_chl_plot = cond_chl_xmin$chl_s +
  theme_ggdist() +
  labs(y = "Lambda (λ)",
       x = "Standardized Chlorophyll")

ggsave(cond_chl_plot, file = "plots/cond_chl_plot.jpg", width = 5.5, height = 4, dpi = 500)

#chlorophyll plot try with posteriors (stuck on type of stat)
posterior_chl = brm_river_chl_update$data %>% 
  distinct(river, site, month, chl_s, xmin, xmax) %>% 
  mutate(counts = 1) %>% 
  add_epred_draws(brm_river_xmin_chl_update, re_formula = NULL)

posterior_chl %>% 
  ggplot(aes(x = chl_s, y = .epred)) +
  stat_intervalh() +
  theme_ggdist()




site_plot <- posterior %>% 
  ggplot(aes(x = site, y = .epred, color = month)) +
  facet_wrap(~river,
             scales = "free",
             labeller = labeller(river = 
                                   c("big" = "Big Sioux",
                                     "james" = "James",
                                     "mis" = "Missouri",
                                     "verm" = "Vermillion"))) +
  stat_pointinterval(position = "dodge") +
  theme_ggdist() +
  labs(y = "Lambda (λ)",
       x = "Site") +
  ylim(-3, -1) +
  scale_color_discrete(name = "Month", labels = c("June", "July")) +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 13,
                                 color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.line=element_line(linewidth = 2),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  #geom_hline(yintercept = -2, linetype = "dashed", color = "darkgray")+
  guides(color = guide_legend(override.aes = list(size = 5)))

