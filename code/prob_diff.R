library(tidyverse)
library(tidybayes)


#organizing posteriors
posts_site = brm_river_xmin_update$data %>% 
  mutate(xmin = min(xmin),
         xmax = max(xmax),
         counts = 1) %>% 
  distinct(river, site, month, xmin, xmax, counts) %>% 
  add_epred_draws(brm_river_xmin_update, re_formula = NA)

#probability difference between big sioux (11, 12, 13) sites
posts_site %>% 
  group_by(site) %>% 
  ungroup %>% 
  select(-.row, -.chain, -.iteration, -river) %>% 
  pivot_wider(names_from = site, values_from = .epred) %>% 
  select(month, .draw, `11`, `12`, `13`) %>% 
  mutate(diff11_12 = `11` - `12`,
         diff11_13 = `11` - `13`,
         diff12_13 = `12` - `13`) %>% 
  group_by(month) %>% 
  reframe(prob11_12 = sum(diff11_12 > 0) / max(.draw),
          prob11_13 = sum(diff11_13 > 0) / max(.draw),
          prob12_13 = sum(diff12_13 > 0) / max(.draw))

#probability difference between vermillion (9, 8, 6) sites
posts_site %>% 
  group_by(site) %>% 
  ungroup %>% 
  select(-.row, -.chain, -.iteration, -river) %>% 
  pivot_wider(names_from = site, values_from = .epred) %>% 
  select(month, .draw, `9`, `8`, `6`) %>% 
  mutate(diff9_8 = `9` - `8`,
         diff9_6 = `9` - `6`,
         diff8_6 = `8` - `6`) %>% 
  group_by(month) %>% 
  reframe(prob9_8 = sum(diff9_8 > 0) / max(.draw),
          prob9_6 = sum(diff9_6 > 0) / max(.draw),
          prob8_6 = sum(diff8_6 > 0) / max(.draw))

#probability difference between james (22, 18, 17, 16, 15) sites 
posts_site %>% 
  group_by(site) %>% 
  ungroup %>% 
  select(-.row, -.chain, -.iteration, -river) %>% 
  pivot_wider(names_from = site, values_from = .epred) %>% 
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
  reframe(prob22_18 = sum(diff22_18 > 0) / max(.draw),
          prob22_17 = sum(diff22_17 > 0) / max(.draw),
          prob22_16 = sum(diff22_16 > 0) / max(.draw),
          prob22_15 = sum(diff22_15 > 0) / max(.draw),
          prob18_17 = sum(diff18_17 > 0) / max(.draw),
          prob18_16 = sum(diff18_16 > 0) / max(.draw),
          prob18_15 = sum(diff18_15 > 0) / max(.draw),
          prob17_16 = sum(diff17_16 > 0) / max(.draw))
          #prob17_15 = sum(diff17_15 > 0) / max(.draw),
          #prob16_15 = sum(diff16_15 > 0) / max(.draw))

#probability difference between missouri (36, 34, 33, 14, 3, 26) sites
#start here
posts_site %>% 
  group_by(site) %>% 
  ungroup %>% 
  select(-.row, -.chain, -.iteration, -river) %>% 
  pivot_wider(names_from = site, values_from = .epred) %>% 
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
  reframe(prob36_34 = sum(diff36_34 > 0) / max(.draw),
    prob36_33 = sum(diff36_33 > 0) / max(.draw),
    prob36_14 = sum(diff36_14 > 0) / max(.draw),
    prob36_3 = sum(diff36_3 > 0) / max(.draw),
    prob36_26 = sum(diff36_26 > 0) / max(.draw),
    prob34_33 = sum(diff34_33 > 0) / max(.draw),
    prob34_14 = sum(diff34_14 > 0) / max(.draw)
    #prob34_3 = sum(diff34_3 > 0) / max(.draw),
    #prob34_26 = sum(diff34_26 > 0) / max(.draw),
    #prob33_14 = sum(diff33_14 > 0) / max(.draw),
    #prob33_3 = sum(diff33_3 > 0) / max(.draw),
    #prob33_26 = sum(diff33_26 > 0) / max(.draw),
    #prob14_3 = sum(diff14_3 > 0) / max(.draw),
    #prob14_26 = sum(diff14_26 > 0) / max(.draw),
    #prob3_26 = sum(diff3_26 > 0) / max(.draw))

#between months
#posteriors
posts_months = brm_river_xmin_update$data %>% 
  mutate(xmin = min(xmin),
         xmax = max(xmax),
         counts = 1) %>% 
  distinct(river, site, month, xmin, xmax, counts) %>% 
  add_epred_draws(brm_river_xmin_update, re_formula = NA)

#probability difference off all month between months
posts_months %>% 
  group_by(river, month, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  ungroup %>% 
  pivot_wider(names_from = month, values_from = .epred) %>% 
  select(river, .draw, `june`, `july`) %>% 
  mutate(diffjuly_june = `july` - `june`) %>% 
  group_by(river) %>% 
  reframe(probjuly_june = sum(diffjuly_june > 0) / max(.draw))

#Between river averages
posts_months %>% 
  group_by(river, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  ungroup %>% 
  pivot_wider(names_from = river, values_from = .epred) %>% 
  #select( .draw, `big`, `james`, `mis`, `verm`) %>% 
  mutate(diffbig_james = `big` - `james`,
         diffbig_mis = `big` - `mis`,
         diffbig_verm = `big` - `verm`,
         diffjames_mis = `james` - `mis`,
         diffjames_verm = `james` - `verm`,
         diffmis_verm = `mis` - `verm`) %>% 
  reframe(probbig_james = sum(diffbig_james < 0) / max(.draw),
          probbig_mis = sum(diffbig_mis < 0) / max(.draw),
          probbig_verm = sum(diffbig_verm < 0) / max(.draw),
          probjames_mis = sum(diffjames_mis < 0) / max(.draw),
          probjames_verm = sum(diffjames_verm < 0) / max(.draw),
          probmis_verm = sum(diffmis_verm < 0) / max(.draw))



