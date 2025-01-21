library(flowCore)
library(tidyverse)

filenames_fcs <- list.files("data/", pattern = "*.fcs", recursive = T)
filenames_fcs = filenames_fcs[grep("CY5", filenames_fcs, invert = TRUE)]

data_fcs = NULL

for(i in seq_along(filenames_fcs)){
  data_fcs[[i]] = read.FCS(paste0("data/", filenames_fcs[i]))
}

get_sizes = function(x){
  tibble(micron = x@exprs[,1]) %>% 
    mutate(file = x@description$GUID)
}

library(dplyr)

size = bind_rows(lapply(data_fcs, FUN = get_sizes))

sizes = size %>% 
  separate(file, c('site', 'river', 'date', 'sample', 'method')) %>% 
  mutate(volume = (1/6*pi*micron^3),
         log10pg_dm = (1.64*log10(volume)^0.82),
         pg_dm = 10^log10pg_dm) %>% 
  subset(method == "BF")

#input of chl data!
library(readxl)
chl <- read_excel("chl.xlsx")

chl$date <- as.character(chl$date)
chl$site <- as.character(chl$site)

#input of sample volume data!
sam_vol <- read_excel("sam_vol.xlsx") %>% 
  select(site, month, river, samvol_ul)

sam_vol$site <- as.character(sam_vol$site)

#joining sizes with chl and sample volume
sizes_full = sizes %>% 
  left_join(chl) %>%
  left_join(sam_vol)

#saving sizes_full
saveRDS(sizes_full, file = "data/sizes_full.rds")

#generic histograms!
library(ggplot2)

ggplot(sizes_full, aes(x = pg_dm)) + 
  geom_histogram() +
  scale_x_log10() +
  facet_wrap(~site + month)

#number of individuals
sizes %>% 
  group_by(site, river) %>% 
  summarise_at(vars("counts"), sum)

ggplot(sizes_full, aes(x = micron))+
  geom_histogram()
