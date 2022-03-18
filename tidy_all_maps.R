# This code prepares S1 S4 from all sources to tables to parametrize Atlantis
# most of this code just reads in tables produced elsewhere and wrangles them to a format that will go to the GOA_Bioparam_Parameter_Table google sheet
# this brings to my attention many minor inconsistencies between the several methods (mostly just naming, output format, etc.)
# TODO: over time try and tidy those up as much as possible and bring them together

# The goal is to have one place that pulls together all maps we produced in many different places
# this is also where we check that the S values add up to 1 exactly
# TODO: move this check to the relevant code for each species

library(tidyverse)
library(data.table)
library(rbgm)
library(sf)

# get all Atlantis groups
atlantis_fg <- read.csv('../data/GOA_Groups.csv')
atlantis_bgm <- read_bgm('../data/GOA_WGS84_V4_final.bgm')
atlantis_box <- atlantis_bgm %>% box_sf()

# Bottom trawl ------------------------------------------------------------

bt_files <- list.files('../output/s/tables/', full.names = TRUE)

bt_groups <- substr(list.files('../output/s/tables/'),
                    7,(nchar(list.files('../output/s/tables/'))-6)) %>% unique()

bt_s <- rbindlist(lapply(bt_files, function(x){
  fg_stg <- substr(x,26,(nchar(x)-4))
  this <- read.csv(x)
  x1 <- this %>% 
    rename(S1 = S) %>%
    mutate(S2 = S1, S3 = S1, S4 = S1, NameStage = fg_stg) %>% # expand to all seasons assuming there are no differences
    pivot_longer(c(S1:S4), names_to = 'S', values_to = 'Prop') %>%
    mutate(Fg_S = paste(NameStage, S, sep='_')) %>%
    select(box_id, Fg_S, Prop) %>%
    mutate(Prop = replace_na(Prop,0))
    arrange(Fg_S)
  x1
}))

# GOAIERP -----------------------------------------------------------------

# here we also need to fill the empty boxes - they are so small that here we just fill them with the smallest value for simplicity,
# instead of using the value of the closest box (save some code)

goaierp_files <- list.files('../../GOAIERP/output/cpue_tables/', full.names = TRUE)

goaierp_groups <- str_match(goaierp_files, 'cpue_tables/(.*?)_')[,2] %>%
  as.data.frame() %>%
  set_names('Code') %>%
  left_join((atlantis_fg %>% select(Code,Name)), by = 'Code') %>%
  drop_na() %>%
  distinct() %>%
  filter(Name %in% setdiff(Name, bt_groups)) %>% # these minus sharks are the ones we want
  pull(Code) %>%
  setdiff('SHP') # drop sharks

to_keep <- list()
for(i in 1:length(goaierp_groups)){to_keep[[i]] <- grep(goaierp_groups[i], goaierp_files)}

goaierp_files <- goaierp_files[unlist(to_keep)]

tidy_goaierp <- function(x){
  
  fg_stg <- str_match(x, 'cpue_tables/(.*?)_GOAIERP')[,2]
  
  # group and stage
  fg <- str_split(fg_stg,'_',simplify = T)[,1]
  fg <- atlantis_fg %>% filter(Code==fg) %>% pull(Name)
  stg <- str_split(fg_stg,'_',simplify = T)[,2]
  stg <- ifelse(stg=='J', 'J', 'A')
  
  this <- read.csv(x)
  min_cpue <- min(this$all_years_kgkm2, na.rm=T)
  
  x1 <- atlantis_box %>%
    st_set_geometry(NULL) %>%
    select(box_id,area,botz,boundary) %>%
    left_join((this %>% select(box_id,all_years_kgkm2)), by = 'box_id') %>%
    mutate(all_years_kgkm2 = replace_na(all_years_kgkm2,0)) %>% # turn NA's to 0
    rowwise() %>%
    mutate(all_years_kgkm2 = ifelse(all_years_kgkm2==0 & botz<0 & boundary==F, min_cpue, all_years_kgkm2)) %>%
    ungroup() %>%
    mutate(biomass = all_years_kgkm2*area,
           prop = biomass/sum(biomass)) %>%
    select(box_id, prop) %>%
    rename(S1 = prop) %>%
    mutate(S2 = S1, S3 = S1, S4 = S1, NameStage = paste(fg,stg,sep='_')) %>% # expand to all seasons assuming there are no differences
    pivot_longer(c(S1:S4), names_to = 'S', values_to = 'Prop') %>%
    mutate(Fg_S = paste(NameStage, S, sep='_')) %>%
    select(box_id, Fg_S, Prop) %>%
    arrange(Fg_S)
  
  return(x1)
  
}

goaierp_s <- rbindlist(lapply(goaierp_files, tidy_goaierp))

# Non-mammal custom maps --------------------------------------------------

# demersal sharks
SHD <- read.csv('../../Sharks/output/sleep_s1_s4.csv') %>% # demersal (sleeper) shark from IPHC
  mutate(S2 = prop, S3 = prop, S4 = prop, NameStage = 'Shark_demersal_A') %>%
  select(.bx0, NameStage, prop, S2:S4) %>%
  rename(S1 = prop, box_id = .bx0) %>%
  pivot_longer(cols=c(S1:S4), names_to = 'S', values_to = 'Prop') %>%
  mutate(Fg_S = paste(NameStage, S, sep = '_')) %>%
  select(box_id, Fg_S, Prop) %>%
  arrange(Fg_S)

SHD <- rbind(SHD, (SHD %>% mutate(Fg_S = str_replace(Fg_S, '_A_', '_J_')))) # assuming juveniles are same as adults

# pelagic sharks
SHP <- cbind(read.csv('../../Sharks/output/salmon_shark_s1.csv'), # pelagic (salmon shark from Weng et al. 2008)
             read.csv('../../Sharks/output/salmon_shark_s2.csv')[,2],
             read.csv('../../Sharks/output/salmon_shark_s3.csv')[,2],
             read.csv('../../Sharks/output/salmon_shark_s4.csv')[,2]) %>%
  set_names(c('box_id','S1','S2','S3','S4')) %>%
  mutate(NameStage = 'Shark_pelagic_A') %>%
  pivot_longer(cols=c(S1:S4), names_to = 'S', values_to = 'Prop') %>%
  mutate(Fg_S = paste(NameStage, S, sep = '_')) %>%
  select(box_id, Fg_S, Prop) %>%
  arrange(Fg_S)

SHP <- rbind(SHP, (SHP %>% mutate(Fg_S = str_replace(Fg_S, '_A_', '_J_')))) # assuming juveniles are same as adults

# capelin
CAP <- read.csv('../../Capelin/output/capelin_s1_s4.csv') %>% # capelin from McGowan et al. (2020)
  mutate(S2 = Prop, S3 = Prop, S4 = Prop, NameStage = 'Capelin_A') %>%
  select(.bx0, NameStage, Prop, S2:S4) %>%
  rename(S1 = Prop, box_id = .bx0) %>%
  pivot_longer(cols=c(S1:S4), names_to = 'S', values_to = 'Prop') %>%
  mutate(Fg_S = paste(NameStage, S, sep = '_')) %>%
  select(box_id, Fg_S, Prop) %>%
  arrange(Fg_S)

CAP <- rbind(CAP, (CAP %>% mutate(Fg_S = str_replace(Fg_S, '_A_', '_J_')))) # assuming juveniles are same as adults

# sandlance
SAN <- read.csv('../../Sandlance/output/sandlance_s1_s4.csv') %>% # capelin from McGowan et al. (2020)
  mutate(S2 = Prop, S3 = Prop, S4 = Prop, NameStage = 'Sandlance_A') %>%
  select(.bx0, NameStage, Prop, S2:S4) %>%
  rename(S1 = Prop, box_id = .bx0) %>%
  pivot_longer(cols=c(S1:S4), names_to = 'S', values_to = 'Prop') %>%
  mutate(Fg_S = paste(NameStage, S, sep = '_')) %>%
  select(box_id, Fg_S, Prop) %>%
  arrange(Fg_S)

SAN <- rbind(SAN, (SAN %>% mutate(Fg_S = str_replace(Fg_S, '_A_', '_J_')))) # assuming juveniles are same as adults

# Marine mammals ----------------------------------------------------------

# whales and dolphins from Rone et al. (2017)
whale_files <- list.files('../../Whales/output/', full.names = T)

whale_s <- rbindlist(lapply(whale_files, function(x){
  this <- read.csv(fg) %>%
    mutate(S2 = Prop, S3 = Prop, S4 = Prop, NameStage = paste(Atlantis_group,'A',sep='_')) %>%
    select(.bx0, NameStage, Prop, S2:S4) %>%
    rename(S1 = Prop, box_id = .bx0) %>%
    pivot_longer(cols=c(S1:S4), names_to = 'S', values_to = 'Prop') %>%
    mutate(Fg_S = paste(NameStage, S, sep = '_')) %>%
    select(box_id, Fg_S, Prop) %>%
    arrange(Fg_S)
  
  this <- rbind(this, (this %>% mutate(Fg_S = str_replace(Fg_S, '_A_', '_J_')))) # assuming juveniles are same as adults
  this
}))

# pinnipeds
pinnipeds_s <- read.csv('../../Pinnipeds/output/s1_s4.csv') %>%
  mutate(S2 = S, S3 = S, S4 = S, NameStage = paste(Name,'A',sep='_')) %>%
  select(box_id, NameStage, S, S2:S4) %>%
  rename(S1 = S) %>%
  pivot_longer(cols=c(S1:S4), names_to = 'S', values_to = 'Prop') %>%
  mutate(Fg_S = paste(NameStage, S, sep = '_')) %>%
  select(box_id, Fg_S, Prop) %>%
  arrange(Fg_S)

pinnipeds_s <- rbind(pinnipeds_s, (pinnipeds_s %>% mutate(Fg_S = str_replace(Fg_S, '_A_', '_J_')))) # assuming juveniles are same as adults

# SSL
ssl_s1 <- read.csv('../../StellerSeaLions/code/ssl_s1_s4.csv') %>%
  mutate(fg_stg = 'Steller_sea_lion_A', S = 'S1',
         Fg_S = paste(fg_stg, S, sep = '_')) %>%
  select(.bx0,Fg_S,prop) %>%
  rename(box_id = .bx0, Prop = prop)

ssl_s1 <- rbind(ssl_s1, (ssl_s1 %>% mutate(Fg_S = str_replace(Fg_S, '_A_', '_J_')))) # assuming juveniles are same as adults)
ssl_s4 <- ssl_s1 %>% mutate(Fg_S = str_replace(Fg_S, 'S1', 'S4'))

ssl_s2 <- read.csv('../../StellerSeaLions/code/ssl_s2_s3.csv') %>%
  mutate(fg_stg = 'Steller_sea_lion_A', S = 'S2',
         Fg_S = paste(fg_stg, S, sep = '_')) %>%
  select(.bx0,Fg_S,prop) %>%
  rename(box_id = .bx0, Prop = prop)

ssl_s2 <- rbind(ssl_s2, (ssl_s2 %>% mutate(Fg_S = str_replace(Fg_S, '_A_', '_J_')))) # assuming juveniles are same as adults)
ssl_s3 <- ssl_s2 %>% mutate(Fg_S = str_replace(Fg_S, 'S2', 'S3'))

ssl_s <- rbind(ssl_s1, ssl_s2, ssl_s3, ssl_s4)

# Seabirds ----------------------------------------------------------------

# winter from pelagic distributions
winter_files <- list.files('../../Birds/outputs/pelagic_winter/s/', full.names = T)
summer_files <- list.files('../../Birds/outputs/coastal_summer/', full.names = T)

key <- data.frame('old'=c('Diving_fish','Diving_plankton','Surface_fish','Surface_plankton'),
                  'new'=c('Seabird_dive_fish','Seabird_dive_inverts','Seabird_surface_fish','Seabird_surface_inverts'))

birds_s_1_4 <- rbindlist(lapply(winter_files, function(x){
  x <- winter_files[1]
  fg <- str_match(x, '/s/(.*?)_winter')[,2]
  fg <- key %>% filter(old==fg) %>% pull(new)
  
  this <- read.csv(x) %>%
    mutate(Fg_S = paste(fg, 'A', 'S1', sep='_')) %>%
    select(.bx0,Fg_S,Prop) %>%
    rename(box_id = .bx0)
  
  this <- rbind(this, (this %>% mutate(Fg_S = str_replace(Fg_S, '_A_', '_J_')))) # assuming juveniles are same as adults)
  this <- rbind(this, (this %>% mutate(Fg_S = str_replace(Fg_S, 'S1', 'S4'))))
  this  
}))

birds_s_2_3 <- rbindlist(lapply(summer_files, function(x){
  x <- summer_files[1]
  fg <- str_match(x, 'summer/(.*?)_s2')[,2]
  fg <- key %>% filter(old==fg) %>% pull(new)
  
  this <- read.csv(x) %>%
    mutate(Fg_S = paste(fg, 'A', 'S2', sep='_')) %>%
    select(.bx0,Fg_S,Prop) %>%
    rename(box_id = .bx0)
  
  this <- rbind(this, (this %>% mutate(Fg_S = str_replace(Fg_S, '_A_', '_J_')))) # assuming juveniles are same as adults)
  this <- rbind(this, (this %>% mutate(Fg_S = str_replace(Fg_S, 'S2', 'S3'))))
  this  
}))

birds_s <- rbind(birds_s_1_4, birds_s_2_3)

# Plankton ----------------------------------------------------------------

plankton_files <- list.files('../../Plankton_and_nutrients/outputs/s1s4/final_for_parameters/',
                             pattern = '.csv', full.names = T)

nuts <- c('Iron','NH3','NO3') # drop the nuts
t <- list()
for(i in 1:length(nuts)){t[[i]] <- grep(nuts[i],plankton_files)}

plankton_files <- setdiff(plankton_files, plankton_files[unlist(t)])

plankton_s <- rbindlist(lapply(plankton_files, function(x){
  this_fg <- str_match(x, 'parameters/(.*?).csv')[,2]
  this <- read.csv(x) %>%
    set_names(c('fg','box_id','S1','S2','S3','S4')) %>%
    pivot_longer(cols=S1:S4, names_to = 'S', values_to = 'Prop') %>%
    mutate(Fg_S = paste(fg,'A',S,sep='_')) %>%
    select(box_id,Fg_S,Prop)
  this
}))

# Infauna -----------------------------------------------------------------

# make infauna based on the substratum

# check what is missing