####-------------------- FUN WITH SISAL v2 ---------------------------####

# load libraries
library(plyr)
library(dplyr)
library(tidyverse)

### load and transform data --------------------------------------###
# load SISAL v2 csv files

load_sisal_data_janica_echam <- function(prefix = "", year_start, year_stop){
  prefix = ""
  #path = "/stacywork/ginnyweasley/02_SISAL/SISAL_v2/"
  path = "/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2/"
  composite_link_entity <- read.csv(paste(path, prefix,'composite_link_entity.csv',sep = ''), header = T,stringsAsFactors = F)
  d13C <- read.csv(paste(path, prefix,'d13c.csv',sep='') ,header = T, stringsAsFactors = F)
  d13C <- plyr::rename(d13C, c("iso_std" = "iso_std_d13C"))
  d18O <- read.csv(paste(path, prefix,'d18o.csv', sep =''),header = T, stringsAsFactors = F)
  d18O <- plyr::rename(d18O, c("iso_std" = "iso_std_d18O"))
  dating_lamina <- read.csv(paste(path, prefix,'dating_lamina.csv', sep = ''), header = T, stringsAsFactors = F)
  dating <- read.csv(paste(path, prefix,'dating.csv',sep = ''), header = T, stringsAsFactors = F)
  entity_link_reference <- read.csv(paste(path, prefix,'entity_link_reference.csv', sep = ''), header =T, stringsAsFactors = F)
  entity <- read.csv(paste(path, prefix,'entity.csv', sep = ''), header = T, stringsAsFactors = F)
  gap <- read.csv(paste(path, prefix,'gap.csv', sep = ''), header = T, stringsAsFactors = F)
  hiatus <- read.csv(paste(path, prefix,'hiatus.csv', sep =''), header = T, stringsAsFactors = F)
  notes <- read.csv(paste(path, prefix,'notes.csv', sep = ''), header = T, stringsAsFactors = F)
  original_chronology <- read.csv(paste(path, prefix,'original_chronology.csv', sep = ''), header = T, stringsAsFactors = F)
  
  reference <- read.csv(paste(path, prefix,'reference.csv', sep = ''), header = T, stringsAsFactors = F)
  sample <- read.csv(paste(path, prefix,'sample.csv', sep = ''), header = T, stringsAsFactors = F)
  sisal_chronology <- read.csv(paste(path, prefix,'sisal_chronology.csv', sep = ''), header = T, stringsAsFactors = F)
  site <- read.csv(paste(path, prefix,'site_countries.csv', sep = ''), header = T, stringsAsFactors = F)


  # build SISAL tables
  site_tb <- left_join(site, entity, by = 'site_id') %>% left_join(., entity_link_reference, by = 'entity_id') %>% 
    left_join(., reference, by = 'ref_id') %>% left_join(., notes, by = 'site_id') %>% mutate_at(vars(site_id, entity_id), as.numeric)
  dating_tb <- left_join(dating, entity) %>% group_by(entity_id) %>%mutate(laminar_dated = if_else((entity_id %in% dating_lamina$entity_id), 'yes', 'no')) %>% 
    mutate_at(vars(dating_id, depth_dating, dating_thickness, X14C_correction, corr_age, corr_age_uncert_pos, corr_age_uncert_neg), as.numeric) %>%ungroup()
  sample_tb <- plyr::join_all(list(sample,hiatus, gap, original_chronology, sisal_chronology, d13C, d18O), by = 'sample_id', type = 'full', match = 'all') %>% 
    mutate_at(vars(entity_id, sample_id, sample_thickness, depth_sample, interp_age, interp_age_uncert_pos, interp_age_uncert_neg, copRa_age,
                   copRa_age_uncert_pos, copRa_age_uncert_neg, lin_interp_age, lin_interp_age_uncert_pos, lin_interp_age_uncert_neg, 
                   lin_reg_age, lin_reg_age_uncert_pos, lin_reg_age_uncert_neg, d13C_measurement,
                   d13C_precision, d18O_measurement, d18O_precision), as.numeric)


  # filter for 'from base' dated entities
  entity_from_base <- site_tb %>% filter(depth_ref == 'from base') %>% distinct(entity_id)
  sample_from_base <- sample_tb %>% filter(entity_id %in% entity_from_base$entity_id) %>% 
    select(entity_id,depth_sample) %>% group_by(entity_id) %>% dplyr::summarise(max = max(depth_sample))

  # transform depths for 'from base' dated entities in dating file
  dating_tb_new <- full_join(dating_tb, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_dating, NA_real_)) %>% 
    mutate(depth_dating = if_else(!is.na(depth_conv), depth_conv, depth_dating)) %>%
    select(-depth_conv) %>% arrange(., depth_dating, .by_group = T)

  #transform depths for 'from base' dated entities in sample file
  sample_tb_new <- full_join(sample_tb, sample_from_base, by = 'entity_id') %>% group_by(entity_id) %>% 
    mutate(depth_conv = if_else(entity_id %in% entity_from_base$entity_id, max-depth_sample, NA_real_)) %>% 
    mutate(depth_sample = if_else(!is.na(depth_conv), depth_conv, depth_sample)) %>%
    select(-depth_conv) %>% arrange(., depth_sample, .by_group = T)

  ### entity objects ---------------------------------------------###
  # In the MCensemble file each object name is designed as follows: entity_id - entity_name
  # Each object is a named list containing information (site information, references, notes), the proxy values, the original chronology,
  # the available dating information, hiatus depths and the chronologies and ensembles for the 6 different AM.

  # load entity object and assign name for, e.g. entity_id = 1
  #eID_1 <- get(load('/stacywork/ginnyweasley/02_SISAL/SISAL_v2_CARLA/1-BT-1'))
  eID_1 <- get(load('/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2_CARLA/1-BT-1'))
  


  ## filter dating table
  # filter if status is current and gets all hiatuses out
  dating_tb_filtered <- dating_tb_new %>% filter(entity_status == 'current') %>% mutate_at(vars(corr_age),as.numeric) %>% 
    filter(date_used == 'yes' & date_type != 'Event; hiatus'& date_type != 'Event; actively forming'& date_type != 'Event; start of laminations'& date_type != 'Event; end of laminations')


  only_h <-  dating %>% group_by(entity_id) %>% filter(all(date_type == 'Event; hiatus')) %>% distinct(entity_id)

  # entities not included in sample table
  no_sample_info <- entity %>% filter(!(entity_id %in% sample_tb$entity_id))

  # entities with more than 3 dates
  nr_dates <- dating_tb %>% filter(date_used == 'yes' & date_used != 'Event; hiatus') %>% dplyr::select(entity_id, corr_age) %>% group_by(entity_id) %>%count() %>% filter(n>=3)

  # entities with missing sample depths
  no_depth_sample <- sample_tb %>% group_by(entity_id) %>% dplyr::summarise(depths = if_else(all(is.na(depth_sample)), FALSE, TRUE)) %>% filter(!depths)

  # entities containing only U/Th dates, enough dates, sample depths; 523
  run <- dating_tb_filtered %>% distinct(entity_id) %>%
    filter(entity_id %in% nr_dates$entity_id) %>%
    #filter(!(entity_id %in% no_depth_sample$entity_id))# %>%
    # filter(!(entity_id %in% not_UTh_dates$entity_id)) %>%
    # filter(!(entity_id %in% no_sample_info$entity_id)) %>%
    filter(!(entity_id %in% only_h$entity_id)) %>%
    arrange(., entity_id)


  dating_tb_filtered_2 <- dating_tb_filtered %>% filter(entity_id %in% run$entity_id)
  sample_tb_new_filtered <- sample_tb_new %>% filter(entity_id %in% run$entity_id)
  site_tb_filtered <- site_tb %>% filter(entity_id %in% run$entity_id)
  
  
  entities_mineralogy <- sample_tb_new_filtered %>% filter(interp_age < (year_start) & interp_age > (year_stop)) %>% distinct(entity_id, mineralogy)
  test <- entities_mineralogy %>% group_by(entity_id) %>% count() %>% filter(n>1)
  entities_mineralogy <- sample_tb_new_filtered %>% filter(interp_age < (year_start) & interp_age > (year_stop)) %>% distinct(entity_id, mineralogy) %>% 
    filter(!entity_id %in% test$entity_id)# %>% group_by(entity_id) %>% count()
  
  dating_join <- dating_tb_filtered_2 %>% filter(entity_id %in% entities_mineralogy$entity_id) %>% distinct(entity_id, date_type) %>% group_by(entity_id)
  
  entities_dating <- sample_tb_new_filtered %>% filter(interp_age<(year_start) & interp_age > (year_stop)) %>% distinct(entity_id, age_model_type) %>% right_join(., dating_join, by = "entity_id") %>% distinct(entity_id, date_type, age_model_type)
  
  # all entities that have at least two measurement the past millennium
  # get rid of all those entities with multiple mineralogy (id = 107, 188, 190, 197, 348, 349, 436)
  entities_past1000 <- sample_tb_new_filtered %>% filter(interp_age < (year_start) & interp_age > (year_stop)) %>% group_by(entity_id) %>% count() %>% 
    filter(!entity_id %in% test$entity_id) %>% right_join(., entities_mineralogy, by = "entity_id") %>% filter(n>2)
  entities_past1000_dating <- entities_dating %>% filter(entity_id %in% entities_past1000$entity_id)
  
  sample_min2 <- sample_tb_new_filtered %>% filter(entity_id %in% entities_past1000$entity_id) %>% distinct(entity_id, interp_age, d18O_measurement) %>% filter(interp_age<(year_start) & interp_age > (year_stop))
  dating_all <- sample_tb_new_filtered %>% filter(entity_id %in% entities_past1000$entity_id) %>% distinct(entity_id, interp_age, lin_interp_age, lin_reg_age, Bchron_age, Bacon_age, OxCal_age, copRa_age, StalAge_age, d18O_measurement) %>% filter(interp_age<(1950-year_start) & interp_age > (1950-year_stop))
  dating_all$Bchron_age <- as.numeric(dating_all$Bchron_age)
  dating_all$Bacon_age <- as.numeric(dating_all$Bacon_age)
  dating_all$OxCal_age <- as.numeric(dating_all$OxCal_age)
  dating_all$StalAge_age <- as.numeric(dating_all$StalAge_age)
  site_min2 <- site_tb_filtered %>% filter(entity_id %in% entities_past1000$entity_id) %>% distinct(site_id, entity_id, distance_entrance, cover_thickness) %>% right_join(., sample_min2, by = "entity_id")
  
  site_period <- site_min2 %>% group_by(entity_id) %>% 
    summarise(min_corr_age = round(min(interp_age, na.rm = T), digits = 2),
              max_corr_age = round(max(interp_age, na.rm = T), digits = 2)) %>% 
    mutate(period = max_corr_age -min_corr_age) %>% filter(period > 4000)
  
  site_to_entity <- site_tb_filtered %>% filter(entity_id %in% entities_past1000$entity_id) %>%distinct(site_id, entity_id, distance_entrance, geology, cover_thickness)
  
  #new_data <- site_to_entity %>% filter(!entity_id %in% DATA_past1000$CAVES$entity_info$entity_id) %>% group_by(entity_id)
  #data <- dating_tb_filtered_2 %>% filter(entity_id == 48)
  
  #filter(entity_id %in% site_period$entity_id) %>% 
  return_data <- entities_past1000 %>% right_join(., site_period, by = "entity_id") %>% right_join(site_to_entity,., by = "entity_id") %>% filter(n>10)
  sample_min2_return <- sample_min2 %>% filter(entity_id %in% return_data$entity_id)
  entities_past1000_dating_return <- entities_past1000_dating %>% filter(entity_id %in% return_data$entity_id)
  site_to_entity_return <- site_to_entity %>% filter(entity_id %in% return_data$entity_id)
  dating_all_return <- dating_all %>% filter(entity_id %in% return_data$entity_id)

return(list(return_data, sample_min2_return, entities_past1000_dating_return, site_to_entity_return, dating_all_return))
  
  # remove(composite_link_entity, d13C, d18O, dating, dating_lamina, dating_tb, dating_tb_filtered, dating_tb_filtered_2, dating_tb_new,
  #       eID_1, entities_mineralogy, entities_past1000, entity, entity_from_base, entity_link_reference, file_name, gap, hiatus, no_depth_sample,
  #       no_sample_info, notes, nr_dates, only_h, original_chronology, reference, return_data, run, sample, sample_from_base, sample_min2, sample_tb,
  #       sample_tb_new, sample_tb_new_filtered, sisal_chronology, site_min2, site_period, site_tb, site_tb_filtered, dating_join,
  #       entities_dating, entities_past1000_dating, path, prefix, site, site_to_entity, test)

}











