library(tidyverse)
library(lubridate)
library(DescTools)
library(arrow)

source(here::here("1_Loading and preparing PI data.R"))

PI_rake_dates <- db_sub_pa %>% 
  mutate(
    SURVEY_START = ymd(SURVEY_START),
    survey_yr = year(SURVEY_START),
    survey_month = month(SURVEY_START)
    , .before = RAKE_MAX)

PI_surveys <- PI_rake_dates %>%
  group_by(DOW, SURVEY_START) %>% 
  slice(1) %>%
  ungroup() %>% 
  dplyr::select(DOW, SURVEY_START, survey_yr, survey_month, RAKE_MAX)

taxonomic_columns <- directory_taxonomic$fieldNames


PI_tax_collapse <- PI_rake_dates %>% 
  group_by(DOW, SURVEY_START) %>% 
  
  mutate(
    station_count = n_distinct(sta_nbr)) %>% #NOTE: is this a good way to determine the # of stations
  
  summarise(
    
    across(
      any_of(taxonomic_columns), 
      ~ sum(.x, na.rm = TRUE)
    ),
    
    across(
      !any_of(taxonomic_columns), 
      ~ first(.x)
    ),
    
    .groups = "drop"
  )

tax_cols <- intersect(taxonomic_columns, names(PI_tax_collapse))

MN_foo <- PI_tax_collapse %>%
  mutate(across(all_of(tax_cols), ~ .x / station_count))
