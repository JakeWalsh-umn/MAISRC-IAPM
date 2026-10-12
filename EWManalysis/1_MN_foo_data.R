library(tidyverse)
library(lubridate)
library(DescTools)
library(arrow)

source(here::here("1_Loading and preparing PI data.R"))

natives <- taxa_xwalk %>%
  filter(ORIGIN == "native") %>%
  pull(fieldNames)

PI_rake_dates <- db_sub_pa %>% 
  mutate(
    SURVEY_START = ymd(SURVEY_START),
    survey_yr = year(SURVEY_START),
    survey_month = month(SURVEY_START),
    depth_ft = as.numeric(depth_ft),
    native.rich = rowSums(across(any_of(natives)) == 1, na.rm = TRUE)
    , .before = RAKE_MAX)

taxonomic_columns <- directory_taxonomic$fieldNames

PI_rake_dates <- PI_rake_dates %>% 
  mutate(SURVEY_START = as.Date(SURVEY_START)) %>% 
  group_by(DOW, SURVEY_START) %>% 
  summarise(
    station_count = n(),
    
    # Average depth across survey stations
    depth_ft = mean(depth_ft, na.rm = TRUE),
    
    avg.rich = mean(native.rich, na.rm = TRUE), 
    
    # Count of survey stations where at least one native plant was found
    nat.ab = sum(native.rich > 0, na.rm = TRUE),
    
    # Sum plant columns only
    across(any_of(taxonomic_columns), ~ sum(.x, na.rm = TRUE)),
    
    # Pass first value for remaining non-taxonomic metadata (excluding native.rich)
    across(
      !any_of(c(taxonomic_columns, "depth_ft", "native.rich")), 
      ~ first(.x)
    ),
    
    .groups = "drop"
  ) %>% 
  mutate(
    # Count unique native species detected anywhere in the lake survey
    total.rich = rowSums(across(any_of(natives)) > 0, na.rm = TRUE),
    #NOTE: add an average native abundance here. 
    .after = nat.ab
  ) 

PI_rake_dates <- PI_rake_dates %>% 
  mutate(
    lake.rich = total.rich/station_count
  )

# I don't think we need this chunk of code
#PI_surveys <- PI_rake_dates %>%
#  group_by(DOW, SURVEY_START) %>% 
#  slice(1) %>%
#  ungroup() %>% 
#  dplyr::select(DOW, SURVEY_START, survey_yr, survey_month, RAKE_MAX)


tax_cols <- intersect(taxonomic_columns, names(PI_rake_dates))

MN_foo <- PI_rake_dates %>%
  mutate(across(all_of(tax_cols), ~ .x /station_count))


MN_foo <- MN_foo %>%
  mutate(DOW_parent = str_replace(DOW, "\\d{2}$", "00"))

# add lake acreage and average depth to the PI surveys
# add lake acreage and average depth to the infested acres dataset
# fix the step in the analysis script that adds uninvaded/invaded and treated/invaded and untreated
## make it so the untreated is looking for lakes that don't have IAPM permits. 
### also make it so the untreated lakes that don't fit the range of depth and acreage is just filtered out. 
#### add some kind of filter for funky outliers. 


