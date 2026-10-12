# last edited 6/17/2026
# this script cleans up the PARs data as downloaded from sharepoint.
# the end result is lumped treatments to the permit-year-treatment-chemical level
# there can be two rows of the same permit-year-day if there are mulitple chemicals

library(tidyverse)
library(lubridate)
library(DescTools)

ifelse(!dir.exists("data"), 
       dir.create("data"), 
       "Folder exists already")
ifelse(!dir.exists("outputdata"), 
       dir.create("outputdata"), 
       "Folder exists already")
ifelse(!dir.exists("draftcode"), 
       dir.create("draftcode"), 
       "Folder exists already")

PARs_raw <- read.csv('data/PARs_raw.csv')

PARs_working <- PARs_raw %>% 
  filter(permit_number != "", permit_number != "NULL")

#check how many rows have a range of rates (and therefore would be hard to combine)
## if the number of rows yielded here is not zero, write code to handle range of rates
PARs_working %>% 
  filter(rate_low != rate_high)


dir_map <- c(
  N = 0,
  NNE = 22.5,
  NE = 45,
  ENE = 67.5,
  E = 90,
  ESE = 112.5,
  SE = 135,
  SSE = 157.5,
  S = 180,
  SSW = 202.5,
  SW = 225,
  WSW = 247.5,
  W = 270,
  WNW = 292.5,
  NW = 315,
  NNW = 337.5
)

PARs_working <- PARs_working %>%
  
  # remove columns we don't care about
  select(-c(fk_PAR_ID, fk_Treatment_ID, temp_nos, PAR_comments, rate_high)) %>% 
  
  mutate(across(c(
    treatment_size,
    average_depth,
    PAR_listed_rate = rate_low,
    amount_total,
    water_temp_low,
    water_temp_high,
    wind_speed_low,
    wind_speed_high,
    air_temp), as.numeric)) %>%
# could fix this to calculate based on the rate itself, but need to convert rates out of ppm from patrick. 
  mutate(
    PAR_vol_rate = (amount_total / average_depth / treatment_size),
    PAR_areal_rate = (amount_total / treatment_size), 
    vol_rate_unit = paste0(rate_units, "/acre-ft"), 
    areal_rate_unit = paste0(rate_units, "/acre"),
    PAR_start_date = as.Date(start_date, format = "%m/%d/%Y"),
    PAR_end_date = as.Date(end_date, format = "%m/%d/%Y"),
    wind_dir = dir_map[wind_direction],
    wind_speed = rowMeans(pick(wind_speed_low, wind_speed_high), na.rm = TRUE))
  
# write function for min and max that avoid inf and -inf 
safe_min <- function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

PARs_treatment_group <- PARs_working %>% 
  
  arrange(permit_number, app_year, species_code, pesticide_trade_name, PAR_start_date) %>%
  
  # create treatment groups based on permit-year-pesticides that have close application dates and can be assumed to be the same treatment.
  ### should I add to this code so treatment group is numbered for easy grouping again
  
  group_by(permit_number, app_year, species_code, pesticide_trade_name) %>%
  mutate(
    day_gap = as.numeric(PAR_start_date - lag(PAR_start_date)),
    new_treatment = if_else(is.na(day_gap) | day_gap > 5, 1, 0), # assumes 5 day treatment window is the max.
    treatment_group = cumsum(new_treatment)
  ) %>%
  
  #summarize based on the treatment groups
  group_by(permit_number, app_year, species_code, pesticide_trade_name, amount_units, treatment_group) %>%
  
  summarise(
    start_date = if(all(is.na(PAR_start_date))) {
      as.Date(NA)
      } else {
        min(PAR_start_date, na.rm = TRUE)
        },
    end_date   = if(all(is.na(PAR_end_date))) {
      as.Date(NA)
    } else {
      max(PAR_end_date, na.rm = TRUE)
    },
    
    PAR_treatment_size = sum(treatment_size, na.rm = TRUE),
    
    avg_depth = if (sum(treatment_size, na.rm = TRUE) > 0)
      weighted.mean(average_depth, treatment_size, na.rm = TRUE)
    else NA_real_,
    
    avg_vol_rate = if (sum(average_depth * treatment_size, na.rm = TRUE) > 0)
      weighted.mean(
        PAR_vol_rate,
        w = average_depth * treatment_size,
        na.rm = TRUE
      )
    else NA_real_,
    
    avg_areal_rate = if (sum(treatment_size, na.rm = TRUE) > 0)
      weighted.mean(PAR_areal_rate, treatment_size, na.rm = TRUE)
    else NA_real_,
    
    amount_total = sum(amount_total, na.rm = TRUE),
    
    water_temp_low  = safe_min(water_temp_low),
    water_temp_high = safe_max(water_temp_high),
    
    wind_speed_low  = safe_min(wind_speed_low),
    wind_speed_high = safe_max(wind_speed_high),
    
    wind_direction = names(sort(table(wind_direction), decreasing = TRUE))[1],
    
    air_temp_low  = safe_min(air_temp),
    air_temp_high = safe_max(air_temp),
    
    u = mean(wind_speed * sin(wind_dir * pi / 180), na.rm = TRUE),
    v = mean(wind_speed * cos(wind_dir * pi / 180), na.rm = TRUE)) %>%
    
  ungroup() %>%
  
  mutate(
    app_year = as.character(app_year),
    avg_speed = sqrt(u^2 + v^2),
    avg_dir = (atan2(u, v) * 180 / pi + 360) %% 360
  )
  
deg_to_dir <- function(deg) {
  dirs <- c(
    "N", "NNE", "NE", "ENE",
    "E", "ESE", "SE", "SSE",
    "S", "SSW", "SW", "WSW",
    "W", "WNW", "NW", "NNW"
  )
  
  # there are 16 different wind direction -> 360/16 = 22.5 degrees
  idx <- floor((deg + 11.25) / 22.5) %% 16 + 1
  dirs[idx]
}

# Apply function
PARs_treatment_group <- PARs_treatment_group %>%
  mutate(wind_dir = deg_to_dir(avg_dir)) %>% 
  select(-c(u, v))

# Cleanup
PARs_treatment_group <- PARs_treatment_group %>% 
  mutate(
    avg_depth = round(avg_depth, digits = 1),
    avg_vol_rate = round(avg_vol_rate, digits = 2),
    avg_areal_rate = round(avg_areal_rate, digits = 2),
    amount_total = round(amount_total, digits = 2),
    avg_speed = round(avg_speed, digits = 1)
  )

