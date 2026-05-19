# last edited 5/18/2026
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

PARs_Raw <- read.csv('data/Pars_Raw.csv')

PARs_working <- PARs_Raw %>% 
  filter(Permit.number != "", Permit.number != "NULL")

#check how many rows have a range of rates (and therefore would be hard to combine)
## if the number of rows yielded here is not zero, write code to handle range of rates
PARs_working %>% 
  filter(PAR_rate.low != PAR_rate.high)


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
  select(-c(fk_PAR.ID, fk_Treatment.ID, PAR_temp.NOS, PAR_comments, PAR_rate.high)) %>% 
  
  mutate(across(c(
    PAR_treatment.size,
    PAR_average.depth,
    PAR_listed_rate = PAR_rate.low,
    PAR_amount.total,
    PAR_water.temp.low,
    PAR_water.temp.high,
    PAR_wind.speed.low,
    PAR_wind.speed.high,
    PAR_air.temp), as.numeric)) %>%
# could fix this to calculate based on the rate itself, but need to convert rates out of ppm from patrick. 
  mutate(
    PAR_vol_rate = (PAR_amount.total / PAR_average.depth / PAR_treatment.size),
    PAR_areal_rate = (PAR_amount.total / PAR_treatment.size), 
    vol_rate_unit = paste0(PAR_rate.units, "/acre-ft"), 
    areal_rate_unit = paste0(PAR_rate.units, "/acre"),
    PAR_start.date = as.Date(PAR_start.date, format = "%m/%d/%Y"),
    PAR_end.date = as.Date(PAR_end.date, format = "%m/%d/%Y"),
    wind_dir = dir_map[PAR_wind.direction],
    wind_speed = rowMeans(pick(PAR_wind.speed.low, PAR_wind.speed.high), na.rm = TRUE))
  
# write function for min and max that avoid inf and -inf 
safe_min <- function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

PARs_treatment_group <- PARs_working %>% 
  
  arrange(Permit.number, App.year, species.code, Pesticide.trade.name, PAR_start.date) %>%
  
  # create treatment groups based on permit-year-pesticides that have close application dates and can be assumed to be the same treatment.
  ### should I add to this code so treatment group is numbered for easy grouping again
  
  group_by(Permit.number, App.year, species.code, Pesticide.trade.name) %>%
  mutate(
    day_gap = as.numeric(PAR_start.date - lag(PAR_start.date)),
    new_treatment = if_else(is.na(day_gap) | day_gap > 5, 1, 0), # assumes 5 day treatment window is the max.
    treatment_group = cumsum(new_treatment)
  ) %>%
  
  #summarize based on the treatment groups
  group_by(Permit.number, App.year, species.code, Pesticide.trade.name, treatment_group) %>%
  
  summarise(
    start_date = if(all(is.na(PAR_start.date))) {
      as.Date(NA)
      } else {
        min(PAR_start.date, na.rm = TRUE)
        },
    end_date   = if(all(is.na(PAR_end.date))) {
      as.Date(NA)
    } else {
      max(PAR_start.date, na.rm = TRUE)
    },
    
    treatment_size = sum(PAR_treatment.size, na.rm = TRUE),
    
    avg_depth = if (sum(PAR_treatment.size, na.rm = TRUE) > 0)
      weighted.mean(PAR_average.depth, PAR_treatment.size, na.rm = TRUE)
    else NA_real_,
    
    avg_vol_rate = if (sum(PAR_average.depth * PAR_treatment.size, na.rm = TRUE) > 0)
      weighted.mean(
        PAR_vol_rate,
        w = PAR_average.depth * PAR_treatment.size,
        na.rm = TRUE
      )
    else NA_real_,
    
    avg_areal_rate = if (sum(PAR_treatment.size, na.rm = TRUE) > 0)
      weighted.mean(PAR_areal_rate, PAR_treatment.size, na.rm = TRUE)
    else NA_real_,
    
    amount_total = sum(PAR_amount.total, na.rm = TRUE),
    
    water_temp_low  = safe_min(PAR_water.temp.low),
    water_temp_high = safe_max(PAR_water.temp.high),
    
    wind_speed_low  = safe_min(PAR_wind.speed.low),
    wind_speed_high = safe_max(PAR_wind.speed.high),
    
    wind_direction = names(sort(table(PAR_wind.direction), decreasing = TRUE))[1],
    
    air_temp_low  = safe_min(PAR_air.temp),
    air_temp_high = safe_max(PAR_air.temp),
    
    u = mean(wind_speed * sin(wind_dir * pi / 180), na.rm = TRUE),
    v = mean(wind_speed * cos(wind_dir * pi / 180), na.rm = TRUE)) %>%
    
     slice(1) %>%
     ungroup() %>%
  
  mutate(
    App.year = as.character(App.year),
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
  mutate(wind_dir = deg_to_dir(avg_dir))

