# last edit was 5/7/2026
# this script is for inspecting the DOW IDs in the APM permits dataset

# load packages
library(tidyverse)
library(purrr)
library (sf)
library(stringdist)


# load the minnesota hydrography data to have a "bank" of all DOW ids
## download the shapefile here: https://gisdata.mn.gov/dataset/water-dnr-hydrography
mn_hydrography <- st_read('data/shp_water_dnr_hydrography/dnr_hydro_features_all.shp')
mn_hydrography <- st_drop_geometry(mn_hydrography)
waterbodies_DOW <- mn_hydrography %>% 
  filter(!is.na(dowlknum))

# this provides list of permits without DOW that can be inspected and addressed 
apm_iapm_permits_dow_missing <- permits_wo_duplicates %>%
  filter(!(water_resource_numbers %in% waterbodies_DOW$dowlknum))

# Standardize DOWs IDs.
apm_iapm_permits <- apm_iapm_permits %>% mutate(water_resource_numbers = as.character(water_resource_numbers))
mn_hydrography <- mn_hydrography %>% mutate(dowlknum = as.character(dowlknum))

# Join by DOW ID and then determine a similarity score for the lake names
validation_results <- apm_iapm_permits %>%
  left_join(
    mn_hydrography %>% select(dowlknum, lake_name_hydro = pw_basin_n, county_hydro = cty_name),
    by = c("water_resource_numbers" = "dowlknum")
  ) %>%
  mutate(
    noise_words = "\\blake\\b|\\blk\\b",
    
    name_p = str_squish(str_remove_all(tolower(water_resource_names), noise_words)),
    name_h = str_squish(str_remove_all(tolower(lake_name_hydro), noise_words)),
    
    # Normalize strings for comparison (lowercase and trimmed)
    cnty_p = str_trim(tolower(counties)),
    cnty_h = str_trim(tolower(county_hydro)),
    
    # Calculate fuzzy similarity
    lv_sim = if_else(
      !is.na(name_p) | !is.na(name_h),
      1 - (stringdist(name_p, name_h, method = "lv") /
             pmax(nchar(name_p), nchar(name_h), na.rm = TRUE)),
      NA_real_
    ),
    
    is_contained = map2_lgl(name_p, name_h, ~{
      str_detect(.y, fixed(.x)) | str_detect(.x, fixed(.y))
    }),
    
    name_sim = if_else(is_contained, 1.0, lv_sim),
    
    # Check if the permit county string is found within the hydro county string
    #fixed() treats the input as a string
    county_contained = map2_lgl(cnty_p, cnty_h, ~str_detect(.y, fixed(.x))),
    
    # assign status based on similarity score. Edit as needed.
    status = case_when(
      is.na(lake_name_hydro) ~ "ID Not Found",
      name_sim > 0.8 & county_contained ~ "Verified",
      name_sim > 0.8 & !county_contained ~ "Name Matches, County Discrepancy",
      name_sim > 0.2 ~ "Potential Name Match - Check Manually",
      TRUE ~ "Review Required"
    )
  ) %>%
  # remove the columns that were used for processing
  select(-name_p, -name_h, -cnty_p, -cnty_h)

# export as csv for easier inspection. 
validation_results %>% 
  st_drop_geometry() %>%
  write_csv("validation_results.csv")
