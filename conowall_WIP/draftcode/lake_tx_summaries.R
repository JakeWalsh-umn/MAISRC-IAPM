# last updated 8/10/26
# this script creates lake-year summaries of treatment activities that are not iapm-chemical
# these datasets will be used for added context for anaylsis and PI charter 
# NOTE: I did not filter based on the type of plants permitted because it was too hard to determine a line to draw
###### fix the script to account for that, if necessary.

library(tidyverse)
library(purrr)

ifelse(!dir.exists("data"), 
       dir.create("data"), 
       "Folder exists already")
ifelse(!dir.exists("outputdata"), 
       dir.create("outputdata"), 
       "Folder exists already")
ifelse(!dir.exists("draftcode"), 
       dir.create("draftcode"), 
       "Folder exists already")

apm_iapm_permits <- read.csv('data/apm_iapm_permit_detail.csv')

apm_iapm_permits <- apm_iapm_permits %>%
  filter(
    !is.na(water_resource_numbers),
    str_detect(str_trim(water_resource_numbers), "^\\d+$")
  ) %>%
  mutate(DOW_parent = str_replace(water_resource_numbers, "\\d{2}$", "00"))

# Add littoral zone data
lake_attribute_df <- read.csv('data/MN_lake_basin_littoral_zone_15_ft_std.csv', colClasses = c(DOWLKNUM = "character"))

# join lit
# all DOWs will be simplified to 6 digit DOW for compatibility 
lz_df <- lake_attribute_df %>%
  
  dplyr::select(DOWLKNUM, LK_LZACRES) %>%
  
  mutate(DOW_parent = str_replace(DOWLKNUM, "\\d{2}$", "00")) %>%
  
  group_by(DOW_parent) %>%
  
  summarise(across(2, sum, na.rm = TRUE))



# lake-year summaries for APM-chemical ----

apm_chem_permits <- apm_iapm_permits %>%
  filter(treatment_method == "Pesticide Control" & permit_type == "Aquatic Plant Management") %>%
  #filter(str_detect(species_groups, "Submersed Plants")) %>%
  mutate(across(everything(), str_trim)) %>%
  dplyr::select(-c(gid, geom_type, the_geom, restoration_type, multi_party, attachment_types)) 

apm_chem_permits <- apm_chem_permits %>%
  pivot_longer(
    cols = starts_with("permit_approved_acres_"),
    names_to = "treatment_year",
    names_prefix = "permit_approved_acres_",
    values_to = "yr_approved_acres",
    values_drop_na = TRUE,
    values_transform = list(
      yr_approved_acres = as.numeric
    )) %>% 
  distinct()

apm_chem_permits <- apm_chem_permits %>%
  group_by(permit_number, treatment_year) %>%
  arrange(desc(yr_approved_acres), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()


apm_chem_lksum <- apm_chem_permits %>%
  dplyr::select(c(lake_name = water_resource_names, DOW_parent, permit_type,
                  treatment_method, species_groups, treatment_year, yr_approved_acres)) %>%
  group_by(DOW_parent, treatment_year) %>%
  summarise(
    species_groups = species_groups %>%
      na.omit() %>%
      str_c(collapse = "; ") %>%          # Combine all group rows together
      str_split_1("[;,/]\\s*") %>%        # Split on semicolons (;), commas (,), or slashes (/)
      str_squish() %>%                   # Clean up extra spaces
      str_to_lower() %>%                 # Convert to lowercase so "Submersed" matches "submersed"
      .[. != ""] %>%                     # Remove empty strings
      unique() %>%                       # Deduplicate individual items
      str_to_sentence() %>%              # Format nicely (e.g., "Filamentous algae")
      str_flatten(collapse = ", ") %>%   # Collapse back into a clean comma list
      na_if(""),
    yr_approved_acres = sum(yr_approved_acres, na.rm = TRUE),
    lake_name = first(lake_name),
    treatment_method = first(treatment_method),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  distinct()

table(apm_chem_lksum$species_groups)

apm_chem_lksum <- apm_chem_lksum %>%
  filter(yr_approved_acres > 0)

apm_chem_lksum <- apm_chem_lksum %>% 
  left_join(lz_df, by = c("DOW_parent")
  )

apm_chem_lksum <- apm_chem_lksum %>%
  mutate(
    prcnt_LZ = round(yr_approved_acres / LK_LZACRES, digits = 5)
    )

# lake-year summaries for APM-mechanical ----

apm_mech_permits <- apm_iapm_permits %>%
  filter((treatment_method == "Mechanical Control" | treatment_method == "Automated Aquatic Plant Control Device (AAPCD)") 
         & permit_type == "Aquatic Plant Management") %>%
  #filter(str_detect(species_groups, "Submersed Plants")) %>%
  mutate(across(everything(), str_trim)) %>%
  dplyr::select(-c(gid, geom_type, the_geom, restoration_type, multi_party, attachment_types)) 

apm_mech_permits <- apm_mech_permits %>%
  pivot_longer(
    cols = starts_with("permit_approved_acres_"),
    names_to = "treatment_year",
    names_prefix = "permit_approved_acres_",
    values_to = "yr_approved_acres",
    values_drop_na = TRUE,
    values_transform = list(
      yr_approved_acres = as.numeric
    )) %>% 
  distinct()

apm_mech_permits <- apm_mech_permits %>%
  group_by(permit_number, treatment_year) %>%
  arrange(desc(yr_approved_acres), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

apm_mech_lksum <- apm_mech_permits %>%
  dplyr::select(c(lake_name = water_resource_names, DOW_parent, permit_type,
                  treatment_method, species_groups, treatment_year, yr_approved_acres)) %>%
  group_by(DOW_parent, treatment_year) %>%
  summarise(
    species_groups = species_groups %>%
      na.omit() %>%
      str_c(collapse = "; ") %>%          # Combine all group rows together
      str_split_1("[;,/]\\s*") %>%        # Split on semicolons (;), commas (,), or slashes (/)
      str_squish() %>%                   # Clean up extra spaces
      str_to_lower() %>%                 # Convert to lowercase so "Submersed" matches "submersed"
      .[. != ""] %>%                     # Remove empty strings
      unique() %>%                       # Deduplicate individual items
      str_to_sentence() %>%              # Format nicely (e.g., "Filamentous algae")
      str_flatten(collapse = ", ") %>%   # Collapse back into a clean comma list
      na_if(""),
    yr_approved_acres = sum(yr_approved_acres, na.rm = TRUE),
    lake_name = first(lake_name),
    treatment_method = first(treatment_method),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  distinct()

apm_mech_lksum <- apm_mech_lksum %>%
  filter(yr_approved_acres > 0)

apm_mech_lksum <- apm_mech_lksum %>% 
  left_join(lz_df, by = c("DOW_parent")
  )

apm_mech_lksum <- apm_mech_lksum %>%
  mutate(
    prcnt_LZ = round(yr_approved_acres / LK_LZACRES, digits = 5)
  )

# lake-year summaries for IAPM-mechanical ----


iapm_mech_permits <- apm_iapm_permits %>%
  filter((treatment_method == "Mechanical Control" | treatment_method == "Automated Aquatic Plant Control Device (AAPCD)") 
         & permit_type == "Invasive Aquatic Plant Management") %>%
  #filter(str_detect(species_groups, "Submersed Plants")) %>% #need to fix this line for species rather than group
  mutate(across(everything(), str_trim)) %>%
  dplyr::select(-c(gid, geom_type, the_geom, restoration_type, multi_party, attachment_types)) 

iapm_mech_permits <- iapm_mech_permits %>%
  pivot_longer(
    cols = starts_with("permit_approved_acres_"),
    names_to = "treatment_year",
    names_prefix = "permit_approved_acres_",
    values_to = "yr_approved_acres",
    values_drop_na = TRUE,
    values_transform = list(
      yr_approved_acres = as.numeric
    )) %>% 
  distinct()

iapm_mech_permits <- iapm_mech_permits %>%
  group_by(permit_number, treatment_year) %>%
  arrange(desc(yr_approved_acres), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()


iapm_mech_lksum <- iapm_mech_permits %>%
  dplyr::select(c(lake_name = water_resource_names, DOW_parent, permit_type,
                  treatment_method, species, treatment_year, yr_approved_acres)) %>%
  group_by(DOW_parent, treatment_year) %>%
  summarise(
    species = species %>%
      na.omit() %>%
      str_c(collapse = "; ") %>%          # Combine all group rows together
      str_split_1("[;,/]\\s*") %>%        # Split on semicolons (;), commas (,), or slashes (/)
      str_squish() %>%                   # Clean up extra spaces
      str_to_lower() %>%                 # Convert to lowercase so "Submersed" matches "submersed"
      .[. != ""] %>%                     # Remove empty strings
      unique() %>%                       # Deduplicate individual items
      str_to_sentence() %>%              # Format nicely (e.g., "Filamentous algae")
      str_flatten(collapse = ", ") %>%   # Collapse back into a clean comma list
      na_if(""),
    yr_approved_acres = sum(yr_approved_acres, na.rm = TRUE),
    lake_name = first(lake_name),
    treatment_method = first(treatment_method),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  distinct()

iapm_mech_lksum <- iapm_mech_lksum %>%
  filter(yr_approved_acres > 0)

iapm_mech_lksum <- iapm_mech_lksum %>% 
  left_join(lz_df, by = c("DOW_parent")
  )

iapm_mech_lksum <- iapm_mech_lksum %>%
  mutate(
    prcnt_LZ = round(yr_approved_acres / LK_LZACRES, digits = 5)
  )

# invasive species reference table where each possible IAPM species has a code
species_reference <- tibble::tribble(
  ~species_name,                ~code,
  "Brittle Naiad",               "BN",
  "Curly-leaf Pondweed",         "CLP",
  "Eurasian Watermilfoil",       "EWM",
  "Flowering Rush",              "FR",
  "Java waterdropwort",          "JW",
  "Purple Loosestrife",          "PL",
  "Starry Stonewort",            "SS",
  "Water Hyacinth",              "WH",
  "Yellow Iris",                 "YI",
  "Zebra Mussels",               "ZM",
  "non-native Phragmites spp.",  "PHR",
  "non-native water lily",       "NWL",
  "NULL",                        "UNKNOWN"
)

ref_species_clean <- species_reference$species_name %>% 
  str_to_lower() %>% 
  str_remove_all("[\\s-]")

iapm_mech_lksum <- iapm_mech_lksum %>%
  mutate(
    # Split on BOTH commas and semicolons (with any surrounding whitespace)
    species_list = str_split(species, "\\s*[,;]\\s*"),
    
    # Compare list to reference table
    species_code = map_chr(species_list, function(x) {
      x_clean <- x %>% 
        str_to_lower() %>% 
        str_remove_all("[\\s-]")
      
      matched_codes <- species_reference$code[match(x_clean, ref_species_clean)]
      matched_codes <- coalesce(matched_codes, "NA")
      
      # Re-combine into a clean comma-separated string
      paste(matched_codes, collapse = ", ")
    })
  ) %>%
  dplyr::select(-c(species_list, species))

# IAPM chem permit history -----

iapm_chem_permits <- apm_iapm_permits %>% 
  filter(treatment_method == "Pesticide Control" & permit_type == "Invasive Aquatic Plant Management") %>%
  #filter(str_detect(species_groups, "Submersed Plants")) %>%
  mutate(across(everything(), str_trim)) %>%
  dplyr::select(-c(gid, geom_type, the_geom, restoration_type, multi_party, attachment_types)) 

iapm_chem_permits <- iapm_chem_permits %>%
  pivot_longer(
    cols = starts_with("permit_approved_acres_"),
    names_to = "treatment_year",
    names_prefix = "permit_approved_acres_",
    values_to = "yr_approved_acres",
    values_drop_na = TRUE,
    values_transform = list(
      yr_approved_acres = as.numeric
    )) %>% 
  distinct()

iapm_chem_permits <- iapm_chem_permits %>%
  group_by(permit_number, treatment_year) %>%
  arrange(desc(yr_approved_acres), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()


iapm_chem_permits <- iapm_chem_permits %>%
  dplyr::select(c(permit_number, lake_name = water_resource_names, DOW_parent, permit_type,
                  treatment_method, species, treatment_year, permit_effective_date,
                  permit_expiration_date, yr_approved_acres))

iapm_chem_permits <- iapm_chem_permits %>%
  filter(yr_approved_acres > 0)

iapm_chem_permits <- iapm_chem_permits %>% 
  left_join(lz_df, by = c("DOW_parent")
  )

iapm_chem_permits <- iapm_chem_permits %>%
  mutate(
    prcnt_LZ = round(yr_approved_acres / LK_LZACRES, digits = 5)
  )


iapm_chem_permits <- iapm_chem_permits %>%
  mutate(
    # Split on BOTH commas and semicolons (with any surrounding whitespace)
    species_list = str_split(species, "\\s*[,;]\\s*"),
    
    # Compare list to reference table
    species_code = map_chr(species_list, function(x) {
      x_clean <- x %>% 
        str_to_lower() %>% 
        str_remove_all("[\\s-]")
      
      matched_codes <- species_reference$code[match(x_clean, ref_species_clean)]
      matched_codes <- coalesce(matched_codes, "NA")
      
      # Re-combine into a clean comma-separated string
      paste(matched_codes, collapse = ", ")
    }),
  ) %>%
  dplyr::select(-c(species_list, species))

apm_lksum <- bind_rows(apm_chem_lksum, apm_mech_lksum)


apm_chem_lksum %>%
  write_csv("outputdata/apm_chem_lksum.csv")

apm_mech_lksum %>%
  write_csv("outputdata/apm_mech_lksum.csv")

iapm_mech_lksum %>%
  write_csv("outputdata/iapm_mech_lksum.csv")

iapm_chem_permits %>% 
  write_csv("outputdata/iapm_chem_permits.csv")

apm_lksum %>% 
  write_csv("outputdata/apm_lksum.csv")
