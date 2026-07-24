# last updated 7/23/2026
# this script mainly updates the all_tx_clean dataset to have complete data for models

# NOTE: eventually replace this with a call to run the script that produces the csv.
MN_trts <- read.csv('outputdata/all_IAPM_treatments.csv')

no_data <- c(
  0, 0.0,            
  NA, NaN,           
  "NA", "NULL",      
  "", " ",           
  "None", "n/a", "."
)

# Add littoral zone data
lake_attribute_df <- read.csv('data/MN_lake_basin_littoral_zone_15_ft_std.csv', colClasses = c(DOWLKNUM = "character"))

# join lit
# all DOWs will be simplified to 6 digit parent DOW for compatibility 
lz_df <- lake_attribute_df %>%

  dplyr::select(DOWLKNUM, LK_LZACRES) %>%

  mutate(DOW_parent = substr(as.character(DOWLKNUM), 1, 6)) %>%
  
  group_by(DOW_parent) %>%
  
  summarise(across(2, sum, na.rm = TRUE))

# this line won't be needed once I update the original code for all_IAPM_treatments
MN_trts <- MN_trts %>% 
  dplyr::select(-c(lake_LZ_acres, prcnt_LZ))

# load permit data so we can grab the permit approved acres
apm_iapm_permits <- read.csv('data/apm_iapm_permit_detail.csv')

iapm_permits <- apm_iapm_permits %>%
  filter (permit_type == "Invasive Aquatic Plant Management") %>%
  arrange(
    permit_number, 
    desc(treatment_method == "Pesticide Control")
  ) %>%
  # Keep only the first row for each permit number (respects the sort order above)
  distinct(permit_number, .keep_all = TRUE)
  
iapm_permits <- iapm_permits %>%
  pivot_longer(
    cols = starts_with("permit_approved_acres_"),
    names_to = "treatment_year",
    names_prefix = "permit_approved_acres_",
    values_to = "yr_approved_acres",
    values_drop_na = TRUE,
    values_transform = list(
      yr_approved_acres = as.numeric
    ))

MN_trts <- MN_trts %>% 
  mutate(
    treatment_year = as.character(treatment_year)
  )

# join the permit approved acres to the treatment data
MN_trts <- MN_trts %>%
  left_join(
    iapm_permits %>%
      # Pre-select to only bring in the join keys and the new column
      dplyr::select(permit_number, treatment_year, yr_approved_acres),
    by = c("permit_number", "treatment_year")
  )

# MN_trts <- MN_trts %>%
#  left_join(
#    lake_attribute_df %>%
#      select(
#        DOWLKNUM,
#        LK_LZACRES
#      ),
#    by = c("DOW" = "DOWLKNUM")
 # )

MN_trts <- MN_trts %>% 
  mutate(p.area = case_when(
    # If either columns value is in your no_data list
    treatment_size %in% no_data | yr_approved_acres %in% no_data ~ NA_real_,
    
    # Otherwise, calculate percent and round to 4 decimal places
    TRUE ~ round(treatment_size / yr_approved_acres, digits = 4)
    )) %>%
    mutate(DOW_parent = substr(as.character(DOW), 1, 6))

MN_trts <- MN_trts %>% 
  left_join(lz_df, by = c("DOW_parent")
  )

# calculate percent littoral zone treated using actual or permitted treated acres
MN_trts <- MN_trts %>%
  mutate(
    prcnt_LZ = case_when(
      # If lake acreage is missing, we can't do any math -> NA
      LK_LZACRES %in% no_data ~ NA_real_,
      
      # If BOTH treatment_size and yr_approved_acres are missing -> NA
      treatment_size %in% no_data & yr_approved_acres %in% no_data ~ NA_real_,
      
      # If treatment_size is missing (but approved acres has data) -> use approved acres
      treatment_size %in% no_data ~ round(yr_approved_acres / LK_LZACRES, digits = 4),
      
      # Otherwise (treatment_size has data) -> use treatment_size
      TRUE ~ round(treatment_size / LK_LZACRES, digits = 4)
    )
  )

# filter out rows that don't have enough data to calculate % littoral treated
MN_trts_acre <- MN_trts %>% 
  filter(!(prcnt_LZ %in% no_data))
