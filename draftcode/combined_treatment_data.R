# This is the second version of combined_treatment_data.R
# last updated 6/17/2026

library(tidyverse)
library(lubridate)
library(DescTools)
library(stringr)
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

# STEP 1: load data 
iapm_annual_survey <- read.csv('data/iapm_annual_survey_all.csv')
apm_iapm_permits <- read.csv('data/apm_iapm_permit_detail.csv')
chem_ref_table <- read.csv('data/chem_ref_table.csv')

# STEP 2
# remove the columns we don't care about
iapm_annual_survey <- iapm_annual_survey %>%
  select(-c(data_date_time, when_permit_expires, starts_with("aapcd"),
            satisified_with_swim_itch_control)) %>% 
  rename(treatment_year = calander_year) %>%
  mutate(
    # clean up insignificant things in survey acreage like extra spaces and periods
    total_treated_area_acres = str_remove_all(total_treated_area_acres, "\\s+"),
    
    # cleans up multiple decimals in a row
    total_treated_area_acres = str_replace_all(total_treated_area_acres, "\\.+", "."),
    
    # cleans up extra decimals between numbers 
    total_treated_area_acres = str_replace(total_treated_area_acres, "^([^.]*\\.[^.]*)\\.", "\\1"),
    
    total_treated_area_acres = as.numeric(total_treated_area_acres)
  )


# STEP 3
# left join by permit number and keep the permit DOW and treatment_method
# NOTE: the treatment method from permit may not be accurate until we get the new data
iapm_annual_survey <- iapm_annual_survey %>%
  left_join(apm_iapm_permits %>% select(DOW = water_resource_numbers, lk_name = water_resource_names, permit_number, 
                                        tx_method_permit = treatment_method), by = "permit_number")

# STEP 4
# figure out if there are surveys without permit details 
surveys_wo_permit_data <- iapm_annual_survey %>% anti_join(apm_iapm_permits, by = "permit_number")

# count how many surveys without permit details
length(unique(surveys_wo_permit_data$permit_number)) # there are 219 permits with survey responses but no permit data. 

#remove the survey responses that don't have corresponding permit details
iapm_annual_survey <- iapm_annual_survey %>% anti_join(surveys_wo_permit_data, by = "apm_annual_survey_id")


# STEP 5
# remove rows where the key treatment details are exact duplicates
# this treats commercial and landowner rows as the same, and is only looking at treatment data.
iapm_annual_survey <- iapm_annual_survey %>%
  arrange(apm_annual_survey_type_id != 2) %>% 
  distinct(
    permit_number, 
    total_treated_area_acres, 
    treatment_year, 
    treatment_dates,
    pick(starts_with("work_done_"), starts_with("hp")), 
    .keep_all = TRUE
  )

# for simplicity throughout process, define what represents "no data"
no_data <- c(
  0, 0.0,            
  NA, NaN,           
  "NA", "NULL",      
  "", " ",           
  "None", "n/a", "."
)


# STEP 6
# add a count of how many months work was done in 
iapm_annual_survey <- iapm_annual_survey %>%
  mutate(
    work_month_count = rowSums(across(
      .cols = starts_with("work_done_in_"),
      .fns = ~ !(.x %in% no_data | is.na(.x))
    ))
  )

# STEP 7
# Add the following columns to the survey data table 

## tx_confirmed (Y/N) 
#### -> general indicator of treatment history

## tx_months (Y/N) 
#### -> broad idea of when treatment happened. Could be used to pair with treatment notifications.

## tx_dates (Y/N) 
#### -> specific date of when treatment happened

## multiple_tx (Y/N) 
#### -> tells us what we can assume about the chemical and acreage data (summed across treatments?)

iapm_annual_survey <- iapm_annual_survey %>%
  mutate(
    
    # define if we think the permit was used
    tx_confirmed = case_when(
      
      str_detect(permit_used, "No") ~ "N",
      
      if_any(
        .cols = c(total_treated_area_acres, starts_with("work_done_in_"), 
                  starts_with("hp_"), treatment_dates),
        .fns = ~ !(.x %in% no_data | is.na(.x))
      ) ~ "Y",
      
      TRUE ~ "N"),
    
    # define if we know the month(s) of treatments
    tx_months = case_when(
      if_any(
        .cols = c(treatment_dates, starts_with("work_done_in")),
        .fns = ~ !(.x %in% no_data | is.na (.x))
      ) ~ "Y",
      
      TRUE ~ "N"),
    
    # define if we have specific dates for a treatment
    tx_dates = if_else(!(treatment_dates %in% no_data | is.na(treatment_dates)), "Y", "N"),
    
    # define if there are multiple treatments for a survey
    multiple_tx = if_else(
      work_month_count > 1 | str_count(treatment_dates, ",") >= 1, "Y", "N")
  )


# STEP 8
# add together all chemical amount data for easy reference of whether there is any chemical data
iapm_annual_survey <- iapm_annual_survey %>%
  mutate(
    sum_chem = rowSums(across(
      .cols = starts_with("hp"),
      .fns = ~ as.numeric(.x) # Ensures values are treated as numbers
    ), na.rm = TRUE)
  )

# STEP 9
iapm_annual_survey <- iapm_annual_survey %>% mutate(treatment_id = paste0("S_", row_number()))

# STEP 10a
# add chemical delimited columns with the original name (unit) and active ingredients
chem_lookup <- setNames(chem_ref_table$active_ing, chem_ref_table$orig_name)

# STEP 10b
# create a delimited list of chemicals(unit)
# use chem_lookup to also make a delimited list of active ingredients. 
iapm_annual_survey <- iapm_annual_survey %>%
  rowwise() %>%
  mutate(
    res = {
      hp_cols <- names(pick(starts_with("hp")))
      hp_vals <- c_across(starts_with("hp"))
      active_cols <- hp_cols[!(hp_vals %in% no_data) & !is.na(hp_vals)]
      
      if (length(active_cols) > 0) {
        # pull out chemical names and units from the hp column names
        chem_names <- str_match(active_cols, "^[^_]+_(.*)_[^_]+$")[, 2]
        unit_names <- str_extract(active_cols, "[^_]+$")
        
        # Look up the active ingredients using chem_lookup
        act_ing_vector <- unique(chem_lookup[chem_names])
        act_ing_vector <- act_ing_vector[!is.na(act_ing_vector)]
        
        # Create the full chem (unit) pairs
        bracket_pairs <- paste0(chem_names, " (", unit_names, ")")
        
        
        tibble(
          chemicals_units = paste(unique(bracket_pairs), collapse = ", "),
          active_ing      = paste(act_ing_vector, collapse = ", ")
        )
      } else {
        tibble(
          chemicals_units = "", 
          active_ing      = ""
        )
      }
    }
  ) %>%
  ungroup() %>%
  unpack(res)


# STEP 11
# set aside subset of the data to build out treatment details dataset
survey_tx_details <- iapm_annual_survey %>% 
  
  select(
    permit_id, 
    treatment_id, 
    treatment_dates, 
    total_treated_area_acres, 
    starts_with("hp"),
    chemicals_units,
    active_ing
  )

# STEP 12
# clean up survey dataset
surveys_narrow <- iapm_annual_survey %>%
  select(-c(starts_with("hp")))

# STEP 13
# pivot the month columns to long form so each row represents a survey-month
surveys_long <- surveys_narrow %>%
  pivot_longer(
    cols = starts_with("work_done_in_"), 
    values_to = "tx_month",           
    values_drop_na = FALSE
  ) %>%
  select(-name)

# STEP 14
# remove duplicate rows within different types of survey responses
# AND
# remove extra rows that were created in the pivot_longer
# work_month_count = 0, tx_confirmed = "Y" keep the row
# work_month_count >= 1 remove rows where treatment_range has no data
# handle rows where there are specific dates of treatment rather than months
# NOTE: parsed_date and collapsed string are redundant - fix 

# # STEP 14a: subset rule1, rule2, and rule3 data
rule1_data <- surveys_long %>%
  filter(work_month_count == 0 & tx_dates == "Y") %>%
  
  group_by(permit_number, treatment_dates) %>%
  
  filter(
    # check for groups with at least one row with valid chemical data
    if (any(!chemicals_units %in% no_data)) {
      # remove any rows without valid chemical data
      !chemicals_units %in% no_data
    } else {
      # If the group contains ONLY no_data in chemical columns, keep them all
      TRUE
    }
  ) %>%
  
  ungroup() %>% 
  
  mutate(
    original_string = treatment_dates,
    string_length = nchar(treatment_dates)) %>%
  
  # split the dates into individual rows so we can evaluate them individually
  separate_longer_delim(treatment_dates, delim = ",") %>%
  mutate(treatment_dates = str_squish(treatment_dates)) %>%
  
  mutate(
    parsed_date = as.Date(parse_date_time(treatment_dates, orders = c("mdy", "ymd", "dmy")))) %>%
  
 # for each treatment_id, only keep rows with distinct dates
  distinct(treatment_id, parsed_date, .keep_all = TRUE) %>%
  
  mutate(
    date_count = str_count(original_string, "\\d+[/-]\\d+[/-]\\d+"),
    date_count = coalesce(date_count, 0)) %>% 
  
  # collapse dates that are part of the same survey AND are within 5 day treatment window
  group_by(treatment_id) %>% 
  arrange(parsed_date, .by_group = TRUE) %>% 
  mutate(
    group_anchor = if(n() > 1) {
      accumulate(parsed_date, ~ if (.y - .x <= 5) .x else .y)} else {
      parsed_date},
    
    day_gap = as.numeric(parsed_date - group_anchor),
    new_treatment = if_else(parsed_date == group_anchor, 1, 0),
    treatment_group = cumsum(new_treatment)) %>% 
  
  select(-group_anchor) %>% 
  
  # group by treatment_groups which represent the 5 day windows 
  # collapse dates into list
  group_by(treatment_id, treatment_group) %>%
  
  mutate(
    collapsed_string = paste(as.character(parsed_date), collapse = ", ")
  ) %>%
  
  # only keep one row, because above steps created as many rows as dates in a treatment group
  filter(row_number() == 1) %>% 
  
  mutate(parsed_date = collapsed_string) %>%
  ungroup() %>%
  
  # arrange by the length of the date strings
  arrange(string_length) %>% 
  
  # keep only one permit-chemical-date combo
  # NOTE: this step could prioritize other variables in the arrange step above
  distinct(permit_number, chemicals_units, parsed_date, .keep_all = TRUE) %>% 
  
  mutate(dup_clean_id = paste0("rule1"))

# handle data where the months of treatment are defined
rule2_data <- surveys_long %>%
  filter(
    work_month_count >= 1,
    !is.na(tx_month),               
    !tx_month %in% no_data) %>%
  
  group_by(permit_number, treatment_year, tx_month) %>%
  filter(
    # Check if any row in the group has a valid chemical recorded
    if (any(!chemicals_units %in% no_data)) {
      # If real chemicals exist, drop any row where the chemical is in the no_data list
      !chemicals_units %in% no_data
    } else {
      # If the group only contains chemical values in your no_data list, keep them all
      TRUE
    }
  ) %>%
  
  ungroup() %>%
  
  arrange(
    # Low work_month_count more likely to have better quality tx details
    work_month_count,
    
    # apm_annual_survey_type_id = 2 (commercial) comes BEFORE = 1
    apm_annual_survey_type_id != 2,
    
    # Highest total_treated_area_acres (descending order)
    desc(total_treated_area_acres)
  ) %>%
  
  # Group by permit-year-month-chemical and slice the first row.
  distinct(permit_number, treatment_year, tx_month, chemicals_units, .keep_all = TRUE) %>%

  
  mutate(dup_clean_id = paste0("rule2"))

# handle data where the survey implies a treatment happened but there is no date or month data
rule3_data <- surveys_long %>%
  filter(work_month_count == 0 & tx_dates == "N") %>%  
  group_by(permit_number, treatment_year) %>%
  
  # sort by indicators of better data and slice to only keep one row per group 
  arrange(
    sum_chem %in% no_data, # FALSE (valid data) comes before TRUE (no-data), 
    apm_annual_survey_type_id != 2, # commercial surveys before landowner
    desc(total_treated_area_acres), # highest acres prioritized
    .by_group = TRUE                # Ensures sorting happens within our groups
  ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(dup_clean_id = paste0("rule3"))

# STEP 14b:add treatment start and end dates using the date strings (rule 1) and months (rule 2)

rule1_data <- rule1_data %>%
  select(-any_of(c("original_string", "string_length", "date_count", "day_gap", "new_treatment"))) %>%
  
  mutate(
    # Use a flexible regex to pull ALL dates out of the string as a clean list
    extracted_dates = str_extract_all(collapsed_string, "\\d{4}[/-]\\d{2}[/-]\\d{2}"),
    
    # Convert the lists of text strings into lists of date objects
    date_objects = map(extracted_dates, ymd),
    
    # extract the true minimum and maximum chronologically from each list
    start_date = map_vec(date_objects, ~ if(length(.x) > 0) min(.x, na.rm = TRUE) else as.Date(NA)),
    end_date   = map_vec(date_objects, ~ if(length(.x) > 0) max(.x, na.rm = TRUE) else as.Date(NA))
    ) %>%
  
  # Clean up the temporary list columns
  select(-extracted_dates, -date_objects)

rule2_data <- rule2_data %>%
  mutate(
    # pull the year from treatment_year, combine with tx_month, 
    # and day "01", then parse it into a date
    clean_date = ymd(str_glue("{treatment_year}-{tx_month}-01")),
    
    # Get the first day of that month
    start_date = floor_date(clean_date, unit = "month"),
    
    # Get the last day of that month (correctly factors in leap years based on the data's year!)
    end_date = ceiling_date(clean_date, unit = "month") - days(1)
  ) %>%
  # Clean up the temporary column
  select(-clean_date)


surveys_long_working <- bind_rows(rule1_data, rule2_data, rule3_data)

# remove subset dataframes to keep environment clean
#rm(rule1_data, rule2_data, rule3_data)

# STEP 14c: build out hierarchical rules that eliminate rule 2 and rule 3 data 
# where rule 1 data exists already 

surveys_dedup <- surveys_long_working %>%
  # Clean up old active_ing column - rebuild at the end
  select(-any_of("active_ing")) %>% 
  mutate(row_id = row_number()) %>%
  
  # Build temporary date intervals for robust timeline mapping
  mutate(
    treatment_interval = interval(
      ymd(start_date), 
      ymd(end_date) + days(1) 
    )
  ) %>% 
  
  # expand the chemical column
  separate_rows(chemicals_units, sep = ",\\s*") %>%
  
  # group by permit-year-chemicals
  group_by(permit_number, treatment_year, chemicals_units) %>%
  filter(
    # prioritize keeping rows from rule1
    dup_clean_id == "rule1" | 
      
      # only keep rule2 rows that don't overlap in dates with rule1 rows
      (dup_clean_id == "rule2" & !any(
        dup_clean_id == "rule1" & 
          int_overlaps(treatment_interval, .data$treatment_interval[dup_clean_id == "rule1"])
      )) |
      
      # only keep rule3 rows when there are no rule1 or rule2 rows
      (dup_clean_id == "rule3" & !any(dup_clean_id %in% c("rule1", "rule2")))
  ) %>%  
  ungroup() %>%  
  
  # re-assemble chemicals row-by-row
  group_by(row_id) %>%
  mutate(
    chemicals_units = paste(unique(chemicals_units), collapse = ", ")
  ) %>%
  ungroup() %>%
  
  distinct(row_id, .keep_all = TRUE) %>%
  select(-row_id, -treatment_interval) %>% 
  
  # rejoin active ingredients from chem_ref_table based on chemicals_units strings
  rowwise() %>% 
    mutate(
      active_ing = {
      # 1. Strip away the brackets and units (e.g., "Weedar 64 (gal)" -> "Weedar 64")
        raw_chems <- str_remove_all(chemicals_units, "\\s\\([^)]+\\)")
        chem_vector <- str_split(raw_chems, ",\\s*")[[1]]
      
        if (length(chem_vector) > 0 && chem_vector[1] != "") {
        # 2. Extract active ingredients, strip out duplicates, drop NAs, collapse
          act_vector <- unique(chem_lookup[chem_vector])
          paste(act_vector[!is.na(act_vector)], collapse = ", ")
      } else {
        ""
      }
    }
  ) %>% 
  ungroup()

# the script at this point mostly sorts duplicates, but there could be some
# inefficient pipelines, some messy naming conventions, columns to clean up
# the next step is to bring in PAR data and sort out duplicates in a similar way 
# with PAR data. 

###############
# Pull in PAR data
###############

# STEP 15
#### run the PAR script "IAPM_PARs_clean.R" ####
source(here::here("draftcode", "IAPM_PARs_clean.R"))

# STEP 16
# add columns to refer to additional details
PARs_treatment_group <- PARs_treatment_group %>% 
  left_join(
    chem_ref_table %>% select(orig_name, active_ing),
    by = c("pesticide_trade_name" = "orig_name")) %>%
  mutate(treatment_id = paste0("P_", row_number())) %>%
  mutate(additional_data = "Y")


# STEP 17
PARs_clean <- PARs_treatment_group %>% 
  mutate(
    chemicals_units = paste0(pesticide_trade_name, " (", amount_units, ")")) %>% 
  select(
    treatment_id, 
    permit_number, 
    treatment_year = app_year, 
    chemicals_units, 
    start_date,
    end_date, 
    treatment_size = PAR_treatment_size, 
    additional_data,
    active_ing
  )

# STEP 18:
# collapse PAR rows by treatment and chemical to match the format of the surveys
PARs_clean <- PARs_clean %>%
  
  group_by(pick(everything(), -c(treatment_id, chemicals_units, active_ing))) %>%
  
  mutate(
    treatment_id = paste(unique(treatment_id), collapse = ", "),
    chemicals_units = paste(unique(chemicals_units), collapse = ", "),
    active_ing = paste(unique(active_ing), collapse = ", "),
  ) %>%
  ungroup()


# STEP 19
surveys_clean <- surveys_dedup %>% 
  mutate(
    treatment_size = if_else(
      multiple_tx == "N", 
      as.numeric(total_treated_area_acres), 
      NA_real_                             
    ), 
    
    additional_data = if_else(
      multiple_tx == "N",
      "Y",
      "N"
    )
  ) %>%
  
  select(
    permit_number,
    treatment_year, 
    start_date,
    end_date,
    treatment_size,
    chemicals_units, 
    active_ing,
    treatment_id,
    additional_data
  )

# STEP 20
all_tx <- bind_rows(
  surveys_clean,
  # fix column(s) with different data types
  PARs_clean %>% mutate(treatment_year = as.integer(treatment_year)))


# STEP 21: sort out overlap between survey treatment data and PAR treatment data
all_tx_clean <- all_tx %>%
  
  # define treament_intervals
  mutate(
    treatment_interval = case_when(
      !is.na(start_date) & !is.na(end_date) ~ interval(ymd(start_date), ymd(end_date) + days(1)),
      TRUE                                  ~ NA
    )
  ) %>% 
  
  # group by permit-year
  group_by(permit_number, treatment_year) %>%
  

  filter(
    # ensure there are valid start and end dates in those being sorted
    (!is.na(start_date) & !is.na(end_date) & (
      
      # if treatment_id starts with "P", always keep it
      str_sub(treatment_id, 1, 1) == "P" |
        
        # determine if there is overlap in each S row with the treatment interval of a P row
        (str_sub(treatment_id, 1, 1) == "S" & !map_lgl(
          treatment_interval, 
          ~ any(int_overlaps(.x, treatment_interval[str_sub(treatment_id, 1, 1) == "P"]), na.rm = TRUE)
        ))
    )) |
      
      # when there is no date data
      # only keep it if it's the only record in its groups
      ((is.na(start_date) | is.na(end_date)) & n() == 1)
  ) %>%  
  
  ungroup() %>% 
  select(-treatment_interval)

# STEP 22: attempt to clean up supplemental data
# ideally we only want survey_tx_details and pars_tx_details to include
# data that we have deemed valuable


# STEP 22a: create a list of treatment_ids that made it through all filters
valid_tx_ids <- all_tx_clean %>%
  
  select(treatment_id) %>% 

  separate_longer_delim(treatment_id, delim = ",") %>% 
  
  mutate(treatment_id = str_squish(treatment_id)) %>% 
  
  # create vector of treatment_ids
  pull(treatment_id) %>% 
  unique()

# create a list of treatment_ids that did make it in the dataset but with limited info
surveys_no_additional_data <- all_tx_clean %>%
  # Isolate only the rows with "N"
  filter(additional_data == "N") %>%
  
  # Keep only the treatment_id column to save memory
  select(treatment_id) %>%
  
  # Unpack any comma-separated strings so we can evaluate every individual ID
  separate_longer_delim(treatment_id, delim = ",") %>%
  mutate(treatment_id = str_squish(treatment_id)) %>%
  
  # Pull the column into a simple vector of unique bad IDs
  pull(treatment_id) %>%
  unique()

# STEP 22b: filter tx_details data based on the IDs in step 22a
survey_tx_details <- survey_tx_details %>%
  # remove treatment_ids that don't appear in treatment dataset
  filter(treatment_id %in% valid_tx_ids) %>%
  # remove tx_details for rows with "no additional data" because limited survey
  filter(!treatment_id %in% surveys_no_additional_data)


pars_tx_details <- PARs_treatment_group %>%
  # remove treatment_ids that don't appear in treatment dataset
  filter(treatment_id %in% valid_tx_ids)

# STEP: 22c clean up par_tx_details columns
# NOTE: this step will be simplified when I fix the crazy column names. 
pars_tx_details <- pars_tx_details %>% 
  select(c(
    permit_number,
    treatment_year = app_year,
    treatment_id,
    target_species = species_code,
    chemical = pesticide_trade_name,
    chem_unit = amount_units,
    start_date,
    end_date,
    treatment_size = PAR_treatment_size,
    avg_depth,
    avg_vol_rate,
    avg_areal_rate,
    amount_total,
    water_temp_low,
    water_temp_high,
    wind_speed_low,
    wind_speed_high,
    air_temp_low,
    air_temp_high,
    avg_speed,
    wind_dir,
    active_ing
  ))

# add lakes and DOW id
# NOTE: this step shouldn't be necessary but I accidentally dropped these data earlier
all_tx_clean <- all_tx_clean %>%
  left_join(
    apm_iapm_permits %>%
      distinct(permit_number, .keep_all = TRUE) %>%
      select(
        permit_number, 
        lake_name = water_resource_names, 
        DOW = water_resource_numbers,
        target_species = species
      ),
    by = "permit_number"
  )

# Add code here that adds littoral zone data to the all_tx_clean
#lake_attribute_df <- read.csv('data/MN_lake_basin_littoral_zone_15_ft_std.csv')


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

# Right now this uses the permit species for all rows, but I could prioritize the PAR species when it exists. 
all_tx_clean <- all_tx_clean %>%
  mutate(
    # 1. Split the text strings into a structured list of individual species
    species_list = str_split(target_species, ";\\s*"),
    
    # 2. Map each list element against the master reference table
    species_code = map_chr(species_list, function(x) {
      # Use match() to find the single-species codes from our reference
      matched_codes <- species_reference$code[match(x, species_reference$species_name)]
      
      # Handle missing translations cleanly
      matched_codes <- coalesce(matched_codes, "NA")
      
      # 3. Stitch them back together into a comma-separated code string
      paste(matched_codes, collapse = "_")
    })
  ) %>%
  # Clean up the helper column
  select(-species_list)

#export data: 
all_tx_clean %>%
  write_csv("outputdata/all_IAPM_treatments.csv")

pars_tx_details %>%
  write_csv("outputdata/IAPM_PAR_details.csv")

survey_tx_details %>%
  write_csv("outputdata/survey_tx_details.csv")



# exploring chemicals when there are repeated treatments
multiple_iapm_tx <- all_tx_clean %>%
  # Group and count simultaneously
 group_by(permit_number, treatment_year) %>%
  # Filter for combinations that appear more than once
  filter(n() > 1) %>%
  arrange(permit_number, treatment_year)
