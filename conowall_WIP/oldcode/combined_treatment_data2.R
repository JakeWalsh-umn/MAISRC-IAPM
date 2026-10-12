# This is the second version of combined_treatment_data.R
# I am trying to improve the process

library(tidyverse)
library(lubridate)
library(DescTools)
library(stringr)

ifelse(!dir.exists("data"), 
       dir.create("data"), 
       "Folder exists already")
ifelse(!dir.exists("outputdata"), 
       dir.create("outputdata"), 
       "Folder exists already")
ifelse(!dir.exists("draftcode"), 
       dir.create("draftcode"), 
       "Folder exists already")

# STEP 1
#### load data #### 
iapm_annual_survey <- read.csv('data/iapm_annual_survey_all.csv')
apm_iapm_permits <- read.csv('data/apm_iapm_permit_detail.csv')

# STEP 8a and 8b
# remove the columns we don't care about
iapm_annual_survey <- iapm_annual_survey %>%
  select(-c(data_date_time, when_permit_expires, starts_with("aapcd"),
            satisified_with_swim_itch_control)) %>% 
  rename(treatment_year = calander_year) %>%
  mutate(
    # Step 1: Remove all whitespace
    total_treated_area_acres = str_remove_all(total_treated_area_acres, "\\s+"),
    
    # Step 2: Replace multiple consecutive periods with a single period
    total_treated_area_acres = str_replace_all(total_treated_area_acres, "\\.+", "."),
    
    # Step 3: Fix scattered decimals (e.g., 45.6.7 -> 45.67)
    # This keeps the first decimal and removes any subsequent ones
    total_treated_area_acres = str_replace(total_treated_area_acres, "^([^.]*\\.[^.]*)\\.", "\\1"),
    
    # Optional Step 4: Convert to numeric if you need to do math on it
    total_treated_area_acres = as.numeric(total_treated_area_acres)
  )


# STEP 2
# left join by permit number and keep the permit DOW and treatment_method
iapm_annual_survey <- iapm_annual_survey %>%
  left_join(apm_iapm_permits %>% select(DOW = water_resource_numbers, lk_name = water_resource_names, permit_number, 
                                        tx_method_permit = treatment_method), by = "permit_number")

# STEP 3
# figure out if there are surveys without permit details 
surveys_wo_permit_data <- iapm_annual_survey %>% anti_join(apm_iapm_permits, by = "permit_number")

# count how many surveys without permit details
length(unique(surveys_wo_permit_data$permit_number)) # there are 219 permits with survey responses but no permit data. 

#remove the survey responses that don't have corresponding permit details
iapm_annual_survey <- iapm_annual_survey %>% anti_join(surveys_wo_permit_data, by = "apm_annual_survey_id")


# STEP 4
# remove rows where the key treatment details are exact duplicates
# this treats commercial and landowner rows as the same, and is only looking at treatment data.
iapm_annual_survey <- iapm_annual_survey %>%
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
  0, 0.0,            # Numeric zeros
  NA, NaN,           # System missing
  "NA", "NULL",      # Common string placeholders
  "", " ",           # Empty strings
  "None", "n/a", "." # Manual entry placeholders
)

# STEP 5a
# try to tease apart what kind of treatments the surveys represent
# first, "chemical" if there is chemical data
iapm_annual_survey <- iapm_annual_survey %>%
  mutate(across(starts_with("hp"), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(treatment_method = if_else(
    if_any(starts_with("hp"), ~ !(.x %in% no_data | is.na(.x))),
    "Chemical", 
    "other"
  ))

# STEP 5b
# second, note when there is chemical and mechanical data in a survey
iapm_annual_survey <- iapm_annual_survey %>%
  mutate(
    treatment_method2 = case_when(
      !(total_cut_area_acres %in% no_data) & 
        (!(total_treated_area_acres %in% no_data) | treatment_method == "Chemical") ~ "Mechanical, Chemical",
      
      # "chemical" is also designated when total_treated_area_acres has data
      !(total_treated_area_acres %in% no_data) | treatment_method == "Chemical" ~ "Chemical",
      
      # mechanical if there is total_cut_area_acres data
      !(total_cut_area_acres %in% no_data) ~ "Mechanical",
      
      TRUE ~ "none"
    )
  )  


# this code may help make sense of how to designate the treatment methods for surveys. 
table(iapm_annual_survey$tx_method_permit, iapm_annual_survey$treatment_method)
table(iapm_annual_survey$tx_method_permit, iapm_annual_survey$treatment_method2)

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

# NOTE: skipped the step here where the month columns are turned into dates rather than the word of month

# STEP 10
# add together all chemical amount data for easy reference of whether there is any chemical data
iapm_annual_survey <- iapm_annual_survey %>%
  mutate(
    sum_chem = rowSums(across(
      .cols = starts_with("hp"),
      .fns = ~ as.numeric(.x) # Ensures values are treated as numbers
    ), na.rm = TRUE)
  )

# STEP 11
iapm_annual_survey <- iapm_annual_survey %>% mutate(treatment_id = paste0("S_", row_number()))


# STEP 17a
# set aside subset of the data to build out treatment details dataset
survey_tx_details <- iapm_annual_survey %>% 
  
  select(
    permit_id, 
    treatment_id, 
    treatment_dates, 
    total_treated_area_acres, 
    starts_with("hp")
  )

iapm_annual_survey <- iapm_annual_survey %>%
  
  rowwise() %>%
  
  mutate(
    chemicals = {
      # Grab the names of all 'hp' columns in this specific row
      hp_cols <- names(pick(starts_with("hp")))
      
      # Grab the values for those columns in this row
      hp_vals <- c_across(starts_with("hp"))
      
      # only keep columns where there is useful data 
      active_cols <- hp_cols[!(hp_vals %in% no_data) & !is.na(hp_vals)]
      
      # If any active columns exist, extract chemical names from between underscores
      if (length(active_cols) > 0) {
        chem_names <- str_match(active_cols, "^[^_]+_(.*)_[^_]+$")[, 2]
        paste(unique(chem_names), collapse = ", ")
      } else {
        "" 
      }
    }
  ) %>%
  
  ungroup()

# STEP 17b
# clean up survey dataset
surveys_narrow <- iapm_annual_survey %>%
  select(-c(starts_with("hp")))

# STEP 12
# pivot the month columns to long form so each row represents a survey-month
surveys_long <- surveys_narrow %>%
  pivot_longer(
    cols = starts_with("work_done_in_"), 
    values_to = "tx_month",           
    values_drop_na = FALSE
  ) %>%
  select(-name)

# STEP 13
# remove duplicate rows within different types of survey responses
# AND
# remove extra rows that were created in the pivot_longer
# work_month_count = 0, tx_confirmed = "Y" keep the row
# work_month_count >= 1 remove rows where treatment_range has no data


# NOTE: potential problem with this step is that there could be cases where 
## a row with dates 5/3/2021, 10/2/2021 may result in only keeping 10/2/2021 
### if 5/3/2021 was entered separately earlier in the season. Is that okay?

# handle data where there is data to retain but didn't show up as months of treatment
rule1_data <- surveys_long %>%
  filter(work_month_count == 0 & tx_dates == "Y") %>%
  mutate(
    original_string = treatment_dates,
    string_length = nchar(treatment_dates)) %>%
  
  # 2. Split the dates into individual rows so we can evaluate them one by one
  separate_longer_delim(treatment_dates, delim = ",") %>%
  mutate(treatment_dates = str_squish(treatment_dates)) %>%
  
  # 3. Parse the dates to ensure exact matches
  mutate(
    parsed_date = as.Date(parse_date_time(treatment_dates, orders = c("mdy", "ymd", "dmy")))) %>%
  
  # -------------------------------------------------------------------------
# FIXED STEP 1: KEEP ONLY UNIQUE DATES PER TREATMENT ID
# -------------------------------------------------------------------------
# This kills duplicate rows and duplicate user entries at the source.
  distinct(treatment_id, parsed_date, .keep_all = TRUE) %>%
  
  mutate(
    date_count = str_count(original_string, "\\d+[/-]\\d+[/-]\\d+"),
    date_count = coalesce(date_count, 0)) %>% 
  
  filter(date_count > 1) %>%
  
  # -------------------------------------------------------------------------
# FIXED STEP 2: GROUP EXPLICITLY BY TREATMENT ID (NOT EVERYTHING)
# -------------------------------------------------------------------------
  group_by(treatment_id) %>% 
  arrange(parsed_date, .by_group = TRUE) %>%
  mutate(
    group_anchor = accumulate(parsed_date, ~ if (.y - .x <= 5) .x else .y),
    day_gap = as.numeric(parsed_date - group_anchor),
    new_treatment = if_else(parsed_date == group_anchor, 1, 0),
    treatment_group = cumsum(new_treatment)) %>% # Added missing closing parenthesis here
  
  # Optional: drop the helper column if you don't need it anymore
  select(-group_anchor) %>% # Added missing pipe here
  
  # -------------------------------------------------------------------------
# FIXED STEP 3: EXPLICIT GROUP COLLAPSE
# -------------------------------------------------------------------------
  group_by(treatment_id, treatment_group) %>%
  
  # Pulling from distinct rows means this string can no longer have duplicates
  mutate(
    collapsed_string = paste(as.character(parsed_date), collapse = ", ")) %>%
  
  filter(row_number() == 1) %>% 
  
  mutate(parsed_date = collapsed_string) %>%
  
  ungroup() %>%
  #select(-c(day_gap, new_treatment, treatment_group, collapsed_string)) %>%
  
  mutate(dup_clean_id = paste0("rule1"))

# handle data where the months of treatment are defined
rule2_data <- surveys_long %>%
  filter(
    work_month_count >= 1,
    !is.na(tx_month),               
    !tx_month %in% no_data) %>%
  mutate(dup_clean_id = paste0("rule2"))

rule3_data <- surveys_long %>%
  filter(work_month_count == 0 & tx_dates == "N") %>%  # 4. Group by the unique combinations you want to keep
  group_by(permit_number, treatment_year) %>%
  
  # 5. For each unique date, keep only the row that came from the shortest original string
  # note to self: should this following arrange step include commercial/landowner
  arrange(
    sum_chem %in% no_data,    # 2nd priority: FALSE (valid data) comes before TRUE (no-data)
    desc(total_treated_area_acres),        # 3rd priority: Highest acres first
    .by_group = TRUE            # Ensures sorting happens within our groups
  ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(dup_clean_id = paste0("rule3"))

# STEP: add treatment start and end dates using the date strings (rule 1) and months (rule 2)
rule1_data <- rule1_data %>%
  select(-c(original_string, string_length, date_count, day_gap, new_treatment, collapsed_string)) %>%
  rowwise() %>%
  mutate(
    # 1. Extract all sequences of text that look like dates (4 digits - 2 digits - 2 digits)
    extracted_dates = list(str_extract_all(parsed_date, "\\d{4}-\\d{2}-\\d{2}")[[1]]),
    
    # 2. Convert those extracted strings into actual Date objects
    date_objects = list(ymd(extracted_dates)),
    
    # 3. Safely calculate min and max row by row
    start_date = if(length(date_objects[[1]]) > 0) min(date_objects[[1]], na.rm = TRUE) else as.Date(NA),
    end_date   = if(length(date_objects[[1]]) > 0) max(date_objects[[1]], na.rm = TRUE) else as.Date(NA)
  ) %>%
  ungroup() %>%
  # Clean up temporary columns
  select(-extracted_dates, -date_objects)

rule2_data <- rule2_data %>%
  mutate(
    # Step 1: Combine the current year (2026), the month name, and day "01" 
    # into a single string, then parse it into a real Date object
    clean_date = ymd(str_glue("2026-{tx_month}-01")),
    
    # Step 2: Get the first day of that month
    start_date = floor_date(clean_date, unit = "month"),
    
    # Step 3: Get the last day of that month (handles Feb 28 vs 29 automatically)
    end_date = ceiling_date(clean_date, unit = "month") - days(1)
  ) %>%
  # Clean up the temporary column
  select(-clean_date)


surveys_long_working <- bind_rows(rule1_data, rule2_data, rule3_data)

#rm(rule1_data, rule2_data, rule3_data)

# STEP: build out hierarchical rules that eliminate rule 2 and rule 3 data 
# where one is already 

surveys_dedup <- surveys_long_working %>%
  mutate(row_id = row_number()) %>%
  separate_rows(chemicals, sep = ",\\s*") %>%
  
  # Step B: Group by permit, year, and the individual chemical string
  group_by(permit_number, treatment_year, chemicals) %>%
  
  # Step C: Run the hierarchical filtering logic
  filter(
    # Clause A: Always keep Rule 1
    dup_clean_id == "rule1" | 
      
      # Clause B: Evaluate Rule 2 rows
      (dup_clean_id == "rule2" & !any(
        dup_clean_id == "rule1" &           # <--- Checked against rule1, not rule2!
          start_date <= end_date & 
          end_date >= start_date
      )) |
      
      (dup_clean_id == "rule3" & !any(dup_clean_id %in% c("rule1", "rule2")))) %>%
    ungroup() %>%
  
  # Step D: Re-assemble the original rows based on row_id so chemicals are bundled back together
  group_by(row_id, permit_number, treatment_year, dup_clean_id, start_date, end_date) %>%
  summarize(
    chemicals = paste(chemicals, collapse = ", "),
    .groups = "drop"
  ) %>%
  
  # Clean up by arranging and removing the temporary row tracking ID
  arrange(permit_number, row_id) %>%
  select(-row_id)





