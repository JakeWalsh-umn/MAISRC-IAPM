# last updated: May 19, 2026
# this script provides the most comprehensive IAPM overview of treatments
# the output is one main dataset for PI charter and analysis with supplemental datasets


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

#load data
iapm_annual_survey <- read.csv('data/iapm_annual_survey_all.csv')
apm_iapm_permits <- read.csv('data/apm_iapm_permit_detail.csv')


# figure out if there are surveys without permit details 
iapm_annual_survey <- iapm_annual_survey %>%
  left_join(apm_iapm_permits %>% select(permit_number, treatment_method), by = "permit_number")

surveys_wo_permit_data <- iapm_annual_survey %>% anti_join(apm_iapm_permits, by = "permit_number")

length(unique(surveys_wo_permit_data$permit_number)) # there are 219 permits with survey responses but no permit data. 

#remove the survey responses that don't have corresponding permit details
iapm_annual_survey <- iapm_annual_survey %>% anti_join(surveys_wo_permit_data, by = "apm_annual_survey_id")

## at this point the survey dataframe should have mechanical/chemical for all

# standardize "hp" column formats
iapm_annual_survey <- iapm_annual_survey %>% 
  mutate(across(starts_with("hp"), ~ suppressWarnings(as.numeric(.x))))

no_data <- c(
  0, 0.0,            # Numeric zeros
  NA, NaN,           # System missing
  "NA", "NULL",      # Common string placeholders
  "", " ",           # Empty strings
  "None", "n/a", "." # Manual entry placeholders
)

# try to tease apart what kind of treatments the surveys represent
# first, "chemical" if there is chemical data
iapm_annual_survey <- iapm_annual_survey %>%
  mutate(treatment_method = if_else(
    if_any(starts_with("hp"), ~ !(.x %in% no_data | is.na(.x))),
    "Chemical", 
    "other"
  ))

# second, note when there is chemical and mechanical data in a survey
iapm_annual_survey <- iapm_annual_survey %>%
  mutate(
    treatment_method2 = case_when(
      !(total_cut_area_acres %in% no_data) & 
        (!(total_treated_area_acres %in% no_data) | treatment_method == "Chemical") ~ "mechanical, chemical",
    
      # "chemical" is also designated when total_treated_area_acres has data
      !(total_treated_area_acres %in% no_data) | treatment_method == "Chemical" ~ "chemical",
      
      # mechanical if there is total_cut_area_acres data
      !(total_cut_area_acres %in% no_data) ~ "mechanical",
      
      TRUE ~ "none"
    )
  )

# third, join the method from the permit df
iapm_annual_survey <- iapm_annual_survey %>%
  left_join(
    # Subset the permit df to ONLY include the key and your target column
    apm_iapm_permits %>% select(permit_number, tx_method_permit = treatment_method), 
    by = "permit_number"
  )


iapm_annual_survey <- iapm_annual_survey %>%
  mutate(
    work_month_count = rowSums(across(
      .cols = starts_with("work_done_in_"),
      .fns = ~ !(.x %in% no_data | is.na(.x))
    ))
  )

#### Add the following columns to the survey data table 

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


# remove the columns we don't care about
iapm_annual_survey <- iapm_annual_survey %>%
  select(-c(data_date_time, when_permit_expires, starts_with("aapcd"),
            satisified_with_swim_itch_control))


# rename the misspelled column
iapm_annual_survey <- iapm_annual_survey %>%
  rename(
    treatment_year = calander_year
  )


# change the month columns into date ranges instead of months
iapm_annual_survey <- iapm_annual_survey %>%
  mutate(across(starts_with("work_done_in_"), ~ {
    
    clean_month = str_trim(tolower(.x))
    
    m_num = match(clean_month, tolower(month.name))
    m_num = coalesce(m_num, match(clean_month, tolower(month.abb)))
    
    start_date <- make_date(year = treatment_year, month = m_num, day = 1)
    end_date <- rollback(start_date + months(1))
    
    ifelse(!is.na(m_num),
           paste0(month(start_date), "/", day(start_date), "/", year(start_date),
                  "-",
                  month(end_date), "/", day(end_date), "/", year(end_date)),
           NA_character_
    )
  }))

# add together all chemical data for easy reference of whether there is any chemical data
iapm_annual_survey <- iapm_annual_survey %>%
  mutate(
    sum_chem = rowSums(across(
      .cols = starts_with("hp"),
      .fns = ~ as.numeric(.x) # Ensures values are treated as numbers
    ), na.rm = TRUE)
  )


 #### code for inspecting duplicate surveys ####
duplicate_survey_idyear <- iapm_annual_survey %>%
  group_by(permit_id, treatment_year) %>%
  filter(n() > 1) %>%
  ungroup()

duplicate_survey_idyear %>%
  count(permit_number, treatment_year) %>%
  pull(n) %>%
  mean()

duplicate_survey_idyear_landowner <- duplicate_survey_idyear %>% 
  filter(apm_annual_survey_type_id == 1)

duplicate_survey_idyear_commercial <- duplicate_survey_idyear %>% 
  filter(apm_annual_survey_type_id == 2)

duplicate_survey_commercial <- duplicate_survey_idyear_commercial %>%
  group_by(permit_id, treatment_year) %>%
  filter(n() > 1) %>%
  ungroup()

duplicate_survey_landowner <- duplicate_survey_idyear_landowner %>%
  group_by(permit_id, treatment_year) %>%
  filter(n() > 1) %>%
  ungroup()

# the above code doesn't actually have an applied use as the script is written now (5/19/2026)

#### next steps ####
# add a unique id to the survey data as it is right now and then pivot the chemical data long by month. 


iapm_annual_survey <- iapm_annual_survey %>% mutate(treatment_id = paste0("S_", row_number()))


# change this code pivoting the surveys longer
surveys_long <- iapm_annual_survey %>%
  pivot_longer(
    cols = starts_with("work_done_in_"), # Selects all your month columns
    values_to = "treatment_range",           # New column will be named 'work_done_in'
    values_drop_na = FALSE                 # Removes the NULL/NA rows automatically
  ) %>%
  select(-name)

# work_month_count = 0, tx_confirmed = "Y" keep the row
# work_month_count >= 1 remove rows where treatment_range has no data

# handle data where there is data to retain but didn't show up as months of treatment
rule1_data <- surveys_long %>%
  filter(work_month_count == 0 & tx_confirmed == "Y") %>%
  distinct(treatment_id, .keep_all = TRUE)

# handle data where the months of treatment are defined
rule2_data <- surveys_long %>%
  filter(work_month_count >= 1 & !is.na(treatment_range))

surveys_long_cleaned <- bind_rows(rule1_data, rule2_data)

# expand comma delimited dates
surveys_long_exdates <- surveys_long_cleaned %>%
  # Create a copy of the original column to be the one we expand
  mutate(individual_treatment_date = treatment_dates) %>%
  # Expand the copy into multiple rows, keeping the original 'treatment_dates' for reference
  separate_rows(individual_treatment_date, sep = ",") %>%
  mutate(individual_treatment_date = trimws(individual_treatment_date))


surveys_long_exdates <- surveys_long_exdates %>%
  mutate(
    # Regular expression matches any sequence of numbers separated by / or -
    date_count = str_count(treatment_dates, "\\d+[/-]\\d+[/-]\\d+"),
    
    # Optional: Turn NA results into 0 instead of leaving them as NA
    date_count = coalesce(date_count, 0)
  )

# group together treatment rows that have different individual dates but are part of the same treatment
# multiple dates become comma delimited
surveys_multi_date <- surveys_long_exdates %>%
  
  filter(date_count > 1) %>%
  
  mutate(individual_treatment_date = lubridate::mdy(individual_treatment_date)) %>%
  
  group_by(pick(everything(), -individual_treatment_date)) %>%
  arrange(treatment_id, individual_treatment_date, .by_group = TRUE) %>%
  
  # define treatment groups by those that happened within 5 days of each other
  mutate(
    day_gap = as.numeric(individual_treatment_date - lag(individual_treatment_date)),
    new_treatment = if_else(is.na(day_gap) | day_gap > 5, 1, 0),
    treatment_group = cumsum(new_treatment)
  ) %>%
  
  group_by(pick(everything(), -individual_treatment_date), treatment_group) %>%
  
  mutate(
    collapsed_string = paste(as.character(individual_treatment_date), collapse = ", ")
  ) %>%
  
  filter(row_number() == 1) %>% 
  
  mutate(individual_treatment_date = collapsed_string) %>%
  
  ungroup() %>%
  select(-c(day_gap, new_treatment, treatment_group, collapsed_string))
  

surveys_single_no_date <- surveys_long_exdates %>%
  filter(date_count < 2) %>%
  mutate(
    individual_treatment_date = format(
      parse_date_time(individual_treatment_date, orders = c("ymd", "mdy"), quiet = TRUE), 
      "%Y-%m-%d"
    )
  )
    
    
# combine the two separately processed datasets into one
surveys_long_clean <- bind_rows(surveys_single_no_date, surveys_multi_date)


# set aside subset of the data to build out treatment details dataset
surveys_tx_details <- surveys_long_clean %>% 
  
  filter(multiple_tx == "N") %>% 
  
  select(
    permit_id, 
    treatment_id, 
    treatment_dates, 
    total_treated_area_acres, 
    starts_with("hp")
  )


# summarize chemical data to be simpler for treatment records dataset
surveys_long_clean <- surveys_long_clean %>%
  
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

#remove the chemical details now that we don't need it. 
surveys_long_clean <- surveys_long_clean %>% 
  select(-starts_with("hp"))

#### run the PAR script "IAPM_PARs_clean.R" ####
source(here::here("draftcode", "IAPM_PARs_clean.R"))

# add columns to refer to additional details
PARs_treatment_group <- PARs_treatment_group %>% mutate(treatment_id = paste0("P_", row_number())) %>%
                                                 mutate(additional_data = "Y")
                                  
PARs_clean <- PARs_treatment_group %>% 
  select(treatment_id, permit_number = Permit.number, treatment_year = App.year, chemicals = Pesticide.trade.name, start_date,
         end_date, treatment_size, additional_data)

# collapse PAR rows by treatment and chemical to match the format of the surveys
PARs_clean <- PARs_clean %>%
  
  group_by(pick(everything(), -c(treatment_id, chemicals))) %>%
  
  summarize(
    treatment_id = paste(unique(treatment_id), collapse = ", "),
    chemicals = paste(unique(chemicals), collapse = ", "),
    .groups = "drop" 
  )



# define a function to help convert multiple date formats into start and end dates
parse_flexible_dates <- function(text_string) {
  if (is.na(text_string) || text_string == "") {
    return(as.Date(NA))
  }
  
  # Split on commas, spaces around dashes, OR the middle dash between two yyyy-mm-dd blocks
  # e.g., converts "2026-05-19-2026-05-22" into "2026-05-19" "2026-05-22"
  raw_fragments <- str_split_1(text_string, ",|\\s+[-–—]\\s+|(?<=\\d)-(?=\\d{4})") %>% str_trim()
  
  # Drop empty strings
  raw_fragments <- raw_fragments[raw_fragments != "" & !is.na(raw_fragments)]
  
  if (length(raw_fragments) == 0) return(as.Date(NA))
  
  parsed_dates <- parse_date_time(raw_fragments, orders = "ymd", quiet = TRUE)
  
  return(as.Date(parsed_dates))
}

# use above function to add start and end dates to the survey data
surveys_long_clean <- surveys_long_clean %>%
  rowwise() %>%
  mutate(
    start_date = {
      if (!is.na(individual_treatment_date) & individual_treatment_date != "") {
        dates <- parse_flexible_dates(individual_treatment_date)
        if (all(is.na(dates))) as.Date(NA) else min(dates, na.rm = TRUE)
        
      } else if (!is.na(treatment_range) & treatment_range != "") {
        dates <- parse_flexible_dates(treatment_range)
        if (all(is.na(dates))) as.Date(NA) else min(dates, na.rm = TRUE)
        
      } else {
        as.Date(NA)
      }
    },
    
    end_date = {
      if (!is.na(individual_treatment_date) & individual_treatment_date != "") {
        dates <- parse_flexible_dates(individual_treatment_date)
        if (all(is.na(dates))) as.Date(NA) else max(dates, na.rm = TRUE)
        
      } else if (!is.na(treatment_range) & treatment_range != "") {
        dates <- parse_flexible_dates(treatment_range)
        if (all(is.na(dates))) as.Date(NA) else max(dates, na.rm = TRUE)
        
      } else {
        as.Date(NA)
      }
    }
  ) %>%
  ungroup()
  

surveys_clean <- surveys_long_clean %>% 
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
    chemicals, 
    treatment_id,
    additional_data
  )

all_tx <- bind_rows(
  surveys_clean,
  # fix column(s) with different data types
  PARs_clean %>% mutate(treatment_year = as.integer(treatment_year)))

#### need to figure out how to handle duplicates 
## can we use the treatment ids to help us? Did we already sort duplicates for each source respectively? 
duplicate_tx_row <- all_tx %>%
  # 1. Create your temporary month column
  mutate(month = lubridate::month(start_date, label = TRUE, abbr = TRUE)) %>%
  
  # 2. Group by your target columns (this keeps all other columns hidden in the background)
  group_by(permit_number, treatment_year, month) %>%
  
  # 3. Filter to keep rows where the group size is greater than 1
  filter(n() > 1) %>%
  
  # 4. Always ungroup at the end so downstream functions work normally
  ungroup()


#### next steps ####
# join the all_iapm_treatments data back to the permit data

#### CAUTION ####
# I think there are duplicate rows between surveys and PARs that still need to be sorted. 


### at this point there are the following datasets
# all_tx which represents the treatments that we know happened. 
# PARs_treatment_group which is the permit-treatment-chemical grouped PAR data
# surveys_tx_details which has the additional details on treatments

#export data: 
all_tx %>%
  write_csv("outputdata/all_IAPM_treatments.csv")

PARs_treatment_group %>%
  write_csv("outputdata/IAPM_PAR_details.csv")

surveys_tx_details %>%
  write_csv("outputdata/survey_tx_details.csv")
  




