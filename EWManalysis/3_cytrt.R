# adapting Melina's treatment history matrix code for my own use. 

#create unique IDs for PI surveys
MN_foo <- MN_foo %>%
  mutate(
    #survey_id = row_number(),
    year = year(SURVEY_START)
  ) 

EWM_trts <- MN_trts_acre %>%
  filter(species_code == "EWM") %>%
  mutate(start_date   = mdy(start_date),
         end_date   = mdy(end_date)) %>%
  filter(!is.na(start_date))

# start building out the treatment history matrix using confirmed treatments ----

treatment_history_matrix <- matrix(0, nrow = nrow(MN_foo), ncol = 15) %>%
  as.data.frame()

colnames(treatment_history_matrix) <- c(
  "survey_id", 
  "lag0_before", 
  "lag0_after", 
  "acres_lag0_before",
  "acres_lag0_after",
  "date_lag0_before",
  "date_lag0_after",
  paste0("lag", 1:8)
)

treatment_history_matrix$date_lag0_before <- as.Date(rep(NA, nrow(MN_foo)))
treatment_history_matrix$date_lag0_after  <- as.Date(rep(NA, nrow(MN_foo)))

treatment_history_matrix$acres_lag0_before <- rep(NA, nrow(MN_foo))
treatment_history_matrix$acres_lag0_after  <- rep(NA, nrow(MN_foo))


treatment_history_matrix$survey_id <- MN_foo$survey_id
treatment_history_matrix$SURVEY_START <- MN_foo$SURVEY_START
treatment_history_matrix$DOW_parent <- MN_foo$DOW_parent
treatment_history_matrix <- treatment_history_matrix %>%
  relocate(DOW_parent, SURVEY_START, .before = 2)

# Tracking vector
target_string <- rep(NA_character_, nrow(MN_foo))


#Populate matrix row-by-row
for (i in 1:nrow(MN_foo)) {
  
  current_dow  <- MN_foo$DOW_parent[i]
  survey_dt    <- MN_foo$SURVEY_START[i]
  survey_yr    <- MN_foo$year[i]
  
  # Check if the lake exists in the CLP treatment database
  if (current_dow %in% EWM_trts$DOW_parent) {
    
    # Subset treatment data for this specific lake
    EWM_trts_i <- EWM_trts %>% 
      filter(DOW_parent == current_dow)
    
    # Track species code if present
    if ("species_code" %in% names(EWM_trts_i)) {
      target_string[i] <- paste(unique(EWM_trts_i$species_code), collapse = ",")
    }
    
    # --- Same-Year Logic (lag0_before vs lag0_after) ---
    same_year_tx <- EWM_trts_i %>% 
      filter(treatment_year == survey_yr) %>% 
      arrange(start_date) # Sort by date to guarantee [1] is the first treatment
    
    if (nrow(same_year_tx) > 0) {
      before_tx <- same_year_tx %>% filter(start_date <= survey_dt)
      after_tx  <- same_year_tx %>% filter(start_date > survey_dt)
      
      if (nrow(before_tx) > 0) {
        treatment_history_matrix[i, "lag0_before"]       <- 1
        treatment_history_matrix[i, "acres_lag0_before"] <- before_tx$treatment_size[1]
        treatment_history_matrix[i, "date_lag0_before"]  <- before_tx$start_date[1]
      }
      
      if (nrow(after_tx) > 0) {
        treatment_history_matrix[i, "lag0_after"]        <- 1
        treatment_history_matrix[i, "acres_lag0_after"]  <- after_tx$treatment_size[1]
        treatment_history_matrix[i, "date_lag0_after"]   <- after_tx$start_date[1]
      }
    }
    
    
    # --- Lag Years 1 through 8 ---
    for (j in 1:8) {
      target_year <- survey_yr - j
      col_name <- paste0("lag", j)
      
      if (target_year %in% EWM_trts_i$treatment_year) {
        treatment_history_matrix[i, col_name] <- 1
      }
    }
    
  } else {
    # No treatment history for this lake DOW
    target_string[i] <- "none"
  }
}


#filter dat_dlnm to surveys that have previous treatments
survey_ids_treated <- treatment_history_matrix %>%
  filter(lag0_before == 1 | if_any(lag1:lag8, ~ .x == 1)) %>%
  pull(survey_id)

#MN_foo_trt_hist <- MN_foo %>% #filter
#  filter(survey_id %in% survey_ids_treated)

# skipped step where Melina removes surveys where we only know the month of treatment.
# add permit data if I want to. 

# skipped loading permit and littoral zone data that Melina loaded here

# NOTE: change this to string detecting Eurasian Watermilfoil
EWM_permits <- iapm_permits %>% 
  filter(species == "Eurasian Watermilfoil") %>%
  mutate(DOW_parent = str_replace(water_resource_numbers, "\\d{2}$", "00"))

#permit treatment history matrix
permit_treatment_history_matrix <- matrix(0, nrow = nrow(MN_foo), ncol = 10) %>%
  as.data.frame()

colnames(permit_treatment_history_matrix) <- c(
  "survey_id", 
  paste0("lag", 0:8)
)

permit_treatment_history_matrix$survey_id <- MN_foo$survey_id
permit_treatment_history_matrix$SURVEY_START <- MN_foo$SURVEY_START
permit_treatment_history_matrix$DOW_parent <- MN_foo$DOW_parent
permit_treatment_history_matrix <- permit_treatment_history_matrix %>%
  relocate(DOW_parent, SURVEY_START, .before = 2)

# Tracking vector
target_string <- rep(NA_character_, nrow(MN_foo))


#Populate matrix row-by-row
for (i in 1:nrow(MN_foo)) {
  
  current_dow  <- MN_foo$DOW_parent[i]
  survey_dt    <- MN_foo$SURVEY_START[i]
  survey_yr    <- MN_foo$year[i]
  
  # Check if the lake exists in the CLP treatment database
  if (current_dow %in% EWM_permits$DOW_parent) {
    
    # Subset treatment data for this specific lake
    EWM_permits_i <- EWM_permits %>% 
      filter(DOW_parent == current_dow)
    
    # Track species code if present
    if ("species" %in% names(EWM_permits_i)) {
      target_string[i] <- paste(unique(EWM_permits_i$species), collapse = ",")
    }
    
    # --- Same-Year Logic (lag0) ---
    same_year_tx <- EWM_permits_i %>% 
      filter(
        treatment_year == survey_yr,
        !is.na(yr_approved_acres),
        yr_approved_acres > 0
      )
    
    permit_treatment_history_matrix[i, "lag0"] <- ifelse(nrow(same_year_tx) > 0, 1, 0)
    
    # --- Lag Years 1 through 8 ---
    for (j in 1:8) {
      target_year <- survey_yr - j
      col_name <- paste0("lag", j)
      
      has_treated_lag <- EWM_permits_i %>% 
        filter(
          treatment_year == target_year,
          !is.na(yr_approved_acres),
          yr_approved_acres > 0
        ) %>% 
        nrow() > 0
      
      permit_treatment_history_matrix[i, col_name] <- ifelse(has_treated_lag, 1, 0)
    }
    
  } else {
    # No treatment history for this lake
    target_string[i] <- "none"
    
    # Explicitly ensure lag columns stay 0
    permit_treatment_history_matrix[i, paste0("lag", 0:8)] <- 0
  }
}

#check the number of surveys that had treatment histories in original matrix
sum(permit_treatment_history_matrix$survey_id %in% survey_ids_treated, na.rm = TRUE)
#check the number of surveys with permit history (at least one 1 in lag matrix)
sum(rowSums(permit_treatment_history_matrix[, 4:12] == 1) > 0)


filtered_permit_matrix <- permit_treatment_history_matrix %>%
  filter(survey_id %in% survey_ids_treated)

filtered_treatment_matrix <- treatment_history_matrix %>%
  filter(survey_id %in% survey_ids_treated)


#check how many 1s there are in each of the filtered matrices
sum(filtered_permit_matrix[, paste0("lag", 0:8)] == 1, na.rm = TRUE)
sum(filtered_treatment_matrix[, 4:13] == 1, na.rm = TRUE)


#check for any 0s in permit matrix that are 1s in treatment matrix (should be o)
# Add temp_col as column 4 temporarily
filtered_permit_matrix <- filtered_permit_matrix %>% 
  mutate(temp_col = 0, .before = 4)

#Count the mismatches
mismatch_count <- sum(
  filtered_permit_matrix[, 6:13] == 0 & filtered_treatment_matrix[, 6:13] == 1,
  na.rm = TRUE
)
mismatch_count

#Create a logical matrix of where mismatches occur (cols 6-13) (lag1-8)
mismatch_matrix <- filtered_permit_matrix[, 6:13] == 0 & filtered_treatment_matrix[, 6:13] == 1
#Total count of mismatched cells
mismatch_count <- sum(mismatch_matrix, na.rm = TRUE)
cat("Total mismatched cells:", mismatch_count, "\n\n")
# 4. Extract rows where at least one mismatch occurred in cols 6-13
mismatch_row_indices <- which(rowSums(mismatch_matrix, na.rm = TRUE) > 0)
# 5. Get unique survey_ids for those rows
mismatched_survey_ids <- filtered_permit_matrix$survey_id[mismatch_row_indices]
# Display the survey_ids
mismatched_survey_ids


#and the other way around. This should be more than 0
mismatch_count2 <- sum(
  filtered_permit_matrix[, 6:13] == 1 & filtered_treatment_matrix[, 6:13] == 0,
  na.rm = TRUE
)


filtered_permit_matrix <- filtered_permit_matrix %>%
  add_column(temp1 = NA, temp2 = NA, temp3 = NA, temp4 = NA, .before = 4)


# 1. Identify logical condition (Permit == 1 & Treatment == 0) for cols 6 to 13
mismatch_mask <- filtered_permit_matrix[, 10:17] == 1 & filtered_treatment_matrix[, 10:17] == 0

# 2. Update filtered_treatment_matrix cells to 2 where condition is TRUE
filtered_treatment_matrix[, 10:17][mismatch_mask] <- 2

#check to make sure there are now 0 mismatched columns
mismatch_count2 <- sum(
  filtered_permit_matrix[, 10:17] == 1 & filtered_treatment_matrix[, 10:17] == 0,
  na.rm = TRUE
)

mismatch_count2


#remove temporary columns
filtered_permit_matrix <- filtered_permit_matrix %>% 
  select(-temp_col, -temp1, -temp2, -temp3, -temp4)

cytrt_matrix <- filtered_treatment_matrix %>%
  rowwise() %>%
  mutate(
    .vals = list(c_across(c(4, 10:17))),
    .first_nz = list(which(!is.na(.vals) & .vals != 0)[1]),
    
    # Calculate consecutive treatment years
    cytrt = ifelse(
      is.na(.first_nz), 
      0, 
      rle(.vals[.first_nz:length(.vals)] != 0 & !is.na(.vals[.first_nz:length(.vals)]))$lengths[1]
    )
  ) %>%
  select(-.vals, -.first_nz)


