library(tidyverse)
library(lubridate)
library(DescTools)
library(arrow)

source(here::here("1_Loading and preparing PI data.R"))

taxonomic_columns <- directory_taxonomic$fieldNames


PI_rake_dates <- PI_rake_dates %>% 
  mutate(
    lake.rich = total.rich/station_count
  )

taxonomic_columns <- directory_taxonomic$fieldNames

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

PI_rake_dates <- PI_rake_dates %>% 
  mutate(SURVEY_START = as.Date(SURVEY_START)) %>% 
  group_by(DOW, SURVEY_START) %>% 
  
  summarise(
    station_count = n(),
    
    avg.rich = mean(native.rich, na.rm = TRUE), 
    
    # Count of survey stations where at least one native plant was found
    nat.ab = sum(native.rich > 0, na.rm = TRUE),
    
    ewm.avg.depth = mean(depth_ft[myriophyllum_spicatum > 0], na.rm = TRUE),
    
    # Sum plant columns only
    across(any_of(taxonomic_columns), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    # Count unique native species detected anywhere in the lake survey
    total.rich = rowSums(across(any_of(natives)) > 0, na.rm = TRUE),
    .after = nat.ab
  ) %>%
  mutate(
    lake.rich = total.rich/station_count,
    .after = total.rich
  ) %>%
  mutate(survey_id = row_number()) %>% 
  select(survey_id, station_count, everything())


tax_cols <- intersect(taxonomic_columns, names(PI_rake_dates))

MN_foo <- PI_rake_dates %>%
  mutate(DOW_parent = str_replace(DOW, "\\d{2}$", "00"))



source(here::here("EWM_Brooks_Model/2_MN_IAPM_trts.R"))

source(here::here("EWM_Brooks_Model/3_cytrt.R"))

source(here::here("EWM_Brooks_Model/4_infested_waters.R"))

source(here::here("EWM_Brooks_Model/5_MN_lake_chem.R"))


# join cytrt matrix and water quality matrix to PI data ----

MN_foo <- MN_foo %>% 
  left_join(
    cytrt_matrix %>% 
      select(1, 4:18),
    by = c("survey_id")
  )

wq_cols <- c("pH", "Phosphorus", "chlorophyll_corrected", "alkalinity_total", 
             "chlorophyll_a", "secchi", "TDS")

MN_foo_plus <- MN_foo %>%
  # Extract survey year and month from survey date
  select(-any_of(c("year", "month"))) %>%
  mutate(
    SURVEY_START = as.Date(SURVEY_START),
    survey_year  = year(SURVEY_START),
    survey_month = month(SURVEY_START)
  ) %>%
  # Join WQ data by lake identifier
  left_join(wq_april_september, by = "DOW_parent", relationship = "many-to-many") %>%
  # Filter: keep prior years (all months) OR survey year (prior months only)
  filter(
    (survey_year - year) <= 5 & 
      (year < survey_year | (year == survey_year & month < survey_month))
  ) %>%
  # Calculate relative year (0 = survey year, 1 = 1 year prior, etc.)
  mutate(relative_year = paste0("year_", survey_year - year)) %>%
  # Group per survey event and water quality year
  group_by(SURVEY_START, DOW_parent, relative_year) %>% #I previously had depth here which may have been a good grouping factor
  # Average all numeric water quality parameters
  summarise(
    across(where(is.numeric) & !c(year, month, survey_year, survey_month), 
           ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = relative_year,
    values_from = any_of(wq_cols),
    names_glue = "{.value}{relative_year}"
  )


########################

# pre and post analysis ----

# filter and format data
MN_foo_plus = MN_foo_plus %>% 
  filter(year(SURVEY_START) >= 2012)

MN_trts_acre <- MN_trts_acre %>%
  filter(
    treatment_year >= 2017,
    !(start_date %in% no_data)
  )

MN_trts_acre <- MN_trts_acre %>% 
  mutate(end_date = as.Date(end_date, format = "%m/%d/%Y")) %>%
  rename(pre_SURVEY_START = SURVEY_START)

MN_foo_plus  <- MN_foo_plus  %>% mutate(SURVEY_START = as.Date(SURVEY_START))

# Match each treatment to the closest survey on or before end_date
MN_trts_acre <- MN_trts_acre %>%
  left_join(
    MN_foo_plus %>% select(DOW_parent, SURVEY_START, survey_id),
    by = join_by(DOW_parent, closest(end_date <= SURVEY_START))
  ) %>% 
  rename(pre_surveyid = survey_id)

# filter foo surveys that have been treated

pi.trt.status <- MN_foo_plus %>% 
  mutate(status = if_else(DOW_parent %in% MN_trts_acre$DOW_parent, "treated", "untreated"))
         
         
pi.trt.status$DOW_parent = as.character(pi.trt.status$DOW_parent)

# merge
df = left_join(pi.trt.status, MN_trts_acre, by = "DOW_parent", relationship = "many-to-many")
df$treatment_year = as.numeric(df$treatment_year)


EWM_iw_lks <- iw %>% 
  filter(tolower(species) == "eurasian watermilfoil") %>% 
  group_by(DOW_parent) %>% 
  arrange(year, .by_group = TRUE) %>% 
  slice(1) %>% 
  ungroup()

df <- df %>%
  left_join(
    EWM_iw_lks %>% 
      # Add dplyr:: to force the correct select function
      dplyr::select(DOW_parent, connected, year_confirm = year),
    by = c("DOW_parent")
  ) 


# calculate time since treatment and time since invasion
df$yrs.post.tx <- year(df$SURVEY_START) - df$treatment_year
df$YSI = df$treatment_year - df$year_confirm

# separate out pre and post treatment surveys
#pre = df %>% filter(yrs.post.tx < 0)
#post = df %>% filter(yrs.post.tx > 0)



EWM_nichelks <- hydrography %>% 
  filter(lk_area < 2000 & lk_avdepth < 20)

# setting up the status for lake surveys ----
# NOTE: adapt the following code to help sort out which lakes should be part of the analysis and which ones shouldn't
# IDEA: add a category that is treated within year, treated ever, invaded not treated, not EWM niche lake 
#### for most analysis we might only care about treated within year and invaded not treated, but others may be interesting.


EWM_nichelks <- EWM_nichelks %>%
  mutate(DOW_parent = str_replace(dowlknum, "\\d{2}$", "00"))

iapm_permits <- iapm_permits %>% 
  mutate(DOW_parent = str_replace(water_resource_numbers, "\\d{2}$", "00"))

niche_vector   <- EWM_nichelks$DOW_parent
permit_vector  <- iapm_permits$DOW_parent  # update column name if different
trt.lakes = unique(MN_trts_acre$DOW_parent)

df_status <- df %>%
  mutate(
    DOW_parent = as.character(DOW_parent),
    
    tx_status = case_when(
      yrs.post.tx == 1                                                        ~ "treated_lastyr",
      DOW_parent %in% trt.lakes     & !is.na(myriophyllum_spicatum) ~ "invaded_treated",
      !DOW_parent %in% permit_vector & !is.na(myriophyllum_spicatum) ~ "invaded_untreated",
      TRUE                                                           ~ NA_character_
    )
  )

pi.dedup = df_status %>% group_by(DOW_parent, SURVEY_START) %>%
  filter(yrs.post.tx > 0 | is.na(yrs.post.tx)) %>%
  arrange(DOW_parent, yrs.post.tx) %>%
  slice_head(n = 1) %>%
  ungroup() 



# separate out pre and post treatment surveys
post = pi.dedup %>% filter(yrs.post.tx > 0)

pre <- df_status %>%
  filter(survey_id %in% post$pre_surveyid)

pre = pre %>% group_by(DOW_parent) %>%
  arrange(DOW_parent, yrs.post.tx) %>%
  slice_tail(n = 1) %>%
  ungroup()

#join pre treatment characteristics to the post data
df <- post %>%
  left_join(
    pre %>%
      mutate(pre.avg.n.ab = nat.ab / station_count) %>% 
      select(
        survey_id,
        pre.avg.n.rich = avg.rich,
        pre.avg.n.ab,
        pre.ewm.avg.depth = ewm.avg.depth,
        preEWM = myriophyllum_spicatum
      ),
    by = c("pre_surveyid" = "survey_id")
  )

df$pre.lag = df$treatment_year - year(df$pre_SURVEY_START)


df$nat.richness = rowSums(dplyr::select(pi.dedup,`bidens_beckii`:`zannichellia_palustris`) > 0, na.rm =T)

df <- df %>% 
  select(-c( year_confirm, YSI))
  
df <- df %>%
  left_join(
    EWM_iw_lks %>% 
      # Add dplyr:: to force the correct select function
      dplyr::select(DOW_parent, year_confirm = year),
    by = c("DOW_parent")
  ) 


# calculate time since treatment and time since invasion
df$YSI = year(df$SURVEY_START)- df$year_confirm

pi.ewmlks <- df %>% 
  filter(YSI > 0 | is.na(YSI))

pi.ewmlks <- pi.ewmlks %>%
  mutate(across(c(cytrt, prcnt_LZ), ~ replace_na(., 0)))

# brm modeling

# Data preparation using raw, unscaled units
pi.ewmlks <- pi.ewmlks %>%
  mutate(
    yr.post.tx_fact = factor(yrs.post.tx), # Post-treatment timing factor (1-5)
    cytrt           = as.numeric(cytrt),  # Raw count: consecutive years treated (1, 2, 3...)
    YSI             = as.numeric(YSI),    # Raw count: years since invasion (1-30)
    prcnt_LZ        = as.numeric(prcnt_LZ)# Raw percentage (0-100 or 0-1)
  )

# Model in raw units
ewm_brm <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
)

ewm_brm_s <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    s(yrs.post.tx, k = 4) + prcnt_LZ + cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  backend = "cmdstanr"
)

ewm_brm_interact <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    yrs.post.tx * prcnt_LZ + 
    prcnt_LZ * cytrt + 
    (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
)

# leave one out comparison to compare models
loo_compare(add_criterion(ewm_brm_s, "loo"), 
            add_criterion(ewm_brm_interact, "loo"))





# just a bunch of junk below here

library(brms)

install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE)
install_cmdstan()

# next step is to add some pre-treatment variables like richness to these models.

ewm_brm <- brm(
  myriophyllum_spicatum | trials(station_count) ~ prcnt_LZ + YSI + lag + cytrt + (1 | DOW),
  family = binomial(link = "logit"),
  data = pi.dedup,
  prior = c(
    prior(normal(0, 1.5), class = "b"),
    prior(student_t(3, 0, 2.5), class = "sd")
  ),
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
)

ewm_brm2 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ lag_cat*(prcnt_LZ + cytrt) + YSI + (1 | DOW_parent),
  family = binomial(link = "logit"),
  data = pi.ewmlks,
  prior = c(
    prior(normal(0, 1.5), class = "b"),
    prior(student_t(3, 0, 2.5), class = "sd")
  ),
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
)



