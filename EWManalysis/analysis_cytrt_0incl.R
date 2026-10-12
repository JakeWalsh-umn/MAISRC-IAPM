library(tidyverse)
library(lubridate)
library(DescTools)
library(arrow)

source(here::here("1_Loading and preparing PI data.R"))

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

# collapse survey rows from a rake to a survey per row 
PI_rake_dates <- PI_rake_dates %>% 
  mutate(SURVEY_START = as.Date(SURVEY_START)) %>% 
  group_by(DOW, SURVEY_START) %>% 
  
  summarise(
    station_count = n(), #add number of stations in a survey
    
    # Calculate average native richness per rake
    avg.rich = mean(native.rich, na.rm = TRUE), 
    
    # Count of survey stations where at least one native plant was found
    nat.ab = sum(native.rich > 0, na.rm = TRUE),
    
    # Calculate the average depth that EWM is growing in lake
    ewm.avg.depth = mean(depth_ft[myriophyllum_spicatum > 0], na.rm = TRUE),
    
    # Sum plant columns only to get # of rakes with each plant
    across(any_of(taxonomic_columns), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    # Count unique native species detected anywhere in the lake survey
    total.rich = rowSums(across(any_of(natives)) > 0, na.rm = TRUE),
    .after = nat.ab
  ) %>%
  mutate(
    lake.rich = total.rich/station_count, #whole lake richness
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

# here add a join with WQ data if it's needed

MN_foo_plus <- MN_foo %>%
  mutate(
    SURVEY_START = as.Date(SURVEY_START),
    survey_year  = year(SURVEY_START),
    survey_month = month(SURVEY_START)
  )


########################

# set up dataframe with pre and post data ----

# PI surveys from 2012 on 
MN_foo_plus = MN_foo_plus %>% 
  filter(year(SURVEY_START) >= 2012)

# treatment records from 2017 on
EWM_trts <- EWM_trts %>%
  filter(
    treatment_year >= 2017,
    !(start_date %in% no_data)
  )

# format dates
EWM_trts <- EWM_trts %>% 
  mutate(end_date = as.Date(end_date, format = "%m/%d/%Y"),
          start_date = as.Date(start_date, format = "%m/%d/%Y"))

MN_foo_plus  <- MN_foo_plus  %>% mutate(SURVEY_START = as.Date(SURVEY_START))

# Match each treatment to the closest PI survey on or before the start date of treatment
#EWM_trts <- EWM_trts %>%
#  left_join(
#    MN_foo_plus %>% select(DOW_parent, SURVEY_START, survey_id),
#    by = join_by(DOW_parent, closest(start_date >= SURVEY_START)) # Retains all 3 surveys
#  ) %>% 
#  rename(
#    pre_surveyid = survey_id,         # define the pre-treatment survey id for later reference
#    pre_SURVEY_START = SURVEY_START   # and the start date
#  ) %>%
#  mutate(effective_end_date = end_date + days(14))

EWM_trts <- EWM_trts %>%
  left_join(
    MN_foo_plus %>% 
      # Exclude spring surveys (keep only surveys on or after June 15)
      filter(month(SURVEY_START) > 6 | (month(SURVEY_START) == 6 & day(SURVEY_START) >= 15)) %>%
      select(DOW_parent, SURVEY_START, survey_id),
    by = join_by(DOW_parent, closest(start_date >= SURVEY_START))
  ) %>% 
  rename(
    pre_surveyid = survey_id,         # define the pre-treatment survey id for later reference
    pre_SURVEY_START = SURVEY_START   # and the start date
  ) %>%
  mutate(effective_end_date = end_date + days(14))
  
EWM_trts <- EWM_trts %>%
  arrange(DOW_parent, end_date, desc(prcnt_LZ)) %>%
  distinct(DOW_parent, end_date, .keep_all = TRUE)

# label surveys as belonging to a treated or untreated lake
# NOTE this isn't fully accurate because it doesn't account for permit only records
pi.trt.status <- MN_foo_plus %>% 
  mutate(status = if_else(DOW_parent %in% EWM_trts$DOW_parent, "treated", "untreated"))


pi.trt.status$DOW_parent = as.character(pi.trt.status$DOW_parent)

# merge treatments to post survey(s)

#df <- pi.trt.status %>%
#  left_join(
#    EWM_trts,
#    by = join_by(DOW_parent, closest(SURVEY_START >= effective_end_date))
#  )

df <- pi.trt.status %>%
  # Exclude spring surveys (keep only surveys occurring on or after June 15)
  filter(month(SURVEY_START) > 6 | (month(SURVEY_START) == 6 & day(SURVEY_START) >= 15)) %>%
  left_join(
    EWM_trts,
    by = join_by(DOW_parent, closest(SURVEY_START >= effective_end_date))
  )


df$treatment_year = as.numeric(df$treatment_year)

# organize milfoil infestation years accounting for DOW standardization
EWM_iw_lks <- iw %>% 
  filter(tolower(species) == "eurasian watermilfoil") %>% 
  group_by(DOW_parent) %>% 
  arrange(year, .by_group = TRUE) %>% 
  slice(1) %>% 
  ungroup()

# join EWM infest data
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
df$EWM_foo = df$myriophyllum_spicatum/df$station_count

df_clean <- df %>%
  # Keep valid, non-negative treatment lags
  filter(!is.na(yrs.post.tx) & yrs.post.tx >= 0)
  

# separate out pre and post treatment surveys
#pre = df %>% filter(yrs.post.tx < 0)
#post = df %>% filter(yrs.post.tx > 0)


# filtered by predeteremined MN EWM infested lake niche
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
trt.lakes = unique(EWM_trts$DOW_parent)

#df_status <- df %>%
#  mutate(
#    DOW_parent = as.character(DOW_parent),
#    
#    tx_status = case_when(
 #     yrs.post.tx == 1                                                        ~ "treated_lastyr",
 #     DOW_parent %in% trt.lakes     & !is.na(myriophyllum_spicatum) ~ "invaded_treated",
#      !DOW_parent %in% permit_vector & !is.na(myriophyllum_spicatum) ~ "invaded_untreated",
#      TRUE                                                           ~ NA_character_
#    )
#  )


pi.tx <- df_clean %>%
  # Filter for valid post-treatment records
  filter(yrs.post.tx > 0 | (yrs.post.tx == 0 & SURVEY_START > end_date))


# separate out pre and post treatment surveys

pre <- df %>%
  filter(survey_id %in% pi.tx$pre_surveyid) 

pre = pre %>% 
  #group_by(DOW_parent) %>%
  #arrange(DOW_parent, yrs.post.tx) %>%
  group_by(survey_id) %>%
  slice_tail(n = 1) %>%
  ungroup()

#join pre treatment characteristics to the post data
paired_tx_surv <- pi.tx %>%
  left_join(
    pre %>%
      mutate(pre.avg.n.ab = nat.ab / station_count) %>% 
      select(
        survey_id,
        pre.avg.n.rich = avg.rich,
        pre.avg.n.ab,
        pre.ewm.avg.depth = ewm.avg.depth,
        preEWM = myriophyllum_spicatum,
        preEWM_foo = EWM_foo,
        pre.richness = total.rich
      ),
    by = c("pre_surveyid" = "survey_id")
  )

# define the number of years before treatment the pre-survey occurred
paired_tx_surv$pre.lag = year(paired_tx_surv$pre_SURVEY_START)- paired_tx_surv$treatment_year

# native richness across survey
paired_tx_surv$nat.richness = rowSums(dplyr::select(pi.tx,`bidens_beckii`:`zannichellia_palustris`) > 0, na.rm =T)

paired_tx_surv <- paired_tx_surv %>% 
  select(-c( year_confirm, YSI))

paired_tx_surv <- paired_tx_surv %>%
  left_join(
    EWM_iw_lks %>% 
      # Add dplyr:: to force the correct select function
      dplyr::select(DOW_parent, year_confirm = year),
    by = c("DOW_parent")
  ) 


# calculate time since treatment and time since invasion
paired_tx_surv$YSI = year(paired_tx_surv$SURVEY_START)- paired_tx_surv$year_confirm

pi.ewmlks <- paired_tx_surv %>% 
  filter(YSI > 0 | is.na(YSI))

pi.ewmlks <- pi.ewmlks %>%
  mutate(across(c(cytrt, prcnt_LZ), ~ replace_na(., 0)))

inspect <- pi.ewmlks %>% select(SURVEY_START, end_date, pre_SURVEY_START, treatment_id, pre.lag, yrs.post.tx)

pi.ewmlks$yrs.post.tx_fac <- factor(pi.ewmlks$yrs.post.tx)


pre.post.tx <- pi.ewmlks %>% 
  filter(EWM_foo > 0 & preEWM_foo > 0) %>%
  mutate(EWMchangeFOO = EWM_foo - preEWM_foo, 
         nRich_change = avg.rich - pre.avg.n.rich ,
         total_richchange = total.rich - pre.richness,
         nAb_change = (nat.ab/station_count) - pre.avg.n.ab)

pre1.post3.tx <- pre.post.tx %>% 
  filter (pre.lag > -2 & yrs.post.tx < 4)

pre1.post.tx <- pre.post.tx %>% 
  filter(pre.lag > -2)
#################

# EWM models: 

pi.ewmlks <- pi.ewmlks %>%
  mutate(
    yrs.post.tx_fac = fct_collapse(
      as.factor(yrs.post.tx_fac),
      "3+" = c("3", "4", "5")
    )
  )

ewm_1 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    prcnt_LZ + cytrt + yrs.post.tx +
    (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
)

ewm_2 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    prcnt_LZ*(cytrt + yrs.post.tx_fac) +
    (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
)


ewm_3 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    cytrt + prcnt_LZ*yrs.post.tx_fac + 
    (1 | DOW_parent/treatment_id),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
)

ewm_4 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    prcnt_LZ + cytrt*yrs.post.tx + 
    (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
)

ewm_4 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + (1 | DOW_parent/treatment_id),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
)

ewm_5 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + cytrt + (1 | DOW_parent/treatment_id),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # prcnt_LZ, yrs.post.tx
   # promising

ewm_6 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    yrs.post.tx_fac + prcnt_LZ + cytrt + (1 | DOW_parent/treatment_id),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) 



# native richness models ----
# might need to rethink this one

nat_rich_1 <- brm(
  total.rich | trials(station_count) ~ 
    prcnt_LZ + yrs.post.tx + cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # cytrt, 1 year post tx

nat_rich_2 <- brm(
  total.rich | trials(station_count) ~ 
    prcnt_LZ*(yrs.post.tx_fac + cytrt) + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # cytrt




# native abundance models ----
nat_ab_1 <- brm(
  nat.ab | trials(station_count) ~ 
    prcnt_LZ*(yrs.post.tx_fac + cytrt) + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # cytrt


nat_ab_2 <- brm(
  nat.ab | trials(station_count) ~ 
    prcnt_LZ + yrs.post.tx_fac + cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # cytrt, 1 year post tx

nat_ab_3 <- brm(
  nat.ab | trials(station_count) ~ 
    prcnt_LZ + yrs.post.tx_fac + cytrt + YSI + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # cytrt, 1 year post tx



# paired treatment plots with lines representing treatments and points representing 
# run a gam with intercept only random effect (1|lake+treatment) and then plot the average values with error bars.

# Here I will do some simple analysis ----




library(ggplot2)

ggplot(pi.ewmlks.1pre.2post, aes(x = factor(yrs.post.tx), y = EWMchangeFOO)) +
  # Box plots show the overall distribution per year
  geom_boxplot(alpha = 0.4, fill = "gray90", color = "gray30", outlier.shape = NA) +
  # Jittered points show individual observations colored by continuous % LZ
  geom_jitter(aes(color = prcnt_LZ), width = 0.2, size = 2.5, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  # Continuous color gradient scale
  scale_color_viridis_c(option = "viridis", name = "% Littoral\nZone Treated") +
  labs(
    title = "Change in EWM Frequency of Occurrence Post-Treatment",
    x = "Number of Years Post-Treatment",
    y = "Change in Frequency of Occurrence (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "right" # Display legend so color values are readable
  )

library(ggplot2)

ggplot(pi.ewmlks.1pre.2post, aes(
  x = factor(yrs.post.tx), 
  y = EWMchangeFOO, 
  fill = cut(prcnt_LZ, breaks = c(0, 25, 50, 75, 100), include.lowest = TRUE)
)) +
  geom_violin(alpha = 0.8, position = position_dodge(width = 0.8), trim = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_viridis_d(name = "% LZ Treated Tier", option = "C") +
  labs(
    title = "Change in EWM Frequency of Occurrence Post-Treatment",
    x = "Number of Years Post-Treatment",
    y = "Change in Frequency of Occurrence (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

library(ggplot2)

# Ensure each row has a unique identifier
pi.ewmlks.1pre.2post$tx_pair_id <- 1:nrow(pi.ewmlks.1pre.2post)

ggplot(pi.ewmlks.1pre.2post) +
  # Paired segment lines connecting Pre and Post
  geom_segment(
    aes(x = "Pre-Tx", xend = "Post-Tx", 
        y = preEWM_foo, yend = EWM_foo, 
        group = tx_pair_id, color = prcnt_LZ),
    alpha = 0.3, linewidth = 0.5
  ) +
  # End points for each line
  geom_point(aes(x = "Pre-Tx", y = preEWM_foo, color = prcnt_LZ), alpha = 0.4, size = 1.5) +
  geom_point(aes(x = "Post-Tx", y = EWM_foo, color = prcnt_LZ), alpha = 0.4, size = 1.5) +
  
  # Split into 4 facets (1 row of 4 plots, or 2x2 grid using ncol = 2)
  facet_wrap(~ factor(yrs.post.tx), ncol = 4, labeller = labeller(`factor(yrs.post.tx)` = function(x) paste("Year", x))) +
  
  # Continuous color gradient for treatment percentage
  scale_color_viridis_c(option = "viridis", name = "% Littoral\nZone Treated") +
  scale_x_discrete(expand = c(0.2, 0.2)) +
  labs(
    title = "Pre vs. Post-Treatment EWM Frequency of Occurrence",
    subtitle = "Faceted across 0–3 years post-treatment",
    x = "Treatment Stage",
    y = "EWM Frequency of Occurrence (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold")
  )


# 1. Reshape data into long format, retaining yrs.post.tx and prcnt_LZ
df_paired <- pi.ewmlks.1pre.2post %>%
  mutate(tx_pair_id = row_number()) %>%
  pivot_longer(
    cols = c(preEWM_foo, EWM_foo),
    names_to = "period",
    values_to = "foo_val"
  ) %>%
  # Set clean labels and correct factor order for the x-axis
  mutate(period = factor(period, levels = c("preEWM_foo", "EWM_foo"), labels = c("Pre-Tx", "Post-Tx")))

# 2. Plot with 4 faceted subplots and per-panel mean trend lines
ggplot(df_paired, aes(x = period, y = foo_val)) +
  # Individual site trajectories
  geom_line(aes(group = tx_pair_id, color = prcnt_LZ), alpha = 0.25, linewidth = 0.5) +
  geom_point(aes(color = prcnt_LZ), alpha = 0.3, size = 1.5) +
  
  # Bold mean trajectory calculated individually for each facet
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 1.3) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", color = "black", size = 3) +
  
  # Facet across years post-treatment (1 row of 4 plots)
  facet_wrap(
    ~ factor(yrs.post.tx), 
    ncol = 4, 
    labeller = labeller(`factor(yrs.post.tx)` = function(x) paste("Year", x))
  ) +
  
  # Scale and layout options
  scale_color_viridis_c(option = "viridis", name = "% Littoral\nZone Treated") +
  scale_x_discrete(expand = c(0.2, 0.2)) +
  labs(
    title = "Pre vs. Post-Treatment EWM Frequency of Occurrence",
    subtitle = "Black lines show the mean trajectory for each post-treatment year",
    x = "Treatment Stage",
    y = "EWM Frequency of Occurrence (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

df_paired <- pi.ewmlks.1pre.2post %>%
  mutate(
    tx_pair_id = row_number(),
    # Bin prcnt_LZ using thresholds: <= 0.05, 0.05–0.15, 0.15–1.0
    lz_tier = cut(
      prcnt_LZ,
      breaks = c(-Inf, 0.15, 1.0),
      labels = c("≤ 15%% Treated", "> 15% Treated"),
      include.lowest = TRUE
    )
  ) %>%
  pivot_longer(
    cols = c(preEWM_foo, EWM_foo),
    names_to = "period",
    values_to = "foo_val"
  ) %>%
  mutate(period = factor(period, levels = c("preEWM_foo", "EWM_foo"), labels = c("Pre-Tx", "Post-Tx")))

# 2. Plot with tier-specific mean trajectory lines across 4 post-treatment years
ggplot(df_paired, aes(x = period, y = foo_val)) +
  # Faint background lines for individual lake trajectories
  geom_line(aes(group = tx_pair_id), alpha = 0.12, color = "gray50", linewidth = 0.4) +
  
  # Distinct summary mean lines calculated per treatment tier
  stat_summary(
    aes(group = lz_tier, color = lz_tier), 
    fun = mean, 
    geom = "line", 
    linewidth = 1.3
  ) +
  stat_summary(
    aes(group = lz_tier, color = lz_tier), 
    fun = mean, 
    geom = "point", 
    size = 3
  ) +
  
  # Facet across 4 post-treatment years
  facet_wrap(
    ~ factor(yrs.post.tx), 
    ncol = 4, 
    labeller = labeller(`factor(yrs.post.tx)` = function(x) paste("Year", x))
  ) +
  
  # Colors and discrete scale
  scale_color_viridis_d(option = "D", name = "% Littoral Zone\nTreated Group") +
  scale_x_discrete(expand = c(0.2, 0.2)) +
  labs(
    title = "Pre vs. Post-Treatment EWM Frequency of Occurrence",
    subtitle = "Mean trajectories categorized by treatment intensity thresholds (0.05, 0.15, 1.0)",
    x = "Treatment Stage",
    y = "EWM Frequency of Occurrence (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )


library(ggplot2)

ggplot(pi.ewmlks.1pre.2post, aes(x = prcnt_LZ, y = EWMchangeFOO)) +
  # Zero-line reference (points below 0 = EWM reduction)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  
  # Individual lake observation points
  geom_point(aes(color = prcnt_LZ), alpha = 0.6, size = 2.5) +
  
  # Trend line showing dose-response relationship per year
  geom_smooth(method = "loess", color = "black", fill = "gray80", alpha = 0.5, se = TRUE) +
  
  # Facet by year post-treatment
  facet_wrap(
    ~ factor(yrs.post.tx), 
    ncol = 4, 
    labeller = labeller(`factor(yrs.post.tx)` = function(x) paste("Year", x))
  ) +
  
  scale_color_viridis_c(option = "viridis", name = "% Littoral\nZone Treated") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + # Formats 0.15 as 15%
  labs(
    title = "EWM Efficacy vs. Treatment Intensity Over Time",
    subtitle = "Dose-response curve (LOESS smooth) across 0–3 years post-treatment",
    x = "% Littoral Zone Treated",
    y = "Change in Frequency of Occurrence (Post - Pre %)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

library(ggplot2)
library(dplyr)

# Categorize treatment intensity into tiers
df_summary <- pi.ewmlks.1pre.2post %>%
  mutate(
    lz_tier = cut(
      prcnt_LZ,
      breaks = c(-Inf, 0.05, 0.15, 1.0),
      labels = c("Low (≤ 5%)", "Medium (5–15%)", "High (> 15%)"),
      include.lowest = TRUE
    )
  )

ggplot(df_summary, aes(x = yrs.post.tx, y = EWMchangeFOO, color = lz_tier, fill = lz_tier)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Smooth mean trajectory over time with standard error shading
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  scale_color_viridis_d(option = "D", name = "% Littoral Zone Treated") +
  scale_fill_viridis_d(option = "D", name = "% Littoral Zone Treated") +
  scale_x_continuous(breaks = 0:3) +
  labs(
    title = "Mean EWM Reduction by Treatment Intensity",
    subtitle = "Ribbons represent ±1 standard error of the mean",
    x = "Years Post-Treatment",
    y = "Mean Change in Frequency of Occurrence"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

library(ggplot2)
library(dplyr)

library(ggplot2)
library(dplyr)

# 1. Categorize cytrt into groups
df_plot <- pi.ewmlks.1pre.2post %>%
  mutate(
    cytrt_group = case_when(
      cytrt == 0 ~ "0",
      cytrt == 1 ~ "1",
      cytrt == 2 ~ "2",
      cytrt >= 3 ~ "3+",
      TRUE ~ as.character(cytrt)
    ),
    cytrt_group = factor(cytrt_group, levels = c("0", "1", "2", "3+"))
  )

# 2. Plot mean trend trajectory across 0–3 years
ggplot(df_plot, aes(x = yrs.post.tx, y = EWMchangeFOO, color = cytrt_group, group = cytrt_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Mean trajectory line and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  # Optional ±1 Standard Error ribbons
  stat_summary(aes(fill = cytrt_group), fun.data = mean_se, geom = "ribbon", alpha = 0.15, color = NA) +
  
  scale_color_viridis_d(option = "D", name = "Treatment Count") +
  scale_fill_viridis_d(option = "D", name = "Treatment Count") +
  scale_x_continuous(breaks = 0:3) +
  labs(
    title = "EWM Trajectory Over Time by Treatment Count Group",
    subtitle = "Lines display group means with ±1 SE shaded ribbons",
    x = "Years Post-Treatment",
    y = "Mean Change in Frequency of Occurrence (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )



library(ggplot2)
library(dplyr)

# 1. Prep data with binned LZ coverage and cytrt tiers
df_trajectories <- pi.ewmlks.1pre.2post %>%
  mutate(
    lz_tier = cut(
      prcnt_LZ, 
      breaks = c(-Inf, 0.05, 0.20, 1.0), 
      labels = c("Low Coverage (≤5%)", "Med Coverage (5-20%)", "High Coverage (>20%)")
    ),
    cytrt_group = case_when(
      cytrt <= 1 ~ "Low Prior Tx (0–1)",
      cytrt >= 2 ~ "High Prior Tx (2+)"
    )
  )

# 2. Longitudinal trend profile
ggplot(df_trajectories, aes(x = yrs.post.tx, y = EWMchangeFOO, color = lz_tier, fill = lz_tier)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Mean trajectories with SE ribbons
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.15, color = NA) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  
  # Facet rows by Prior Treatment History
  facet_wrap(~ cytrt_group, ncol = 2) +
  
  scale_color_viridis_d(option = "D", name = "Treatment Coverage") +
  scale_fill_viridis_d(option = "D", name = "Treatment Coverage") +
  scale_x_continuous(breaks = 0:3) +
  labs(
    title = "EWM Rebound Rate Across Treatment Tiers",
    subtitle = "Tracking change in FOO over 0–3 years post-treatment",
    x = "Years Post-Treatment",
    y = "Mean Change in FOO (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

library(ggplot2)

ggplot(pi.ewmlks.1pre.2post, aes(x = prcnt_LZ, y = EWMchangeFOO, color = factor(yrs.post.tx), fill = factor(yrs.post.tx))) +
  # Zero-line reference (values below 0 indicate EWM reduction)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  
  # Scatter points for individual observations colored by year
  geom_point(alpha = 0.4, size = 2) +
  
  # Separate trend line for each year post-treatment
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.2) +
  
  # Discrete color palette for the 4 years
  scale_color_viridis_d(option = "viridis", name = "Years Post-Tx") +
  scale_fill_viridis_d(option = "viridis", name = "Years Post-Tx") +
  
  # Formats 0.15 as 15% (remove labels function if prcnt_LZ is stored as whole numbers like 15)
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "EWM Efficacy vs. Treatment Coverage by Post-Treatment Year",
    subtitle = "Unique trend lines show how the dose-response relationship shifts from Year 0 to Year 3",
    x = "% Littoral Zone Treated",
    y = "Change in Frequency of Occurrence (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

library(ggplot2)

ggplot(pi.ewmlks.1pre.2post, aes(x = prcnt_LZ, y = EWMchangeFOO)) +
  # Reference line: points below 0 indicate EWM reduction
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  
  # Raw data points
  geom_point(aes(color = prcnt_LZ), alpha = 0.5, size = 2) +
  
  # Empirical smooth curve (LOESS) showing non-linear trend per year
  geom_smooth(method = "loess", color = "black", fill = "gray80", alpha = 0.5, se = TRUE) +
  
  # Facet across post-treatment years
  facet_wrap(
    ~ factor(yrs.post.tx), 
    ncol = 4, 
    labeller = labeller(`factor(yrs.post.tx)` = function(x) paste("Year", x))
  ) +
  
  scale_color_viridis_c(option = "viridis", name = "% Littoral\nZone Treated") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "EWM Response vs. Treatment Coverage Over Time",
    subtitle = "Points below zero show EWM reduction; black curves show non-linear empirical trend (LOESS)",
    x = "% Littoral Zone Treated",
    y = "Change in FOO (% Post - % Pre)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

library(ggplot2)
library(dplyr)

# 1. Categorize cytrt into clear historical groups
df_summary <- pre1.post3.tx %>%
  mutate(
    cytrt_group = case_when(
      cytrt == 1 ~ "0 Prior Tx",
      cytrt == 2 ~ "1 Prior Tx",
      cytrt >= 3 ~ "2+ Prior Tx"
    ),
    cytrt_group = factor(cytrt_group, levels = c("0 Prior Tx", "1 Prior Tx", "2+ Prior Tx"))
  )

# 2. Plot longitudinal summary trajectory
ggplot(df_summary, aes(x = yrs.post.tx, y = EWMchangeFOO, color = cytrt_group, fill = cytrt_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  # Standard Error ribbons representing raw data variance
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  scale_color_brewer(palette = "Set1", name = "Treatment History") +
  scale_fill_brewer(palette = "Set1", name = "Treatment History") +
  scale_x_continuous(breaks = 0:3) +
  labs(
    title = "EWM Rebound Trajectory by Treatment History",
    subtitle = "Lines show mean change in FOO; shaded ribbons display ±1 Standard Error",
    x = "Years Post-Treatment",
    y = "Mean Change in FOO (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
library(ggplot2)
library(dplyr)

library(ggplot2)
library(dplyr)

# 1. Categorize cytrt AND group yrs.post.tx_fac into 0, 1, 2, 3+
df_summary <- pre1.post3.tx %>%
  mutate(
    # Treatment history groups
    cytrt_group = case_when(
      cytrt == 1 ~ "0 Prior Tx",
      cytrt == 2 ~ "1 Prior Tx",
      cytrt >= 3 ~ "2+ Prior Tx"
    ),
    cytrt_group = factor(cytrt_group, levels = c("0 Prior Tx", "1 Prior Tx", "2+ Prior Tx")),
    
    # Categorize yrs.post.tx_fac safely whether it's factor, character, or numeric
    yrs_cat = case_when(
      yrs.post.tx == 0 ~ "0",
      yrs.post.tx == 1 ~ "1",
      yrs.post.tx == 2 ~ "2",
      yrs.post.tx >= 3 ~ "3+"
    ),
    yrs_cat = factor(yrs_cat, levels = c("0", "1", "2", "3+"))
  ) %>%
  filter(!is.na(yrs_cat))

# 2. Plot longitudinal summary trajectory
ggplot(df_summary, aes(x = yrs_cat, y = EWMchangeFOO, color = cytrt_group, fill = cytrt_group, group = cytrt_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  # Standard Error ribbons representing raw data variance
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  scale_color_brewer(palette = "Set1", name = "Treatment History") +
  scale_fill_brewer(palette = "Set1", name = "Treatment History") +
  labs(
    title = "EWM Rebound Trajectory by Treatment History",
    subtitle = "Lines show mean change in FOO; shaded ribbons display ±1 Standard Error",
    x = "Years Post-Treatment",
    y = "Mean Change in FOO (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

library(ggplot2)
library(dplyr)

# 1. Categorize prcnt_LZ into coverage groups
library(ggplot2)
library(dplyr)

# 1. Categorize coverage and time tiers safely
df_summary <- pre.post.tx %>%
  filter(!is.na(prcnt_LZ)) %>%
  mutate(
    # Corrected string labels so case_when and factor levels match exactly
    lz_group = case_when(
      prcnt_LZ < 0.15 ~ "< 15% LZ",
      #prcnt_LZ >= 0.08 & prcnt_LZ <= 0.15 ~ "8–15% LZ",
      prcnt_LZ >= 0.15 ~ ">= 15% LZ"
    ),
    lz_group = factor(lz_group, levels = c("< 15% LZ", ">= 15% LZ")),
    
    # Safe conversion of yrs.post.tx in case it is stored as factor or character
    yrs_num = as.numeric(as.character(yrs.post.tx)),
    yrs_cat = case_when(
      yrs_num == 0 ~ "0",
      yrs_num == 1 ~ "1",
      yrs_num == 2 ~ "2",
      yrs_num >= 3 ~ "3+"
    ),
    yrs_cat = factor(yrs_cat, levels = c("0", "1", "2", "3+"))
  ) %>%
  filter(!is.na(lz_group), !is.na(yrs_cat))


# 2. Plot longitudinal summary trajectory
ggplot(df_summary, aes(x = yrs_cat, y = EWMchangeFOO, color = lz_group, fill = lz_group, group = lz_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Standard Error ribbons representing raw data variance
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  scale_color_brewer(palette = "Set1", name = "% LZ Treated") +
  scale_fill_brewer(palette = "Set1", name = "% LZ Treated") +
  labs(
    title = "EWM Rebound Trajectory by Treatment Coverage",
    subtitle = "Lines show mean change in FOO; shaded ribbons display ±1 Standard Error",
    x = "Years Post-Treatment",
    y = "Mean Change in FOO (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# 1. Categorize prcnt_LZ into coverage groups


library(ggplot2)
library(dplyr)



df_summary <- pi.ewmlks.1pre.2post %>%
  mutate(
    cytrt_group = case_when(
      cytrt == 0 ~ "0 Prior Tx",
      cytrt == 1 ~ "1 Prior Tx",
      cytrt >= 2 ~ "2+ Prior Tx"
    ),
    cytrt_group = factor(cytrt_group, levels = c("0 Prior Tx", "1 Prior Tx", "2+ Prior Tx"))
  )

ggplot(df_summary, aes(x = yrs.post.tx, y = change.n.ab, color = cytrt_group, fill = cytrt_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  # Standard Error ribbons representing raw data variance
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  scale_color_brewer(palette = "Set1", name = "Treatment History") +
  scale_fill_brewer(palette = "Set1", name = "Treatment History") +
  scale_x_continuous(breaks = 0:3) +
  labs(
    title = "Native abundance changes by EWM Treatment History",
    subtitle = "Lines show mean change in abundance; shaded ribbons display ±1 Standard Error",
    x = "Years Post-Treatment",
    y = "Mean Change in Abundance (change.n.ab)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggplot(df_summary, aes(x = yrs.post.tx, y = change.n.rich, color = cytrt_group, fill = cytrt_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  # Standard Error ribbons representing raw data variance
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  scale_color_brewer(palette = "Set1", name = "Treatment History") +
  scale_fill_brewer(palette = "Set1", name = "Treatment History") +
  scale_x_continuous(breaks = 0:3) +
  labs(
    title = "Native richness changes by EWM Treatment History",
    subtitle = "Lines show mean change in abundance; shaded ribbons display ±1 Standard Error",
    x = "Years Post-Treatment",
    y = "Mean Change in Abundance (change.n.ab)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

all_PI <- PI_rake_dates %>%
  mutate(DOW_parent = str_replace(DOW, "\\d{2}$", "00")) %>%
  add_count(DOW_parent, name = "n_surveys")

all_PI_data %>%
  write_csv("outputdata/EWM_PI_data.csv")

EWM_trts %>% 
  write_csv("outputdata/EWM_trts.csv")

pi.ewmlks.1pre.2post %>% 
  write_csv("outputdata/EWM_pi_tx.csv")

n_distinct(pi.ewmlks$DOW_parent)
