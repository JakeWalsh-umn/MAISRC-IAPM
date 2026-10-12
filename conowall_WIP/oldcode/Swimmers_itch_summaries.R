# swimmer's itch summaries

apm_iapm_permits <- read.csv('data/apm_iapm_permit_detail.csv')

apm_permits_filtered <- apm_iapm_permits %>%
  filter(treatment_method == "Pesticide Control" & permit_type == "Aquatic Plant Management") %>%
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

apm_permits_filtered <- apm_permits_filtered %>%
    group_by(permit_number) %>%
    arrange(
      desc(treatment_year),
      desc(coalesce(str_detect(species_groups, "Snails"), FALSE)),
      desc(yr_approved_acres),
      .by_group = TRUE
   ) %>%
    slice(1) %>%
    ungroup()

swimitch_permit_summary_yrs <- apm_permits_filtered %>%
  group_by(treatment_year) %>%
  summarize(
    total_permits = n_distinct(permit_number, na.rm = TRUE),
    snail_permits = n_distinct(permit_number[str_detect(species_groups, "Snails")], na.rm = TRUE),
    
    total_approved_acres = round(sum(yr_approved_acres, na.rm = TRUE), 2),
    snail_approved_acres = round(sum(yr_approved_acres[str_detect(species_groups, "Snails")], na.rm = TRUE), 2),
    
    .groups = "drop"
  ) %>%
  mutate(
    pct_snail_permits = round((snail_permits / total_permits) * 100, 2),
    pct_snail_approved_acres = round((snail_approved_acres / total_approved_acres) * 100, 2)
  )

swimitch_permit_summary_2017pres <- apm_permits_filtered %>%
  filter(treatment_year >= 2017) %>%
  summarize(
    total_permits = n_distinct(permit_number, na.rm = TRUE),
    snail_permits = n_distinct(permit_number[str_detect(species_groups, "Snails")], na.rm = TRUE),
    total_approved_acres = round(sum(yr_approved_acres, na.rm = TRUE), 2),
    snail_approved_acres = round((sum(yr_approved_acres[str_detect(species_groups, "Snails")], na.rm = TRUE)), 2)
  ) %>%
  mutate(
    pct_snail_permits = round((snail_permits / total_permits) * 100, 2),
    pct_snail_approved_acres = round((snail_approved_acres / total_approved_acres) * 100, 2)
  )

apm_permits_filtered %>%
  filter(treatment_year >= 2017) %>%
  filter(str_detect(species_groups, "Snail")) %>%
  summarize(unique_snail_permits = n_distinct(permit_number, na.rm = TRUE))

swimitch_permit_summary_yrs %>%
  write_csv("outputdata/swimitch_permit_summary_yrs.csv")

swimitch_permit_summary_2017pres %>%
  write_csv("outputdata/swimitch_permit_summary_2017pres.csv")
