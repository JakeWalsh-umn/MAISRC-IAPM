# this file organizes the water quality data for lakes in MN 
# the main output data are wq_april_september 

wq_data <- list.files(
  path = "data/wq_data", 
  pattern = "\\.csv$", 
  full.names = TRUE
) %>%
  read_csv(id = "source_file")

# select only the fields that we need for our analysis

wq_data_clean <- wq_data %>%
  dplyr::select("ActivityStartDate", "CharacteristicName", "ActivityDepthHeightMeasure/MeasureValue", 
         "ActivityDepthHeightMeasure/MeasureUnitCode", "MonitoringLocationIdentifier",
         "MonitoringLocationName", "ResultSampleFractionText", "ResultMeasureValue", 
         "ResultMeasure/MeasureUnitCode")

wq_data_clean <- wq_data_clean %>%
  rename(test_date = ActivityStartDate,
         parameter = CharacteristicName,
         depth_measure = "ActivityDepthHeightMeasure/MeasureValue",
         depth_unit = "ActivityDepthHeightMeasure/MeasureUnitCode",
         location = MonitoringLocationIdentifier,
         lk_name = MonitoringLocationName,
         resultsamplefraction = ResultSampleFractionText, 
         result = ResultMeasureValue,
         unit = "ResultMeasure/MeasureUnitCode")

wq_data_clean <- wq_data_clean %>%
  mutate(DOW = str_sub(str_remove_all(location, "\\D"), 1, 8)) %>%
  mutate(DOW_parent = str_replace(DOW, "\\d{2}$", "00")) 

wq_data_clean <- wq_data_clean %>%
  distinct()

wq_data_clean <- wq_data_clean %>%
  filter(!(result %in% no_data))

wq_month_yr <- wq_data_clean %>%
  mutate(
    year = year(test_date),
    month = month(test_date)
  ) %>%
  group_by(year, month, DOW_parent, lk_name, parameter, unit) %>%
  summarise(
    result = mean(result, na.rm = TRUE),
    .groups = "drop"
  )

wq_month_yr <- wq_month_yr %>% 
  mutate(parameter = case_match(
    parameter, 
    "Alkalinity, total"                       ~ "alkalinity_total", 
    "Chlorophyll a"                           ~ "chlorophyll_a", 
    "Chlorophyll a, corrected for pheophytin" ~ "chlorophyll_corrected", 
    "Depth, Secchi disk depth"                ~ "secchi", 
    "Total dissolved solids"                  ~ "TDS", 
    .default = parameter                      # Keeps any unmapped categories as-is 
  )) %>% 
  pivot_wider(
    names_from = parameter, 
    values_from = c(result, unit), 
    names_glue = "{parameter}_{.value}",
    names_vary = "slowest",
    # FIX: This flattens multiple entries into a single average value per cell
    values_fn = list(result = \(x) mean(x, na.rm = TRUE), unit = first) 
  ) %>% 
  rename_with(~ str_remove(., "_result$")) # Changed from _value$ to match your result column

wq_month_yr <- wq_month_yr %>%
  mutate(
    Phosphorus = if_else(Phosphorus_unit %in% c("ug/L"), Phosphorus / 1000, Phosphorus),
    # Update unit string to mg/l
    Phosphorus_unit = if_else(Phosphorus_unit %in% c("ug/L"), "mg/L", Phosphorus_unit)
  ) %>%
  mutate(
    # Multiply by 1000 if mg/L; leave mg/m3 and existing ug/L unchanged
    chlorophyll_a = if_else(
      chlorophyll_a_unit %in% c("mg/L"), 
      chlorophyll_a * 1000, 
      chlorophyll_a
    ),
    # Update unit labels to standard ug/L
    chlorophyll_a_unit = if_else(
      chlorophyll_a_unit %in% c("mg/L", "mg/m3"), 
      "ug/L", 
      chlorophyll_a_unit
    )
  ) %>%
  mutate(
    chlorophyll_corrected = if_else(chlorophyll_corrected_unit %in% c("mg/L"), chlorophyll_corrected * 1000, chlorophyll_corrected),
    # Update unit string to ug/l
    chlorophyll_corrected_unit = if_else(chlorophyll_corrected_unit %in% c("mg/L"), "ug/L", chlorophyll_corrected_unit)
  )



wq_april_september <- wq_month_yr %>%
  filter(month < 10 & month > 3)

# example graph to make sure the data structure is what I think it is ---- 
# NOTE: this makes me realize I should 

secchi_cedar <- wq_april_september %>%
  filter(DOW_parent == 27003900) %>%
  # If your data is in long format, uncomment the line below:
  # filter(parameter == "Secchi") %>% 
  
  # Create a unified date column if you only have year and month
  mutate(date = make_date(year, month, 1)) %>%
  drop_na(secchi)

secchi_cedar <- secchi_cedar %>%
  filter(DOW_parent == 27003900) %>%
  # Expand list elements into individual numeric rows
  unnest(secchi) %>% 
  mutate(
    secchi = as.numeric(secchi),
    date = make_date(year, month, 1)
  ) %>%
  drop_na(secchi)

# 2. Generate plot
ggplot(secchi_cedar, aes(x = date, y = secchi)) +
  # Points for individual samples
  geom_point(color = "#005a32", alpha = 0.6, size = 2) +
  
  # Faint line connecting observations over time
  geom_line(color = "#41ab5d", alpha = 0.4, linewidth = 0.6) +
  
  # Smooth trend line (LOESS curve to capture non-linear trends)
  geom_smooth(method = "loess", color = "#00441b", se = TRUE, fill = "#c7e9c0") +
  
  # Invert Y-axis: deeper Secchi readings (clearer water) extend further down
  scale_y_reverse() +
  
  # Format X-axis for multi-year time series
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  
  labs(
    title = paste("Secchi Transparency Over Time, Cedar Lake"),
    subtitle = "Inverted Y-axis: values further down represent deeper light penetration (clearer water)",
    x = "Year",
    y = "Secchi Depth (meters)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

