library(sf)
hydro <- st_read("C:/Users/conow001/Downloads/Dnr_Hydro_Features_All/Dnr_Hydro_Features_All.shp")

# 1. Inspect sample values and string lengths
head(EWM_PI_data$DOW_parent)
head(hydro$DOW_parent) # or whatever the DOW column is named in the shapefile

# 2. Check for unmatched lakes
unmatched_dows <- setdiff(EWM_PI_data$DOW_parent, hydro$DOW_parent)
unmatched_dows <- setdiff(EWM_trts$DOW_parent, hydro$DOW_parent)
unmatched_dows <- setdiff(pi.ewmlks.1pre.2post$DOW_parent, hydro$DOW_parent)

hydro_unique <- hydro %>%
  distinct(DOW_parent, .keep_all = TRUE)

n_lakes <- MN_foo %>%
  filter(myriophyllum_spicatum > 0) %>% 
  distinct(DOW_parent) %>% 
  nrow()


pi_postsurv <- pi.trt.status %>%
  left_join(
    EWM_trts,
    by = join_by(DOW_parent, closest(SURVEY_START >= effective_end_date))
  ) %>% 
  filter()

hydro_joined <- hydro_unique %>%
  mutate(
    PI_data    = if_else(DOW_parent %in% all_PI$DOW_parent, 1, 0),
    EWM_PI     = if_else(DOW_parent %in% EWM_PI$DOW_parent, 1, 0),
    treated    = if_else(DOW_parent %in% EWM_trts$DOW_parent, 1, 0),
    pi_post    = if_else(DOW_parent %in% pi.tx$DOW_parent, 1, 0),
    pi_pre_post = if_else(DOW_parent %in% pi.ewmlks.1pre.2post$DOW_parent, 1, 0)
  )

hydro_joined %>%
  st_drop_geometry() %>%
  summarise(
    total_unique_lakes = n(),
    total_PI_data      = sum(PI_data),
    total_EWM_PI       = sum(EWM_PI),
    total_treated      = sum(treated),
    total_post         = sum(pi_post),
    total_pi_pre_post   = sum(pi_pre_post)
  )

st_write(
  hydro_joined, 
  "C:/Users/conow001/Downloads/Dnr_Hydro_Features_All/hydro_joined.shp", 
  delete_layer = TRUE
)

library(ggplot2)
library(dplyr)

# Create dataset with both metrics
df_attrition <- data.frame(
  stage = factor(
    c("1. Total EWM Infested Pool", 
      "2. PI Surveyed Lakes", 
      "3. Recorded Treatments", 
      "4. Treatment + Post-Survey", 
      "5. Pre + Treatment + Post-Survey"),
    levels = rev(c("1. Total EWM Infested Pool", 
                   "2. PI Surveyed Lakes", 
                   "3. Recorded Treatments", 
                   "4. Treatment + Post-Survey", 
                   "5. Pre + Treatment + Post-Survey"))
  ),
  lakes = c(400, 205, 164, 43, 25),
  treatments = c(NA, NA, 840, 195, 151)
) %>%
  mutate(
    # Create combined label for display
    label = ifelse(
      is.na(treatments),
      paste0(lakes, " Lakes"),
      paste0(lakes, " Lakes | ", treatments, " Treatments")
    )
  )

# Build plot
ggplot(df_attrition, aes(x = stage, y = lakes)) +
  geom_col(fill = "#2b5c8f", width = 0.55) +
  geom_text(aes(label = label), hjust = -0.1, color = "black", size = 3.8, fontface = "bold") +
  scale_y_continuous(limits = c(0, 480), expand = c(0, 0)) +
  coord_flip() +
  labs(
    title = "Spatiotemporal Pairing Data Attrition",
    subtitle = "Severe data reduction when enforcing temporal sequence constraints",
    x = NULL,
    y = "Unique Lakes (N)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

library(plotly)

# Node labels with treatment counts embedded in parentheses
nodes <- c(
  "Study Pool\n(400 Lakes)",                  # 0
  "No Treatment Record\n(236 Lakes)",         # 1
  "Recorded Treatments\n(164 Lakes | 840 Tx)",# 2
  "Missing Post-Survey\n(121 Lakes)",         # 3
  "Treated + Post\n(43 Lakes | 195 Tx)",      # 4
  "Missing Pre-Survey\n(18 Lakes)",           # 5
  "Complete Pair: Pre + Tx + Post\n(25 Lakes | 151 Tx)" # 6
)

df_filtered <- pi.ewmlks.1pre.2post %>%
  # Filter for cytrt = 0, 1, 2, 3
  filter(cytrt %in% c(0, 1, 2, 3)) %>%
  # Exclude missing values and categorize LZ coverage threshold
  filter(!is.na(prcnt_LZ)) %>%
  mutate(
    lz_group = case_when(
      prcnt_LZ < 0.15 ~ "< 15% LZ Treated",
      prcnt_LZ >= 0.15 ~ "> 15% LZ Treated"
    ),
    lz_group = factor(lz_group, levels = c("< 15% LZ Treated", "> 15% LZ Treated"))
  )

library(ggplot2)
library(dplyr)
library(ggnewscale) # Allows multiple color scales in one ggplot

# 1. Filter and split into two datasets
df_filtered <- pi.ewmlks.1pre.2post %>%
  filter(cytrt %in% c(0, 1, 2, 3), !is.na(prcnt_LZ)) %>%
  mutate(
    lz_group = case_when(
      prcnt_LZ < 0.15 ~ "< 15% LZ Treated",
      prcnt_LZ >= 0.15 ~ "> 15% LZ Treated"
    ),
    lz_group = factor(lz_group, levels = c("< 15% LZ Treated", "> 15% LZ Treated"))
  )

df_low  <- filter(df_filtered, lz_group == "< 15% LZ Treated")
df_high <- filter(df_filtered, lz_group == "> 15% LZ Treated")

# 2. Build plot with distinct color scales per panel
ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_boxplot(data = df_filtered, aes(x = factor(yrs.post.tx), y = EWMchangeFOO),
               fill = "gray95", color = "gray50", outlier.shape = NA, alpha = 0.5) +
  
  # Layer 1: Low coverage group (< 15%) using a Blue scale
  geom_jitter(data = df_low, aes(x = factor(yrs.post.tx), y = EWMchangeFOO, color = prcnt_LZ),
              width = 0.18, alpha = 0.8, size = 2.5) +
  scale_color_gradient(low = "#9ecae1", high = "#08306b", 
                       name = "% LZ (< 15%)", labels = scales::percent_format(accuracy = 1)) +
  
  # Reset color scale for second group
  new_scale_color() +
  
  # Layer 2: High coverage group (> 15%) using an Orange/Red scale
  geom_jitter(data = df_high, aes(x = factor(yrs.post.tx), y = EWMchangeFOO, color = prcnt_LZ),
              width = 0.18, alpha = 0.8, size = 2.5) +
  scale_color_gradient(low = "#fdd0a2", high = "#7f2704", 
                       name = "% LZ (> 15%)", labels = scales::percent_format(accuracy = 1)) +
  
  facet_wrap(~ lz_group, ncol = 2) +
  labs(
    title = "EWM Change in FOO Before vs. Post-Treatment",
    subtitle = "Separate color scales highlight variation within low vs. high coverage tiers",
    x = "Years Post-Treatment",
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

# 1. Calculate relative coverage within each group
df_rel <- pi.ewmlks.1pre.2post %>%
  filter(cytrt %in% c(0, 1, 2, 3), !is.na(prcnt_LZ)) %>%
  mutate(
    lz_group = case_when(
      prcnt_LZ < 0.15 ~ "< 15% LZ Treated",
      prcnt_LZ >= 0.15 ~ "> 15% LZ Treated"
    ),
    lz_group = factor(lz_group, levels = c("< 15% LZ Treated", "> 15% LZ Treated"))
  ) %>%
  group_by(lz_group) %>%
  mutate(
    # Rescales prcnt_LZ from 0 to 1 within each group
    lz_relative = (prcnt_LZ - min(prcnt_LZ)) / (max(prcnt_LZ) - min(prcnt_LZ))
  ) %>%
  ungroup()

# 2. Plot using relative color mapping
ggplot(df_rel, aes(x = factor(yrs.post.tx), y = EWMchangeFOO)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_boxplot(fill = "gray95", color = "gray50", outlier.shape = NA, alpha = 0.5) +
  
  # Color mapped to intra-group relative intensity
  geom_jitter(aes(color = lz_relative), width = 0.18, alpha = 0.8, size = 2.5) +
  
  facet_wrap(~ lz_group, ncol = 2) +
  scale_color_viridis_c(
    option = "viridis",
    name = "Relative Intensity\nWithin Group",
    labels = scales::percent
  ) +
  labs(
    title = "EWM Change in FOO Before vs. Post-Treatment",
    subtitle = "Color scale stretched across the min/max range within each treatment tier",
    x = "Years Post-Treatment",
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
library(ggnewscale) # Enables distinct color scales per group in a single plot

# 1. Calculate groups with updated thresholds (<8%, 8–15%, >15%)
df_rel <- pi.ewmlks.1pre.2post %>%
  filter(cytrt %in% c(0, 1, 2, 3), !is.na(prcnt_LZ)) %>%
  mutate(
    lz_group = case_when(
      prcnt_LZ < 0.08 ~ "< 8% LZ Treated",
      prcnt_LZ >= 0.08 & prcnt_LZ <= 0.15 ~ "8–15% LZ Treated",
      prcnt_LZ > 0.15 ~ "> 15% LZ Treated"
    ),
    lz_group = factor(
      lz_group, 
      levels = c("< 8% LZ Treated", "8–15% LZ Treated", "> 15% LZ Treated")
    )
  )

# Split data by group for separate color scale layers
df_g1 <- filter(df_rel, lz_group == "< 8% LZ Treated")
df_g2 <- filter(df_rel, lz_group == "8–15% LZ Treated")
df_g3 <- filter(df_rel, lz_group == "> 15% LZ Treated")

# 2. Plot with dedicated color scale per group
ggplot(df_rel, aes(x = factor(yrs.post.tx), y = EWMchangeFOO)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_boxplot(fill = "gray95", color = "gray50", outlier.shape = NA, alpha = 0.5) +
  
  # Group 1 (< 8% LZ): Blue gradient
  geom_jitter(data = df_g1, aes(color = prcnt_LZ), width = 0.18, alpha = 0.8, size = 2.5) +
  scale_color_distiller(palette = "Blues", direction = 1, name = "< 8% LZ\nCoverage", labels = scales::percent) +
  
  new_scale_color() +
  
  # Group 2 (8–15% LZ): Green gradient
  geom_jitter(data = df_g2, aes(color = prcnt_LZ), width = 0.18, alpha = 0.8, size = 2.5) +
  scale_color_distiller(palette = "Greens", direction = 1, name = "8–15% LZ\nCoverage", labels = scales::percent) +
  
  new_scale_color() +
  
  # Group 3 (> 15% LZ): Purple gradient
  geom_jitter(data = df_g3, aes(color = prcnt_LZ), width = 0.18, alpha = 0.8, size = 2.5) +
  scale_color_distiller(palette = "Purples", direction = 1, name = "> 15% LZ\nCoverage", labels = scales::percent) +
  
  facet_wrap(~ lz_group, ncol = 3) +
  labs(
    title = "EWM Change in FOO Before vs. Post-Treatment",
    subtitle = "Independent color gradients mapped to littoral zone coverage within each tier",
    x = "Years Post-Treatment",
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

df_1pre.2post.no0 <- pi.ewmlks.1pre.2post %>% 
  filter(EWM_foo > 0 & preEWM_foo > 0)

library(ggplot2)
library(dplyr)
library(ggnewscale) # Enables distinct color scales per group in a single plot

# 1. Calculate groups with updated thresholds (<8%, 8–15%, >15%)
df_rel <- df_1pre.2post.no0 %>%
  filter(cytrt %in% c(0, 1, 2, 3), !is.na(prcnt_LZ)) %>%
  mutate(
    lz_group = case_when(
      prcnt_LZ < 0.08 ~ "< 8% LZ Treated",
      prcnt_LZ >= 0.08 & prcnt_LZ <= 0.15 ~ "8–15% LZ Treated",
      prcnt_LZ > 0.15 ~ "> 15% LZ Treated"
    ),
    lz_group = factor(
      lz_group, 
      levels = c("< 8% LZ Treated", "8–15% LZ Treated", "> 15% LZ Treated")
    )
  )

# Split data by group for separate color scale layers
df_g1 <- filter(df_rel, lz_group == "< 8% LZ Treated")
df_g2 <- filter(df_rel, lz_group == "8–15% LZ Treated")
df_g3 <- filter(df_rel, lz_group == "> 15% LZ Treated")

# 2. Plot with dedicated color scale per group
ggplot(df_rel, aes(x = factor(yrs.post.tx), y = EWMchangeFOO)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_boxplot(fill = "gray95", color = "gray50", outlier.shape = NA, alpha = 0.5) +
  
  # Group 1 (< 8% LZ): Blue gradient
  geom_jitter(data = df_g1, aes(color = prcnt_LZ), width = 0.18, alpha = 0.8, size = 2.5) +
  scale_color_distiller(palette = "Blues", direction = 1, name = "< 8% LZ\nCoverage", labels = scales::percent) +
  
  new_scale_color() +
  
  # Group 2 (8–15% LZ): Green gradient
  geom_jitter(data = df_g2, aes(color = prcnt_LZ), width = 0.18, alpha = 0.8, size = 2.5) +
  scale_color_distiller(palette = "Greens", direction = 1, name = "8–15% LZ\nCoverage", labels = scales::percent) +
  
  new_scale_color() +
  
  # Group 3 (> 15% LZ): Purple gradient
  geom_jitter(data = df_g3, aes(color = prcnt_LZ), width = 0.18, alpha = 0.8, size = 2.5) +
  scale_color_distiller(palette = "Purples", direction = 1, name = "> 15% LZ\nCoverage", labels = scales::percent) +
  
  facet_wrap(~ lz_group, ncol = 3) +
  labs(
    title = "EWM Change in FOO Before vs. Post-Treatment",
    subtitle = "Independent color gradients mapped to littoral zone coverage within each tier",
    x = "Years Post-Treatment",
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
library(ggnewscale) # Allows multiple color scales in one ggplot

# 1. Filter and split into two datasets
df_filtered <- pre1.post3.tx %>%
  filter(cytrt %in% c(0, 1, 2, 3), !is.na(prcnt_LZ)) %>%
  mutate(
    lz_group = case_when(
      prcnt_LZ < 0.15 ~ "< 15% LZ Treated",
      prcnt_LZ >= 0.15 ~ "> 15% LZ Treated"
    ),
    lz_group = factor(lz_group, levels = c("< 15% LZ Treated", "> 15% LZ Treated"))
  )

df_low  <- filter(df_filtered, lz_group == "< 15% LZ Treated")
df_high <- filter(df_filtered, lz_group == "> 15% LZ Treated")

# 2. Build plot with distinct color scales per panel
ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_boxplot(data = df_filtered, aes(x = factor(yrs.post.tx), y = EWMchangeFOO),
               fill = "gray95", color = "gray50", outlier.shape = NA, alpha = 0.5) +
  
  # Layer 1: Low coverage group (< 15%) using a Blue scale
  geom_jitter(data = df_low, aes(x = factor(yrs.post.tx), y = EWMchangeFOO, color = prcnt_LZ),
              width = 0.18, alpha = 0.8, size = 2.5) +
  scale_color_gradient(low = "#9ecae1", high = "#08306b", 
                       name = "% LZ (< 15%)", labels = scales::percent_format(accuracy = 1)) +
  
  # Reset color scale for second group
  new_scale_color() +
  
  # Layer 2: High coverage group (> 15%) using an Orange/Red scale
  geom_jitter(data = df_high, aes(x = factor(yrs.post.tx), y = EWMchangeFOO, color = prcnt_LZ),
              width = 0.18, alpha = 0.8, size = 2.5) +
  scale_color_gradient(low = "#fdd0a2", high = "#7f2704", 
                       name = "% LZ (> 15%)", labels = scales::percent_format(accuracy = 1)) +
  
  facet_wrap(~ lz_group, ncol = 2) +
  labs(
    title = "EWM Change in FOO Before vs. Post-Treatment",
    subtitle = "Separate color scales highlight variation within low vs. high coverage tiers",
    x = "Years Post-Treatment",
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

# EWM plot
# 1. Categorize coverage and time tiers safely
df_summary <- pre1.post.tx %>%
  filter(!is.na(prcnt_LZ)) %>%
  mutate(
    # Corrected string labels so case_when and factor levels match exactly
    lz_group = case_when(
      prcnt_LZ < 0.08 ~ "< 8% LZ",
      prcnt_LZ >= 0.08 & prcnt_LZ <= 0.15 ~ "8-15% LZ",
      prcnt_LZ > 0.15 ~ "> 15% LZ"
    ),
    lz_group = factor(lz_group, levels = c("< 8% LZ", "8-15% LZ", "> 15% LZ")),
    
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

pd <- position_dodge(width = 0.25)

ggplot(df_summary, aes(x = yrs_cat, y = EWMchangeFOO, color = lz_group, fill = lz_group, group = lz_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Standard Error ribbons
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA, position = pd) +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, position = pd) +
  stat_summary(fun = mean, geom = "point", size = 3, position = pd) +
  
  # Sample size (n) text labels above each point
  stat_summary(
    fun.data = function(y) {
      data.frame(
        y = mean(y),
        label = paste0("n=", length(y))
      )
    },
    geom = "text",
    vjust = -1.3,           # Offsets text slightly above the point
    size = 3.5,
    fontface = "bold",
    position = pd,
    show.legend = FALSE
  ) +
  
  scale_color_brewer(palette = "Set1", name = "% LZ Treated") +
  scale_fill_brewer(palette = "Set1", name = "% LZ Treated") +
  # Expand y-axis limits slightly to ensure top labels aren't clipped
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
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

ggplot(pre.post.tx, aes(x = prcnt_LZ, y = EWMchangeFOO)) +
  # Zero reference line (points below 0 = EWM reduction)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  
  # Scatterplot of raw data points
  geom_point(alpha = 0.5, color = "#2b83ba", size = 2.5) +
  
  # Overall fitted linear trendline with confidence band
  # (Change method = "lm" to method = "loess" for a smooth non-linear curve)
  geom_smooth(method = "lm", color = "firebrick", fill = "firebrick", alpha = 0.2, linewidth = 1.2) +
  
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "EWM Change in FOO vs. Treatment Coverage",
    subtitle = "Fitted linear trend across all post-treatment observations",
    x = "% Littoral Zone Treated (prcnt_LZ)",
    y = "Change in FOO (% Post - % Pre)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

df_plot <- pre.post.tx %>%
  filter(!is.na(prcnt_LZ), !is.na(EWMchangeFOO)) %>%
  mutate(
    yrs_num = as.numeric(as.character(yrs.post.tx)),
    yrs_cat = case_when(
      yrs_num == 0 ~ "0",
      yrs_num == 1 ~ "1",
      yrs_num == 2 ~ "2",
      yrs_num >= 3 ~ "3+"
    ),
    yrs_cat = factor(yrs_cat, levels = c("0", "1", "2", "3+"))
  ) %>%
  filter(!is.na(yrs_cat))

# 2. Build scatterplot with separate fitted trendlines per yrs_cat
ggplot(df_plot, aes(x = prcnt_LZ, y = EWMchangeFOO, color = yrs_cat, fill = yrs_cat)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  
  # Raw data points
  geom_point(alpha = 0.5, size = 2.5) +
  
  # Fitted linear regression line with confidence intervals per yrs_cat
  # (Change method = "lm" to method = "loess" if you prefer non-linear curves)
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.2) +
  
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_brewer(palette = "Set1", name = "Years Post-Tx") +
  scale_fill_brewer(palette = "Set1", name = "Years Post-Tx") +
  labs(
    title = "EWM Change in FOO vs. Treatment Coverage",
    subtitle = "Fitted linear trends evaluated across post-treatment year categories (0, 1, 2, 3+)",
    x = "% Littoral Zone Treated (prcnt_LZ)",
    y = "Change in FOO (% Post - % Pre)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

#native richness plot
# 1. Categorize coverage and time tiers safely
df_summary <- pre1.post.tx %>%
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
ggplot(df_summary, aes(x = yrs_cat, y = nRich_change, color = lz_group, fill = lz_group, group = lz_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Standard Error ribbons representing raw data variance
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  scale_color_brewer(palette = "Set1", name = "% LZ Treated") +
  scale_fill_brewer(palette = "Set1", name = "% LZ Treated") +
  labs(
    title = "Native Richness Changes by Treatment Coverage",
    subtitle = "Lines show mean change in FOO; shaded ribbons display ±1 Standard Error",
    x = "Years Post-Treatment",
    y = "Mean Change in average rake richness (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# 2. Plot longitudinal summary trajectory
ggplot(df_summary, aes(x = yrs_cat, y = total_richchange, color = lz_group, fill = lz_group, group = lz_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Standard Error ribbons representing raw data variance
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  scale_color_brewer(palette = "Set1", name = "% LZ Treated") +
  scale_fill_brewer(palette = "Set1", name = "% LZ Treated") +
  labs(
    title = "Native Richness Changes by Treatment Coverage",
    subtitle = "Lines show mean change in FOO; shaded ribbons display ±1 Standard Error",
    x = "Years Post-Treatment",
    y = "Mean Change in average rake richness (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggplot(pre.post.tx, aes(x = prcnt_LZ, y = total_richchange)) +
  # Zero reference line (points below 0 = EWM reduction)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  
  # Scatterplot of raw data points
  geom_point(alpha = 0.5, color = "#2b83ba", size = 2.5) +
  
  # Overall fitted linear trendline with confidence band
  # (Change method = "lm" to method = "loess" for a smooth non-linear curve)
  geom_smooth(method = "lm", color = "firebrick", fill = "firebrick", alpha = 0.2, linewidth = 1.2) +
  
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "EWM Change in FOO vs. Treatment Coverage",
    subtitle = "Fitted linear trend across all post-treatment observations",
    x = "% Littoral Zone Treated (prcnt_LZ)",
    y = "Change in FOO (% Post - % Pre)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

#native abundance plot

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
ggplot(df_summary, aes(x = yrs_cat, y = nAb_change, color = lz_group, fill = lz_group, group = lz_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Standard Error ribbons representing raw data variance
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  scale_color_brewer(palette = "Set1", name = "% LZ Treated") +
  scale_fill_brewer(palette = "Set1", name = "% LZ Treated") +
  labs(
    title = "Native Abundance Changes by Treatment Coverage",
    subtitle = "Lines show mean change in FOO; shaded ribbons display ±1 Standard Error",
    x = "Years Post-Treatment",
    y = "Mean Change in average rake richness (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggplot(pre.post.tx, aes(x = prcnt_LZ, y = nAb_change)) +
  # Zero reference line (points below 0 = EWM reduction)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  
  # Scatterplot of raw data points
  geom_point(alpha = 0.5, color = "#2b83ba", size = 2.5) +
  
  # Overall fitted linear trendline with confidence band
  # (Change method = "lm" to method = "loess" for a smooth non-linear curve)
  geom_smooth(method = "lm", color = "firebrick", fill = "firebrick", alpha = 0.2, linewidth = 1.2) +
  
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "EWM Change in FOO vs. Treatment Coverage",
    subtitle = "Fitted linear trend across all post-treatment observations",
    x = "% Littoral Zone Treated (prcnt_LZ)",
    y = "Change in FOO (% Post - % Pre)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
C:/Users/conow001/Documents/GitHub/MAISRC-IAPM/figures/figures.R