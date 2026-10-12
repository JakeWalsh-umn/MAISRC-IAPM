# trajectory plots using % littoral zone, years post treatment and change in foo or total richness 

pre1.post.tx <- pre1.post.tx %>%
  mutate(nAb_change100 = nAb_change*100,
         EWMchange100 = EWMchangeFOO*100, 
         prcntchangeEWM = (EWM_foo/preEWM_foo)*100,
         relprcntchangeEWM = (EWM_foo-preEWM_foo)/preEWM_foo*100,
         relprcntchangeNAB = nAb_change/pre.avg.n.ab*100)


model <- lm(nAb_change100 ~ EWMchangeFOO, data = pre1.post.tx)


# setup ---------------------------------------------
df_summary <- pre1.post.tx %>%
  filter(!is.na(prcnt_LZ)) %>%
  mutate(
    # Corrected string labels so case_when and factor levels match exactly
   # lz_group = case_when(
   #   prcnt_LZ < 0.15 ~ "< 15%",
   #   prcnt_LZ >= 0.15 ~ ">= 15%"
   # ),
   # lz_group = factor(lz_group, levels = c("< 15%", ">= 15%")),
    
    lz_group = case_when(
      prcnt_LZ < 0.05 ~ "< 5%",
      prcnt_LZ >= 0.05 & prcnt_LZ <= 0.15 ~ "5-15%",
      prcnt_LZ > 0.15 ~ "> 15%"
   ),
   lz_group = factor(lz_group, levels = c("< 5%", "5-15%", "> 15%")),
    
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

is_not_outlier <- function(x, coef = 3) {
  if (all(is.na(x))) return(rep(TRUE, length(x)))
  q25 <- quantile(x, 0.25, na.rm = TRUE)
  q75 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q75 - q25
  lower <- q25 - (coef * iqr)
  upper <- q75 + (coef * iqr)
  is.na(x) | (x >= lower & x <= upper)
}

# Target variables for outlier removal
target_vars <- c("EWMchangeFOO", "relprcntchangeEWM", "total_richchange", "nAb_change", "relprcntchangeNAB")

df_summary_clean <- df_summary %>%
  # Apply 3x IQR filter within each post-treatment year group to preserve group-level variance
  group_by(yrs_cat) %>%
  filter(if_all(all_of(target_vars), is_not_outlier)) %>%
  ungroup()

# EWM plot --------------------------------------------------------

pd <- position_dodge(width = 0.25)

ggplot(df_summary_clean, aes(x = yrs_cat, y = EWMchangeFOO, color = lz_group, fill = lz_group, group = lz_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Standard Error ribbons
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
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
  
  scale_color_brewer(palette = "Set1", direction = -1, name = "% Treated") +
  scale_fill_brewer(palette = "Set1", direction = -1, name = "% Treated") +
  # Expand y-axis limits slightly to ensure top labels aren't clipped
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title = "Change in EWM frequency of occurrence by % littoral zone treated",
    x = "Years Post-Treatment",
    y = "Mean Change in FOO (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggplot(df_summary_clean, aes(x = yrs_cat, y = relprcntchangeEWM, color = lz_group, fill = lz_group, group = lz_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Standard Error ribbons
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.20, color = NA) +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  # Sample size (n) text labels above each point
  stat_summary(
    fun.data = function(y) {
      data.frame(
        y = mean(y),
        label = paste0("n=", length(y))
      )
    },
    geom = "text",
    vjust = -1.3,
    size = 4,
    fontface = "bold",
    position = pd,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c("< 5%" = "#B07AA1", "5-15%" = "#F28E2B", "> 15%" =  "#174D79"), 
    name = "% Littoral \nZone Treated"
  ) +
  scale_fill_manual(
    values = c("< 5%" = "#B07AA1", "5-15%" = "#F28E2B", "> 15%" =  "#174D79"), 
    name = "% Littoral \nZone Treated"
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title = "Change in EWM Frequency of Occurrence \nafter EWM treatment",
    x = "# Years Post-Treatment",
    y = "Mean Change in EWM FOO (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # 1. Bigger Title (add `size = ...`)
    plot.title = element_text(size = 24, face = "bold"),
    
    # 2. Bigger Axis Titles ("# Years Post-Treatment" & "Mean % Change...")
    axis.title = element_text(size = 20),
    
    # 3. Bigger Axis Tick Labels (the numbers/categories on x & y axes)
    axis.text = element_text(size = 16),
    
    # 4. Bigger Legend Title & Text
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"), # Optional: enlarges legend symbol boxes
    
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


# native richness plot -----------------------------------------------
ggplot(df_summary_clean, aes(x = yrs_cat, y = total_richchange, color = lz_group, fill = lz_group, group = lz_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Standard Error ribbons
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
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
    size = 4,
    fontface = "bold",
    position = pd,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("< 5%" = "#B07AA1", "5-15%" = "#F28E2B", "> 15%" =  "#174D79"), 
    name = "% Littoral \nZone Treated"
  ) +
  scale_fill_manual(
    values = c("< 5%" = "#B07AA1", "5-15%" = "#F28E2B", "> 15%" =  "#174D79"), 
    name = "% Littoral \nZone Treated"
  ) +
  # Expand y-axis limits slightly to ensure top labels aren't clipped
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title = "Change in Total Native Species Richness \nafter EWM treatment",
    x = "# Years Post-Treatment",
    y = "Mean Change in Total Native Species Richness"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # 1. Bigger Title (add `size = ...`)
    plot.title = element_text(size = 24, face = "bold"),
    
    # 2. Bigger Axis Titles ("# Years Post-Treatment" & "Mean % Change...")
    axis.title = element_text(size = 20),
    
    # 3. Bigger Axis Tick Labels (the numbers/categories on x & y axes)
    axis.text = element_text(size = 16),
    
    # 4. Bigger Legend Title & Text
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"), # Optional: enlarges legend symbol boxes
    
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# native abundance plot -----------------------------------------------
ggplot(df_summary_clean, aes(x = yrs_cat, y = nAb_change, color = lz_group, fill = lz_group, group = lz_group)) +
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
  
  scale_color_brewer(palette = "Set1", direction = -1, name = "% LZ Treated") +
  scale_fill_brewer(palette = "Set1", direction = -1, name = "% LZ Treated") +
  # Expand y-axis limits slightly to ensure top labels aren't clipped
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title = "Change in overall native species frequency of occurrence after EWM treatment by % littoral zone treated",
    x = "# Years Post-Treatment",
    y = "Mean Change in Overall Native Species FOO (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggplot(df_summary_clean, aes(x = yrs_cat, y = relprcntchangeNAB, color = lz_group, fill = lz_group, group = lz_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Standard Error ribbons
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.18, color = NA) +
  
  # Mean trajectory lines and points
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
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
    size = 4,
    fontface = "bold",
    position = pd,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c("< 5%" = "#B07AA1", "5-15%" = "#F28E2B", "> 15%" =  "#174D79"), 
    name = "% Littoral \nZone Treated"
  ) +
  scale_fill_manual(
    values = c("< 5%" = "#B07AA1", "5-15%" = "#F28E2B", "> 15%" =  "#174D79"), 
    name = "% Littoral \nZone Treated"
  ) +
  # Expand y-axis limits slightly to ensure top labels aren't clipped
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title = "Change in Native Species Frequency of \nOccurrence after EWM treatment",
    x = "# Years Post-Treatment",
    y = "Mean Change in Overall Native Species FOO (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # 1. Bigger Title (add `size = ...`)
    plot.title = element_text(size = 24, face = "bold"),
    
    # 2. Bigger Axis Titles ("# Years Post-Treatment" & "Mean % Change...")
    axis.title = element_text(size = 20),
    
    # 3. Bigger Axis Tick Labels (the numbers/categories on x & y axes)
    axis.text = element_text(size = 16),
    
    # 4. Bigger Legend Title & Text
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"), # Optional: enlarges legend symbol boxes
    
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# uncomment to see color palette
#"#4E79A7" "#F28E2B" "#59A14F" "#E15759" "#B07AA1" "#76B7B2"


