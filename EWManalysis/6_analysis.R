
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
  group_by(SURVEY_START, DOW_parent, relative_year, depth_ft) %>%
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

# filter foo surveys that have been treated
pi.trt = MN_foo_plus %>% filter(DOW_parent %in% MN_trts_acre$DOW_parent)
pi.trt$DOW_parent = as.character(pi.trt$DOW_parent)

# merge
df = left_join(pi.trt, MN_trts_acre, by = "DOW_parent", relationship = "many-to-many")
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
df$lag <- year(df$SURVEY_START) - df$treatment_year
df$YSI = df$treatment_year - df$year_confirm

# separate out pre and post treatment surveys
pre = df %>% filter(lag < 0)
post = df %>% filter(lag > 0)

post = post %>% group_by(DOW_parent, SURVEY_START) %>%
  arrange(DOW_parent, lag) %>%
  slice_head(n = 1) %>%
  ungroup() %>% filter(DOW_parent %in% pre$DOW_parent)

pre = pre %>% group_by(DOW_parent) %>%
  arrange(DOW_parent, lag) %>%
  slice_tail(n = 1) %>%
  ungroup() %>% filter(DOW_parent %in% post$DOW_parent) %>% 
  left_join(dplyr::select(post, DOW_parent), by = "DOW_parent") %>%
  mutate(lag = post$lag)


# extract just the native species

natives <- taxa_xwalk %>%
  filter(ORIGIN == "native") %>%
  pull(fieldNames)

pre.n = pre %>%
  dplyr::select(any_of(natives))

post.n = post %>%
  dplyr::select(any_of(natives))

# remove missing species
post.n = post.n[,-(which(colSums(pre.n[ ,1:ncol(pre.n)], na.rm = T) == 0))]
pre.n = pre.n[,-(which(colSums(pre.n[ ,1:ncol(pre.n)], na.rm = T) == 0))]


# get pre-post difference in EWM, native abundance, and native richness
ewm_result <- dplyr::select(post,`myriophyllum_spicatum`) - dplyr::select(pre,`myriophyllum_spicatum`)
pre$EWM.post = post$myriophyllum_spicatum
pre$trt.effect = ifelse(ewm_result$`myriophyllum_spicatum` < 0, 1, 0)

nab_result <- post.n - pre.n
pre$nab.change = rowSums(nab_result, na.rm = T)
pre$nab.effect = ifelse(rowSums(nab_result, na.rm = T) < 0, 1, 0)

com_result <- rowSums(post.n > 0, na.rm = T) - rowSums(pre.n > 0)
pre$com.change = com_result
pre$com.effect = ifelse(com_result < 0, 1, 0)


# pre-treatment metrics
pre$preEWM = pre$`myriophyllum_spicatum`
pre$pre.richness = rowSums(dplyr::select(pre.n,`bidens_beckii`:`zannichellia_palustris`) > 0, na.rm =T)
pre$pre.nab = rowSums(dplyr::select(pre.n,`bidens_beckii`:`zannichellia_palustris`), na.rm =T)

# note: fix this later 
# there are some lakes where the littoral zone appears to be bigger than the lake area itself
pre <- pre %>% 
  mutate(
    p.area = (treatment_size / lk_area))

# remove duplicates from the paired treatment/surveys

pre <- pre %>%
  distinct(DOW_parent, survey_id, lag, .keep_all = TRUE) 

pre <- pre %>% 
  mutate(cytrt = replace(cytrt, is.nan(cytrt), 0))

pre <- pre %>%
  mutate(cytrt = replace(cytrt, is.nan(cytrt), 0))

# modeling ----

## trt.effect models ----

# Come back to this:
# Filter for post surveys that are only 1 year after treatment

m1 = glm(formula = trt.effect ~ pre.richness + prcnt_LZ , family = binomial(link = "logit"), data = pre)         

m1 = glm(formula = trt.effect ~ pre.richness + p.area + cytrt + offset(preEWM), family = binomial(link = "logit"), data = pre)

m1 = glm(formula = trt.effect ~ pre.richness + prcnt_LZ + cytrt + lag + offset(preEWM), family = binomial(link = "logit"), data = pre)

m1 = glm(formula = trt.effect ~ pre.richness + YSI + pre.nab, family = binomial(link = "logit"), data = pre)

m1 = glm(formula = trt.effect ~ pre.richness + cytrt + prcnt_LZ + pre.nab, family = binomial(link = "logit"), data = pre)

m1 = glm(formula = trt.effect ~ pre.richness + cytrt + prcnt_LZ, family = binomial(link = "logit"), data = pre)

m1 = glm(formula = trt.effect ~ pre.richness + pre.nab + YSI + prcnt_LZ, family = binomial(link = "logit"), data = pre)

m1 = glm(formula = trt.effect ~ pre.richness + prcnt_LZ + cytrt + lag + YSI, family = binomial(link = "logit"), data = pre)

m1 = glm(formula = trt.effect ~ pre.richness + prcnt_LZ, family = binomial(link = "logit"), data = pre)

m1.1 = glm(formula = trt.effect ~ pre.richness, family = binomial(link = "logit"), data = pre) 

m1.2 = glm(formula = trt.effect ~ pre.richness*pre.nab, family = binomial(link = "logit"), data = pre) 

summary(m1.1)

summary(m1.2)

plot_models(m1.1, m1.2)

### trt.effect model viz ----

install.packages(c("ggplot2", "ggeffects"))
library(ggplot2)
library(ggeffects)

pred_m1.1 <- predict_response(m1.1, terms = "pre.richness")

# Render plot
plot(pred_m1.1) +
  labs(
    x = "Native richness before treatment",
    y = "Probability of EWM decline post-treatment",
    title = NULL
  ) +
  theme_minimal()

plot_model(m1.2, type = "pred", terms = c("pre.richness[0:22]", "pre.nab[0.3, 1, 2, 3]")) +
  labs(
    x = "native species richness before treatment",
    y = "Probability of EWM decline post-treatment",
    title = NULL,
    color = "pre.nab"
  ) +
  theme_minimal()


## nab.effect models ----

m2.1 = glm(formula = nab.effect ~ preEWM + prcnt_LZ, family = binomial(link = "logit"), data = pre) 

m2 = glm(formula = nab.effect ~ cytrt + prcnt_LZ, family = binomial(link = "logit"), data = pre)

m2.2 = glm(formula = nab.effect ~ preEWM + prcnt_LZ + cytrt, family = binomial(link = "logit"), data = pre) 

m2.3 = glm(formula = nab.effect ~ prcnt_LZ, family = binomial(link = "logit"), data = pre) 

m2 = glm(formula = nab.effect ~ prcnt_LZ + pre.richness, family = binomial(link = "logit"), data = pre) 

summary(m2.1)
summary(m2.2)
summary(m2.3)

### nab.effect model viz ----


#visualize all 3 models

plot_models(
  m2.1, m2.2, m2.3)


# visualize m2.1

plot_model(m2.1, type = "pred", terms = c("prcnt_LZ[0:1]", "preEWM[0.1, 0.3, 0.5]")) +
  labs(
    x = "% littoral zone treated",
    y = "Probability of reduced native abundance post-treatment",
    title = NULL,
    color = "preEWM"
  ) +
  theme_minimal()



## com.effect models ----

m3 = glm(formula = com.effect ~ cytrt + prcnt_LZ + preEWM + YSI + lag, family = binomial(link = "logit"), data = pre) 

m3.1 = glm(formula = com.effect ~ prcnt_LZ, family = binomial(link = "logit"), data = pre) 

summary(m3.1)

### com.effect model viz ----

pred_m3.1 <- predict_response(m3.1, terms = "prcnt_LZ")

# Render plot
plot(pred_m3.1) +
  labs(
    x = "% littoral zone treated",
    y = "probability of a negative change in native plant community",
    title = NULL
  ) +
  theme_minimal()

# rerunning models with lag = 1 ----

# rerunning models with only lag = 1 data
pre_lag1 <- pre %>% 
  filter(lag == 1) # n=71

m1.1.2 = glm(formula = trt.effect ~ pre.richness, family = binomial(link = "logit"), data = pre_lag1) 

m1.2.2 = glm(formula = trt.effect ~ pre.richness*pre.nab, family = binomial(link = "logit"), data = pre_lag1) 
# better than m1.2

summary(m1.1.2)

summary(m1.2.2)

summary(m1.3.2)

pred_m1.1 <- predict_response(m1.1.2, terms = "pre.richness")

# Render plot
plot(pred_m1.1) +
  labs(
    x = "Native richness before treatment",
    y = "Probability of EWM decline post-treatment",
    title = NULL
  ) +
  theme_minimal()

plot_model(m1.2.2, type = "pred", terms = c("pre.richness[0:22]", "pre.nab[0.3, 1, 2, 3]")) +
  labs(
    x = "native species richness before treatment",
    y = "Probability of EWM decline post-treatment",
    title = NULL,
    color = "pre.nab"
  ) +
  theme_minimal()


library(ggplot2)

# 1. Fit your model (if not already stored)

# 2. Create a grid of values across predictor ranges
grid <- expand.grid(
  pre.richness = seq(min(pre_lag1$pre.richness, na.rm = TRUE), 
                     max(pre_lag1$pre.richness, na.rm = TRUE), 
                     length.out = 100),
  pre.nab = seq(min(pre_lag1$pre.nab, na.rm = TRUE), 
                max(pre_lag1$pre.nab, na.rm = TRUE), 
                length.out = 100)
)

# 3. Predict probability P(trt.effect = 1)
grid$prob <- predict(m1.2.2, newdata = grid, type = "response")

# 4. Plot heatmap with contours
ggplot(grid, aes(x = pre.richness, y = pre.nab, fill = prob)) +
  geom_tile() +
  geom_contour(aes(z = prob), color = "black", alpha = 0.4, linetype = "dashed") +
  scale_fill_viridis_c(option = "YlGnBu", name = "P(trt.effect = 1)") +
  labs(
    title = "Predicted probability of EWM reduction post-treatment",
    x = "Pre-treatment native species richness",
    y = "Pre-treatment native species abundance"
  ) +
  theme_minimal(base_size = 14)

m2.1.2 = glm(formula = nab.effect ~ preEWM + prcnt_LZ, family = binomial(link = "logit"), data = pre_lag1)

m2.2.2 = glm(formula = nab.effect ~ preEWM + prcnt_LZ + cytrt, family = binomial(link = "logit"), data = pre_lag1) 

m2.3.2 = glm(formula = nab.effect ~ prcnt_LZ, family = binomial(link = "logit"), data = pre_lag1) 

summary(m2.1.2)
summary(m2.2.2)
summary(m2.3.2)

# visualize all 3 models
plot_models(
  m2.1.2, m2.2.2, m2.3.2)


# visualize m2.1

plot_model(m2.1.2, type = "pred", terms = c("prcnt_LZ[0:1]", "preEWM[0.1, 0.3, 0.5]")) +
  labs(
    x = "% littoral zone treated",
    y = "Probability of reduced native abundance post-treatment",
    title = NULL,
    color = "preEWM"
  ) +
  theme_minimal()


m3.1.2 = glm(formula = com.effect ~ prcnt_LZ, family = binomial(link = "logit"), data = pre_lag1) 

m3.2.2 = glm(formula = com.effect ~ prcnt_LZ + preEWM, family = binomial(link = "logit"), data = pre_lag1)

summary(m3.1.2)
summary(m3.2.2)

pred_m3.1.2 <- predict_response(m3.1.2, terms = "prcnt_LZ")

# Render plot
plot(pred_m3.1.2) +
  labs(
    x = "% littoral zone treated",
    y = "probability of a negative change in native plant community",
    title = NULL
  ) +
  theme_minimal()
