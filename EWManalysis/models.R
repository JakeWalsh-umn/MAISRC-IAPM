# This script shows the models as of 9/21/2026

# ordered factor models

# EWM models ----
df_forModel <- pi.ewmlks %>% 
  filter(EWM_foo > 0 & preEWM_foo > 0) %>%
  mutate(EWMchangeFOO = EWM_foo - preEWM_foo, 
         nRich_change = avg.rich - pre.avg.n.rich ,
         nAb_change = (nat.ab/station_count) - pre.avg.n.ab)

# Additive Models
ewm_1 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    yrs.post.tx_fac + prcnt_LZ + cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = df_forModel,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # yrs.post.tx and prcnt_LZ 


ewm_brm_additive3 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # yrs.post.tx and prcnt_LZ
   #BEST MODEL

ewm_brm_s <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    s(yrs.post.tx, k = 4) + prcnt_LZ + cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  backend = "cmdstanr"
) # prcnt_LZ 

ewm_brm_additive4 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + pre.ewm.avg.depth + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # yrs.post.tx and prcnt_LZ

ewm_brm_interact <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    yrs.post.tx*prcnt_LZ*cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) 

ewm_2 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    prcnt_LZ*(yrs.post.tx + cytrt) + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = df_forModel,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) 

ewm_2 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    prcnt_LZ*(yrs.post.tx_fac + cytrt) + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
)

ewm_3 <- brm(
  myriophyllum_spicatum | trials(station_count) ~ 
    prcnt_LZ*yrs.post.tx_fac + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) 

library(brms)
library(loo)

# 1. Add LOO criteria to your fitted models
ewm_brm_additive <- add_criterion(ewm_brm_additive, "loo")
ewm_brm_additive2 <- add_criterion(ewm_brm_additive2, "loo")
ewm_brm_additive3 <- add_criterion(ewm_brm_additive3, "loo")
ewm_brm_additive4 <- add_criterion(ewm_brm_additive4, "loo")
ewm_brm_interact <- add_criterion(ewm_brm_interact, "loo")

# 2. Compare the models directly
loo_compare(ewm_brm_additive, ewm_brm_additive3, ewm_brm_interact)




# Native richness models ----

nat_rich_brm_1 <- brm(
  total.rich | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # yrs.post.tx and cytrt

nat_rich_brm_2 <- brm(
  total.rich | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + cytrt + preEWM + pre.avg.n.ab +(1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # yrs.post.tx, prcnt_LZ, cytrt, pre.avg.n.ab (preEWM is zeros across the stats!)

nat_rich_brm_3 <- brm(
  total.rich | trials(station_count) ~ 
    prcnt_LZ*(yrs.post.tx + cytrt) + preEWM + pre.avg.n.ab +(1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # pre.avg.n.ab

nat_rich_brm_4 <- brm(
  total.rich | trials(station_count) ~ 
    prcnt_LZ*(yrs.post.tx + cytrt) + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # prcnt_LZ and cytrt

nat_rich_brm_5 <- brm(
  total.rich | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + cytrt + pre.ewm.avg.depth +(1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # bad

nat_rich_interact <- brm(
  total.rich | trials(station_count) ~ 
    yrs.post.tx*prcnt_LZ*cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # prcnt_LZ, yrs.post.tx:prcnt_LZ, prcnt_LZ:cytrt, yrs.post.tx:prcnt_LZ:cytrt
  # BEST MODEL

# 1. Add LOO criteria to your fitted models
nat_rich_brm_1 <- add_criterion(nat_rich_brm_1 , "loo")
nat_rich_brm_2 <- add_criterion(nat_rich_brm_2 , "loo")
nat_rich_brm_3 <- add_criterion(nat_rich_brm_3, "loo")
nat_rich_brm_4 <- add_criterion(nat_rich_brm_4, "loo")
nat_rich_interact <- add_criterion(nat_rich_interact, "loo")

# 2. Compare the models directly
loo_compare(nat_rich_brm_1, nat_rich_brm_4, nat_rich_interact)

# native abundance models ----

nat_ab_brm_1 <- brm(
  nat.ab | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # yrs.post.tx and cytrt

nat_ab_brm_2 <- brm(
  nat.ab | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + cytrt + preEWM +(1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # bad

nat_ab_brm_3 <- brm(
  nat.ab | trials(station_count) ~ 
    yrs.post.tx + prcnt_LZ + cytrt + pre.avg.n.rich +(1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # bad

nat_ab_brm_4 <- brm(
  nat.ab | trials(station_count) ~ 
    prcnt_LZ*(yrs.post.tx + cytrt) + preEWM + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) # bad

nat_ab_interact <- brm(
  nat.ab | trials(station_count) ~ 
    yrs.post.tx*prcnt_LZ*cytrt + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) #prcnt_LZ, yrs.post.tx:prcnt_LZ, prcnt_LZ:cytrt, yrs.post.tx:prcnt_LZ:cytrt
  # BEST MODEL

nat_ab_interact <- brm(
  nat.ab | trials(station_count) ~ 
    prcnt_LZ*(yrs.post.tx + cytrt) + (1 | DOW_parent),
  family = beta_binomial(link = "logit"),
  data = pi.ewmlks,
  chains = 4, cores = 4, iter = 2000,
  backend = "cmdstanr"
) #prcnt_LZ, yrs.post.tx:prcnt_LZ, prcnt_LZ:cytrt, yrs.post.tx:prcnt_LZ:cytrt
# BEST MODEL

# 1. Add LOO criteria to your fitted models
nat_ab_interact <- add_criterion(nat_ab_interact, "loo")
nat_ab_brm_1 <- add_criterion(nat_ab_brm_1, "loo")

# 2. Compare the models directly
loo_compare(nat_ab_interact, nat_ab_brm_1)

