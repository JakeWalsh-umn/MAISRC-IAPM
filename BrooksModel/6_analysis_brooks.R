library(readr)
library(readxl)
library(robustbase)
library(ggeffects)
library(ggplot2)
library(ggdist)
library(ggpubr)
library(ggrepel)
library(ggforce)
library(RColorBrewer)
library(tidyverse)
library(MASS)
library(rstatix)
library(lubridate)
library(rasterdiv)
library(glmtoolbox)
library(emmeans)
library(multcomp)
library(MuMIn)
library(report)
library(qpcR)
library(viridis)
library(dplyr)


# IDs of treated lakes
trt.lakes = unique(MN_trts_acre$DOW_parent)

#################
# 1. AOR analysis
#################

####################################
# extract and summarize species data

spp_dat = MN_foo %>%
  dplyr::select(DOW_parent, all_of(tax_cols))

spp_dat = spp_dat %>% 
  group_by(DOW_parent) %>%
  summarise(across(everything(), mean)) 


spp_dat[spp_dat == 0] = NA
spp_dat = as.data.frame(spp_dat)


# NOTE: add the step here that filters out non-eurasian niche lakes and permits
hydrography <- read.csv('data/Lake_Basin_Morphology.csv')

hydrography <- hydrography %>% 
  filter(island == "N") %>%
  mutate(dowlknum = as.character(dowlknum))

EWM_lakes <- iw %>% 
  filter(species == "Eurasian watermilfoil") 

filtered_hydro <- hydrography %>%
  filter(dowlknum %in% EWM_lakes$DOW_parent)

summary(filtered_hydro$lk_area)

summary(filtered_hydro$lk_avdepth)

untreated_lks <- EWM_lakes %>%
  mutate(DOW6 = str_sub(as.character(DOW_parent), 1, 6)) %>%
  anti_join(
    iapm_permits %>% mutate(DOW6 = str_sub(as.character(water_resource_numbers), 1, 6)),
    by = "DOW6"
  ) %>%
  dplyr::select(-DOW6)

untreated_lks <- untreated_lks %>% 
  inner_join(
    hydrography %>%
      dplyr::select(dowlknum, lk_avdepth, lk_area),
    by = c("DOW_parent" = "dowlknum")
  )

EWM_nichelks <- hydrography %>% 
  filter(lk_area < 2000 & lk_avdepth < 20)

##################
# invasion status

inv.status = vector()


for(j in 1:length(spp_dat[,1])){
  
  if((!spp_dat$DOW_parent[j] %in% EWM_nichelks & is.na(spp_dat$'myriophyllum_spicatum'[j])))
    inv.status[j] = "non-candidate"
  else
    if((spp_dat$DOW_parent[j] %in% EWM_nichelks & is.na(spp_dat$`myriophyllum_spicatum`[j]))) 
    inv.status[j] = "uninvaded" 
  else 
    if((spp_dat$DOW_parent[j] %in% trt.lakes & !is.na(spp_dat$`myriophyllum_spicatum`[j]))) 
      inv.status[j] = "invaded_treated" 
    else 
      if((!spp_dat$DOW_parent[j] %in% iapm_permits & !is.na(spp_dat$`myriophyllum_spicatum`[j])))  #NOTE: fix this to refer to permits
        inv.status[j] = "invaded_untreated" 
    else
      NA
}

 EWM_nichelks <- EWM_nichelks %>%
  mutate(DOW_parent = str_replace(dowlknum, "\\d{2}$", "00"))

iapm_permits <- iapm_permits %>% 
  mutate(DOW_parent = str_replace(water_resource_numbers, "\\d{2}$", "00"))

niche_vector   <- EWM_nichelks$DOW_parent
permit_vector  <- iapm_permits$DOW_parent  # update column name if different

spp_dat <- spp_dat %>%
  mutate(
    # Optional: ensure DOW_parent is character across all datasets to avoid type mismatches
    DOW_parent = as.character(DOW_parent),
    
    inv_status = case_when(
      !DOW_parent %in% niche_vector  & is.na(myriophyllum_spicatum)  ~ "non-candidate",
      DOW_parent %in% niche_vector  & is.na(myriophyllum_spicatum)  ~ "uninvaded",
      DOW_parent %in% trt.lakes    & !is.na(myriophyllum_spicatum) ~ "invaded_treated",
      !DOW_parent %in% permit_vector & !is.na(myriophyllum_spicatum) ~ "invaded_untreated",
      TRUE                                                           ~ NA_character_
    )
  )

# remove DOW for easier processing across all columns at once.
spp_dat = spp_dat %>%
  dplyr::select(-DOW_parent)

spp_dat = spp_dat %>%
  filter(!is.na(inv_status) & inv_status != "non-candidate")

inv.status <- spp_dat$inv_status

spp_dat <- spp_dat%>%
  dplyr::select(-c(inv_status))
  


####################
# evenness by status

natives <- taxa_xwalk %>%
  filter(ORIGIN == "native") %>%
  pull(fieldNames)

n.com = spp_dat %>%
  dplyr::select(any_of(natives))

n.com[is.na(n.com)] = 0

pielous.evenness <- function(species_counts) {
  # 2. Calculate Species Richness (S) - number of species with counts > 0
  S <- sum(species_counts > 0)
  # 3. Calculate Relative Abundances (p_i)
  p_i <- species_counts[species_counts > 0] / sum(species_counts)
  # 4. Calculate Shannon Diversity (H')
  H_prime <- -sum(p_i * log(p_i))
  # 5. Calculate Pielou's Evenness (J')
  J_prime <- H_prime / log(S)
  
  print(J_prime)
  
}

p.e = apply(n.com, 1, function(x) pielous.evenness(x))

df = data.frame(inv.status, p.e = as.numeric(p.e))

df %>% group_by(inv.status) %>%
  summarise(mean = mean(p.e, na.rm = T), sd = sd(p.e, na.rm = T), n = sum(!is.na(p.e)))


##########################
# AOR regression by status

spp_dat = split(spp_dat, f = inv.status)

aor_dat = data.frame(occupancy = lapply(spp_dat, FUN = function(x) apply(x, 2, FUN = function(y) sum(!is.na(y)))),
                     abundance = lapply(spp_dat, FUN = function(x) apply(x, 2, FUN = function(y) mean(y, na.rm = T))))



aor_dat$spp = row.names(aor_dat)

abs = aor_dat %>% dplyr::select(spp, abundance.invaded_treated, abundance.invaded_untreated, abundance.uninvaded) %>% 
  gather(key = "invaded", value = "abundance", -spp)

occ = aor_dat %>% dplyr::select(spp, occupancy.invaded_treated, occupancy.invaded_untreated, occupancy.uninvaded) %>% 
  gather(key = "invaded", value = "occupancy", -spp) %>% dplyr::select(occupancy)

abs$occupancy = occ$occupancy

abs[is.na(abs)] = 0


abs$n.sites = NA

for(j in 1:length(abs$spp)){
  
  if((abs$invaded[j] == "abundance.uninvaded")) 
    abs$n.sites[j] = length(spp_dat$uninvaded$aquatic_moss)
  else 
    if((abs$invaded[j] == "abundance.invaded_treated"))
      abs$n.sites[j] = length(spp_dat$invaded_treated$aquatic_moss)
    else
      if((abs$invaded[j] == "abundance.invaded_untreated"))
        abs$n.sites[j] = length(spp_dat$invaded_untreated$aquatic_moss)
      
}

taxa_xwalk_clean <- taxa_xwalk %>%
  distinct(fieldNames, .keep_all = TRUE)

abs = merge(abs, taxa_xwalk_clean, by.x = "spp", by.y = "fieldNames")

abs = abs %>% 
  filter(ORIGIN == "native")

abs$invaded = as.factor(abs$invaded)
levels(abs$invaded) = c("Treated", "Untreated", "Uninvaded")
abs$invaded = fct_rev(abs$invaded)
abs = abs[abs$occupancy > 0, ]

m.all <- glmrob(
  cbind(occupancy, n.sites-occupancy) ~ abundance * invaded, 
  family = binomial, 
  data = abs, 
  maxit = 1000)

summary(m.all)

#abs$abundance_std <- scale(abs$abundance)

#library(glmmTMB)

#m_fixed <- glmmTMB(
#  cbind(occupancy, n.sites - occupancy) ~ abundance_std * invaded + (1 | spp),
#  family = betabinomial,
#  data = abs
#)

summary(m_fixed)


library(ggeffects)
install.packages("ggthemes")
library(ggthemes)
plot(ggpredict(m.all, terms = c("abundance [all]", "invaded"))) + 
 # theme_Publication() +
  theme_classic() +
  ggtitle("") + 
  xlab("Abundance") + xlim(0, 0.4) +
  ylab("Occupancy") +
  geom_point(data = abs, aes(x = abundance, y = occupancy / n.sites, fill = invaded, col = invaded)) +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2")+
  theme(legend.title = element_blank())



######################
# 2. pre-post analysis
######################


########################
# filter and format data
MN_foo = MN_foo %>% 
  filter(year(SURVEY_START) >= 2012)

MN_trts_acre = MN_trts_acre %>%
  filter(treatment_year >= 2017)

# filter foo surveys that have been treated
pi.trt = MN_foo %>% filter(DOW_parent %in% MN_trts_acre$DOW_parent)
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
pre.n = pre %>%
  dplyr::select(any_of(natives))

post.n = post %>%
  dplyr::select(any_of(natives))

# remove missing species
post.n = post.n[,-(which(colSums(pre.n[ ,1:ncol(pre.n)], na.rm = T) == 0))]
pre.n = pre.n[,-(which(colSums(pre.n[ ,1:ncol(pre.n)], na.rm = T) == 0))]


# get pre-post difference in EWM, native abundance, and native richness
ewm_result <- dplyr::select(post,`myriophyllum_spicatum`) - dplyr::select(pre,`myriophyllum_spicatum`)
pre$trt.effect = ifelse(ewm_result$`myriophyllum_spicatum` < 0, 1, 0)

nab_result <- post.n - pre.n
pre$nab.effect = ifelse(rowSums(nab_result, na.rm = T) < 0, 1, 0)

com_result <- rowSums(post.n > 0, na.rm = T) - rowSums(pre.n > 0)
pre$com.effect = ifelse(com_result < 0, 1, 0)


# pre-treatment metrics
pre$preEWM = pre$`myriophyllum_spicatum`
pre$pre.richness = rowSums(dplyr::select(pre.n,`bidens_beckii`:`zannichellia_palustris`) > 0, na.rm =T)
pre$pre.nab = rowSums(dplyr::select(pre.n,`bidens_beckii`:`zannichellia_palustris`), na.rm =T)


############################################
# ecological predictors of treatment success

# going to try running this without all of the same fields ----

#lchem_wide = pivot_wider(lchem, names_from = DNR_PARAMETER_DESCRIPTION, values_from = RESULT_AMT, values_fn = function(x) mean(x, na.rm = T))

#df = left_join(pre, lsize, by = "WBIC")
#df = left_join(df, lchem_wide, by = "WBIC")
#df$chem = as.factor(df$`Standardized Name`)
pre = pre %>% drop_na(active_ing)
pre <- pre %>% 
  
  
  mutate(
    active_ing = fct_collapse(active_ing,
                              "24D"                  = c("2,4-D, Triclopyr", "2,4-D"),
                              "DiquatEndothall"      = c("Endothall", "Diquat Dibromide", "Endothall, Diquat Dibromide"),
                              "FlorpyrauxifenBenzyl" = "Florpyrauxifen-benzyl"
    ),
    active_ing = fct_drop(active_ing)
  ) %>%
  
  # 2. Filter for only the recoded target ingredients AFTER recoding
  filter(active_ing == "24D" | active_ing ==  "FlorpyrauxifenBenzyl" | active_ing ==  "DiquatEndothall") 

#####################
# logisitc regression

# top models for each response variable
# NOTE: swapping p.area for prcnt_LZ. We can use p.area (% permitted, treated)
      # but the swap would mean we lose a lot of data 
# NOTE: remove the active_ing from the models. 
# NOTE: try to return the models to the form that Brooks had these in once I fix my station numbers. 

# here are two different version 

# 8.11 notes from Jake
  # don't worry about YSI because we might not have the same data to inform this.
  # consider looking at specific species and how they are affected by certain variables. 
  # consider looking at APM and working it into model
  # could consider water quality parameters, but since George didn't use them, it may not be a priority. 
  # calculate the proportion of lake treated and look at LZ and lake percentages. 


#NOTE: avoid highly correlated variables in the same model, but in the case of lag here, 
# the lag:YSI is the relationship between the variables which is different than the variables themselves. 
mtop = glm(formula = trt.effect ~ prcnt_LZ + active_ing + lag + YSI +  
             lag:YSI + offset(preEWM), 
           family = binomial(link = "logit"), data = pre)


mtop2 = glm(
  formula = nab.effect ~ preEWM * prcnt_LZ + offset(pre.nab), 
  family = binomial(link = "logit"), 
  data = pre
)

mtop3 = glm(formula = com.effect ~ lag * prcnt_LZ + active_ing, family = binomial(link = "logit"), data = pre)

# NOTE: add another model that looks at specific species and how they're affected by treatment, etc.

# At this point look at the testmodels.R script where I play around with different models

# playing around with different models that remove active_ing ----
mtop = glm(formula = trt.effect ~ lag + pre.richness + prcnt_LZ +
             offset(preEWM) + YSI, 
           family = binomial(link = "logit"), data = pre)

summary(mtop)

mtop2 = glm(formula = nab.effect ~ lag + pre.richness + prcnt_LZ +
             + YSI, 
            family = binomial(link = "logit"), data = pre)

summary(mtop2)

#####################
# model validation

install.packages("Rdocumentation")
library(Rdocumentation)

report(mtop)
report(mtop2)
report(mtop3)

hltest(mtop)
hltest(mtop2)
hltest(mtop3)

library(emmeans)
em_means <- emmeans(mtop, "active_ing")
pairs(em_means, adjust = "holm")

em_means <- emmeans(mtop3, "active_ing")
pairs(em_means, adjust = "holm")

  
# plotting
# model predictions

# relationship between years since treatment, time since invasion, and EWM decline
plot(ggpredict(mtop, terms = c("lag [1:3]", "YSI [30, 1]", "prcnt_LZ [0.0001, 0.3]"), ci_level = 0.95)) + 
  theme_classic() + 
  ggtitle("") + 
  xlab("Years since treatment") +
  ylab("Probability that EWM abundance\n is lower post-treatment") + ylim(0, 1) + 
  theme(legend.position = "top") +
  labs(color = "Years since invasion") +
  guides(colour = guide_legend(title.position = "top")) + 
  scale_fill_viridis_d(option = "plasma") + 
  scale_color_viridis_d(option = "plasma")

plot(ggpredict(mtop, terms = c("active_ing"), ci_level = 0.95)) + 
  theme_classic() + 
  ggtitle("") + 
  xlab("Chemical used") +
  ylab("Probability that EWM abundance\n is lower post-treatment") + ylim(0, 1) 

plot(ggpredict(mtop3, terms = c("active_ing"), ci_level = 0.95)) + 
  theme_classic() + 
  ggtitle("") + 
  xlab("Chemical used") +
  ylab("Probability that native richness\n is lower post-treatment") + ylim(0, 1) 



########################
# YSI ~ littoral zone treated contours
# I think the following plots show the probability of EWM decline ----
# predictions for contours
p1 = ggpredict(mtop, terms = c("YSI [0:50]", "prcnt_LZ [0.1:1 by=0.05]", "active_ing [FlorpyrauxifenBenzyl]"), ci_level = NA)
p1$group = as.numeric(as.character(p1$group))
p1 = p1 %>% filter(group >= 0.1 & group <= 1)

p1d = ggpredict(mtop, terms = c("YSI [0:50]", "prcnt_LZ [0.1:1 by=0.05]", "active_ing [24D]"), ci_level = NA)
p1d$group = as.numeric(as.character(p1d$group))
p1d = p1d %>% filter(group >= 0.1 & group <= 1)

p1e = ggpredict(mtop, terms = c("YSI [0:50]", "prcnt_LZ [0.1:1 by=0.05]", "active_ing [DiquatEndothall]"), ci_level = NA)
p1e$group = as.numeric(as.character(p1e$group))
p1e = p1e %>% filter(group >= 0.1 & group <= 1)


p2 = ggpredict(mtop2, terms = c("prcnt_LZ [0.1:1 by=0.05]"), ci_level = NA)
p2 = p2 %>% filter(x >= 0.1 & x <= 1)
p2 = data.frame(t.area = rep(p2$x, length(unique(p1$x))), p = rep(p2$predicted, length(unique(p1$x))))


# EWM plot - FlorpyrauxifenBenzyl
pred.dat = data.frame(cb = p1$predicted, t.area = p1$group, YSI = p1$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)


filled.contour(unique(p1$group), unique(p1$x), pred.dat, zlim = c(0,1), plot.title = title(main = "", xlab = "Proportion of littoral zone treated", ylab = "Years since invasion"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, drawlabels = T)
  
}
)

# EWM plot - 24D
pred.dat = data.frame(cb = p1d$predicted, t.area = p1d$group, YSI = p1d$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)

filled.contour(unique(p1d$group), unique(p1d$x), pred.dat, zlim = c(0,1), plot.title = title(main = "", xlab = "Proportion of littoral zone treated", ylab = "Years since invasion"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, drawlabels = T)
  
}
)

# EWM plot - Diquot/Endothall
pred.dat = data.frame(cb = p1e$predicted, t.area = p1e$group, YSI = p1e$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)


filled.contour(unique(p1e$group), unique(p1e$x), pred.dat, zlim = c(0,1), plot.title = title(main = "", xlab = "Proportion of littoral zone treated", ylab = "Years since invasion"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1e$group), unique(p1e$x), pred.dat, add = TRUE, lwd = 2, drawlabels = T)
  
}
)



# native plot
pred.dat = data.frame(cb = p2$p, t.area = p1$group, YSI = p1$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)


filled.contour(unique(p1$group), unique(p1$x), pred.dat, zlim = c(0,1), plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Years since invasion"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, drawlabels = T)
  
}
)


##########
#SOS plots

# SOS - 24D  
pred.dat = data.frame(cb = p1d$predicted - p2$p, t.area = p1d$group, YSI = p1d$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)

pred.dat2 = data.frame(cb = p2$p, t.area = p1$group, YSI = p1$x)
pred.dat2 = pred.dat2 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat2  =as.matrix(pred.dat2)

pred.dat3 = data.frame(cb = p1d$predicted, t.area = p1d$group, YSI = p1d$x)
pred.dat3 = pred.dat3 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat3  =as.matrix(pred.dat3)

filled.contour(unique(p1$group), unique(p1$x), pred.dat, zlim = c(-0.6,0.6), color.palette = magma, plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Years since invasion"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, levels = 0, drawlabels = F)
  points(pre$prcnt_LZ, pre$YSI, col = "gray", pch = 16, cex = 0.5)
  contour(unique(p1$group), unique(p1$x), pred.dat2, add = TRUE, lwd = 2, levels = 0.25, drawlabels = F, lty = "dashed")
  contour(unique(p1d$group), unique(p1d$x), pred.dat3, add = TRUE, lwd = 2, levels = 0.5, drawlabels = F, lty = "dashed")
  
}
)


# SOS - FlorpyrauxifenBenzyl
pred.dat = data.frame(cb = p1$predicted - p2$p, t.area = p1$group, YSI = p1$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)

pred.dat2 = data.frame(cb = p2$p, t.area = p1$group, YSI = p1$x)
pred.dat2 = pred.dat2 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat2  =as.matrix(pred.dat2)

pred.dat3 = data.frame(cb = p1$predicted, t.area = p1$group, YSI = p1$x)
pred.dat3 = pred.dat3 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat3  =as.matrix(pred.dat3)


filled.contour(unique(p1$group), unique(p1$x), pred.dat, zlim = c(-0.6,0.6), color.palette = magma, plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Years since invasion"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, levels = 0, drawlabels = F)
  points(pre$prcnt_LZ, pre$YSI, col = "gray", pch = 16, cex = 0.5)
  contour(unique(p1$group), unique(p1$x), pred.dat2, add = TRUE, lwd = 2, levels = 0.25, drawlabels = F, lty = "dashed")
  contour(unique(p1d$group), unique(p1d$x), pred.dat3, add = TRUE, lwd = 2, levels = 0.5, drawlabels = F, lty = "dashed")
  
}
)


# SOS - Diquot/Endothall
pred.dat = data.frame(cb = p1e$predicted - p2$p, t.area = p1e$group, YSI = p1e$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)

pred.dat2 = data.frame(cb = p2$p, t.area = p1e$group, YSI = p1e$x)
pred.dat2 = pred.dat2 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat2  =as.matrix(pred.dat2)

pred.dat3 = data.frame(cb = p1e$predicted, t.area = p1e$group, YSI = p1e$x)
pred.dat3 = pred.dat3 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat3  =as.matrix(pred.dat3)


filled.contour(unique(p1$group), unique(p1$x), pred.dat, zlim = c(-0.6,0.6), color.palette = magma, plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Years since invasion"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, levels = 0, drawlabels = F)
  points(pre$prcnt_LZ, pre$YSI, col = "gray", pch = 16, cex = 0.5)
  contour(unique(p1$group), unique(p1$x), pred.dat2, add = TRUE, lwd = 2, levels = 0.25, drawlabels = F, lty = "dashed")
  contour(unique(p1d$group), unique(p1d$x), pred.dat3, add = TRUE, lwd = 2, levels = 0.3, drawlabels = F, lty = "dashed")
  
}
)



#################
# preEWM ~ t.area

# predictions for contours
p1 = ggpredict(mtop, terms = c("preEWM [all]", "prcnt_LZ [0.1:1 by=0.05]", "active_ing [FlorpyrauxifenBenzyl]"), ci_level = NA)
p1$group = as.numeric(as.character(p1$group))
p1 = p1 %>% filter(group >= 0.1 & group <= 1)

p1d = ggpredict(mtop, terms = c("preEWM [all]", "prcnt_LZ [0.1:1 by=0.05]", "active_ing [24D]"), ci_level = NA)
p1d$group = as.numeric(as.character(p1d$group))
p1d = p1d %>% filter(group >= 0.1 & group <= 1)

p1e = ggpredict(mtop, terms = c("preEWM [all]", "prcnt_LZ [0.1:1 by=0.05]", "active_ing [DiquatEndothall]"), ci_level = NA)
p1e$group = as.numeric(as.character(p1e$group))
p1e = p1e %>% filter(group >= 0.1 & group <= 1)

p2 = ggpredict(mtop2, terms = c("preEWM [all]", "prcnt_LZ [0.1:1 by=0.05]"), ci_level = NA)
p2$group = as.numeric(as.character(p2$group))
p2 = p2 %>% filter(group >= 0.1 & group <= 1)



# EWM plot - 24D
pred.dat = data.frame(cb = p1$predicted, t.area = p1$group, EWM = p1$x)
pred.dat = pred.dat %>% spread(EWM, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)

filled.contour(unique(p1$group), unique(p1$x), pred.dat, zlim = c(0,1), plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Invader abundance"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, drawlabels = T)
  
}
)

# EWM plot - FlorpyrauxifenBenzyl
pred.dat = data.frame(cb = p1d$predicted, t.area = p1d$group, YSI = p1d$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)

filled.contour(unique(p1d$group), unique(p1d$x), pred.dat, zlim = c(0,1), plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Invader abundance"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, drawlabels = T)
  
}
)


# EWM plot - Diquot/Endothall
pred.dat = data.frame(cb = p1e$predicted, t.area = p1e$group, YSI = p1e$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)


filled.contour(unique(p1e$group), unique(p1e$x), pred.dat, zlim = c(0,1), plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Invader abundance"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1e$group), unique(p1e$x), pred.dat, add = TRUE, lwd = 2, drawlabels = T)
  
}
)



# native plot
pred.dat = data.frame(cb = p2$p, t.area = p1$group, YSI = p1$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)


filled.contour(unique(p1$group), unique(p1$x), pred.dat, zlim = c(0,1), plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Invader abundance"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, drawlabels = T)
  
}
)


## SOS - 24D
pred.dat = data.frame(cb = p1d$predicted - p2$p, t.area = p1d$group, YSI = p1d$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)

pred.dat2 = data.frame(cb = p2$p, t.area = p1$group, YSI = p1$x)
pred.dat2 = pred.dat2 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat2  =as.matrix(pred.dat2)

pred.dat3 = data.frame(cb = p1d$predicted, t.area = p1d$group, YSI = p1d$x)
pred.dat3 = pred.dat3 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat3  =as.matrix(pred.dat3)


filled.contour(unique(p1$group), unique(p1$x), pred.dat, #ylim = c(0, 0.45), 
               zlim = c(-0.6,0.6), color.palette = magma, plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Invader abundance"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, levels = 0, drawlabels = F)
  points(pre$prcnt_LZ, pre$preEWM, col = "gray", pch = 16, cex = 0.5)
  contour(unique(p1$group), unique(p1$x), pred.dat2, add = TRUE, lwd = 2, levels = 0.25, drawlabels = F, lty = "dashed")
  contour(unique(p1d$group), unique(p1d$x), pred.dat3, add = TRUE, lwd = 2, levels = 0.5, drawlabels = F, lty = "dashed")
  
}
)


# SOS - FlorpyrauxifenBenzyl
pred.dat = data.frame(cb = p1$predicted - p2$p, t.area = p1$group, YSI = p1$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)

pred.dat2 = data.frame(cb = p2$p, t.area = p1$group, YSI = p1$x)
pred.dat2 = pred.dat2 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat2  =as.matrix(pred.dat2)

pred.dat3 = data.frame(cb = p1$predicted, t.area = p1$group, YSI = p1$x)
pred.dat3 = pred.dat3 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat3  =as.matrix(pred.dat3)


filled.contour(unique(p1$group), unique(p1$x), pred.dat, #ylim = c(0, 0.45), 
               zlim = c(-0.6,0.6), color.palette = magma, plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Invader abundance"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, levels = 0, drawlabels = F)
  points(pre$prcnt_LZ, pre$YSI, col = "gray", pch = 16, cex = 0.5)
  contour(unique(p1$group), unique(p1$x), pred.dat2, add = TRUE, lwd = 2, levels = 0.25, drawlabels = F, lty = "dashed")
  contour(unique(p1d$group), unique(p1d$x), pred.dat3, add = TRUE, lwd = 2, levels = 0.5, drawlabels = F, lty = "dashed")
  
}
)


# SOS - Diquot/Endothall
pred.dat = data.frame(cb = p1e$predicted - p2$p, t.area = p1e$group, YSI = p1e$x)
pred.dat = pred.dat %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat  =as.matrix(pred.dat)

pred.dat2 = data.frame(cb = p2$p, t.area = p1e$group, YSI = p1e$x)
pred.dat2 = pred.dat2 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat2  =as.matrix(pred.dat2)

pred.dat3 = data.frame(cb = p1e$predicted, t.area = p1e$group, YSI = p1e$x)
pred.dat3 = pred.dat3 %>% spread(YSI, cb) %>% dplyr::select(-t.area)
pred.dat3  =as.matrix(pred.dat3)


filled.contour(unique(p1$group), unique(p1$x), pred.dat, #ylim = c(0, 0.45), 
               zlim = c(-0.6,0.6), color.palette = magma, plot.title = title(main = "", xlab = "Proportion of lake treated", ylab = "Invader abundance"), plot.axes = {
  axis(1)
  axis(2)
  contour(unique(p1$group), unique(p1$x), pred.dat, add = TRUE, lwd = 2, levels = 0, drawlabels = F)
  points(pre$prcnt_LZ, pre$YSI, col = "gray", pch = 16, cex = 0.5)
  contour(unique(p1$group), unique(p1$x), pred.dat2, add = TRUE, lwd = 2, levels = 0.25, drawlabels = F, lty = "dashed")
  contour(unique(p1d$group), unique(p1d$x), pred.dat3, add = TRUE, lwd = 2, levels = 0.3, drawlabels = F, lty = "dashed")
  
}
)

# NOTE: move this to the proper part of the code: 

#pre <- pre %>%
#  left_join(
#    post %>% select(survey_id, post.EWM = myriophyllum_spicatum),
#    by = "survey_id"  # Replace 'lake_id' with your matching key column
#  )

#pre$cytrt[is.nan(pre$cytrt)] <- 0
