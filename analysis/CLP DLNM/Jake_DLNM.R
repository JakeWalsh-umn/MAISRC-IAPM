# Dynamic lag nonlinear modeling
# Still not annotated but slightly less buggy and slightly more organized


# Packages ----
library(tidyverse)
library(arrow)
library(readxl)
library(lubridate)
library(sf)
library(tigris)
library(tidybayes)

library(mgcv)
library(gratia)
library(dlnm)
library(brms)

library(ggplot2)
library(viridis)
library(scico)
library(gridExtra)
library(cowplot)
library(ggdist)
library(patchwork)
library(ggrepel)
library(plotly)
library(htmlwidgets)

library(dtwclust)
library(clusterCrit)

# Data set up ----

## IAPM Permits ----
permits <- permits <- read.csv("JakeOutput/apm_iapm_permit_detail.csv",
                               na.strings=c("NA", "NULL"))

permits <- permits %>%
  dplyr::select(permit_number, permit_type, 
                water_resource_names, water_resource_numbers,
                treatment_method, species,
                permit_approved_acres_1997:permit_approved_acres_2025)

permits[, 7:35][is.na(permits[, 7:35])] <- 0

iapm_permits_long <- permits

colnames(iapm_permits_long)[7:35] <- 1997:2025

iapm_permits_long <- iapm_permits_long %>%
  filter(permit_type=="Invasive Aquatic Plant Management" & 
           treatment_method%in%c("Pesticide Control", "Mechanical Control")) %>%
  pivot_longer(cols=`1997`:`2025`, 
               names_to="year",
               values_to="approved_acres") %>%
  mutate(year=as.numeric(year))

iapm_permits_long2 <- permits

colnames(iapm_permits_long2)[7:35] <- 1997:2025

iapm_permits_long2 <- iapm_permits_long2 %>%
  filter(permit_type=="Invasive Aquatic Plant Management") %>%
  pivot_longer(cols=`1997`:`2025`, 
               names_to="year",
               values_to="approved_acres") %>%
  mutate(year=as.numeric(year))


## Working with "db_sub_pa" from "1_Loading and preparing PI data.R" ----
dat_dlnm <- db_sub_pa %>%
  filter(!is.na(depth_ft) & depth_ft <= 15) %>%
  mutate(SURVEY_START=as.Date(SURVEY_START),
         year=year(SURVEY_START),
         month=month(SURVEY_START),
         doy=yday(SURVEY_START),
         taxa_rich=rowSums(across(ceratophyllum_demersum:vallisneria_sp)),
         vegetated=ifelse(taxa_rich==0, 0, 1),
         EWM=apply(across(c(myriophyllum_spicatum, myriophyllum_spicatum_x_sibiricum)), MARGIN=1, FUN=max),
         native_rich=rowSums(across(c(ceratophyllum_demersum, 
                                      elodea_canadensis:eleocharis_acicularis,
                                      myriophyllum_sibiricum:potamogeton_robbinsii, 
                                      nitella_sp:potamogeton_strictifolius,
                                      fontinalis_sp:drepanocladus_sp,
                                      nitella_furcata:vallisneria_sp))),
         native_vegetated=ifelse(taxa_rich==0, 0, 1)) %>%
  filter(year>=2017) %>%
  group_by(DOW, SURVEY_START, year, month, doy) %>%
  mutate(n_points=n()) %>%
  summarize_at(vars(ceratophyllum_demersum:vallisneria_sp, taxa_rich:n_points), mean) %>%
  ungroup() %>%
  mutate(survey_taxa_rich=rowSums(across(ceratophyllum_demersum:vallisneria_sp)>0),
         survey_native_rich=rowSums(across(c(ceratophyllum_demersum, 
                                             elodea_canadensis:eleocharis_acicularis,
                                             myriophyllum_sibiricum:potamogeton_robbinsii, 
                                             nitella_sp:potamogeton_strictifolius,
                                             fontinalis_sp:drepanocladus_sp,
                                             nitella_furcata:vallisneria_sp))>0))

dat_dlnm <- left_join(dat_dlnm,
                      lake_basin_morphology_df %>%
                        filter(ISLAND=="N" & !is.na(LK_LZACR_1)) %>%
                        dplyr::select(DOWLKNUM,
                                      LK_LZACR_1),
                      by=join_by(DOW==DOWLKNUM))

dat_dlnm <- left_join(dat_dlnm,
                      apm_permits_long %>%
                        group_by(water_resource_numbers, permit_number, year) %>%
                        summarize(APM_approved_acres=mean(approved_acres)) %>% 
                        ungroup() %>%
                        group_by(water_resource_numbers, year) %>%
                        summarize(APM_approved_acres=sum(APM_approved_acres)) %>%
                        ungroup(),
                      by=join_by(DOW==water_resource_numbers, year==year))

dat_dlnm$APM_approved_acres[is.na(dat_dlnm$APM_approved_acres)] <- 0

dat_dlnm <- dat_dlnm %>%
  mutate(DOW=factor(DOW),
         fyear=factor(year)) %>% 
  mutate(APM_prop=APM_approved_acres/LK_LZACR_1)

dat_dlnm_ewm_lakes <- as.character((dat_dlnm %>%
                                      group_by(DOW) %>%
                                      summarize(myriophyllum_spicatum=max(myriophyllum_spicatum)) %>%
                                      filter(myriophyllum_spicatum > 0))$DOW)

dat_dlnm_clp_lakes <- as.character((dat_dlnm %>%
                                      group_by(DOW) %>%
                                      summarize(potamogeton_crispus=max(potamogeton_crispus)) %>%
                                      filter(potamogeton_crispus > 0))$DOW)

# Treatment history: Annual total (APM + IAPM) treatment area
# as a fraction of littoral zone (<=15') area in a lagged matrix for DLNM.
# rows = surveys, columns = lags
# lag0 = treatment in year of survey 
# (includes surveys before and after treatment in the same year, for now)
# lag1 = treatment year before survey
# lag2 = treatment 2 years before survey ...
# lag8 = treatment 8 years before survey ...
# max possible lags is 8 (a 2017 treatment for a 2025 survey?), 
# so that's 9 columns
treatment_history_matrix <- matrix(nrow=dim(dat_dlnm)[1], ncol=9) %>%
  as.data.frame()
colnames(treatment_history_matrix) <- paste0("lag", 0:8)

target_string <- rep(NA, dim(dat_dlnm)[1])

# populate the matrix
for(i in 1:dim(dat_dlnm)[1]){
  
  if(dat_dlnm$DOW[i]%in%iapm_permits_long2$water_resource_numbers){
    
    permits_i <- iapm_permits_long2 %>%
      filter(water_resource_numbers==dat_dlnm$DOW[i])
    
    target_string[i] <- paste(unique(permits_i$species), collapse=",")
    
    year_i <- dat_dlnm$year[i]
    
    treatment_history_matrix[i, 1] <- ifelse(year_i%in%permits_i$year,
                                             permits_i$approved_acres[permits_i$year==year_i]/dat_dlnm$LK_LZACR_1[i],
                                             0)
    for(j in 2:9){
      treatment_history_matrix[i, j] <- ifelse((year_i+1-j)%in%permits_i$year,
                                               permits_i$approved_acres[permits_i$year==year_i+1-j]/dat_dlnm$LK_LZACR_1[i],
                                               0)
    }
    
  } else {
    
    treatment_history_matrix[i, ] <- rep(0, 9)
    target_string[i] <- "none"
    
  }
  
  print(i)
  
}

treatment_history_matrix[treatment_history_matrix>1] <- 1

treatment_history_matrix_long <- treatment_history_matrix %>%
  pivot_longer(lag0:lag8, names_to="Lag", values_to="% LZ Permitted")


active_treatments <- treatment_history_matrix[treatment_history_matrix>0]
treat_knots  <- quantile(active_treatments, probs = c(0.25, 0.50, 0.75), na.rm=T)

treatment_cb <- crossbasis(treatment_history_matrix[, 1:9],
                           lag=c(0, 8),
                           argvar=list(fun="ns"),
                           arglag=list(fun="ns", df=4))

treatment_cb_ewm <- crossbasis(treatment_history_matrix[dat_dlnm$DOW%in%dat_dlnm_ewm_lakes, 1:9],
                               lag=c(0, 8),
                               argvar=list(fun="ns"),
                               arglag=list(fun="ns", df=4))

treatment_cb_clp <- crossbasis(treatment_history_matrix[dat_dlnm$DOW%in%dat_dlnm_clp_lakes, 1:9],
                               lag=c(0, 8),
                               argvar=list(fun="ns"),
                               arglag=list(fun="ns", df=4))

dat_dlnm$YearsSince <- NA
dat_dlnm$CumulativeTrt <- NA
dat_dlnm$CumulativeTrtW <- NA
dat_dlnm$Consecutive <- NA

for(i in 1:dim(dat_dlnm)[1]){
  
  thm_i <- treatment_history_matrix[i, ]
  
  dat_dlnm$YearsSince[i] <- ifelse(sum(thm_i)==0, 0, min((0:8)[thm_i>0]))
  dat_dlnm$CumulativeTrt[i] <- sum(thm_i)
  dat_dlnm$CumulativeTrtW[i] <- sum(thm_i/(1:9))
  
  treated <- thm_i > 0
  runs <- rle(c(treated))
  
  dat_dlnm$Consecutive[i] <- ifelse(runs$values[1]==TRUE, runs$lengths[1], 0)
  
  print(i)
}


# Modeling ----

## Native richness ----
gam_dlnm_native <- gam(survey_native_rich ~ treatment_cb + 
                         s(APM_prop) + 
                         s(doy) + s(DOW, bs='re'),
                       family="nb",
                       data=dat_dlnm)

summary(gam_dlnm_native)

pred_native <- crosspred(treatment_cb, gam_dlnm_native,
                         at=seq(0, 1, 0.05),
                         cen=0)

plot(pred_native, ptype="contour", 
     xlab="Fraction Permitted", 
     ylab="Lag (Years)",
     key.title=title(main="Rel. Richness", cex.main=0.8))


pred_native_long <- pred_native$matRRfit %>%
  as.data.frame() %>%
  mutate(FractionPermitted=rownames(pred_native$matRRfit)) %>%
  pivot_longer(lag0:lag8, names_to="Lag", values_to="Estimate")

pred_native_long <- left_join(pred_native_long,
                              pred_native$matRRlow %>%
                                as.data.frame() %>%
                                mutate(FractionPermitted=rownames(pred_native$matRRfit)) %>%
                                pivot_longer(lag0:lag8, names_to="Lag", values_to="Lower"),
                              by=c("Lag", "FractionPermitted"))

pred_native_long <- left_join(pred_native_long,
                              pred_native$matRRhigh %>%
                                as.data.frame() %>%
                                mutate(FractionPermitted=rownames(pred_native$matRRfit)) %>%
                                pivot_longer(lag0:lag8, names_to="Lag", values_to="Upper"),
                              by=c("Lag", "FractionPermitted"))

pred_native_long <- pred_native_long %>%
  mutate(LagNum=as.numeric(str_remove(Lag, "lag")))


ggplot(pred_native_long %>%
         filter(FractionPermitted%in%c(0.05, 0.15, 0.75)),
       aes(x=LagNum, y=Estimate, 
           lty=factor(FractionPermitted))) + 
  geom_hline(yintercept=1) + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2) + 
  geom_line() +
  theme_classic(12)

## Native occurrence----
gam_dlnm_nativeFOO_binom <- gam(
  cbind(pres, abs) ~ treatment_cb + 
    s(APM_prop) + 
    s(doy) + s(DOW, bs='re'),
  family=binomial(),
  data=dat_dlnm %>%
    mutate(pres=round(native_vegetated*n_points),
           abs=n_points-pres)
)

summary(gam_dlnm_nativeFOO_binom)

pred_nativeFOO_binom <- crosspred(treatment_cb, gam_dlnm_nativeFOO_binom,
                                  at=seq(0, 1, 0.05),
                                  cen=0)

plot(pred_nativeFOO_binom, ptype="contour", 
     xlab="Fraction Permitted", 
     ylab="Lag (Years)")


pred_nativeFOO_binom_long <- pred_nativeFOO_binom$matRRfit %>%
  as.data.frame() %>%
  mutate(FractionPermitted=rownames(pred_nativeFOO_binom$matRRfit)) %>%
  pivot_longer(lag0:lag5, names_to="Lag", values_to="Estimate")

pred_nativeFOO_binom_long <- left_join(pred_nativeFOO_binom_long,
                                       pred_nativeFOO_binom$matRRlow %>%
                                         as.data.frame() %>%
                                         mutate(FractionPermitted=rownames(pred_nativeFOO_binom$matRRfit)) %>%
                                         pivot_longer(lag0:lag5, names_to="Lag", values_to="Lower"),
                                       by=c("Lag", "FractionPermitted"))

pred_nativeFOO_binom_long <- left_join(pred_nativeFOO_binom_long,
                                       pred_nativeFOO_binom$matRRhigh %>%
                                         as.data.frame() %>%
                                         mutate(FractionPermitted=rownames(pred_nativeFOO_binom$matRRfit)) %>%
                                         pivot_longer(lag0:lag5, names_to="Lag", values_to="Upper"),
                                       by=c("Lag", "FractionPermitted"))

pred_nativeFOO_binom_long <- pred_nativeFOO_binom_long %>%
  mutate(LagNum=as.numeric(str_remove(Lag, "lag")))


ggplot(pred_nativeFOO_binom_long %>%
         filter(FractionPermitted%in%c(0.05, 0.15, 0.75)),
       aes(x=LagNum, y=Estimate, 
           lty=factor(FractionPermitted))) + 
  geom_hline(yintercept=1) + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2) + 
  geom_line() +
  theme_classic(12)


## EWM effectiveness ----

gam_dlnm_ewm_binom <- gam(
  cbind(pres, abs) ~ treatment_cb_ewm + 
    s(APM_prop) + 
    s(doy) + s(DOW, bs='re'),
  family=quasibinomial(),
  data=dat_dlnm %>%
    filter(DOW%in%dat_dlnm_ewm_lakes) %>%
    mutate(pres=round(EWM*n_points),
           abs=n_points-pres))

summary(gam_dlnm_ewm_binom)

pred_ewm_binom <- crosspred(treatment_cb_ewm, gam_dlnm_ewm_binom,
                            at=seq(0, 1, 0.05),
                            cen=0)

plot(pred_ewm_binom, ptype="contour", 
     xlab="Fraction Permitted", 
     ylab="Lag (Years)")


pred_ewm_binom_long <- pred_ewm_binom$matRRfit %>%
  as.data.frame() %>%
  mutate(FractionPermitted=rownames(pred_ewm_binom$matRRfit)) %>%
  pivot_longer(lag0:lag5, names_to="Lag", values_to="Estimate")

pred_ewm_binom_long <- left_join(pred_ewm_binom_long,
                                 pred_ewm_binom$matRRlow %>%
                                   as.data.frame() %>%
                                   mutate(FractionPermitted=rownames(pred_ewm_binom$matRRfit)) %>%
                                   pivot_longer(lag0:lag5, names_to="Lag", values_to="Lower"),
                                 by=c("Lag", "FractionPermitted"))

pred_ewm_binom_long <- left_join(pred_ewm_binom_long,
                                 pred_ewm_binom$matRRhigh %>%
                                   as.data.frame() %>%
                                   mutate(FractionPermitted=rownames(pred_ewm_binom$matRRfit)) %>%
                                   pivot_longer(lag0:lag5, names_to="Lag", values_to="Upper"),
                                 by=c("Lag", "FractionPermitted"))

pred_ewm_binom_long <- pred_ewm_binom_long %>%
  mutate(LagNum=as.numeric(str_remove(Lag, "lag")))


ggplot(pred_ewm_binom_long %>%
         filter(FractionPermitted%in%c(0.05, 0.15, 0.75)),
       aes(x=LagNum, y=Estimate, 
           lty=factor(FractionPermitted))) + 
  geom_hline(yintercept=1) + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2) + 
  geom_line() +
  theme_classic(12)

## CLP effectiveness ----
gam_dlnm_clp_binom <- gam(
  cbind(pres, abs) ~ treatment_cb_clp + 
    s(APM_prop) + 
    s(doy) + s(DOW, bs='re'),
  family=binomial(),
  data=dat_dlnm %>%
    filter(DOW%in%dat_dlnm_clp_lakes) %>%
    mutate(pres=round(potamogeton_crispus*n_points),
           abs=n_points-pres)
)

summary(gam_dlnm_clp_binom)

pred_clp_binom <- crosspred(treatment_cb_clp, gam_dlnm_clp_binom,
                            at=seq(0, 1, 0.05),
                            cen=0)

plot(pred_clp_binom, ptype="contour", 
     xlab="Fraction Permitted", 
     ylab="Lag (Years)")

pred_clp_binom_long <- pred_clp_binom$matRRfit %>%
  as.data.frame() %>%
  mutate(FractionPermitted=rownames(pred_clp_binom$matRRfit)) %>%
  pivot_longer(lag0:lag5, names_to="Lag", values_to="Estimate")

pred_clp_binom_long <- left_join(pred_clp_binom_long,
                                 pred_clp_binom$matRRlow %>%
                                   as.data.frame() %>%
                                   mutate(FractionPermitted=rownames(pred_clp_binom$matRRfit)) %>%
                                   pivot_longer(lag0:lag5, names_to="Lag", values_to="Lower"),
                                 by=c("Lag", "FractionPermitted"))

pred_clp_binom_long <- left_join(pred_clp_binom_long,
                                 pred_clp_binom$matRRhigh %>%
                                   as.data.frame() %>%
                                   mutate(FractionPermitted=rownames(pred_clp_binom$matRRfit)) %>%
                                   pivot_longer(lag0:lag5, names_to="Lag", values_to="Upper"),
                                 by=c("Lag", "FractionPermitted"))

pred_clp_binom_long <- pred_clp_binom_long %>%
  mutate(LagNum=as.numeric(str_remove(Lag, "lag")))

ggplot(pred_clp_binom_long %>%
         filter(FractionPermitted%in%c(0.05, 0.15, 0.75)),
       aes(x=LagNum, y=Estimate, 
           lty=factor(FractionPermitted))) + 
  geom_hline(yintercept=1) + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2) + 
  geom_line() +
  theme_classic(12)


## Accounting for streaks ----
# Recency: Years since last treatment
# Intensity: Cumulative treatment history
# Streak: consecutive years of treatment

### Models ----
gam_dlnm_native2 <- gam(survey_native_rich ~ treatment_cb + 
                          s(Consecutive, k=4, bs='ts') + 
                          s(APM_prop, k=4, bs='ts') + 
                          s(doy, k=4) + s(DOW, bs='re') + s(fyear, bs='re'),
                        family="nb",
                        data=dat_dlnm)

summary(gam_dlnm_native2)
draw(gam_dlnm_native2)

gam_dlnm_nativeFOO2_binom <- gam(
  cbind(pres, abs) ~ treatment_cb + 
    s(Consecutive, k=4, bs='ts') + 
    s(APM_prop, k=4, bs='ts') + 
    s(doy, k=4) + s(DOW, bs='re') + s(fyear, bs='re'),
  family=quasibinomial(),
  data=dat_dlnm %>%
    mutate(pres=round(native_vegetated*n_points),
           abs=n_points-pres)
)

summary(gam_dlnm_nativeFOO2_binom)
draw(gam_dlnm_nativeFOO2_binom)

gam_dlnm_ewm2_binom <- gam(
  cbind(pres, abs) ~ treatment_cb_ewm + 
    s(Consecutive, k=4, bs='ts') + 
    s(APM_prop, k=4, bs='ts') + 
    s(doy, k=4) + s(DOW, bs='re') + s(fyear, bs='re'),
  family=quasibinomial(),
  data=dat_dlnm %>%
    filter(DOW%in%dat_dlnm_ewm_lakes) %>%
    mutate(pres=round(EWM*n_points),
           abs=n_points-pres)
)

summary(gam_dlnm_ewm2_binom)
draw(gam_dlnm_ewm2_binom)

gam_dlnm_clp2_binom <- gam(
  cbind(pres, abs) ~ treatment_cb_clp + 
    s(Consecutive, k=4, bs='ts') + 
    s(APM_prop, k=4, bs='ts') + 
    s(doy, k=4) + s(DOW, bs='re') + s(fyear, bs='re'),
  family=quasibinomial(),
  data=dat_dlnm %>%
    filter(DOW%in%dat_dlnm_clp_lakes) %>%
    mutate(pres=round(potamogeton_crispus*n_points),
           abs=n_points-pres)
)

summary(gam_dlnm_clp2_binom)
draw(gam_dlnm_clp2_binom)


### Predictions ----

# Native FOO
pred_native2 <- crosspred(treatment_cb, gam_dlnm_native2,
                          at=seq(0, 1, 0.05),
                          cen=0)

pred_native2_long <- pred_native2$matRRfit %>%
  as.data.frame() %>%
  mutate(FractionPermitted=rownames(pred_native2$matRRfit)) %>%
  pivot_longer(lag0:lag5, names_to="Lag", values_to="Estimate")

pred_native2_long <- left_join(pred_native2_long,
                               pred_native2$matRRlow %>%
                                 as.data.frame() %>%
                                 mutate(FractionPermitted=rownames(pred_native2$matRRfit)) %>%
                                 pivot_longer(lag0:lag5, names_to="Lag", values_to="Lower"),
                               by=c("Lag", "FractionPermitted"))

pred_native2_long <- left_join(pred_native2_long,
                               pred_native2$matRRhigh %>%
                                 as.data.frame() %>%
                                 mutate(FractionPermitted=rownames(pred_native2$matRRfit)) %>%
                                 pivot_longer(lag0:lag5, names_to="Lag", values_to="Upper"),
                               by=c("Lag", "FractionPermitted"))

pred_native2_long <- pred_native2_long %>%
  mutate(LagNum=as.numeric(str_remove(Lag, "lag")))

# Native FOO
pred_nativeFOO2_binom <- crosspred(treatment_cb, gam_dlnm_nativeFOO2_binom,
                                   at=seq(0, 1, 0.05),
                                   cen=0)

pred_nativeFOO2_binom_long <- pred_nativeFOO2_binom$matRRfit %>%
  as.data.frame() %>%
  mutate(FractionPermitted=rownames(pred_nativeFOO2_binom$matRRfit)) %>%
  pivot_longer(lag0:lag5, names_to="Lag", values_to="Estimate")

pred_nativeFOO2_binom_long <- left_join(pred_nativeFOO2_binom_long,
                                        pred_nativeFOO2_binom$matRRlow %>%
                                          as.data.frame() %>%
                                          mutate(FractionPermitted=rownames(pred_nativeFOO2_binom$matRRfit)) %>%
                                          pivot_longer(lag0:lag5, names_to="Lag", values_to="Lower"),
                                        by=c("Lag", "FractionPermitted"))

pred_nativeFOO2_binom_long <- left_join(pred_nativeFOO2_binom_long,
                                        pred_nativeFOO2_binom$matRRhigh %>%
                                          as.data.frame() %>%
                                          mutate(FractionPermitted=rownames(pred_nativeFOO2_binom$matRRfit)) %>%
                                          pivot_longer(lag0:lag5, names_to="Lag", values_to="Upper"),
                                        by=c("Lag", "FractionPermitted"))

pred_nativeFOO2_binom_long <- pred_nativeFOO2_binom_long %>%
  mutate(LagNum=as.numeric(str_remove(Lag, "lag")))


# EWM
pred_ewm2_binom <- crosspred(treatment_cb_ewm, gam_dlnm_ewm2_binom,
                             at=seq(0, 1, 0.05),
                             cen=0)

pred_ewm2_binom_long <- pred_ewm2_binom$matRRfit %>%
  as.data.frame() %>%
  mutate(FractionPermitted=rownames(pred_ewm2_binom$matRRfit)) %>%
  pivot_longer(lag0:lag5, names_to="Lag", values_to="Estimate")

pred_ewm2_binom_long <- left_join(pred_ewm2_binom_long,
                                  pred_ewm2_binom$matRRlow %>%
                                    as.data.frame() %>%
                                    mutate(FractionPermitted=rownames(pred_ewm2_binom$matRRfit)) %>%
                                    pivot_longer(lag0:lag5, names_to="Lag", values_to="Lower"),
                                  by=c("Lag", "FractionPermitted"))

pred_ewm2_binom_long <- left_join(pred_ewm2_binom_long,
                                  pred_ewm2_binom$matRRhigh %>%
                                    as.data.frame() %>%
                                    mutate(FractionPermitted=rownames(pred_ewm2_binom$matRRfit)) %>%
                                    pivot_longer(lag0:lag5, names_to="Lag", values_to="Upper"),
                                  by=c("Lag", "FractionPermitted"))

pred_ewm2_binom_long <- pred_ewm2_binom_long %>%
  mutate(LagNum=as.numeric(str_remove(Lag, "lag")))


# CLP
pred_clp2_binom <- crosspred(treatment_cb_clp, gam_dlnm_clp2_binom,
                             at=seq(0, 1, 0.05),
                             cen=0)

pred_clp2_binom_long <- pred_clp2_binom$matRRfit %>%
  as.data.frame() %>%
  mutate(FractionPermitted=rownames(pred_clp2_binom$matRRfit)) %>%
  pivot_longer(lag0:lag5, names_to="Lag", values_to="Estimate")

pred_clp2_binom_long <- left_join(pred_clp2_binom_long,
                                  pred_clp2_binom$matRRlow %>%
                                    as.data.frame() %>%
                                    mutate(FractionPermitted=rownames(pred_clp2_binom$matRRfit)) %>%
                                    pivot_longer(lag0:lag5, names_to="Lag", values_to="Lower"),
                                  by=c("Lag", "FractionPermitted"))

pred_clp2_binom_long <- left_join(pred_clp2_binom_long,
                                  pred_clp2_binom$matRRhigh %>%
                                    as.data.frame() %>%
                                    mutate(FractionPermitted=rownames(pred_clp2_binom$matRRfit)) %>%
                                    pivot_longer(lag0:lag5, names_to="Lag", values_to="Upper"),
                                  by=c("Lag", "FractionPermitted"))

pred_clp2_binom_long <- pred_clp2_binom_long %>%
  mutate(LagNum=as.numeric(str_remove(Lag, "lag")))

# Visualization ----

## GAM summary plots, minus crossbasis ----

png("JakeOutput/DLNM_gamsummaries.png", width=7, height=9, units='in', res=300)
plot_grid(
  
  draw(gam_dlnm_native2, nrow=1) & 
    plot_annotation(title="a) Native taxonomic richness") & 
    theme_minimal(9) & 
    theme(plot.title=element_text(face="bold"),
          plot.title.position='plot'),
  
  draw(gam_dlnm_nativeFOO2_binom, nrow=1) & 
    plot_annotation(title="b) Native occurrence") & 
    theme_minimal(9) & 
    theme(plot.title=element_text(face="bold"),
          plot.title.position='plot'),
  
  draw(gam_dlnm_ewm2_binom, nrow=1) & 
    plot_annotation(title="c) EWM occurrence") & 
    theme_minimal(9) & 
    theme(plot.title=element_text(face="bold"),
          plot.title.position='plot'),
  
  draw(gam_dlnm_clp2_binom, nrow=1) & 
    plot_annotation(title="d) CLP occurrence") & 
    theme_minimal(9) & 
    theme(plot.title=element_text(face="bold"),
          plot.title.position='plot'),
  
  nrow=4)
dev.off()


## Cross-basis contour plots ----
p_native2_contour <- ggplot(as.data.frame(as.table(pred_native2$matfit)) %>%
                              rename(Fraction=Var1, Lag=Var2, Fit=Freq) %>%
                              mutate(Fraction=as.numeric(as.character(Fraction)),
                                     Lag=as.numeric(str_remove(as.character(Lag), "lag")),
                                     SE=as.vector(pred_native2$matse),
                                     Low=Fit - 1.96 * SE,
                                     High=Fit + 1.96 * SE,
                                     is_sig=(Low > 0 & High > 0) | (Low < 0 & High < 0)), 
                            aes(x=Lag, y=Fraction, fill=Fit)) +
  geom_raster(interpolate=TRUE) +
  scale_fill_scico(palette="vik", name="Predicted\nEffect", 
                   midpoint=0) +
  stat_contour(aes(z = as.numeric(is_sig)),
               breaks = 0.5,
               color = "black",
               linewidth = 0.8,
               linetype = "dashed") +
  geom_vline(aes(xintercept=6.5), lty=3, col='orange', lwd=1.25) + 
  labs(y="Fraction Permitted",
       x="Lag (Years)",
       title="a) Effect on native taxonomic richness (all NS)") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw(10) +
  theme(plot.title=element_text(size=11, face="bold", hjust=0),
        panel.grid=element_blank(),
        plot.title.position='plot')


p_nativeFOO2_contour <- ggplot(as.data.frame(as.table(pred_nativeFOO2_binom$matfit)) %>%
                                 rename(Fraction=Var1, Lag=Var2, Fit=Freq) %>%
                                 mutate(Fraction=as.numeric(as.character(Fraction)),
                                        Lag=as.numeric(str_remove(as.character(Lag), "lag")),
                                        SE=as.vector(pred_nativeFOO2_binom$matse),
                                        Low=Fit - 1.96 * SE,
                                        High=Fit + 1.96 * SE,
                                        is_sig=(Low > 0 & High > 0) | (Low < 0 & High < 0)),  
                               aes(x=Lag, y=Fraction, fill=Fit)) +
  geom_raster(interpolate=TRUE) +
  stat_contour(aes(z = as.numeric(is_sig)),
               breaks = 0.5,
               color = "black",
               linewidth = 0.8,
               linetype = "dashed") +
  geom_vline(aes(xintercept=6.5), lty=3, col='orange', lwd=1.25) + 
  scale_fill_scico(palette="vik", name="Predicted\nEffect\n(Log-OR)", midpoint=0) +
  labs(y="Fraction Permitted",
       x="Lag (Years)",
       title="b) Effect on native occurrence") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw(10) +
  theme(plot.title=element_text(size=11, face="bold", hjust=0),
        panel.grid=element_blank(),
        plot.title.position='plot')

p_ewm2_contour <- ggplot(as.data.frame(as.table(pred_ewm2_binom$matfit)) %>%
                           rename(Fraction=Var1, Lag=Var2, Fit=Freq) %>%
                           mutate(Fraction=as.numeric(as.character(Fraction)),
                                  Lag=as.numeric(str_remove(as.character(Lag), "lag")),
                                  SE=as.vector(pred_ewm2_binom$matse),
                                  Low=Fit - 1.96 * SE,
                                  High=Fit + 1.96 * SE,
                                  is_sig=(Low > 0 & High > 0) | (Low < 0 & High < 0)),  
                         aes(x=Lag, y=Fraction, fill=Fit)) +
  geom_raster(interpolate=TRUE) +
  stat_contour(aes(z = as.numeric(is_sig)),
               breaks = 0.5,
               color = "black",
               linewidth = 0.8,
               linetype = "dashed") +
  geom_vline(aes(xintercept=6.5), lty=3, col='orange', lwd=1.25) + 
  scale_fill_scico(palette="vik", name="Predicted\nEffect\n(Log-OR)", midpoint=0) +
  labs(y="Fraction Permitted",
       x="Lag (Years)",
       title="c) Effect on EWM occurrence") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw(10) +
  theme(plot.title=element_text(size=11, face="bold", hjust=0),
        panel.grid=element_blank(),
        plot.title.position='plot')

p_clp2_contour <- ggplot(as.data.frame(as.table(pred_clp2_binom$matfit)) %>%
                           rename(Fraction=Var1, Lag=Var2, Fit=Freq) %>%
                           mutate(Fraction=as.numeric(as.character(Fraction)),
                                  Lag=as.numeric(str_remove(as.character(Lag), "lag")),
                                  SE=as.vector(pred_clp2_binom$matse),
                                  Low=Fit - 1.96 * SE,
                                  High=Fit + 1.96 * SE,
                                  is_sig=(Low > 0 & High > 0) | (Low < 0 & High < 0)), 
                         aes(x=Lag, y=Fraction, fill=Fit)) +
  geom_raster(interpolate=TRUE) +
  stat_contour(aes(z = as.numeric(is_sig)),
               breaks = 0.5,
               color = "black",
               linewidth = 0.8,
               linetype = "dashed") +
  geom_vline(aes(xintercept=6.5), lty=3, col='orange', lwd=1.25) + 
  scale_fill_scico(palette="vik", name="Predicted\nEffect\n(Log-OR)", midpoint=0) +
  labs(y="Fraction Permitted",
       x="Lag (Years)",
       title="d) Effect on CLP occurrence") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw(10) +
  theme(plot.title=element_text(size=11, face="bold", hjust=0),
        panel.grid=element_blank(),
        plot.title.position='plot')

p_contours <- plot_grid(p_native2_contour, p_nativeFOO2_contour,
                        p_ewm2_contour, p_clp2_contour,
                        nrow=2)

png("JakeOutput/DLNM_CrossBasisContours.png", width=7, height=5, units='in', res=300)
p_contours
dev.off()

## Cross-basis contour plots - Odds Ratio Scale ----
p_native2_contour_OR <- ggplot(as.data.frame(as.table(pred_native2$matfit)) %>%
                                 rename(Fraction=Var1, Lag=Var2, Fit=Freq) %>%
                                 mutate(Fraction=as.numeric(as.character(Fraction)),
                                        Lag=as.numeric(str_remove(as.character(Lag), "lag")),
                                        SE=as.vector(pred_native2$matse),
                                        Low=Fit - 1.96 * SE,
                                        High=Fit + 1.96 * SE,
                                        is_sig=(Low > 0 & High > 0) | (Low < 0 & High < 0)) %>%
                                 filter(Lag<=6), 
                               aes(x=Lag, y=Fraction, fill=Fit)) +
  geom_raster(interpolate=TRUE) +
  scale_fill_scico(palette="vik", name="Predicted\nEffect", 
                   midpoint=0) +
  stat_contour(aes(z = as.numeric(is_sig)),
               breaks = 0.5,
               color = "black",
               linewidth = 0.8,
               linetype = "dashed") +
  labs(y="Fraction Permitted",
       x="Lag (Years)",
       title="a) Effect on native taxonomic richness (all NS)") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw(10) +
  theme(plot.title=element_text(size=11, face="bold", hjust=0),
        panel.grid=element_blank(),
        plot.title.position='plot')

p_nativeFOO2_contour_OR <- ggplot(as.data.frame(as.table(pred_nativeFOO2_binom$matfit)) %>%
                                    rename(Fraction=Var1, Lag=Var2, Fit=Freq) %>%
                                    mutate(Fraction=as.numeric(as.character(Fraction)),
                                           Lag=as.numeric(str_remove(as.character(Lag), "lag")),
                                           SE=as.vector(pred_nativeFOO2_binom$matse),
                                           Low=Fit - 1.96 * SE,
                                           High=Fit + 1.96 * SE,
                                           is_sig=(Low > 0 & High > 0) | (Low < 0 & High < 0)) %>%
                                    filter(Lag<=6),  
                                  aes(x=Lag, y=Fraction, fill=Fit)) +
  geom_raster(interpolate=TRUE) +
  stat_contour(aes(z = as.numeric(is_sig)),
               breaks = 0.5,
               color = "black",
               linewidth = 0.8,
               linetype = "dashed") +
  scale_fill_scico(palette="vik", name="Predicted\nEffect (OR)", midpoint=0,
                   breaks=log(c(0.36, 0.5, 0.75, 1, 2, 3.7)),
                   labels=c("0.36", "0.50", "0.75", "1.0", "2.0", "3.7"),
                   limits=log(c(0.36, 3.8))) +
  labs(y="Fraction Permitted",
       x="Lag (Years)",
       title="b) Effect on native occurrence") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw(10) +
  theme(plot.title=element_text(size=11, face="bold", hjust=0),
        panel.grid=element_blank(),
        plot.title.position='plot')

p_ewm2_contour_OR <- ggplot(as.data.frame(as.table(pred_ewm2_binom$matfit)) %>%
                              rename(Fraction=Var1, Lag=Var2, Fit=Freq) %>%
                              mutate(Fraction=as.numeric(as.character(Fraction)),
                                     Lag=as.numeric(str_remove(as.character(Lag), "lag")),
                                     SE=as.vector(pred_ewm2_binom$matse),
                                     Low=Fit - 1.96 * SE,
                                     High=Fit + 1.96 * SE,
                                     is_sig=(Low > 0 & High > 0) | (Low < 0 & High < 0)) %>%
                              filter(Lag<=6),  
                            aes(x=Lag, y=Fraction, fill=Fit)) +
  geom_raster(interpolate=TRUE) +
  stat_contour(aes(z = as.numeric(is_sig)),
               breaks = 0.5,
               color = "black",
               linewidth = 0.8,
               linetype = "dashed") +
  scale_fill_scico(palette="vik", name="Predicted\nEffect (OR)", midpoint=0,
                   breaks=log(c(0.02, 0.05, 0.1, 0.25, 0.5, 1, 5, 12)),
                   labels=c("0.02", "0.05", "0.10", "0.25", "0.50", "1.0", "5.0", "12"),
                   limit=log(c(0.02, 12.3))) +
  labs(y="Fraction Permitted",
       x="Lag (Years)",
       title="c) Effect on EWM occurrence") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw(10) +
  theme(plot.title=element_text(size=11, face="bold", hjust=0),
        panel.grid=element_blank(),
        plot.title.position='plot')

p_clp2_contour_OR <- ggplot(as.data.frame(as.table(pred_clp2_binom$matfit)) %>%
                              rename(Fraction=Var1, Lag=Var2, Fit=Freq) %>%
                              mutate(Fraction=as.numeric(as.character(Fraction)),
                                     Lag=as.numeric(str_remove(as.character(Lag), "lag")),
                                     SE=as.vector(pred_clp2_binom$matse),
                                     Low=Fit - 1.96 * SE,
                                     High=Fit + 1.96 * SE,
                                     is_sig=(Low > 0 & High > 0) | (Low < 0 & High < 0)) %>%
                              filter(Lag<=6), 
                            aes(x=Lag, y=Fraction, fill=Fit)) +
  geom_raster(interpolate=TRUE) +
  stat_contour(aes(z = as.numeric(is_sig)),
               breaks = 0.5,
               color = "black",
               linewidth = 0.8,
               linetype = "dashed") +
  scale_fill_scico(palette="vik", name="Predicted\nEffect (OR)", midpoint=0,
                   breaks=log(c(0.48, 0.6, 0.8, 1, 1.2)),
                   labels=c("0.48", "0.60", "0.80", "1.0", "1.2"),
                   limits=log(c(0.47, 1.2))) +
  labs(y="Fraction Permitted",
       x="Lag (Years)",
       title="d) Effect on CLP occurrence") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw(10) +
  theme(plot.title=element_text(size=11, face="bold", hjust=0),
        panel.grid=element_blank(),
        plot.title.position='plot')

p_contours_OR <- plot_grid(p_native2_contour_OR, p_nativeFOO2_contour_OR,
                           p_ewm2_contour_OR, p_clp2_contour_OR,
                           nrow=2)

png("JakeOutput/DLNM_CrossBasisContours_OR.png", width=7, height=5, units='in', res=300)
p_contours_OR
dev.off()


## Discrete predictions ----
ggplot(pred_native2_long %>%
         filter(FractionPermitted%in%c(0.05, 0.15, 0.75)),
       aes(x=LagNum, y=Estimate, 
           lty=factor(FractionPermitted))) + 
  geom_hline(yintercept=1) + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2) + 
  geom_line() +
  theme_classic(12)

ggplot(pred_nativeFOO2_binom_long %>%
         filter(FractionPermitted%in%c(0.05, 0.15, 0.75)),
       aes(x=LagNum, y=Estimate, 
           lty=factor(FractionPermitted))) + 
  geom_hline(yintercept=1) + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2) + 
  geom_line() +
  theme_classic(12)

ggplot(pred_ewm2_binom_long %>%
         filter(FractionPermitted%in%c(0.05, 0.15, 0.75)),
       aes(x=LagNum, y=Estimate, 
           lty=factor(FractionPermitted))) + 
  geom_hline(yintercept=1) + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2) + 
  geom_line() +
  theme_classic(12)

ggplot(pred_clp2_binom_long %>%
         filter(FractionPermitted%in%c(0.05, 0.15, 0.75)),
       aes(x=LagNum, y=Estimate, 
           lty=factor(FractionPermitted))) + 
  geom_hline(yintercept=1) + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2) + 
  geom_line() +
  theme_classic(12)

## Which lakes are driving the lag 7/8 extremes for EWM and CLP? ----

late_lag_lakes <- unique((treatment_history_matrix_long %>%
                            mutate(DOW=rep(dat_dlnm$DOW, each=9),
                                   year=rep(dat_dlnm$year, each=9)) %>%
                            filter(Lag=="lag7"|Lag=="lag8") %>%
                            group_by(DOW, year, Lag) %>%
                            summarize(`% LZ Permitted`=max(`% LZ Permitted`)) %>%
                            filter(`% LZ Permitted`>0))$DOW)

late_lag_lakes_trt_years <- treatment_history_matrix_long %>%
  mutate(DOW=rep(dat_dlnm$DOW, each=9),
         year=rep(dat_dlnm$year, each=9),
         trt=ifelse(`% LZ Permitted`>0, 1, 0),
         trt_year=year-as.numeric(str_remove(Lag, "lag"))) %>%
  filter(DOW%in%late_lag_lakes & trt==1) %>%
  dplyr::select(DOW, `% LZ Permitted`, trt, trt_year) %>%
  distinct()

# 9 lakes isn't nothing...

p_lag7or8 <- ggplot() + 
  facet_wrap(~DOW, nrow=4) + 
  geom_vline(data=late_lag_lakes_trt_years,
             aes(xintercept=trt_year)) +
  geom_point(data=dat_dlnm %>%
               filter(DOW%in%late_lag_lakes) %>%
               mutate(Taxa="EWM"),
             aes(x=year, y=EWM, pch=Taxa, col=Taxa), alpha=0.6) + 
  geom_point(data=dat_dlnm %>%
               filter(DOW%in%late_lag_lakes) %>%
               mutate(Taxa="CLP"),
             aes(x=year, y=potamogeton_crispus, pch=Taxa, col=Taxa), alpha=0.6) + 
  geom_point(data=dat_dlnm %>%
               filter(DOW%in%late_lag_lakes) %>%
               mutate(Taxa="Native"),
             aes(x=year, y=native_vegetated, pch=Taxa, col=Taxa), alpha=0.6) + 
  scale_y_continuous("Frequency of Occurrence",
                     breaks=seq(0, 1, 0.25), minor_breaks=seq(0, 1, 0.25)) + 
  scale_x_continuous(breaks=seq(2017, 2025, 4), minor_breaks=2017:2025) + 
  scale_color_manual("Taxa", values=viridis(4)[1:3]) + 
  scale_linetype_manual("Taxa") + 
  ggtitle("Lakes with permits at lag 7 or lag 8 of at least one survey\n(i.e., above the orange dotted line in cross-basis plots)") + 
  theme_minimal(11) + theme(legend.position='top',
                            legend.title=element_blank(),
                            axis.title.x=element_blank(),
                            plot.title.position="plot")

png("JakeOutput/DLNM_Lag 7 or 8 lakes.png", width=5.5, height=4, units='in', res=300)
p_lag7or8
dev.off()

## Which lakes are driving the lag 4 and up rebound in natives? ----

late4_lag_lakes <- unique((treatment_history_matrix_long %>%
                             mutate(DOW=rep(dat_dlnm$DOW, each=9),
                                    year=rep(dat_dlnm$year, each=9)) %>%
                             filter(Lag=="lag4"|Lag=="lag5"|Lag=="lag6"|
                                      Lag=="lag7"|Lag=="lag8") %>%
                             group_by(DOW, year, Lag) %>%
                             summarize(`% LZ Permitted`=max(`% LZ Permitted`)) %>%
                             filter(`% LZ Permitted`>0))$DOW)

late4_lag_lakes_trt_years <- treatment_history_matrix_long %>%
  mutate(DOW=rep(dat_dlnm$DOW, each=9),
         year=rep(dat_dlnm$year, each=9),
         trt=ifelse(`% LZ Permitted`>0, 1, 0),
         trt_year=year-as.numeric(str_remove(Lag, "lag"))) %>%
  filter(DOW%in%late4_lag_lakes & trt==1) %>%
  dplyr::select(DOW, `% LZ Permitted`, trt, trt_year) %>%
  distinct()

p_lag4up <- ggplot() + 
  facet_wrap(~DOW, nrow=8) + 
  geom_vline(data=late4_lag_lakes_trt_years,
             aes(xintercept=trt_year)) +
  geom_point(data=dat_dlnm %>%
               filter(DOW%in%late4_lag_lakes) %>%
               mutate(Taxa="EWM"),
             aes(x=year, y=EWM, pch=Taxa, col=Taxa), alpha=0.6) + 
  geom_point(data=dat_dlnm %>%
               filter(DOW%in%late4_lag_lakes) %>%
               mutate(Taxa="CLP"),
             aes(x=year, y=potamogeton_crispus, pch=Taxa, col=Taxa), alpha=0.6) + 
  geom_point(data=dat_dlnm %>%
               filter(DOW%in%late4_lag_lakes) %>%
               mutate(Taxa="Native"),
             aes(x=year, y=native_vegetated, pch=Taxa, col=Taxa), alpha=0.6) + 
  scale_y_continuous("Frequency of Occurrence",
                     breaks=seq(0, 1, 0.25), minor_breaks=seq(0, 1, 0.25)) + 
  scale_x_continuous(breaks=seq(2013, 2025, 4), minor_breaks=2013:2025, limits=c(2013, 2025)) + 
  scale_color_manual("Taxa", values=viridis(4)[1:3]) + 
  scale_linetype_manual("Taxa") + 
  ggtitle("Lakes with permits lagged 4 or higher for at least one survey\n(e.g., the native occurrence rebound in the contour plots)") + 
  theme_minimal(9) + theme(legend.position='top',
                           legend.title=element_blank(),
                           axis.title.x=element_blank(),
                           plot.title.position="plot")

png("JakeOutput/DLNM_Lag 4 and up lakes.png", width=7, height=9, units='in', res=300)
p_lag4up
dev.off()


# Extra exploration: what do treatment trajectories/schedules look like? ----
# Looking to account for massive decline in EWM and CLP in larger treatments
# after increase from lag ~4-6.

yearof_trts <- data.frame(DOW=dat_dlnm$DOW,
                          year=dat_dlnm$year,
                          approved_acres=treatment_history_matrix[, 1],
                          ewm=dat_dlnm$EWM,
                          clp=dat_dlnm$potamogeton_crispus) %>%
  mutate(trt=ifelse(approved_acres>0, 1, 0),
         lag_ewm=lag(ewm, 1),
         lag_clp=lag(clp, 1))

yearof_trts <- yearof_trts %>%
  arrange(DOW, year) %>%
  group_by(DOW) %>%
  mutate(trt_start=if_else(trt == 1 & lag(trt, default=0) == 0, 1, 0),
         episode=cumsum(trt_start),
         episode_start_year=if_else(trt_start == 1, year, NA_integer_)) %>%
  fill(episode_start_year, .direction="downup") %>%
  mutate(trt_time=year - episode_start_year) %>%
  ungroup()

treatment_gaps <- yearof_trts %>%
  filter(trt > 0) %>%  
  arrange(DOW, year) %>%
  group_by(DOW) %>%
  mutate(prev_treatment_year=lag(year),
         years_between=year - prev_treatment_year) %>%
  filter(!is.na(years_between))

treatment_gaps %>%
  ungroup() %>%
  summarize(median_gap=median(years_between),
            mean_gap=round(mean(years_between), 2),
            iqr_lower=quantile(years_between, 0.25),
            iqr_upper=quantile(years_between, 0.75),
            total_retreatments=n())

yearof_trts_model <- yearof_trts %>%
  filter(trt_time >= 1, !(is.na(lag_ewm)|is.na(lag_clp)))

gam_retrt_ewm <- gam(trt ~ trt_time*lag_ewm,
                     family=quasibinomial(link="logit"),
                     data=yearof_trts_model %>%
                       filter(DOW%in%dat_dlnm_ewm_lakes))

gam_retrt_clp <- gam(trt ~ trt_time*lag_clp,
                     family=quasibinomial(link="logit"),
                     data=yearof_trts_model %>%
                       filter(DOW%in%dat_dlnm_clp_lakes))
summary(gam_retrt_ewm)
summary(gam_retrt_clp)

# Extra exploration: Subset to lakes with larger programs ----

# How do we quantify a larger treatment program?
# Maybe lakes with at least one treatment greater than 10%?

max_trt_bylake <- apply(treatment_history_matrix, MARGIN=1, FUN=max)

hist(max_trt_bylake)

higher_intensity_lakes <- unique(as.character(dat_dlnm$DOW[max_trt_bylake>=0.1]))

dat_dlnm_lg <- dat_dlnm[dat_dlnm$DOW%in%higher_intensity_lakes,]

treatment_cb_lg <- crossbasis(treatment_history_matrix[dat_dlnm$DOW%in%higher_intensity_lakes, 1:6],
                              lag=c(0, 5),
                              argvar=list(fun="lin"),
                              arglag=list(fun="ns", df=4))

gam_dlnm_native2_lg <- gam(survey_native_rich ~ treatment_cb_lg + 
                             s(Consecutive, k=4, bs='ts') + 
                             s(APM_prop, k=4, bs='ts') + 
                             s(doy, k=4) + s(DOW, bs='re') + s(fyear, bs='re'),
                           family="nb",
                           data=dat_dlnm_lg)

summary(gam_dlnm_native2_lg)

pred_native2_lg <- crosspred(treatment_cb_lg, gam_dlnm_native2_lg,
                             at=seq(0, 1, 0.05),
                             cen=0)

plot(pred_native2_lg, ptype="contour", 
     xlab="Fraction Permitted", 
     ylab="Lag (Years)",
     key.title=title(main="Rel. Richness", cex.main=0.8))


pred_native2_lg_long <- pred_native2_lg$matRRfit %>%
  as.data.frame() %>%
  mutate(FractionPermitted=rownames(pred_native2_lg$matRRfit)) %>%
  pivot_longer(lag0:lag5, names_to="Lag", values_to="Estimate")

pred_native2_lg_long <- left_join(pred_native2_lg_long,
                                  pred_native2_lg$matRRlow %>%
                                    as.data.frame() %>%
                                    mutate(FractionPermitted=rownames(pred_native2_lg$matRRfit)) %>%
                                    pivot_longer(lag0:lag5, names_to="Lag", values_to="Lower"),
                                  by=c("Lag", "FractionPermitted"))

pred_native2_lg_long <- left_join(pred_native2_lg_long,
                                  pred_native2_lg$matRRhigh %>%
                                    as.data.frame() %>%
                                    mutate(FractionPermitted=rownames(pred_native2_lg$matRRfit)) %>%
                                    pivot_longer(lag0:lag5, names_to="Lag", values_to="Upper"),
                                  by=c("Lag", "FractionPermitted"))

pred_native2_lg_long <- pred_native2_lg_long %>%
  mutate(LagNum=as.numeric(str_remove(Lag, "lag")))


ggplot(pred_native2_lg_long %>%
         filter(FractionPermitted%in%c(0.05, 0.15, 0.75)),
       aes(x=LagNum, y=Estimate, 
           lty=factor(FractionPermitted))) + 
  geom_hline(yintercept=1) + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2) + 
  geom_line() +
  theme_classic(12)