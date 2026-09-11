#'---
#' title: "Supplemental Code for Verhoeven et al."
#' author: "Mike Verhoeven"
#' output: 
#'    html_document:
#'       toc: true
#'       theme: default
#'       toc_depth: 2
#'       toc_float:
#'           collapsed: false
#'---


#' This code will load in the MN CLP dataset and analyze and visualize the
#' material presented in Verhoeven et al. 2019. Additonal dataset assembly and 
#' processing code are available upon request from the authors.

# header ------------------------------------------------------------------



#' # Tools
#' Load libraries:

  #+warning=FALSE, message=FALSE 
  library(ggplot2)  
  library(lme4)
  library(dplyr)
  library(data.table)
  library(numDeriv)
  library(sf)
  library(mapdata)
  library(gridExtra)
  library(rlang)
  library(devtools)
  devtools::install_github('oswaldosantos/ggsn')
  library(ggsn)
  library(raster)
  library(ggthemes)
  library(mapproj)
  library(maps)
  library(GISTools) 
  library(ggmap)

  #' Load custom functions:
  
  # diagnostics
  mv_diag <- function(model_to_diag, plots = c(1:4)) {
    par(mfrow = c(2,2))
    plot(model_to_diag, which = plots) 
    par(mfrow = c(1,1))
  }
  
    # ben bolkers CI calc fn (https://github.com/bbolker/asaglmm/blob/master/papers/bolker_chap.rmd)
  easyPredCI <- function(model,newdata,alpha=0.05) {
    ## baseline prediction, on the linear predictor (log-odds) scale:
    pred0 <- predict(model,re.form=NA,newdata=newdata)
    ## fixed-effects model matrix for new data
    X <- model.matrix(formula(model,fixed.only=TRUE)[-2],
                      newdata)
    beta <- fixef(model) ## fixed-effects coefficients
    V <- vcov(model)     ## variance-covariance matrix of beta
    pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
    ## inverse-link (logistic) function: could also use plogis()
    plogis <- model@resp$family$plogis
    ## construct 95% Normal CIs on the link scale and
    ##  transform back to the response (probability) scale:
    crit <- -qnorm(alpha/2)
    plogis(cbind(lwr=pred0-crit*pred.se,
                 upr=pred0+crit*pred.se))
  }

# load in data ------------------------------------------------------------
#' # Load Data

#' Load in clp dataset from the p_crispus dataset_prep script:
  pcri.e <- fread("data/output_data/pcri.e_akDNRmv_combined.csv", drop = 1:2)
  pcri.p <- fread("data/output_data/pcri.p_akDNRmv_combined.csv", drop = 1:2)

#' Visualize the distributions of the responses we are interested in:
  #foc: shown with my calcs(gray) and Adam's calcs (red)
  ggplot(pcri.e, aes(foc_clp))+
      geom_histogram(alpha = .5)
  ggplot(pcri.p, aes(foc_clp))+
      geom_histogram(alpha = .5)

  #rake density
  ggplot(pcri.e, aes(rd_pcri))+
      geom_histogram(alpha = .5 )
  ggplot(pcri.p, aes(rd_pcri))+
      geom_histogram(alpha = .5)


# cleaning data -----------------------------------------------------------
#' # Data Cleaning   
#' Are our datasets complete for the parameters we want to estimate on? Note 
#' that the rake densities will be pres for a smaller proportion of the data 
#' (see descriptive stats for these values).
  #check 6 num summaries
  summary(pcri.e[,c("foc_clp","rd_pcri","nsamp","perc.cheme","cytrt","Snow","IceCover","PrevAUGSecchi")])
  summary(pcri.p[,c("foc_clp","rd_pcri","nsamp","perc.lit.chem","cytrt","Snow","IceCover","PrevAUGSecchi")])

  #how many cases?
  pcri.p[,unique(.(Lake, Year))] # 145
  pcri.e[,unique(.(Lake, Year))] # 101
 
  #which lakes?
  pcri.p[,unique(Lake)] # 45
  pcri.e[,unique(Lake)] # 35
  
  #Trt as factor:
  str(pcri.e)
  pcri.e[ , Trt := as.factor(Trt)]
  pcri.e[ , Trte := as.factor(Trte)]
  pcri.p[ , Trt := as.factor(Trt)]
    
  #Lake as a factor:
  pcri.e[ , Lake := as.factor(Lake)]
  pcri.p[ , Lake := as.factor(Lake)]
  
  #convert inches to centimeters in the Snow column:
  pcri.e[,Snow:= Snow*2.54]
  pcri.p[,Snow:= Snow*2.54]
  
  #whats in each dataset?
  names(pcri.e)
  names(pcri.p)
  

# distances to nearest NWS Station ----------------------------------------

nws_dist <- fread(file = "data/raw/SnowDepth/Lakes_NWS.txt")
  lk_dist <- nws_dist[ , .(dowlknum, NEAR_DIST)]
  
  
e_dist <-   merge(pcri.e, lk_dist, all.x = TRUE, by.x = "Lake.ID", by.y = "dowlknum")
p_dist <-   merge(pcri.p, lk_dist, all.x = TRUE, by.x = "Lake.ID", by.y = "dowlknum")

all_dist <- rbind(e_dist, p_dist, fill = TRUE)

distances <- all_dist[, max(NEAR_DIST), .(Lake.ID) ]

summary(distances$V1)

# descriptive stats -----------------------------------------------------------
#' # Summary statistics
  
  # some simple summary stats first two data frames are tables 1 & 2 in the text.
  data.summary <- list(
      "Table_1_dat" = data.frame("data_type" = rep(c("freqoccurrence", "meanrakedens"), each = 3),
                                 "group" = rep(c("total", "treated", "untreated"), 2),
                                 "n lakes" = c(nrow(unique(rbind(pcri.p[,c("Lake","Year")],pcri.e[,c("Lake","Year")])[,1])), # total foc
                                               nrow(unique(subset(rbind(pcri.p[,c("Lake","Year","Trt")],pcri.e[,c("Lake","Year","Trt")]), Trt == 1)[,1])), #treated foc nlakes
                                               nrow(unique(subset(rbind(pcri.p[,c("Lake","Year","Trt")],pcri.e[,c("Lake","Year","Trt")]), Trt == 0)[,1])), #untrt foc nlakes
                                               nrow(unique(rbind(subset(pcri.e, is.na(pcri.e$rd_pcri) == F)[,c("Lake","Trt")], subset(pcri.p, is.na(pcri.p$rd_pcri) == F)[,c("Lake","Trt")])[,1])), # total rake density nlakes
                                               nrow(unique(subset(rbind(subset(pcri.e, is.na(pcri.e$rd_pcri) == F)[,c("Lake","Trt")], subset(pcri.p, is.na(pcri.p$rd_pcri) == F)[,c("Lake","Trt")]), Trt == 1))), # treated rake density nlakes
                                               nrow(unique(subset(rbind(subset(pcri.e, is.na(pcri.e$rd_pcri) == F)[,c("Lake","Trt")], subset(pcri.p, is.na(pcri.p$rd_pcri) == F)[,c("Lake","Trt")]), Trt == 0))) # untrt rake dens nlakes
                                               ),
                                 "n lake-years" = c(nrow(unique(rbind(pcri.p[,c("Lake","Year","Trt")],pcri.e[,c("Lake","Year","Trt")]))[,1]),# foc total n lakeyr
                                                    nrow(unique(subset(rbind(pcri.p[,c("Lake","Year","Trt")],pcri.e[,c("Lake","Year","Trt")]), Trt == 1))[,1]),  # trt foc n lakeyr
                                                    nrow(unique(subset(rbind(pcri.p[,c("Lake","Year","Trt")],pcri.e[,c("Lake","Year","Trt")]), Trt == 0))[,1]), # untrt foc n lakeyr
                                                    nrow(unique(rbind(subset(pcri.e, is.na(pcri.e$rd_pcri) == F)[,c("Lake","Year")], subset(pcri.p, is.na(pcri.p$rd_pcri) == F)[,c("Lake","Year")]))),# rd total n lakeyr
                                                    nrow(unique(rbind(subset(pcri.e, is.na(pcri.e$rd_pcri) == F & pcri.e$Trt == 1)[,c("Lake","Year")], subset(pcri.p, is.na(pcri.p$rd_pcri) == F & pcri.p$Trt == 1)[,c("Lake","Year")]))),# rd trt n lakeyr
                                                    nrow(unique(rbind(subset(pcri.e, is.na(pcri.e$rd_pcri) == F & pcri.e$Trt == 0)[,c("Lake","Year")], subset(pcri.p, is.na(pcri.p$rd_pcri) == F & pcri.p$Trt == 0)[,c("Lake","Year")])))# rd untrt n lakeyr
                                                    ),
                                 "n late surveys" = c(nrow(unique(pcri.p[,c("Lake","Year")])),# foc total nsurveys
                                                      nrow(unique(pcri.p[Trt == 1,c("Lake","Year")])),# foc trt nsurveys
                                                      nrow(unique(pcri.p[Trt == 0,c("Lake","Year")])),# foc untrt nsurveys
                                                      nrow(pcri.p[is.na(pcri.p$rd_pcri) == F,c("Lake","Year")]),# rd total nsurveys
                                                      nrow(pcri.p[is.na(pcri.p$rd_pcri) == F & Trt == 1,c("Lake","Year")]),# rd trt nsurveys
                                                      nrow(pcri.p[is.na(pcri.p$rd_pcri) == F & Trt == 0,c("Lake","Year")])# rd untrt nsurveys
                                                      ),
                                 "n early surveys" = c(nrow(unique(pcri.e[,c("Lake","Year")])[,1]),# foc total nsurveys
                                                       nrow(unique(pcri.e[Trte == 1,c("Lake","Year")])[,1]),# foc trt nsurveys
                                                       nrow(unique(pcri.e[Trte == 0,c("Lake","Year")])[,1]),# foc untrt nsurveys
                                                       nrow(pcri.e[is.na(pcri.e$rd_pcri) == F ,c("Lake","Year")]),# rd total nsurveys
                                                       nrow(pcri.e[is.na(pcri.e$rd_pcri) == F & Trte == 1,c("Lake","Year")]),# rd trt nsurveys
                                                       nrow(pcri.e[is.na(pcri.e$rd_pcri) == F & Trte == 0,c("Lake","Year")]))# rd untrt nsurveys
      ),
      "variables_summary" = data.frame(
        "group" = c(rep(c("trt", "untrt"), each = 2)),
        "season" = c(rep(c("early", "late"),2)),
        "foc_mean" = c(mean(subset(pcri.e, Trte == 1)[,foc_clp]),
                       mean(subset(pcri.p, Trt == 1)[,foc_clp]),
                       mean(subset(pcri.e, Trte == 0)[,foc_clp]),
                       mean(subset(pcri.p, Trt == 0)[,foc_clp])),
        "foc_sd" = c(sd(subset(pcri.e, Trte == 1)[,foc_clp]),
                     sd(subset(pcri.p, Trt == 1)[,foc_clp]),
                     sd(subset(pcri.e, Trte == 0)[,foc_clp]),
                     sd(subset(pcri.p, Trt == 0)[,foc_clp])),
        "rd_mean" = c(mean(subset(pcri.e, is.na(rd_pcri) == F & Trte == 1)[,rd_pcri]),
                      mean(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,rd_pcri]),
                      mean(subset(pcri.e, is.na(rd_pcri) == F & Trte == 0)[,rd_pcri]),
                      mean(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,rd_pcri])),
        "rd_sd" = c(sd(subset(pcri.e, is.na(rd_pcri) == F & Trte == 1)[,rd_pcri]),
                    sd(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,rd_pcri]),
                    sd(subset(pcri.e, is.na(rd_pcri) == F & Trte == 0)[,rd_pcri]),
                    sd(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,rd_pcri])),
        "cytrt_mean" = c(mean(subset(pcri.e, Trte == 1)[,cytrt]),
                         mean(subset(pcri.p, Trt == 1)[,cytrt]),
                         NA,
                         NA),
        "cytrt_sd" = c(sd(subset(pcri.e, Trte == 1)[,cytrt]),
                         sd(subset(pcri.p, Trt == 1)[,cytrt]),
                         NA,
                         NA),
        "Secchi_mean_meter" = c(mean(subset(pcri.e, Trte == 1)[,PrevAUGSecchi]),
                          mean(subset(pcri.p, Trt == 1)[,PrevAUGSecchi]),
                          mean(subset(pcri.e, Trte == 0)[,PrevAUGSecchi]),
                          mean(subset(pcri.p, Trt == 0)[,PrevAUGSecchi])),
        "Secchi_sd_meter" = c(sd(subset(pcri.e, Trte == 1)[,PrevAUGSecchi]),
                        sd(subset(pcri.p, Trt == 1)[,PrevAUGSecchi]),
                        sd(subset(pcri.e, Trte == 0)[,PrevAUGSecchi]),
                        sd(subset(pcri.p, Trt == 0)[,PrevAUGSecchi])),
        "Snow_mean_in" = c(mean(subset(pcri.e, Trte == 1)[,Snow]),
                        mean(subset(pcri.p, Trt == 1)[,Snow]),
                        mean(subset(pcri.e, Trte == 0)[,Snow]),
                        mean(subset(pcri.p, Trt == 0)[,Snow])),
        "Snow_sd_in" = c(sd(subset(pcri.e, Trte == 1)[,Snow]),
                      sd(subset(pcri.p, Trt == 1)[,Snow]),
                      sd(subset(pcri.e, Trte == 0)[,Snow]),
                      sd(subset(pcri.p, Trt == 0)[,Snow])),
        "ice_mean_days" = c(mean(subset(pcri.e, Trte == 1)[,IceCover]),
                       mean(subset(pcri.p, Trt == 1)[,IceCover]),
                       mean(subset(pcri.e, Trte == 0)[,IceCover]),
                       mean(subset(pcri.p, Trt == 0)[,IceCover])),
        "ice_sd_days" = c(sd(subset(pcri.e, Trte == 1)[,IceCover]),
                       sd(subset(pcri.p, Trt == 1)[,IceCover]),
                       sd(subset(pcri.e, Trte == 0)[,IceCover]),
                       sd(subset(pcri.p, Trt == 0)[,IceCover]))
      ),
      "Descriptive_stats" = data.frame(
        "Group" = rep(c("treated", "untreated"), each = 4),
        "Metric" = rep(c("early rd", "peak rd", "early foc", "peak foc"), 2),
        "Mean" = c(
          mean(subset(pcri.e, is.na(rd_pcri) == F & Trte == 1)[,rd_pcri]),
          mean(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,rd_pcri]),
          mean(subset(pcri.e, Trte == 1)[,foc_clp]),
          mean(subset(pcri.p, Trt == 1)[,foc_clp]),
          mean(subset(pcri.e, is.na(rd_pcri) == F & Trte == 0)[,rd_pcri]),
          mean(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,rd_pcri]),
          mean(subset(pcri.e, Trte == 0)[,foc_clp]),
          mean(subset(pcri.p, Trt == 0)[,foc_clp])
        ),
        "stdev" = c(
          sd(subset(pcri.e, is.na(rd_pcri) == F & Trte == 1)[,rd_pcri]),
          sd(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,rd_pcri]),
          sd(subset(pcri.e, Trte == 1)[,foc_clp]),
          sd(subset(pcri.p, Trt == 1)[,foc_clp]),
          sd(subset(pcri.e, is.na(rd_pcri) == F & Trte == 0)[,rd_pcri]),
          sd(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,rd_pcri]),
          sd(subset(pcri.e, Trte == 0)[,foc_clp]),
          sd(subset(pcri.p, Trt == 0)[,foc_clp])
        ),
        "Median" = c(
          median(subset(pcri.e, is.na(rd_pcri) == F & Trte == 1)[,rd_pcri]),
          median(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,rd_pcri]),
          median(subset(pcri.e, Trte == 1)[,foc_clp]),
          median(subset(pcri.p, Trt == 1)[,foc_clp]),
          median(subset(pcri.e, is.na(rd_pcri) == F & Trte == 0)[,rd_pcri]),
          median(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,rd_pcri]),
          median(subset(pcri.e, Trte == 0)[,foc_clp]),
          median(subset(pcri.p, Trt == 0)[,foc_clp])
        ),
        "Quartile1st" = c(
          quantile(subset(pcri.e, is.na(rd_pcri) == F & Trte == 1)[,rd_pcri])[2],
          quantile(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,rd_pcri])[2],
          quantile(subset(pcri.e, Trte == 1)[,foc_clp])[2],
          quantile(subset(pcri.p, Trt == 1)[,foc_clp])[2],
          quantile(subset(pcri.e, is.na(rd_pcri) == F & Trte == 0)[,rd_pcri])[2],
          quantile(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,rd_pcri])[2],
          quantile(subset(pcri.e, Trte == 0)[,foc_clp])[2],
          quantile(subset(pcri.p, Trt == 0)[,foc_clp])[2]
        ),
        "Quartile3rd" = c(
          quantile(subset(pcri.e, is.na(rd_pcri) == F & Trte == 1)[,rd_pcri])[4],
          quantile(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,rd_pcri])[4],
          quantile(subset(pcri.e, Trte == 1)[,foc_clp])[4],
          quantile(subset(pcri.p, Trt == 1)[,foc_clp])[4],
          quantile(subset(pcri.e, is.na(rd_pcri) == F & Trte == 0)[,rd_pcri])[4],
          quantile(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,rd_pcri])[4],
          quantile(subset(pcri.e, Trte == 0)[,foc_clp])[4],
          quantile(subset(pcri.p, Trt == 0)[,foc_clp])[4]
        ),
        "min" = c(
          min(subset(pcri.e, is.na(rd_pcri) == F & Trte == 1)[,rd_pcri]),
          min(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,rd_pcri]),
          min(subset(pcri.e, Trte == 1)[,foc_clp]),
          min(subset(pcri.p, Trt == 1)[,foc_clp]),
          min(subset(pcri.e, is.na(rd_pcri) == F & Trte == 0)[,rd_pcri]),
          min(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,rd_pcri]),
          min(subset(pcri.e, Trte == 0)[,foc_clp]),
          min(subset(pcri.p, Trt == 0)[,foc_clp])
        ),
        "max" = c(
          max(subset(pcri.e, is.na(rd_pcri) == F & Trte == 1)[,rd_pcri]),
          max(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,rd_pcri]),
          max(subset(pcri.e, Trte == 1)[,foc_clp]),
          max(subset(pcri.p, Trt == 1)[,foc_clp]),
          max(subset(pcri.e, is.na(rd_pcri) == F & Trte == 0)[,rd_pcri]),
          max(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,rd_pcri]),
          max(subset(pcri.e, Trte == 0)[,foc_clp]),
          max(subset(pcri.p, Trt == 0)[,foc_clp])
        )
      ),
      "nsamples_per_survey" = list(
        "foc" = data.frame( "grp" = c("trt", "untrt"), 
          "late mean" = c(mean(subset(pcri.p, Trt == 1)[,nsamp]),
                        mean(subset(pcri.p, Trt == 0)[,nsamp])), 
          "late sd" = c(sd(subset(pcri.p, Trt == 1)[,nsamp]),
                        sd(subset(pcri.p, Trt == 0)[,nsamp])),
          "early mean" = c(mean(subset(pcri.e, Trte == 1)[,nsamp]),
                          mean(subset(pcri.e, Trte == 0)[,nsamp])), 
          "early sd" = c(sd(subset(pcri.e, Trt == 1)[,nsamp]),
                        sd(subset(pcri.e, Trt == 0)[,nsamp]))
        ),
        "rd"= data.frame( "grp" = c("trt", "untrt"), 
          "late mean" = c(mean(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,nocc_pcri]),
                          mean(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,nocc_pcri])), 
          "late sd" = c(sd(subset(pcri.p, is.na(rd_pcri) == F & Trt == 1)[,nocc_pcri]),
                        sd(subset(pcri.p, is.na(rd_pcri) == F & Trt == 0)[,nocc_pcri])),
          "early mean" = c(mean(subset(pcri.e, is.na(rd_pcri) == F & Trte == 1)[,nocc_pcri]),
                           mean(subset(pcri.e, is.na(rd_pcri) == F & Trte == 0)[,nocc_pcri])), 
          "early sd" = c(sd(subset(pcri.e, is.na(rd_pcri) == F & Trt == 1)[,nocc_pcri]),
                         sd(subset(pcri.e, is.na(rd_pcri) == F & Trt == 0)[,nocc_pcri]))
        ))
        
    )
  data.summary

#' ## Export Tables 1 & 2
  # export data for Tables in manuscript
  # write.csv(data.summary$variables_summary, file = "data/output_data/summary_stats_Tbl2.csv")
  # write.csv(data.summary$Table_1_dat, file = "data/output_data/data_summary_Tbl1.csv")

# mixed eff models --------------------------------------------------------
#' # Mixed effects models

#' ## No interactions
#' Here we will start our modeling without considering any interactionsusing a 
#' glmme eval. First, we'll need to add an observation level random effect to
#' account for overdispersion (Harrsion, 2015). The following section uses 
#' centered & scaled  predictor variables (done using function scale()) in order
#' to converge during fitting.
#' glmme for frequency of occurrence data:
  # add obs level key to dataset
  pcri.e[,obs := seq(nrow(pcri.e))]
  pcri.p[,obs := seq(nrow(pcri.p))]
  
#' ### Early Survey Freq
  glme.foc.e <- glmer(foc_clp ~ Trte + cytrt + PrevAUGSecchi.sc + Snow.sc + IceCover.sc + (1 | Lake) + (1 | obs), family = binomial(logit), data = pcri.e, weights = nsamp)
  summary(glme.foc.e)

  plot(glme.foc.e)  # not too bad
  plot(glme.foc.e, resid(.)~ fitted(.) |Lake) # by group
  plot(glme.foc.e, resid(.)~ Snow)  # residuals versus Snow
  plot(glme.foc.e, resid(.)~ PrevAUGSecchi) #secchi
  plot(glme.foc.e, resid(.)~ IceCover) # icecover
  plot(glme.foc.e, resid(.)~ cytrt) # cytrt - def more var at low cytrt. That is also where most of our data lie

#' ### Late Survey Freq
  glme.foc.p <- glmer(foc_clp~ Trt + cytrt + PrevAUGSecchi.sc + Snow.sc + IceCover.sc +
                      (1 | Lake)+ (1 | obs), family = binomial(logit), data = pcri.p, weights = nsamp)
  summary(glme.foc.p)
  plot(glme.foc.p) # this also loks o.k. but note the tails again are a bit wonky 
  plot(glme.foc.p, resid(.)~ fitted(.) |Lake) # by group
  plot(glme.foc.p, resid(.)~ Snow)  # residuals versus Snow
  plot(glme.foc.p, resid(.)~ PrevAUGSecchi) #secchi- here we can see that the secchi data have relatively poor coverage past ~3.5m
  plot(glme.foc.p, resid(.)~ IceCover) # icecover
  plot(glme.foc.p, resid(.)~ cytrt) # cytrt

#' ### Early Survey Rake Density:
  glme.rd.e <- glmer(rd_pcri~ Trte + cytrt + PrevAUGSecchi.sc + Snow.sc + IceCover.sc +
                      (1 | Lake), family = Gamma(inverse), data = pcri.e)
  summary(glme.rd.e)
  plot(glme.rd.e) # this looks o.k. but note the lack of data at large fitted vals 

#' ### Late Survey Rake Densities
  glme.rd.p <- glmer(rd_pcri~ Trt + cytrt + PrevAUGSecchi.sc + Snow.sc + IceCover.sc +
                     (1 | Lake), family = Gamma(inverse), data = pcri.p)
  summary(glme.rd.p)
  plot(glme.rd.p) # probably the best so far, but there is a weird lower bound on this one

# mixed models with interactions by trt ------------------------------------------------
#' ## Interactions (trt X env. var.)

#' In order to evaluate whether or not chemical treatements can "harness" the 
#' natural fluctuations of CLP, we need to test for interactions between 
#' treatment and our environmental variables. Our approach is thus to test for
#' interactions, but retain only significant (p = .05) interactions for 
#' estimation of fixed effects.
#' 

#' ### Late Surveys Freq

  #peak foc
  glme.foc.p <- glmer(foc_clp~ Trt + cytrt + Trt*PrevAUGSecchi.sc + Trt*Snow.sc + Trt*IceCover.sc +
                        (1 | Lake)+ (1 | obs), family = binomial(logit), data = pcri.p, weights = nsamp)
  summary(glme.foc.p)

#' Why the warnings? There are a few weird values in there, thus "fit with non-integer data". These are likely the lakes where I used Adam's calculated values. 

  pcri.p[ , nsamp*foc_clp] # note non-integer values

#' Why non-convergence? Explore a few reasons (adapted from https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html):

  # a singularity?
  tt <- getME(glme.foc.p,"theta")
  ll <- getME(glme.foc.p,"lower")
  min(tt[ll==0]) # likely not because of a singularity

  # a gradient calculation error?
  derivs1 <- glme.foc.p@optinfo$derivs
  sc_grad1 <- with(derivs1,solve(Hessian,gradient))
  max(abs(sc_grad1))
  max(pmin(abs(sc_grad1),abs(derivs1$gradient)))
  dd <- update(glme.foc.p,devFunOnly=TRUE)
  pars <- unlist(getME(glme.foc.p,c("theta","fixef")))
    grad2 <- grad(dd,pars)
  hess2 <- hessian(dd,pars)
  sc_grad2 <- solve(hess2,grad2)
  max(pmin(abs(sc_grad2),abs(grad2))) # not lack of gradient (usually set to require > 0.001)
  
  # Restart from previous fit:
  ss <- getME(glme.foc.p,c("theta","fixef"))
  glme.foc.p2 <- update(glme.foc.p,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
  # that worked
  summary(glme.foc.p2)
  
  #then retain significant interactions (PAUGSecchi & Snow):
  glme.foc.p <- glmer(foc_clp~ Trt + cytrt + PrevAUGSecchi.sc+ Snow.sc + IceCover.sc  + Trt:PrevAUGSecchi.sc + Trt:Snow.sc+
                        (1 | Lake)+ (1 | obs), family = binomial(logit), data = pcri.p, weights = nsamp)
  #model.matrix(glme.foc.p)
  summary(glme.foc.p)

#' So we can see that the secchi and snow slopes are different in Trt v Untrt lakes, but but now we want to ask if the values for the Trt interactions are significantly non-zero (vs. sig diff from UNtrt). To do this we'll flip our intercept value, and estimate the slopes for treated lakes first:
  str(pcri.p)
  pcri.p$TrtN <- as.numeric(pcri.p$Trt == 0)
  pcri.p[,.(Trt,TrtN)]
  glme.foc.p.a <- glmer(foc_clp ~ 
                         TrtN + cytrt +
                         + IceCover.sc +
                         TrtN*PrevAUGSecchi.sc + TrtN*Snow.sc+ 
                      (1 | Lake)+ (1 | obs), family = binomial(logit), data = pcri.p, weights = nsamp)
  # model.matrix(glme.foc.p.a)

  # Restart from previous fit:
  ss <- getME(glme.foc.p.a ,c("theta","fixef"))
  glme.foc.p.a2 <- update(glme.foc.p.a,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
  summary(glme.foc.p.a2) # ests for treated slopes as intercepts
  summary(glme.foc.p) # ests for untreated slopes as intercepts
#' So secchi is non-significant in the treated lakes, whereas snow has a 
#' significantly positive effect on clp growth in the treated subset.From these
#' models we can say that the influence on peak clp freq are: trt(-), cytrt (-),
#' Prev Secchi (- in untrt lakes, non-sig in trt lakes), and Snow (non-sig in 
#' untrt lakes, but + in Trt lakes), ice cover (non-sig).
#' 
#' ### Early Surveys Freq

  #early foc
  glme.foc.e <- glmer(foc_clp~ Trte + cytrt + PrevAUGSecchi.sc*Trte + Snow.sc*Trte + IceCover.sc*Trte +
                        (1 | Lake)+ (1 | obs), family = binomial(logit), data = pcri.e, weights = nsamp)
  summary(glme.foc.e)
  
  #no interactions are significant so include none:
  glme.foc.e <- glmer(foc_clp~ Trte + cytrt +  PrevAUGSecchi.sc + Snow.sc + IceCover.sc +
                        (1 | Lake)+ (1 | obs), family = binomial(logit), data = pcri.e, weights = nsamp)
  # model.matrix(glme.foc.e)
  summary(glme.foc.e)

#' Early freq affected by cumulative yrs of treatment (-), productivity (+), and
#' snow cover (-), whereas ice had no effect.
#' 
#' ### Early Rake Density

  glme.rd.e <- glmer(rd_pcri~ Trte + cytrt + PrevAUGSecchi.sc*Trte + Snow.sc*Trte + IceCover.sc*Trte +
                       (1 | Lake), family = Gamma(link = "log"), data = pcri.e)
  #restart at last fit:
  ss <- getME(glme.rd.e,c("theta","fixef"))
  glme.rd.e2 <- update(glme.rd.e,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
  # worked again
  summary(glme.rd.e2)
  
  # retain no interactions
  glme.rd.e <- glmer(rd_pcri~ Trte + cytrt + PrevAUGSecchi.sc + Snow.sc + IceCover.sc + (1 | Lake), family = Gamma(link = "log"), data = pcri.e)
  summary(glme.rd.e)

#'Only Trt last year (-) and productivity (+) are affecting rake densities
#'

#' ### Late Rake Density

  # peak rake density
  glme.rd.p <- glmer(rd_pcri~ Trt + cytrt + PrevAUGSecchi.sc*Trt + Snow.sc*Trt + IceCover.sc*Trt +
                       (1 | Lake), family = Gamma(link = "log"), data = pcri.p)
  # #restart at last fit:
  # ss <- getME(glme.rd.p,c("theta","fixef"))
  # glme.rd.p2 <- update(glme.rd.p,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
  # # worked again
  summary(glme.rd.p)
  
  #retain no interactions:
  glme.rd.p <- glmer(rd_pcri~ Trt + cytrt + PrevAUGSecchi.sc + Snow.sc + IceCover.sc +
                       (1 | Lake), family = Gamma(link = "log"), data = pcri.p)
  summary(glme.rd.p)

#' Only Trt reduced rd, but snow was marginally significant (-)
#' 

#' ### Final models:
  summary(glme.foc.p)
  summary(glme.foc.e)
  summary(glme.rd.e)
  summary(glme.rd.p)

# unscaled models for use in visualization ---------------------------------------
#' Here we will try to re-construct the same models with unscaled data. This
#' will make sure that our visualizations are on the actual scales of the data.
#' 
#' ## Models on unscaled predictors
#' ### Late & Early Survey Freq 
   
  # peak foc
  glme.foc.pu <- glmer(foc_clp~ Trt + cytrt + PrevAUGSecchi+ Snow + IceCover  + Trt:PrevAUGSecchi + Trt:Snow+
                        (1 | Lake)+ (1 | obs), family = binomial(logit), data = pcri.p, weights = nsamp)
  # restart at last fit:
  ss <- getME(glme.foc.pu,c("theta","fixef"))
  glme.foc.pu2 <- update(glme.foc.pu,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
  summary(glme.foc.pu2)# worked again

  # early foc
  glme.foc.eu <- glmer(foc_clp~ Trte + cytrt +  PrevAUGSecchi + Snow + IceCover +
                        (1 | Lake)+ (1 | obs), family = binomial(logit), data = pcri.e, weights = nsamp)
  summary(glme.foc.eu)
  
#' ### Late & Early Rake Densities

  # early rd
  glme.rd.eu <- glmer(rd_pcri~ Trte + cytrt + PrevAUGSecchi + Snow + IceCover + (1 | Lake), family = Gamma(inverse), data = pcri.e)
  # the scaling issue was induced when we unscaled the predictor vars. Because they lie on radically dif scales, it doesn't like them. Try restart at last fit:
  ss <- getME(glme.rd.eu,c("theta","fixef"))
  glme.rd.eu2 <- update(glme.rd.eu,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))# no worky
  #peak rd?
  glme.rd.pu <- glmer(rd_pcri~ Trt + cytrt + PrevAUGSecchi + Snow + IceCover +
                       (1 | Lake), family = Gamma(inverse), data = pcri.p)
  # Again, scaling issue induced when we unscaled the predictors. Try restart at last fit:
  ss <- getME(glme.rd.pu,c("theta","fixef"))
  glme.rd.pu2 <- update(glme.rd.pu,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))# no worky

#' We have actual scale models for frequency, but not rake density. This is
#' likley because the gamma models do not do well with the variation among the
#' predictor variables. 

# visualizations  ------------------------------------

#' # Visualizing results


# secchi ----------------------------------------------------
#' ## Secchi

  # Early Surveys
  #Prev Secchi
  length(pcri.e$cytrt)
  newdat1<-data.frame(cytrt =rep(0,length(pcri.e$cytrt)),
                      Trte = factor(x = rep(1, length(pcri.e$cytrt)), levels = c(0,1)),
                      Snow = rep(mean(pcri.e$Snow),length(pcri.e$Snow)), 
                      PrevAUGSecchi = seq(min(pcri.e$PrevAUGSecchi), max(pcri.e$PrevAUGSecchi), length.out = length(pcri.e$cytrt)),
                      IceCover = rep(mean(pcri.e$IceCover),length(pcri.e$cytrt)),
                      Lake = pcri.e$Lake,
                      obs = pcri.e$obs)
  #fitted values
  phat1<-predict(glme.foc.eu, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.eu,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  Trt.pred <- newdat1
  newdat1<-data.frame(cytrt =rep(0,length(pcri.e$cytrt)),
                      Trte = factor(x = rep(0, length(pcri.e$cytrt)), levels = c(0,1)),
                      Snow = rep(mean(pcri.e$Snow),length(pcri.e$Snow)), 
                      PrevAUGSecchi = seq(min(pcri.e$PrevAUGSecchi), max(pcri.e$PrevAUGSecchi), length.out = length(pcri.e$cytrt)),
                      IceCover = rep(mean(pcri.e$IceCover),length(pcri.e$cytrt)),
                      Lake = pcri.e$Lake,
                      obs = pcri.e$obs)
  #fitted values
  phat1<-predict(glme.foc.eu, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.eu,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  UnTrt.pred <- newdat1
  
  predtype <- list("Untreated", "Treated")
  pred <- rbindlist(list("Untrt" = UnTrt.pred,"Trt" = Trt.pred), idcol = TRUE)
  
  secchi_early <- ggplot(pred, aes(PrevAUGSecchi,fitted))+
    geom_line(size = 2, aes(linetype = .id))+ 
    geom_ribbon(aes(ymin=low, ymax=up, group = .id), color = NA , fill = "black", alpha = .15)+
    ylab("Frequency of Occurrence")+
    xlab("Previous August Secchi (m)")  +
    scale_linetype_manual("Previous-year", values = c("Trt" = "solid", "Untrt" = "dashed" ), labels = c("Treated","Untreated"))+
    theme(text = element_text(size=20), legend.position = c(.60, .85), legend.key.width = unit(1.25, "in"), legend.background = element_blank())+
    ggtitle("Early Survey Period")+
    theme(plot.title = element_text(hjust = 0.5))
  
  secchi_early 
  

  # mean change in foc per meter of secchi
  lm(fitted~PrevAUGSecchi, data = rbind(UnTrt.pred,Trt.pred))
  
  #' ON average, the early season freq declines 9.2% for each added meter of
#' Secchi clarity. 
  

  # LATE
  #note that we cut off our predicitons at 3m Secchi because beyond that we don't 
  # have enough data to inform about trend/pattern
  #Prev Secchi
  length(pcri.p$cytrt)
  #trtd
  newdat1<-data.frame(cytrt =rep(0,length(pcri.p$cytrt)),
                      Trt = factor(x = rep(1, length(pcri.p$cytrt)), levels = c(0,1)),
                      Snow = rep(mean(pcri.p$Snow),length(pcri.p$Snow)), 
                      PrevAUGSecchi = seq(min(pcri.p$PrevAUGSecchi), 3, length.out = length(pcri.p$cytrt)),
                      IceCover = rep(mean(pcri.p$IceCover),length(pcri.p$cytrt)),
                      Lake = pcri.p$Lake,
                      obs = pcri.p$obs)
  #fitted values
  phat1<-predict(glme.foc.pu2, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.pu2,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  Trt.pred <- newdat1
  #untrtd
  newdat1<-data.frame(cytrt =rep(0,length(pcri.p$cytrt)),
                      Trt = factor(x = rep(0, length(pcri.p$cytrt)), levels = c(0,1)),
                      Snow = rep(mean(pcri.p$Snow),length(pcri.p$Snow)), 
                      PrevAUGSecchi = seq(min(pcri.p$PrevAUGSecchi), 3, length.out = length(pcri.p$cytrt)),
                      IceCover = rep(mean(pcri.p$IceCover),length(pcri.p$cytrt)),
                      Lake = pcri.p$Lake,
                      obs = pcri.p$obs)
  #fitted values
  phat1<-predict(glme.foc.pu2, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.pu2,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  UnTrt.pred <- newdat1

  predtype <- list("Untreated", "Treated")
  pred <- rbindlist(list("Untrt" = UnTrt.pred,"Trt" = Trt.pred), idcol = TRUE)
  
  secchi_late <- ggplot(pred, aes(PrevAUGSecchi,fitted))+
    geom_line(size = 2, aes(linetype = .id))+ 
    geom_ribbon(aes(ymin=low, ymax=up, group = .id), color = NA , fill = "black", alpha = .15)+
    ylab("Frequency of Occurrence")+
    xlab("Previous August Secchi (m)")  +
    scale_linetype_manual("Within-year", values = c("Trt" = "solid", "Untrt" = "dashed" ), labels = c("Treated","Untreated"))+
    theme(text = element_text(size=20), legend.position = c(.60, .85), legend.key.width = unit(1.25, "in"), legend.background = element_blank())+
    ggtitle("Late Survey Period")+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  secchi_late
  
  # secchi_late+
  # geom_point(aes(x = pcri.p$PrevAUGSecchi, y = pcri.p$foc_clp, shape = pcri.p$Trt ), size = 5, alpha = .4)+
  # geom_line(aes(x = pcri.p$PrevAUGSecchi, y=pcri.p$foc_clp ,group = pcri.p$Lake), alpha = .3)

  
  # mean change in foc per meter of secchi
  lm(fitted~PrevAUGSecchi, data = rbind(UnTrt.pred,Trt.pred))
#' Again, a 9.2% decline in frequency for each additonal meter of Secchi 

  # secchi_late+
  # geom_point(aes(x = pcri.p$PrevAUGSecchi, y = pcri.p$foc_clp, shape = pcri.p$Trt ), size = 5, alpha = .4)+
  # geom_line(aes(x = pcri.p$PrevAUGSecchi, y=pcri.p$foc_clp ,group = pcri.p$Lake), alpha = .3)

  
  # mean change in foc per meter of secchi
  lm(fitted~PrevAUGSecchi, data = UnTrt.pred)
#' Again, a 9.4% decline in frequency for each additonal meter of Secchi 

#' clarity.

# snow ------------------------------------------------------

#' ## Snow

  # EARLY Surveys
  length(pcri.e$cytrt)
  newdat1<-data.frame(cytrt =rep(0,length(pcri.e$cytrt)),
                      Trte = factor(x = rep(1, length(pcri.e$cytrt)), levels = c(0,1)),
                      Snow = seq(min(pcri.e$Snow), max(pcri.e$Snow), length.out = length(pcri.e$cytrt)), 
                      PrevAUGSecchi = rep(mean(pcri.e$PrevAUGSecchi),length(pcri.e$cytrt)),
                      IceCover = rep(mean(pcri.e$IceCover),length(pcri.e$cytrt)),
                      Lake = pcri.e$Lake,
                      obs = pcri.e$obs)
  #fitted values
  phat1<-predict(glme.foc.eu, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.eu,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  Trt.pred <- newdat1
  newdat1<-data.frame(cytrt =rep(0,length(pcri.e$cytrt)),
                      Trte = factor(x = rep(0, length(pcri.e$cytrt)), levels = c(0,1)),
                      Snow = seq(min(pcri.e$Snow), max(pcri.e$Snow), length.out = length(pcri.e$cytrt)), 
                      PrevAUGSecchi = rep(mean(pcri.e$PrevAUGSecchi),length(pcri.e$cytrt)),
                      IceCover = rep(mean(pcri.e$IceCover),length(pcri.e$cytrt)),
                      Lake = pcri.e$Lake,
                      obs = pcri.e$obs)
  #fitted values
  phat1<-predict(glme.foc.eu, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.eu,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  UnTrt.pred <- newdat1

  head(Trt.pred)
  tail(Trt.pred)
  
  snow_early <- ggplot(Trt.pred, aes(Snow,fitted))+
    geom_line( size = 2)+ 
    geom_ribbon(aes(ymin=low, ymax=up), alpha=0.2)+
    geom_line(data = UnTrt.pred,aes(Snow,fitted),  size = 2, linetype = 2)+
    geom_ribbon(data = UnTrt.pred,aes(ymin=low, ymax=up), alpha=0.2)+
    theme(text = element_text(size=20), legend.position = "none")+
    ylab("Frequency of Occurrence")+
    xlab("Mean Winter Snow Cover (cm)")
    
    
  snow_early
  
  # show data and average value line
  snow_early+
    geom_line(aes(x = pcri.e$Snow, y=pcri.e$foc_clp ,group = pcri.e$Lake), alpha = .3)+
    geom_point(aes(x = pcri.e$Snow, y = pcri.e$foc_clp, shape = pcri.e$Trte ), size = 5, alpha = .4)+
    geom_smooth( method = "lm", color = "red")
  
  # mean change in foc per centimeter of snow
  bound <- bind_rows(UnTrt.pred,Trt.pred)
  lm(fitted~Snow, data = bound)
  lm(fitted~Snow, data = bound)$'coefficients'[2]*2.54 #convert to inches

  #' This comes out to approximately a 1.1% change in foc per addt'l inch of snow

  #LATE Surveys
  length(pcri.p$cytrt)
  newdat1<-data.frame(cytrt =rep(0,length(pcri.p$cytrt)),
                      Trt = factor(x = rep(1, length(pcri.p$cytrt)), levels = c(0,1)),
                      Snow = seq(min(pcri.p$Snow), max(pcri.p$Snow), length.out = length(pcri.p$cytrt)), 
                      PrevAUGSecchi = rep(mean(pcri.p$PrevAUGSecchi),length(pcri.p$cytrt)),
                      IceCover = rep(mean(pcri.p$IceCover),length(pcri.p$cytrt)),
                      Lake = pcri.p$Lake,
                      obs = pcri.p$obs)
  #fitted values
  phat1<-predict(glme.foc.pu2, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.pu2,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  Trt.pred <- newdat1
  newdat1<-data.frame(cytrt =rep(0,length(pcri.p$cytrt)),
                      Trt = factor(x = rep(0, length(pcri.p$cytrt)), levels = c(0,1)),
                      Snow = seq(min(pcri.p$Snow), max(pcri.p$Snow), length.out = length(pcri.p$cytrt)), 
                      PrevAUGSecchi = rep(mean(pcri.p$PrevAUGSecchi),length(pcri.p$cytrt)),
                      IceCover = rep(mean(pcri.p$IceCover),length(pcri.p$cytrt)),
                      Lake = pcri.p$Lake,
                      obs = pcri.p$obs)
  #fitted values
  phat1<-predict(glme.foc.pu2, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.pu2,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  UnTrt.pred <- newdat1
  
  head(Trt.pred)
  tail(Trt.pred)
  
  snow_late <- ggplot(Trt.pred, aes(Snow,fitted))+
    geom_line(size = 2)+ 
    geom_ribbon(aes(ymin=low, ymax=up), alpha=0.2)+
    geom_line(data = UnTrt.pred,aes(Snow,fitted), size = 2, linetype = 2)+
    geom_ribbon(data = UnTrt.pred,aes(ymin=low, ymax=up), alpha=0.2)+
    theme(text = element_text(size=20), legend.position = "none")+
    ylab("Frequency of Occurrence")+
    xlab("Mean Winter Snow Cover (cm)")
  snow_late
  
  # uncomment these lines to add data in background:
  snow_late+  
  geom_point(aes(x = pcri.p$Snow, y = pcri.p$foc_clp, shape = pcri.p$Trt ), size = 5, alpha = .4)+
    geom_line(aes(x = pcri.p$Snow, y=pcri.p$foc_clp ,group = pcri.p$Lake), alpha = .1)

  # mean change in foc per inch of snow
  bound <- bind_rows(Trt.pred)
  lm(fitted~Snow, data = bound)
  lm(fitted~Snow, data = bound)$'coefficients'[2]*2.54 #convert to inches
  # summary(trt.snow.pred)

#' Approximately a 0.7% increase in foc per addt'l inch of snow in the treated lakes


# ice -------------------------------------------------------
#' ## Ice
#' 

  # Early Surveys
  length(pcri.e$cytrt)
  newdat1<-data.frame(cytrt =rep(0,length(pcri.e$cytrt)),
                      Trte = factor(x = rep(1, length(pcri.e$cytrt)), levels = c(0,1)),
                      Snow = rep(mean(pcri.e$Snow),length(pcri.e$Snow)), 
                      PrevAUGSecchi = rep(mean(pcri.e$PrevAUGSecchi),length(pcri.e$cytrt)),
                      IceCover = seq(min(pcri.e$IceCover), max(pcri.e$IceCover), length.out = length(pcri.e$cytrt)),
                      Lake = pcri.e$Lake,
                      obs = pcri.e$obs)
  #fitted values
  phat1<-predict(glme.foc.eu, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.eu,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  Trt.pred <- newdat1
  newdat1<-data.frame(cytrt =rep(0,length(pcri.e$cytrt)),
                      Trte = factor(x = rep(0, length(pcri.e$cytrt)), levels = c(0,1)),
                      Snow = rep(mean(pcri.e$Snow),length(pcri.e$Snow)), 
                      PrevAUGSecchi = rep(mean(pcri.e$PrevAUGSecchi),length(pcri.e$cytrt)),
                      IceCover = seq(min(pcri.e$IceCover), max(pcri.e$IceCover), length.out = length(pcri.e$cytrt)),
                      obs = pcri.e$obs)
  #fitted values
  phat1<-predict(glme.foc.eu, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.eu,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  UnTrt.pred <- newdat1
  
  ice_early <- ggplot(Trt.pred, aes(IceCover,fitted))+
    geom_line( size = 2)+ 
    geom_ribbon(aes(ymin=low, ymax=up), alpha=0.2)+
    geom_line(data = UnTrt.pred,aes(IceCover, fitted), size = 2, linetype = 2)+
    geom_ribbon(data = UnTrt.pred,aes(ymin=low, ymax=up), alpha=0.2)+
    theme(text = element_text(size=20), legend.position = "none")+
    ylab("Frequency of Occurrence")+
    xlab("Winter Ice Cover (days)")  +
    scale_color_manual(breaks = c(0,1), values = c("blue", "red"))
  ice_early
  
  
  # plots data in background:
    ice_early+ 
      geom_point(aes(x = pcri.e$IceCover, y = pcri.e$foc_clp, shape = pcri.e$Trte ), size = 5, alpha = .4)+
      geom_line(aes(x = pcri.e$IceCover, y=pcri.e$foc_clp ,group = pcri.e$Lake), alpha = .3)

  # Late Surveys
  length(pcri.p$cytrt)
  newdat1<-data.frame(cytrt =rep(0,length(pcri.p$cytrt)),
                      Trt = factor(x = rep(1, length(pcri.p$cytrt)), levels = c(0,1)),
                      Snow = rep(mean(pcri.p$Snow),length(pcri.p$Snow)), 
                      PrevAUGSecchi = rep(mean(pcri.p$PrevAUGSecchi),length(pcri.p$cytrt)),
                      IceCover = seq(min(pcri.p$IceCover), max(pcri.p$IceCover), length.out = length(pcri.p$cytrt)),
                      Lake = pcri.p$Lake,
                      obs = pcri.p$obs)
  #fitted values
  phat1<-predict(glme.foc.pu2, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.pu2,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  Trt.pred <- newdat1
  newdat1<-data.frame(cytrt =rep(0,length(pcri.p$cytrt)),
                      Trt = factor(x = rep(0, length(pcri.p$cytrt)), levels = c(0,1)),
                      Snow = rep(mean(pcri.p$Snow),length(pcri.p$Snow)), 
                      PrevAUGSecchi = rep(mean(pcri.p$PrevAUGSecchi),length(pcri.p$cytrt)),
                      IceCover = seq(min(pcri.p$IceCover), max(pcri.p$IceCover), length.out = length(pcri.p$cytrt)),
                      Lake = pcri.p$Lake,
                      obs = pcri.p$obs)
  #fitted values
  phat1<-predict(glme.foc.pu2, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.pu2,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  UnTrt.pred <- newdat1

  ice_late <- ggplot(Trt.pred, aes(IceCover,fitted))+
    geom_line(size = 2)+ 
    geom_ribbon(aes(ymin=low, ymax=up), alpha=0.2)+
    geom_line(data = UnTrt.pred,aes(IceCover, fitted),size = 2, linetype = 2 )+
    geom_ribbon(data = UnTrt.pred,aes(ymin=low, ymax=up), alpha=0.2)+
    theme(text = element_text(size=20), legend.position = "none")+
    ylab("Frequency of Occurrence")+
    xlab("Winter Ice Cover (days)")
  ice_late
  
  ice_late+
    geom_point(aes(x = pcri.p$IceCover, y = pcri.p$foc_clp, shape = pcri.p$Trt ), size = 5, alpha = .4)+
    geom_line(aes(x = pcri.p$IceCover, y=pcri.p$foc_clp ,group = pcri.p$Lake), alpha = .3)+
    scale_color_manual(breaks = c(0,1), values = c("blue", "red"))

#' these are much nicer, but they do not explain how the effect that is observed is after accounting for other drivers (background data are "raw")


# cytrt -------------------------------------------------------------------
#' ## Consecutive Years of Treatment

  # EARLY Surveys
  length(pcri.e$cytrt)
  newdat1<-data.frame(cytrt = seq(min(pcri.e$cytrt), 
                                  max(pcri.e$cytrt),  length.out = length(pcri.e$cytrt)),
                      Trte = factor(x = rep(1, length(pcri.e$cytrt)), levels = c(0,1)),
                      Snow = rep(mean(pcri.e$Snow),length(pcri.e$cytrt)), 
                      PrevAUGSecchi = rep(mean(pcri.e$PrevAUGSecchi),length(pcri.e$cytrt)),
                      IceCover = rep(mean(pcri.e$IceCover),length(pcri.e$cytrt)),
                      Lake = pcri.e$Lake,
                      obs = pcri.e$obs)
  #fitted values
  phat1<-predict(glme.foc.eu, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.eu,newdata = newdat1) # see fn documentation in load fns section
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  Trt.pred <- newdat1

#' The line we'll show is smoothed, but here are the declines for each
#' added, consecutive year of treatment:
  Trt.pred[Trt.pred$cytrt == 0,"fitted"] - Trt.pred[Trt.pred$cytrt == 1,"fitted"] # This means going from a single year of treating to two years (1 addt'l consecutive year)
  Trt.pred[Trt.pred$cytrt == 1,"fitted"] - Trt.pred[Trt.pred$cytrt == 2,"fitted"]
  Trt.pred[Trt.pred$cytrt == 2,"fitted"] - Trt.pred[Trt.pred$cytrt == 3,"fitted"]
  Trt.pred[Trt.pred$cytrt == 3,"fitted"] - Trt.pred[Trt.pred$cytrt == 4,"fitted"]
  Trt.pred[Trt.pred$cytrt == 4,"fitted"] - Trt.pred[Trt.pred$cytrt == 5,"fitted"]

#' If we strike a line through the model we find mean change in freq per
#' additonal year of treatment. 
  
  bound <- bind_rows(Trt.pred)
  lm(fitted~cytrt, data = bound)
  
  # note that I have shifted the axis labels so that the first year of trt is not zero
cytrt_early <-   ggplot(Trt.pred, aes(cytrt,fitted))+
    geom_line(size = 1.5)+ 
    geom_ribbon(aes(ymin=low, ymax=up), alpha=0.4)+
    theme(text = element_text(size=20), legend.position = "none")+
    ylab("Frequency of Occurrence")+
    xlab("Consecutive Years Treated")  +
    scale_x_continuous( breaks = seq(0,5, by = 1),
                        labels = seq(1,6, by = 1),
                        limits = c(0,5))
  cytrt_early  

  # plot data in backgroud
cytrt_early + 
    geom_point(x = pcri.e$cytrt, y = pcri.e$foc_clp)+
    geom_smooth(method = "lm")
  
  # LATE Surveys
  length(pcri.p$cytrt)
  newdat1<-data.frame(cytrt = seq(min(pcri.p$cytrt), 
                                  max(pcri.p$cytrt),  length.out = length(pcri.p$cytrt)),
                      Trt = factor(x = rep(1, length(pcri.p$cytrt)), levels = c(0,1)),
                      Snow = rep(mean(pcri.p$Snow),length(pcri.p$cytrt)), 
                      PrevAUGSecchi = rep(mean(pcri.p$PrevAUGSecchi),length(pcri.p$cytrt)),
                      IceCover = rep(mean(pcri.p$IceCover),length(pcri.p$cytrt)),
                      Lake = pcri.p$Lake,
                      obs = pcri.p$obs)
  #fitted values
  phat1<-predict(glme.foc.pu2, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.pu2,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  Trt.pred <- newdat1
  
  # declines for each added, consecutive year of treatment:
  Trt.pred[Trt.pred$cytrt == 0,"fitted"] - Trt.pred[Trt.pred$cytrt == 1,"fitted"] # This means going from a single year of treating to two years (1 addt'l consecutive year)
  Trt.pred[Trt.pred$cytrt == 1,"fitted"] - Trt.pred[Trt.pred$cytrt == 2,"fitted"]
  Trt.pred[Trt.pred$cytrt == 2,"fitted"] - Trt.pred[Trt.pred$cytrt == 3,"fitted"]
  Trt.pred[Trt.pred$cytrt == 3,"fitted"] - Trt.pred[Trt.pred$cytrt == 4,"fitted"]
  Trt.pred[Trt.pred$cytrt == 4,"fitted"] - Trt.pred[Trt.pred$cytrt == 5,"fitted"]
  Trt.pred[Trt.pred$cytrt == 5,"fitted"] - Trt.pred[Trt.pred$cytrt == 6,"fitted"]
  Trt.pred[Trt.pred$cytrt == 6,"fitted"] - Trt.pred[Trt.pred$cytrt == 7,"fitted"]
  Trt.pred[Trt.pred$cytrt == 7,"fitted"] - Trt.pred[Trt.pred$cytrt == 8,"fitted"]
  
  
  # avg change per year modeled
  bound <- bind_rows(Trt.pred)
  lm(fitted~cytrt, data = bound)
  
  
  cytrt_late <- ggplot(Trt.pred, aes(cytrt,fitted))+
    geom_line( size = 1.5)+ 
    geom_ribbon(aes(ymin=low, ymax=up), alpha=0.4)+
    theme(text = element_text(size=20), legend.position = "none")+
    ylab("Frequency of Occurrence")+
    xlab("Consecutive Years Treated")  +
    scale_x_continuous( breaks = seq(0,8, by = 1),
                        labels = seq(1,9, by = 1),
                        limits = c(0,9))
  cytrt_late
  
  # mean change in foc per yr trt
cytrt_late+
    geom_point(x = pcri.p$cytrt, y = pcri.p$foc_clp, alpha = .5)+
    geom_smooth(method = "lm")


# treatment effects -------------------------------------------------------

#' ## Treatments
#' 
    #Trte switchover for EARLY freq:
  length(pcri.e$cytrt)
  newdat1<-data.frame(Trte = factor(x = c(rep(0,50),rep(1,51)), levels = c(0,1)),
                      IceCover = rep(mean(pcri.e$IceCover),length(pcri.e$IceCover)),
                      Snow = rep(mean(pcri.e$Snow),length(pcri.e$Snow)), 
                      cytrt = rep(0,length(pcri.e$cytrt)),
                      PrevAUGSecchi = rep(mean(pcri.e$PrevAUGSecchi),length(pcri.e$PrevAUGSecchi)))
  
  phat1<-predict(glme.foc.eu, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.eu,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  UnTrt.pred <- newdat1
  
  modvals = data.frame(Trte= c(0,1), modfoc=c(max(UnTrt.pred$fitted),min(UnTrt.pred$fitted)), lwr= c(.75, 1.75), upr = c(1.25, 2.25) )
  

  # ggplot(pcri.e, aes(x = pcri.e$foc_clp))+
  #   geom_density(aes(group = Trte, lty = Trte))+
  #   geom_segment(data = modvals,aes(y = 0, yend = 2, x = modfoc, xend = modfoc, lty = as.factor(Trte)))+
  #   theme(text = element_text(size=20), legend.position = "none")+
  #   ylab("Density")+xlab("Frequency of Occurrence")+
  #   scale_fill_manual( values = c("blue", "red"))+
  #   theme_bw()
  # 

  trt_foc_early <- ggplot(pcri.e, aes(x=foc_clp)) +
    geom_density( aes(linetype = Trte), alpha=0.25, position="identity", lwd = 1.5)+
    scale_linetype_manual("Previous-year", values = c("0" = "dotted", "1" = "solid"), labels = c("Untreated", "Treated"))+
    geom_vline(xintercept = modvals$modfoc, size = 1.5, lty = c("0" = "dotted", "1" = "solid"))+
    theme(legend.position = c(.75, .75), legend.background = element_blank())+
    ylab("Kernel Density")+xlab("Frequency of Occurrence")+
    geom_text( aes(.42,.75, label = "a"), size = 7)+
    geom_text( aes(.30,.75, label = "a"), size = 7)+
    theme(text = element_text(size=20))+
    ggtitle("Early Survey Period")+
    theme(plot.title = element_text(hjust = 0.5))
  
  trt_foc_early
  
  
#Trt switchover for Late freq:
  length(pcri.p$cytrt)
  newdat1<-data.frame(Trt = factor(x = c(rep(0,72),rep(1,73)), levels = c(0,1)),
                      IceCover = rep(mean(pcri.p$IceCover),length(pcri.p$IceCover)),
                      Snow = rep(mean(pcri.p$Snow),length(pcri.p$Snow)), 
                      cytrt = rep(0,length(pcri.p$cytrt)),
                      PrevAUGSecchi = rep(mean(pcri.p$PrevAUGSecchi),length(pcri.p$PrevAUGSecchi)))
  
  phat1<-predict(glme.foc.pu, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-plogis(phat1)
  
  cpred1.CI <- easyPredCI(glme.foc.pu,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  UnTrt.pred <- newdat1
  
  modvals = data.frame(Trt= c(0,1), modfoc=c(max(UnTrt.pred$fitted),min(UnTrt.pred$fitted)), lwr= c(.75, 1.75), upr = c(1.25, 2.25) )

  # ggplot(pcri.p, aes(x = pcri.p$foc_clp))+
  #   geom_density(aes(group = Trt, lty = Trt))+
  #   theme(text = element_text(size=20), legend.position = "none")+
  #   ylab("Density")+xlab("Frequency of Occurrence")+
  #   geom_segment(data = modvals,aes(y = 0, yend = 3, x = modfoc, xend = modfoc, lty = as.factor(Trt)))+
  #   theme_bw()

  
  trt_foc_late <- ggplot(pcri.p, aes(x=foc_clp)) +
    geom_density(aes(linetype = Trt), alpha=0.25, position="identity", size = 1.5)+
    theme(text = element_text(size=20))+
    scale_linetype_manual("Within-year", values = c("0" = "dotted", "1" = "solid"), labels = c("Untreated","Treated"))+
    geom_vline(xintercept = modvals$modfoc, size = 1.5,  lty = c("0" = "dotted", "1" = "solid"))+
    ylab("Kernel Density")+xlab("Frequency of Occurrence")+
    geom_text( aes(.07,.75, label = "a"), size = 7)+
    geom_text( aes(.30,.75, label = "b"), size = 7)+
    theme(legend.position = c(.75, .75), legend.background = element_blank())+
    ggtitle("Late Survey Period")+
    theme(plot.title = element_text(hjust = 0.5))
    
  
  trt_foc_late

  # Trt switchover for early rd:
  
  length(pcri.e$cytrt)
  newdat1<-data.frame(Trte = factor(x = c(rep(0,50),rep(1,51)), levels = c(0,1)),
                      IceCover.sc = rep(mean(pcri.e$IceCover.sc),length(pcri.e$IceCover.sc)),
                      Snow.sc = rep(mean(pcri.e$Snow.sc),length(pcri.e$Snow.sc)), 
                      cytrt = rep(0,length(pcri.e$cytrt)),
                      PrevAUGSecchi.sc = rep(mean(pcri.e$PrevAUGSecchi.sc),length(pcri.e$PrevAUGSecchi.sc)))
  
  phat1<-predict(glme.rd.e, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-exp(phat1)
  
  cpred1.CI <- easyPredCI(glme.rd.e,newdata = newdat1)
  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  UnTrt.pred <- newdat1
  
  modvals = data.frame(Trte= c(0,1), modfoc=c(max(UnTrt.pred$fitted),min(UnTrt.pred$fitted)), lwr= c(.75, 1.75), upr = c(1.25, 2.25) )
  
  trt_rd_early <- ggplot(pcri.e, aes(x=rd_pcri)) +
    geom_density(aes(linetype=Trte), alpha=0.25, position="identity", size = 1.5)+
    xlim(c(0,4))+
    theme(text = element_text(size=20))+
    ylab("Kernel Density")+xlab("Mean Rake Density")+
    scale_linetype_manual("Previous-year", values = c("0" = "dotted", "1" = "solid"), labels = c("Untreated", "Treated"))+
    geom_vline(xintercept = modvals$modfoc, size = 1.5, lty = c("0" = "dotted", "1" = "solid"))+
    geom_text( aes(1.25,.75, label = "a"), size = 7)+
    geom_text( aes(1.86,.75, label = "b"), size = 7)+
    theme(legend.position = c(.75, .75), legend.background = element_blank())
  
  trt_rd_early
  
  #trt switchover for peak rd
  
  length(pcri.p$cytrt)
  newdat1<-data.frame(Trt = factor(x = c(rep(0,72),rep(1,73)), levels = c(0,1)),
                      IceCover.sc = rep(mean(pcri.p$IceCover.sc),length(pcri.p$IceCover.sc)),

                      Snow.sc = rep(mean(pcri.p$Snow.sc),length(pcri.p$Snow.sc)), 
                      cytrt = rep(0,length(pcri.p$cytrt)),
                      PrevAUGSecchi.sc = rep(mean(pcri.p$PrevAUGSecchi.sc),length(pcri.p$PrevAUGSecchi.sc)))
  
  phat1<-predict(glme.rd.p, newdat=newdat1, re.form=NA, type="link")
  newdat1$fitted<-exp(phat1)
  
  cpred1.CI <- easyPredCI(glme.rd.p,newdata = newdat1)

  newdat1$up<-cpred1.CI[,2]
  newdat1$low<-cpred1.CI[,1]
  UnTrt.pred <- newdat1
  
  modvals = data.frame(Trt= c(0,1), modfoc=c(max(UnTrt.pred$fitted),min(UnTrt.pred$fitted)), lwr= c(.75, 1.75), upr = c(1.25, 2.25) )

  trt_rd_late <- ggplot(pcri.p, aes(x=rd_pcri)) +
    geom_density(aes( linetype = Trt), alpha=0.25, position="identity", size = 1.5)+
    xlim(c(0,4))+
    theme(text = element_text(size=20))+
    ylab("Kernel Density")+xlab("Mean Rake Density")+
    scale_linetype_manual("Within-year", values = c("0" = "dotted", "1" = "solid"), labels = c("Untreated", "Treated"))+
    geom_vline(xintercept = modvals$modfoc, size = 1.5, lty = c("0" = "dotted", "1" = "solid"))+
    geom_text( aes(1.1,.75, label = "a"), size = 7)+
    geom_text( aes(2.05,.75, label = "b"), size = 7)+
    theme(legend.position = c(.75, .75), legend.background = element_blank())
  
  trt_rd_late


# area treated declines -----------------------------------------------
#' # Area Managed
#' Evidence for decreased intervention thru time?

# does acreage decline with increased cum yrs trt (results in hectares)

ggplot(subset(pcri.p, Trt == 1), aes(jitter(cytrt), Acres_Chem*0.404686))+
  geom_point(size = 4, alpha = .2)+
  geom_smooth(method = "lm", color = "black")+
  theme_bw()

area_trt <- summary(lmer(Acres_Chem ~ cytrt + ( 1 | Lake ), data = subset(pcri.p, Trt == 1)))
summary(area_trt)

#convert to ha
area_trt$coefficients*0.404686

#' so the acreage treated declines by 5.5 hectares +- 1.96*2.23 for each additonal year of treatment
#' and avearge lakes treated 45.1 ha in first treatment


# treating based on spring populations? -----------------------------------

#' # Bias in Treated Late Surveys
#' Do you get more treatment in years with high early season foc?

ggplot(subset(pcri.e, Trte == 1), aes(foc_clp,perc.lit.chem))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(subset(pcri.e, Trte == 1), aes(Snow,perc.lit.chem))+
  geom_point()+
  geom_smooth(method = "lm")

summary(lm(data = subset(pcri.e, Trte == 1), perc.lit.chem ~ foc_clp))
summary(lm(data = subset(pcri.e, Trte == 1), perc.lit.chem ~ Snow))


plot(pcri.p$perc.lit.chem ~ pcri.p$foc_clp)

plot(pcri.p$foc_clp~ pcri.p$perc.lit.chem)

ggplot(pcri.p, aes(cytrt, Acres_Chem))+
  geom_point()+
  geom_path(aes(group = Lake))

# export centered scaled models for reporting -----------------------------

#' # Review and export modeled estimates
#' Reported parameter estimates should be on a centered scaled model, thus allowing intuitive inference of sign.

summary(glme.foc.p)# June_FOC_split_secchi(notsigdifzero)&snow(nonzero)   (centered scaled)
summary(glme.foc.e)# April_FOC_split_none (centered scaled)
summary(glme.rd.e)# April_DENS_split_none (centered scaled)
summary(glme.rd.p)# June_DENS_split_none (centered scaled)

summary(glme.foc.p)[["coefficients"]]

#' ## Build and export Table 3

early_mod_foc <- summary(glme.foc.e)[["coefficients"]]
peak_mod_foc <-  summary(glme.foc.p)[["coefficients"]]
early_mod_rd <- summary(glme.rd.e)[["coefficients"]]
peak_mod_rd  <-  summary(glme.rd.p)[["coefficients"]]

final_model_coefs <- cbind(c(rep("early_foc",6),rep("peak_foc",8),rep("early_rd",6),rep("peak_rd",6)), rbind(early_mod_foc,peak_mod_foc,early_mod_rd,peak_mod_rd))

# write.csv(final_model_coefs, file = "data/output_data/final_model_coefs.csv")




# mapping data ------------------------------------------------------------

#' # Map of Study

states <- map_data("state")
mn_df <- subset(states, region == "minnesota")

mn_lakes <- st_read("G:/My Drive/Documents/UMN/Grad School/Larkin Lab/R_projects/curlyleaf_mgmt/data/spatial_data/MN_lakesonly_trimmed_strytrk.shp")

study_lakes <- st_read("G:/My Drive/Documents/UMN/Grad School/Larkin Lab/R_projects/curlyleaf_mgmt/data/spatial_data/lake_loc_fromRN.shp")

clp_lakes <- st_read("G:/My Drive/Documents/UMN/Grad School/Larkin Lab/R_projects/curlyleaf_mgmt/data/spatial_data/clppoints.shp")
clp_lakes <- clp_lakes %>% st_set_crs(2027) 

study_map <- ggplot(study_lakes)+
  geom_sf(data = clp_lakes, size = 3 ,shape = 1, alpha = 0.5)+
  geom_sf(size = 8, shape = "+", color= "black")+
  geom_polygon(data = mn_df,aes(x = long, y = lat), color = "black", alpha = .05)+
  scale_shape_discrete(solid = FALSE)+
  theme(text = element_text(size=20), legend.position = )+
  ylab("Longitude")+
  xlab("Latitude")

study_map

#add scalebar & North arrow 

study_map+
  north(study_lakes, location = "topright", anchor = c("x"=-90, "y"= 49), scale = 0.2)


#(https://stackoverflow.com/questions/39067838/parsimonious-way-to-add-north-arrow-and-scale-bar-to-ggmap):

scalebar = function(x,y,w,n,d, units="Degrees"){
  # x,y = lower left coordinate of bar
  # w = width of bar
  # n = number of divisions on bar
  # d = distance along each division
  
  bar = data.frame( 
    xmin = seq(0.0, n*d, by=d) + x,
    xmax = seq(0.0, n*d, by=d) + x + d,
    ymin = y,
    ymax = y+w,
    z = rep(c(1,0),n)[1:(n+1)],
    fill.col = rep(c("black","white"),n)[1:(n+1)])
  
  labs = data.frame(
    xlab = c(seq(0.0, (n+1)*d, by=d) + x, x), 
    ylab = c(rep(y-w*1.5, n+2), y-3*w),
    text = c(as.character(seq(0.0, ((n+1)*d), by=d)), units)
  )
  list(bar, labs)
}

sb = scalebar( - 92.5, 49.1, .1, 2, .5 )

# Plot map


study_map+
  north(study_lakes, location = "topright", anchor = c("x"=-90, "y"= 49), scale = 0.2)+
  geom_rect(data=sb[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = sb[[1]]$fill.col) +
  geom_text(data=sb[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) 


# add inset:

study_map <- ggplot(study_lakes)+
  geom_sf(data = clp_lakes, size = 3 ,shape = 1, alpha = 0.5)+
  geom_sf(size = 8, shape = "+", color= "black")+
  geom_polygon(data = mn_df,aes(x = long, y = lat), color = "black", alpha = .05)+
  scale_shape_discrete(solid = FALSE)+
  theme(text = element_text(size=20), legend.position = )+
  ylab("Longitude")+
  xlab("Latitude")

study_map

#inset
foo <- map_data("state")

g2 <- ggplotGrob(
  ggplot() +
    geom_polygon(data = foo,
                 aes(x = long, y = lat, group = group),
                 fill = NA, color = "black", alpha = 0.5) +
    geom_polygon(data = mn_df,aes(x = long, y = lat), color = "black", alpha = .95)+
    theme(panel.background = element_rect(fill = NULL))+
    coord_map("polyconic")+
    theme_bw()+
    theme_inset()
    
)     

g3 <- study_map +
  annotation_custom(grob = g2, xmin = -92.5, xmax = -89.1,
                    ymin = 44.4, ymax = 46.1)


sb = scalebar( - 91.7, 46.5, .1, 2, .5 )



study_map <- ggplot(study_lakes)+
  geom_sf(data = clp_lakes, size = 3 ,shape = 1, alpha = 0.5)+
  geom_sf(size = 8, shape = "+", color= "black")+
  geom_polygon(data = mn_df,aes(x = long, y = lat), color = "black", alpha = .0)+
  scale_shape_discrete(solid = FALSE)+

  annotation_custom(grob = g2, xmin = -92.5, xmax = -89.1,
                    ymin = 44.4, ymax = 46.1)+
  north(study_lakes, location = "topright", anchor = c("x"=-89.3, "y"= 46.7), scale = 0.2)+
  geom_rect(data=sb[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = sb[[1]]$fill.col) +
  geom_text(data=sb[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F)+
  theme_bw()+
  theme(text = element_text(size=20), legend.position = )+
  ylab("Longitude")+
  xlab("Latitude")




# figure construction -----------------------------------------------------
#Figure 1
tiff("Fig1.tiff", res = 600, width = 9, height = 12, units = "in", compression = "lzw")
plot(study_map) # Make plot
dev.off()

# Figure 2
fig2 <- grid.arrange(trt_foc_early,trt_foc_late,trt_rd_early,trt_rd_late, cytrt_early, cytrt_late )
tiff("Fig2.tiff", res = 600, width = 9, height = 12, units = "in", compression = "lzw")
plot(fig2) # Make plot
dev.off()

#Figure 3
fig3 <- grid.arrange( secchi_early, secchi_late, snow_early, snow_late, ice_early, ice_late )
tiff("Fig3.tiff", res = 600, width = 9, height = 12, units = "in", compression = "lzw")
plot(fig3) # Make plot
dev.off()




# footer ------------------------------------------------------------------


#' # Session Information
sessionInfo()
#' ## Cite packages:
citation("ggplot2")
citation("lme4")
citation("dplyr")
citation("data.table")
citation("numDeriv")
citation("sf")
citation("mapdata")
citation("ggmap")

