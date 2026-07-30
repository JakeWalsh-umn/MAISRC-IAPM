# header ------------------------------------------------------------------


#' # P crispus dataset prep
#' ### Author: Mike Verhoeven
#' ### Date: 27 Sep 2018
#' 
#' 
#' This code has is the final cleanup script in the clp dataset prep. It 
#' combines the A.Kautza, dnr, and point level collation data into a final
#' dataset. it has script to create summary plots of the data such as the
#' treatment timelines.
#' 
#' On 12 August 2019 I revised this script to run with the "multiple 
#' littoral zone cutoff" dataset (surveylaketaxastat) that we ran to 
#' meet the requests of reviewers. 

#' # Preamble
#' Load libraries
  #+warning=FALSE, message=FALSE 
  library(knitr)
  library(nlme)
  library(ggplot2)  
  library(lme4)
  library(ezknitr)
  library(stringr)
  library(GGally)
  library(ResourceSelection)
  library(LaplacesDemon)
  library(effects)
  library(ggrepel)
  library(dplyr)
  library(data.table)



#' Remove anything in memory
  # rm(list = ls())

#' Load custom functions
  # diagnostics
  mv_diag <- function(model_to_diag, plots = c(1:4)) {
    par(mfrow = c(2,2))
    plot(model_to_diag, which = plots) 
    par(mfrow = c(1,1))
  }
#' Ensure working directory is set to project location
  


# load in data ------------------------------------------------------------

#' Load in clp dataset :
  pcri <- fread("data/output_data/pcrispus_wDNRtrtdat.csv")
#' Load in the pcri_foc dataset (from my statewide PI compilation)
  foc <- fread("data/output_data/surveylktaxastat.csv")


# visualize ---------------------------------------------------------------



  foc$juliandate <- as.POSIXlt(foc$datemv, format = "%Y-%m-%d")$yday
  ggplot(foc[taxon== "potamogeton.crispus"], aes(juliandate, litfoc))+  
    geom_point(aes(juliandate, litfoc),size = 3, alpha = .25)+
    theme(text = element_text(size=20))+
    ylab("Proportion of Littorral")+
    xlab("Julian Day")+
    scale_x_continuous(breaks = c(100,150,200,250,300), limits = c(50,300))+
    geom_smooth(se=T, colour="red")
  ggplot(foc[taxon== "potamogeton.crispus"], aes(juliandate, rd))+  
    geom_point(aes(juliandate, rd),size = 3, alpha = .25)+
    theme(text = element_text(size=20))+
    ylab("Mean Relative Rake Dens")+
    xlab("Julian Day")+
    scale_x_continuous(breaks = c(100,150,200,250,300), limits = c(50,300))+
    geom_smooth(se=T, colour="red")
  
  # how many lake surveys do we have?
  foc[,.N, c("lknamemv", "datasourcemv", "datemv")]
  summarydat <- foc[]


# clean up data -----------------------------------------------------------
#' Clean and prep the two datasets:
  summary(foc)

  # need to parse the date in order to match to p.cri sampling timeframes
  foc[, "year" := word(foc$datemv, start = 1,  sep = "-"), ]
  foc[, "month" := word(foc$datemv, start = 2,  sep = "-"), ]
  foc[, "day" := word(foc$datemv, start = 3,  sep = "-"), ]


  
  ggplot(foc[taxon=="potamogeton.crispus",], aes(month, lslfoc))+  
    geom_point(aes(month, lslfoc), alpha = .25)
  ggplot(foc[taxon=="potamogeton.crispus",], aes(month, rd))+  
    geom_point(aes(month, rd), alpha = .25)
  
  head(foc)
  head(pcri)
  
  # foc <- data.frame(foc)
  # pcri <- data.frame(pcri)
  
  #clean up lake names
  names(pcri)[4] <- "Lake"
  sort(unique(pcri[,Lake]))
  sort(unique(foc[,lknamemv]))
  
  sort(unique(pcri[,Lake.ID]))
  sort(unique(foc[,Lake.ID]))
  
    #check results
  cbind (sort(unique(pcri[,Lake])),sort(unique(foc[,lknamemv])))
  cbind (sort(unique(pcri[,Lake.ID])),sort(unique(foc[,Lake.ID])))
  #which lakes in pcri (Adams data) are missing foc matches?
  unique(pcri[,Lake])[unique(pcri[,Lake])%in%unique(foc[,lknamemv]) == F]
  unique(pcri[,Lake.ID])[unique(pcri[,Lake.ID])%in%unique(foc[,Lake.ID]) == F]
 

# create centered scaled vars --------------------------------------------

#' Create centered variables (used in glms):
  #all except Yrs Trt
  pcri$Snow.sc <- c(scale(pcri$Snow))
  pcri$PrevAUGSecchi.sc <- c(scale(pcri$PrevAUGSecchi))
  pcri$IceCover.sc <- c(scale(pcri$IceCover))
  
# add in native species richness values -----------------------------------
  
  #calculate a total and native species richness value or each survey 
  foc[, trichness := .N, c("datemv","datasourcemv","lknamemv","Lake.ID") ]  
  
  foc[taxon != "potamogeton.crispus" &
        taxon!="myriophyllum.spicatum", nrichness := .N, c("datemv","datasourcemv","lknamemv","Lake.ID") ]


# create early season dataset ---------------------------------------------
  #only pcri data
  foc <- subset(foc, foc$taxon == "potamogeton.crispus")
  
  
#' split the early (.e) and peak (.p) data from p.cri and create datasets

  #drop the peak data from early set:
  pcri.e <- pcri
  names(pcri.e)
  pcri.e <- pcri.e[, - c(17:24)]
  # grab out a subset of the pcri foc data that include only the APR/MAY surveys
  foc.e <- subset(foc, foc$month == "04" |
                    foc$month == "05")
  foc.e[,1]
  
  #test a match based on lake and year  
  unique(pcri.e[,c("Lake", "Year")])
  unique(foc.e[,c("lknamemv", "year")])
  sort(unique(pcri.e$Lake))
  sort(unique(foc.e$lknamemv))
  
  #combine lake and year to give a lk-yr field
  pcri.e$lake_year <- paste(pcri.e$Lake, pcri.e$Year, sep = "")
  foc.e$lkyr <- paste(foc.e$lknamemv, foc.e$year, sep = "")
  
  # commbine new calcs into datasets 
  pcri.e$tnsamp <- foc.e[match(pcri.e$lake_year, foc.e$lkyr),"tnsamp"]
  pcri.e$lnsamp <- foc.e[match(pcri.e$lake_year, foc.e$lkyr),"lnsamp"]
  pcri.e$lslsamp <- foc.e[match(pcri.e$lake_year, foc.e$lkyr),"lslsamp"]
  
  pcri.e$totfoc <- foc.e[match(pcri.e$lake_year, foc.e$lkyr),"totfoc"]
  pcri.e$litfoc <- foc.e[match(pcri.e$lake_year, foc.e$lkyr),"litfoc"]
  pcri.e$lslfoc <- foc.e[match(pcri.e$lake_year, foc.e$lkyr),"lslfoc"]
  
  pcri.e$maxdep <- foc.e[match(pcri.e$lake_year, foc.e$lkyr), "maxdep"]
  
  hist(na.omit(foc.e[match(pcri.e$lake_year, foc.e$lkyr),rd]), )
  pcri.e$rd_pcri <- foc.e[match(pcri.e$lake_year, foc.e$lkyr),"rd"]
  pcri.e$datasource <- foc.e[match(pcri.e$lake_year, foc.e$lkyr), "datasourcemv"]
  
  # visualize new stats

  #plot new v old foc calcs:
    # check lake specific littoral values against littoral cutoff
    ggplot(pcri.e, aes(litfoc,lslfoc))+
      geom_point()+
      geom_abline(slope = 1, intercept = 0)  
    ggplot(pcri.e, aes(totfoc,litfoc))+
      geom_point()+
      geom_abline(slope = 1, intercept = 0)



# create peak season dataset ----------------------------------------------
  #drop earlyseason stuff from peak dataset
  pcri.p <- pcri

  # grab out a subset of the pcri foc data that include only the June surveys
  foc.p <- subset(foc, foc$month == "06")
  foc.p[,1]

  #' make sure the names and years will match up
  #test a match based on lake and year  
  unique(pcri.p[,c("Lake", "Year")])
  unique(foc.p[,c("lknamemv", "year")])
  sort(unique(pcri.p$Lake))
  sort(unique(foc.p$lknamemv))

#' combine lake and year to give a lk-yr field
  pcri.p$lake_year <- paste(pcri.p$Lake, pcri.p$Year, sep = "")
  foc.p$lkyr <- paste(foc.p$lknamemv, foc.p$year, sep = "")

  #combine:
  pcri.p$tnsamp <- foc.p[match(pcri.p$lake_year, foc.p$lkyr),"tnsamp"]
  pcri.p$lnsamp <- foc.p[match(pcri.p$lake_year, foc.p$lkyr),"lnsamp"]
  pcri.p$lslsamp <- foc.p[match(pcri.p$lake_year, foc.p$lkyr),"lslsamp"]
  
  pcri.p$totfoc <- foc.p[match(pcri.p$lake_year, foc.p$lkyr),"totfoc"]
  pcri.p$litfoc <- foc.p[match(pcri.p$lake_year, foc.p$lkyr),"litfoc"]
  pcri.p$lslfoc <- foc.p[match(pcri.p$lake_year, foc.p$lkyr),"lslfoc"]
  
  pcri.p$maxdep <- foc.p[match(pcri.p$lake_year, foc.p$lkyr), "maxdep"]
  
  hist(foc.p[match(pcri.p$lake_year, foc.p$lkyr),rd])
  pcri.p$rd_pcri <- foc.p[match(pcri.p$lake_year, foc.p$lkyr),"rd"]
  pcri.p$datasource <- foc.p[match(pcri.p$lake_year, foc.p$lkyr), "datasourcemv"]
  
  
  
  #plot new v old foc calcs:
  names(pcri.p)
  ggplot(pcri.p, aes(lslfoc,litfoc))+
    geom_point()+
    geom_abline(slope = 1, intercept = 0)


# calculating consecutive years of treatment ------------------------------

  
  
#' The YrsTrt measures in June surveys is correct 1; = Treated this year, However
#' the spring surveys should be asking about wheter they were treated last spring. Thne
#' we need to calculate cumulative years of trt as years in additon to year one. 
#' This separates immediate and cumulative effects.
#' Option 2: Reaasign a value for consecutive years treated or untreated --
#' 
  
  # convert data back to df... only because I didn't know data.table when I wrote this script
  pcri.p <- data.frame(pcri.p)
  pcri.e <- data.frame(pcri.e)
  
  pcri.p$cytrt <- NA
  pcri.p$yutrt <- NA
  for (j in (unique(pcri.p$Lake) )) {
    # j = "silver"
    print( pcri.p[pcri.p[,"Lake"] == j , c( "Lake" , "Year", "Trt", "YrsTrt") ])
    a = pcri.p[pcri.p[,"Lake"] == j , c( "Lake" , "Year", "Trt", "YrsTrt", "cytrt", "yutrt") ]
    
    #calc yrs trt
    for (i in (1:nrow(a))) {
      if (i==1) a$cytrt[i] <- 0 else 
        ifelse(a$Trt[i] == "U", a$cytrt[i] <- 0, 
               ifelse(a$Trt[i-1]=="U", a$cytrt[i] <-  0 ,
                      ifelse( a$Year[i] - a$Year[i-1] == 1, a$cytrt[i] <- a$cytrt[i-1]+1 , 
                              a$cytrt[i] <- 0 
                      )
               )
        )          
    }
    #yrs un trt calculation for peak data
    for (i in (1:nrow(a))) {
      #i <- 6
      if (i==1) a$yutrt[i] <- 0 else
        ifelse(a$Trt[i]== "T", a$yutrt[i] <- 0 , 
               ifelse(a$Trt[i-1] == "T", ifelse(a$Year[i] - a$Year[i-1] == 1, a$yutrt[i] <- (a$yutrt[i-1]+1), a$yutrt[i] <- (a$Year[i] - a$Year[i-1])),
                      ifelse(a$yutrt[i-1] == 0, a$yutrt[i] <- 0,  a$yutrt[i] <- (a$yutrt[i-1]+1))
               )
        )
    }
    
    #add into the dataset
    pcri.p[pcri.p[,"Lake"] == j , c( "Lake" , "Year", "Trt", "YrsTrt", "cytrt", "yutrt") ] <- a
    print(a)
  }
  
  plot(lslfoc ~ yutrt, data = pcri.p) 
  ggplot(pcri.p, aes(yutrt, lslfoc))+
    geom_path(aes(yutrt, lslfoc, group = Lake))
  
  plot(lslfoc ~ cytrt, data = pcri.p)
  ggplot(pcri.p, aes(cytrt, lslfoc))+
    geom_path(aes(cytrt, lslfoc, group = Lake))
  
  # and now for early data:
  pcri.e$cytrt <- NA
  pcri.e$yutrt <- NA
  pcri.e$Trte <- NA
  pcri.e$perc.cheme <- NA
  for (j in (unique(pcri.e$Lake) )) {
    # j = "reshanau"
    print( pcri.e[pcri.e[,"Lake"] == j , c( "Lake" , "Year", "Trt", "YrsTrt") ])
    a = pcri.e[pcri.e[,"Lake"] == j , c( "Lake" , "Year", "Trt","Trte", "YrsTrt", "cytrt", "yutrt", "perc.lit.chem") ]
    
    # yrs trt calculation
    for (i in (1:nrow(a))) {
      # assign appropriate treatment value (U/T)
      # i <- 6
      if (i==1) a$Trte[i] <- "U" else
        ifelse( a$Year[i] - a$Year[i-1] != 1, a$Trte[i] <- "U",
                ifelse(a$Trt[i-1] == "U", a$Trte[i] <- "U",
                       a$Trte[i] <- "T"
                )
        )
      # assign appropriate treatment acreage (last year's)
      # i <- 6
      if (i==1) a$perc.cheme[i] <- 0 else
        ifelse( a$Year[i] - a$Year[i-1] != 1, a$perc.cheme[i] <- 0,
                ifelse(a$Trt[i-1] == "U", a$perc.cheme[i] <- 0,
                       a$perc.cheme[i] <- a$perc.lit.chem[i-1]
                )
        )
      #assign cum yrs trt value
      if (i==1) a$cytrt[i] <- 0 else 
        ifelse( a$Year[i] - a$Year[i-1] != 1, a$cytrt[i] <- 0,
                ifelse(a$Trte[i] == "U", a$cytrt[i] <-  0, 
                  ifelse(a$Trte[i-1] == "T", a$cytrt[i] <- a$cytrt[i-1] + 1 , a$cytrt[i] <-  0
                  )
                )
        )          
    }
    
    # yrs un-trt calculation
    for (i in (1:nrow(a))) {
      #i <- 6
      if (i==1) a$yutrt[i] <- 0 else
        ifelse(a$Trte[i]== "T", a$yutrt[i] <- 0 , 
               ifelse(a$Trte[i-1] == "U", a$yutrt <- 0, a$yutrt[i] <- (a$Year[i] - a$Year[i-1])
               )
        )
    }
    #add into the dataset
    pcri.e[pcri.e[,"Lake"] == j , c( "Lake" , "Year", "Trt","Trte", "YrsTrt", "cytrt", "yutrt", "perc.lit.chem", "perc.cheme") ] <- a
    print(a)
  }
  
  pcri.e <- pcri.e[order(pcri.e$Lake, pcri.e$Year),]
  summary(pcri.e$yutrt)
  plot(lslfoc~yutrt, data = pcri.e) 
  # pcri.e$trtcolplot <- ifelse(pcri.e$Trte)
  ggplot(pcri.e, aes(yutrt, APRMAYPcri))+
    geom_path(aes(jitter(yutrt), APRMAYPcri, group = Lake, colour = Trte), arrow = arrow(ends = "last", length = unit(.1, "inches"), type = "closed"))
  
  plot(lslfoc~cytrt, data = pcri.e)
  ggplot(pcri.e, aes(cytrt, lslfoc))+
    geom_path(aes(cytrt, lslfoc, group = Lake))
  ggplot(pcri.e, aes(cytrt, perc.cheme))+
    geom_path(aes(cytrt, perc.cheme, group = Lake))

  # make trt a 0/1:
  str(pcri.e)
  
  pcri.e$Trt <- as.character(pcri.e$Trt)
  pcri.e[pcri.e[, "Trt"] == "T", "Trt"] <- "1"
  pcri.e[pcri.e[, "Trt"] == "U", "Trt"] <- "0"
  pcri.e$Trt <- as.factor(pcri.e$Trt)
  
  pcri.e$Trte <- as.character(pcri.e$Trte)
  pcri.e[pcri.e[, "Trte"] == "T", "Trte"] <- "1"
  pcri.e[pcri.e[, "Trte"] == "U", "Trte"] <- "0"
  pcri.e$Trte <- as.factor(pcri.e$Trte)
  
  pcri.p$Trt <- as.character(pcri.p$Trt)
  pcri.p[pcri.p[, "Trt"] == "T", "Trt"] <- "1"
  pcri.p[pcri.p[, "Trt"] == "U", "Trt"] <- "0"
  pcri.p$Trt <- as.factor(pcri.p$Trt)

# check for missing data --------------------------------------------------
  
  #' what else is missing?
  summary(pcri.e[,c("Snow", "IceCover", "PrevAUGSecchi", "YrsTrt", "tnsamp", "lnsamp", "lslsamp", "totfoc", "litfoc","lslfoc")])
  summary(pcri.p[,c("Snow", "IceCover", "PrevAUGSecchi", "YrsTrt", "tnsamp", "lnsamp", "lslsamp", "totfoc", "litfoc","lslfoc")])
  #' remove rows without secchi or occurrance (submitted at point level) data
  #early
  pcri.e <- subset(pcri.e, is.na(pcri.e$tnsamp) == FALSE )
  summary(pcri.e[,"tnsamp"])
  pcri.e <- subset(pcri.e, is.na(pcri.e$PrevAUGSecchi) == FALSE )
  summary(pcri.e[,"PrevAUGSecchi"])
  #peak
  pcri.p<- subset(pcri.p, is.na(pcri.p$tnsamp) == FALSE )
  summary(pcri.p[,"tnsamp"])  
  pcri.p<- subset(pcri.p, is.na(pcri.p$PrevAUGSecchi) == FALSE )  
  summary(pcri.p[,"PrevAUGSecchi"])
  

# rake density cleaning ---------------------------------------------------
  #' - Attempt to squeeze all rake densities into a 1-3 (1, middle, max), especially sblood, brasch, and crwd

  #check for completeness of data used in 
  summary(pcri.e[,c("rd_pcri","Snow.sc","IceCover.sc", "PrevAUGSecchi.sc","cytrt","Trte","Lake" )]) 
  
  # corrections for rake density data:
  #' We need to check who contributed these surveys because not everyone runs 1-4 scales
  #' 
  pcri.e[which(is.na(pcri.e[,"rd_pcri"]) == F), c("datasource", "Lake", "Year", "datasource") ]  
  
  # definitely delete entries from fisheries data: they are all bad (collected 0/1).
  pcri.e[which(is.na(pcri.e[,"rd_pcri"]) == F & 
                 pcri.e[,"datasource"] == "dustin"), "rd_pcri"] <- NA
  # sblood data are 1-5, delete for now
  pcri.e[which(is.na(pcri.e[,"rd_pcri"]) == F & 
                 pcri.e[,"datasource"] == "sblood"), "rd_pcri"] <- NA
  # crwd data are 1-5, delete for now
  pcri.e[which(is.na(pcri.e[,"rd_pcri"]) == F & 
                 pcri.e[,"datasource"] == "crwd"), "rd_pcri"] <- NA
  # brasch data are 1-5, delete for now
  pcri.e[which(is.na(pcri.e[,"rd_pcri"]) == F & 
                 pcri.e[,"datasource"] == "brasch"), "rd_pcri"] <- NA
  
  #' Again, but with Peak season data:
  #' 
  pcri.p[which(is.na(pcri.e[,"rd_pcri"]) == F), c("datasource", "Lake", "Year", "datasource") ]  
  
  # definitely delete entries from fisheries data: they are all bad (collected 0/1).
  pcri.p[which(is.na(pcri.p[,"rd_pcri"]) == F & 
                 pcri.p[,"datasource"] == "dustin"), "rd_pcri"] <- NA
  # sblood data are 1-5, delete for now
  pcri.p[which(is.na(pcri.p[,"rd_pcri"]) == F & 
                 pcri.p[,"datasource"] == "sblood"), "rd_pcri"] <- NA
  # crwd data are 1-5, delete for now
  pcri.p[which(is.na(pcri.p[,"rd_pcri"]) == F & 
                 pcri.p[,"datasource"] == "crwd"), "rd_pcri"] <- NA
  # brasch data are 1-5, delete for now
  pcri.p[which(is.na(pcri.p[,"rd_pcri"]) == F & 
                 pcri.p[,"datasource"] == "brasch"), "rd_pcri"] <- NA
  

# delete old calculated values  -------------------------------------------
  pcri.e <- data.table(pcri.e)
  pcri.p <- data.table(pcri.p)
  # drop unneeded columns (calculated without method tracking by prev collab)
  names(pcri.e)
  pcri.e[ , c("APRMAYPcri", "APRMAY3.7", "APRMAYDIFF", "APRMAYRelDens", "APRMAYRelDensall",
              "Latitude.y", "Lake.y", "lake_year", "County") := NULL ]
  names(pcri.p)
  pcri.p[ , c("APRMAYPcri", "APRMAY3.7", "APRMAYDIFF", "APRMAYRelDens", "APRMAYRelDensall",
              "Latitude.y", "Lake.y", "lake_year", "JUNPcri", "JUN3.7", "JUNDIFF",
              "JUNRelDens", "JUNRelDensall", "DIFF.4.6", "DIFF.3.7", "DIFFRelDens", "County" ) := NULL ]
  names(pcri.e)
  names(pcri.p)
  
# export data -------------------------------------------------------------
  
  
  # write.csv(pcri.e, file = "data/output_data/pcri.e_laklvllit.csv")
  # write.csv(pcri.p, file = "data/output_data/pcri.p_laklvllit.csv")
  
  

# trajectories through time, gridded by lake ----------------------------------------------------------
  
  
  # management trajectories through time, gridded by lake.
  pcri.p <- pcri.p[order(pcri.p$Year, pcri.p$Lake),]
  pcri.p$Year <- as.integer(pcri.p$Year)
  
  #reorder by first year of sampling per lake
  pcri.e$Lake <- factor(pcri.e$Lake, levels = unique(pcri.e$Lake[order(pcri.e$Year)]))
  summary(pcri.e[,"tnsamp"])
  #peak
  summary(pcri.p[,"tnsamp"])
  pcri.p$Lake <- factor(pcri.p$Lake, levels = unique(pcri.p$Lake[order(pcri.p$Year)]))
  #x$name <- factor(x$name, levels = x$name[order(x$val)])
  
  # ggplot(pcri.p, aes( x = Year, y = foc_clp))+
  #   geom_line()+
  #   geom_point(aes(shape = Trt))+
  #   geom_line(data = pcri.e, aes(x = Year, y = foc_clp), linetype = 2)+
  #   geom_point(data = pcri.e, aes(x = Year, y = foc_clp, shape = Trte))+
  #   facet_wrap( ~ Lake, ncol = 5, dir = "v")+
  #   scale_shape_manual(values = c(21,16))+
  #   scale_fill_manual(values = c("white", "black"))+
  #   scale_y_continuous(breaks = c(0,1), minor_breaks = c(), labels = c(0,1))+
  #   scale_x_continuous(breaks = c(2005,2010,2015), minor_breaks = seq(2005,2015,1))
  # 
  
  # ggplot(pcri.e, aes( x = Year, y = foc_clp))+
  #   geom_line()+
  #   geom_point(aes(shape = Trte))+
  #   facet_wrap( ~ Lake, nrow = 10, dir = "v")+
  #   scale_shape_manual(values = c(21,16))+
  #   scale_fill_manual(values = c("white", "black"))
  # 
  # pcri$Lake <- factor(pcri$Lake, levels = unique(pcri$Lake[order(pcri$Year)]))
  # 
  # ggplot(pcri, aes( x = as.numeric(Year), y = Latitude, group = Lake))+
  #   geom_line(aes( x = as.numeric(Year), y = APRMAYPcri, group = Lake), linetype = 2)+
  #   geom_point(aes(y = 1, shape = Trt))+
  #   geom_point(aes(y = 1, shape = Trt))+
  #   geom_line(aes( x = as.numeric(Year), y = JUNPcri, group = Lake))+
  #   facet_wrap( ~ Lake, ncol = 6, dir = "v")+
  #   scale_y_continuous(breaks = c(0,1), minor_breaks = c(), labels = c(0,1))+
  #   scale_shape_manual(values = c(16,21))+
  #   #scale_color_manual(values = c("black", "white"))+
  #   scale_x_continuous(breaks = c(2006,2010,2015), minor_breaks = seq(2006,2015,1), labels = c("'06","'10","'15"))+
  #   theme(legend.position = "none")+
  #   ylab("Frequency of Occurrence")
  # 
  # hill <- subset(pcri, pcri$Lake == "hill")
  # 
  # ggplot(hill, aes(x = as.numeric(Year), y = Latitude))+
  #   scale_y_continuous(breaks = c(0,1), minor_breaks = c(), labels = c(0,1))+
  #   geom_point(aes(y = APRMAYPcri))+
  #   geom_point(aes(y = JUNPcri), shape = 21)+
  #   geom_line(aes(x = as.numeric(Year), y = as.numeric(APRMAYPcri)), na.rm = T)
  # 
  
  plotdat.e <- data.frame(lake = pcri.e$Lake, trt = pcri.e$Trte ,
                          year = pcri.e$Year , period = rep("early", length(pcri.e$Lake)) ,
                          foc = pcri.e$litfoc )  
  plotdat.p <- data.frame(lake = pcri.p$Lake, trt = pcri.p$Trt ,
                          year = pcri.p$Year , period = rep("peak", length(pcri.p$Lake)) ,
                          foc = pcri.p$litfoc )  
  plotdat <- rbind(plotdat.e,plotdat.p)
  plotdat$lake <- factor(plotdat$lake, levels = unique(plotdat$lake[order(plotdat$year)]))
  plotdat$period <- factor(plotdat$period, levels = c("peak", "early"))
  
  ggplot(plotdat, aes( x = year, y = foc))+
    geom_line(aes( x = year, y = foc, linetype = period))+
    geom_point(aes( y = foc, shape = trt, colour = period))+
    facet_wrap( ~ lake, ncol = 5, dir = "v")+
    scale_y_continuous(breaks = c(0,1), minor_breaks = c(), labels = c(0,1))+
    scale_shape_manual(values = c(21,16))+
    scale_color_manual(values = c("black", "red"))+
    scale_x_continuous(breaks = c(2006,2010,2015), minor_breaks = seq(2006,2015,1), labels = c("'06","'10","'15"))+
    theme(legend.position = "none")+
    ylab("Frequency of Occurrence")+
    xlab("Year")
  
  
# Note that 10 lakes have no data derived from the point-based dataset 
  length(unique(plotdat$lake))
  length(unique(pcri$Lake))
  sort(unique(plotdat$lake))
  sort(unique(pcri$Lake))
  
# footer ------------------------------------------------------------------


#' ## Document footer 
#' 
#' Document spun with: ezspin("f_p_crispus_dataset_prep.R", out_dir = "coding_html_outputs/p_crispus_dataset_prep", fig_dir = "figures", keep_md=FALSE)
#' 
#' Session Information:
#+ sessionInfo
sessionInfo()





























