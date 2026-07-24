# Visualizing data from the MN DNR Infested Waters List
# https://www.dnr.state.mn.us/invasives/ais/infested.html

# Packages ----
library(tidyverse)
library(readxl)

library(ggplot2)
library(cowplot)
library(scales)
library(viridis)

library(mgcv)
library(gratia)
library(strucchange)

# Data (RUN SECTION TO UPDATE INFESTED WATERS DATA) ----
# The code below will automatically download the newest infested-waters.xlsx document from our webpage
# as long as the url for the download hasn't changed. The code will also get the data into 
# a workable/summarized formats for the visualizations below. 

# Click the down arrow by the line number next to "Data ..." to hide/show
# the code for downloading and updating the infested waters data used here.

# Infested Waters URL
url_iw <- "https://files.dnr.state.mn.us/eco/invasives/infested-waters.xlsx"

# Path to where I want data stored
main_dir <- getwd()
data_path <- "Infested Waters Data/"

# Create a folder to store infested waters data downloads
if(!file.exists(data_path)){
  dir.create(file.path(main_dir, data_path))
} else {
  print("Infested Waters Data folder already exists.")
}

# Date data was downloaded
# This will make it so that a new file is generated in my Data folder every time I run this code
date_downloaded <- Sys.Date()

# Name I want to save infested waters under
data_name <- paste("infested-waters_", date_downloaded, ".xlsx", sep="")

# Paste together data_path and data_name to tell R where to save downloaded infested waters list
destfile <- paste(data_path, data_name, sep="")

# Download infested waters list
download.file(url_iw, destfile, mode='wb')

# Read in infested waters list
iw <- read_excel(path=destfile, skip=1)

# Check to make sure cleaning code is still relevant/useful.
head(iw)
str(iw)
summary(iw)

iw <- iw[, -dim(iw)[2]]

# Change column names
# I normally don't do this, but the default are pretty bulky to code with
colnames(iw) <- c("waterbody_name", "county", "species", "year", "year_conf", "dowlknum")

iw$connected <- factor(ifelse(grepl("connect", iw$year_conf) | grepl("Connect", iw$year_conf) | grepl("conect", iw$year_conf),
                              "connected", "confirmed"))

unique(iw$species)

iw$species[iw$species=="Eurasian Watermilfoil"] <- "Eurasian watermilfoil"

# Fixing some dowlknum entries, hyphenating these really helps here
# This code fixes dowlknum issues up to 4/26/2023
# CC-LLLL-BB, where C=county, L=lake or parent dow, B=basin, is a nice convention to standardize to
unique(iw$dowlknum[!grepl("-", iw$dowlknum)])

iw$dowlknum[iw$dowlknum=="NA"] <- NA
iw$dowlknum[iw$dowlknum=="na"] <- NA
iw$dowlknum[iw$dowlknum=="none"] <- NA
iw$dowlknum[iw$dowlknum=="NONE"] <- NA

iw$dowlknum[iw$dowlknum=="none, part of Winnibigoshish"] <- "11-0147"



# There are some dowlknum that use hyphens, but not according to the convention above
# These all look good: CC-LLLL
unique(iw$dowlknum[grepl("-", iw$dowlknum) & nchar(iw$dowlknum)==7])
# None here:
unique(iw$dowlknum[grepl("-", iw$dowlknum) & nchar(iw$dowlknum)==8])
# Here's a mistake where there isn't a second hyphen
unique(iw$dowlknum[grepl("-", iw$dowlknum) & nchar(iw$dowlknum)==9])
# fix second hyphen
iw$dowlknum[iw$dowlknum=="18-012601"] <- "18-0126-01"
# These all look good; CC-LLLL-BB
unique(iw$dowlknum[grepl("-", iw$dowlknum) & nchar(iw$dowlknum)==10])

iw$parentdow <- substr(iw$dowlknum, 1, 7)

iw <- iw %>%
  mutate(
    parentdow = str_remove(parentdow, "-")
  )

# End of script use for me - MC, 072426 ----
