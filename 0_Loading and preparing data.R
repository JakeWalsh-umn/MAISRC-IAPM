# MAISRC IAPM Phase I
# 0 - Loading and preparing data

# This script will load in and do minor cleaning for the following core datasets:
## PI Charter Database
## IAPM permitting data
## AIS Control Grant PARs


# Packages ----
library(tidyverse)
library(arrow)
library(readxl)

# Set up repository folders ----
# This will store data I pull from the source locations

ifelse(!dir.exists("data"), 
       dir.create("data"), 
       "Folder exists already")

# You can save figures, tables, and output data products in these folders
# These are ignored in .gitignore, so they won't be pushed as you generate
# them.
ifelse(!dir.exists("figures"), 
       dir.create("figures"), 
       "Folder exists already")
ifelse(!dir.exists("tables"), 
       dir.create("tables"), 
       "Folder exists already")
ifelse(!dir.exists("exploration"), 
       dir.create("exploration"), 
       "Folder exists already")
ifelse(!dir.exists("outputdata"), 
       dir.create("outputdata"), 
       "Folder exists already")
ifelse(!dir.exists("draftcode"), 
       dir.create("draftcode"), 
       "Folder exists already")

# Reference tables ----

# Column descriptions for db_unified.parquet
directory <- read.csv("data/PICharterDirectory.csv")

# PI Charter common and scientific name look up table
taxa <- read.csv("data/commonsci_name_lookup.csv")

# MN DNR Taxa List Export
taxa_xwalk <- read_xlsx("data/MNDNR_AquaticPlant_TaxaListExport.xlsx")

# PI Charter Database ----

db <- read_parquet("data/db_unified_2026-02-12.parquet")

# reorganize db

db <- cbind(db[, colnames(db)%in%directory$fieldNames[directory$Taxonomic=="N"]],
                db[, colnames(db)%in%directory$fieldNames[directory$Taxonomic=="Y"]])

## Sorting out cell values for taxonomic fieldNames ----
# These are not all rake ratings, or even all numbers, so we'll need to figure out
# what's in each column and what we need to do with them in order to 
# work with them. As a first pass, we can just figure out what's a 0 and what's a 1
# then we can try to sort out rake ratings.



## Subsetting submersed species ----

taxa_xwalk$TAXON <- taxa_xwalk$SCIENTIFIC_NAME
taxa_xwalk$TAXON <- tolower(taxa_xwalk$TAXON)
taxa_xwalk$TAXON <- gsub(x=taxa_xwalk$TAXON, pattern=", ", replacement="_")
taxa_xwalk$TAXON <- gsub(x=taxa_xwalk$TAXON, pattern="\\. ", replacement="_")
taxa_xwalk$TAXON <- gsub(x=taxa_xwalk$TAXON, pattern=",", replacement="_")
taxa_xwalk$TAXON <- gsub(x=taxa_xwalk$TAXON, pattern="\\.", replacement="_")
taxa_xwalk$TAXON <- gsub(x=taxa_xwalk$TAXON, pattern=" ", replacement="_")

taxa_xwalk$TAXON[grepl("genus", taxa_xwalk$COMMON_NAME)] <- paste0(taxa_xwalk$TAXON[grepl("genus", taxa_xwalk$COMMON_NAME)], "_sp")

# Joining taxa_xwalk to directory

directory_taxonomic <- left_join(directory %>% filter(Taxonomic=="Y"),
                       taxa_xwalk %>% dplyr::select(TAXON, VEGETATION_ID, ORIGIN),
                       by=c("fieldNames"="TAXON"))

# Export directory_taxonomic for filling in missing VEGETATION_ID and ORIGIN

write.csv(x=directory_taxonomic, file="outputdata/directory_taxonomic.csv")


## COMMENTED OUT UNTIL WE HAVE THE XLSX FILE FILLED ----
# # Creating "submersed" object to filter "db" to just submersed taxa
# 
# directory_taxonomic_updated <- read_xlsx("outputdata/directory_taxonomic_updated.xlsx")
# 
# submersed <- directory_taxonomic_updated$fieldNames[grepl("S", directory_taxonomic_updated$VEGETATION_ID)]
# 
# # Creating "db_sub"
# 
# db_sub <- cbind(db[, colnames(db)%in%directory$fieldNames[directory$Taxonomic=="N"]],
#                 db[, colnames(db)%in%submersed])
# 
