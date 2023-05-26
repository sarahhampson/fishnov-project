# Title: Novelty Detection Setup
# Author: Sarah Hampson


# 1.  Clear environment and install packages/functions --------------------

# Clear environment
rm(list = ls())
graphics.off()

# Set working directory
#setwd("YOUR OWN CHOICE HERE")
setwd("~/Documents/GitHub/sarahfishproject")

# Obtain required functions from 'functions' subfolder
sapply(list.files("./Functions", pattern="\\.R", full.names=TRUE), source)

# Other random functions

# A little function to add dates to output files
date.wrap <- function(string, ext){
  paste0(string, " ", Sys.Date(), ext)
}

# Small function to add cat.after to results matrices
add.cat.after <- function(x){
  x$cat.after <- "None"
  for(i in 1:(nrow(x))){
    if(i != nrow(x)){
      x$cat.after[i] = x$cat[i+1]
    }else{
      x$cat.after[i] = "None"
      return(x)
    }}}

#### Load packages
install.packages(c("mgcv", "vegan", "lme4", "nlme",
                   "DHARMa", "merTools", "shape",
                   "multcomp", "maptools", "sp", 
                   "divDyn", "plotrix", "raster",
                   "rgeos", "fun", "analogue",
                   "brms", "data.table", "dplyr",
                   "tidyverse", "funrar", "rfishbase",
                   "data.table", "svMisc", "data.table",
                   "tidyr"))

# Libraries required
library(tidyverse)
library(dplyr)

# 2.  Load time series data -----------------------------------------------

# RivFishTIME data
time_series_data <-read.csv("inputs/raw_data/1873_2_RivFishTIME_TimeseriesTable.csv")
survey_data <- read.csv("inputs/raw_data/1873_2_RivFishTIME_SurveyTable.csv")

# 3.  Load fishmorph data -------------------------------------------------

# FishMORPH data plus change species variable name
fishmorph_data <- read.csv("inputs/raw_data/FISHMORPH_Database.csv", sep = ";")
fishmorph_data <- fishmorph_data %>% dplyr::rename("Species" = "Genus.species")

# Align all species names across fishmorph with fishbase
fishmorph_data <- standardise.species.names(fishmorph_data)

# Save to inputs
write.csv(fishmorph_data, "./inputs/data/fishmorph_data.csv")

# 4.  Select time series data ---------------------------------------------

# Select time series data using function
survey_data <- select.timeseries(survey_data, time_series_data)

# Standardise species names
survey_data <- standardise.species.names(survey_data)

# Filter out survey data again with less than 10 years of data 
# (needed again after standardising names)
survey_data <- survey_data %>% 
  group_by(sqID) %>% 
  filter(n_distinct(Year)>=10)

# Modify time series to only include time series from new ts_table
ts_data <- time_series_data[time_series_data$TimeSeriesID %in%
                              survey_data$TimeSeriesID,]

# Look at number of sites, time series, and species
n_distinct(survey_data$TimeSeriesID)
n_distinct(survey_data$sqID)
n_distinct(survey_data$Species)

# 5.  Subset data for species in fishmorph --------------------------------

# Create a list of species from the selected TS
species_list <- survey_data %>% dplyr::summarise(Species = unique(Species))

# Create list of species missing from fishmorph
missing_species <- species_list %>%
  filter(!Species %in% fishmorph_data$Species)

# Identify TS which have these missing species
survey_data_missing_sp <- survey_data %>%
  filter(Species %in% missing_species$Species) %>% 
  dplyr::summarise(sqID = unique(sqID))

# Remove TS which do not have data from survey data
survey_data <- survey_data %>% 
  filter(!sqID %in% survey_data_missing_sp$sqID)

# Remove sites which do not appear in the new survey data
ts_data <- ts_data %>% 
  filter(TimeSeriesID %in% survey_data$TimeSeriesID)

# Make new species list
species_list <- survey_data %>% dplyr::summarise(Species = unique(Species))

# 6.  Save subset survey data and ts data ---------------------------------

write.csv(survey_data, "./inputs/data/survey_data.csv")
write.csv(ts_data, "./inputs/data/ts_data.csv")

# 7.  Cleanup -------------------------------------------------------------

# If necessary
rm(fishmorph_data, species_list, survey_data, time_series_data, ts_data,
   missing_species, survey_data_missing_sp)
