################################################################################
############# SCRIPT FOR PRE-PROCESSING INDICATORS BEFORE ANALYSIS #############
################################################################################

# Author: 
# Date created: 28-07-2025
# Last edited: 28-07-2025

# In this script, I will extract and pre-process indicators from the NLSS survey 
# ready for analyses.

# INSTALL AND LOAD PACKAGES:

rq_packages <- c("readr", "tidyverse")

installed_packages <- rq_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(rq_packages[!installed_packages])
}

lapply(rq_packages, require, character.only = T)

rm(list= c("rq_packages", "installed_packages"))

#-------------------------------------------------------------------------------

# Read in pre-processed indicators:
# source("src/02indicator_preprocessing.R")
# source("src/helper_functions.R")

# Read in household locations that are matched to shapefiles: 
household_locations <- read_rds("shapefiles/rwa_hh_adm.Rds")|> rename(adm1_name = adm1, adm2_name = adm2)

# Get further household info:
hh_info <- read_csv("processed_data/rwa_eicv2324_hh_info.csv") |> rename(adm1_code = adm1, amd2_code = adm2)

hh_info <- hh_info |> 
  left_join(household_locations, by = 'hhid') 

# sub-sampling to make it admin1 \times Urban-rural representative


