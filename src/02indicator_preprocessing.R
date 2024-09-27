################################################################################
############# SCRIPT FOR PRE-PROCESSING INDICATORS BEFORE ANALYSIS #############
################################################################################

# Author: Mo Osman
# Date created: 19-09-2024
# Last edited: 

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

# Read in functions required for compiling base model: 
source("src/01base_model_functions.R")

#-------------------------------------------------------------------------------

# Get base case apparent intake for each MN: 
nga_base_ai <- apparent_intake("nga_lss1819")

# Tidy environment: 
rm(list = setdiff(ls(), c("nga_base_ai", "allen_ear")))

#-------------------------------------------------------------------------------

# Sense check base apparent intake by plotting histogram of energy intake: 
# hist(nga_base_ai$energy_kcal, main = "Histogram of energy intake") 
# Distribution of energy intake appears plausible, therefore continue with analyses. 

# Binarise vitamin B12 inadequacy: 
nga_base_ai$vb12_inadequate <- ifelse(nga_base_ai$vitb12_mcg < allen_ear$ear_value[allen_ear$nutrient == "vitb12_mcg"], 1, 0)

# Keep only required indicators: 
nga_base_ai <- nga_base_ai %>% 
  dplyr::select(hhid, vb12_inadequate)

rm(allen_ear)

#-------------------------------------------------------------------------------

# SUMMARISE DATA AVAILABILITY AND SAMPLE SIZES: 

hh_locations <- read_csv("shapefiles/household_locations.csv")

nga_base_ai <- nga_base_ai %>% 
  left_join(hh_locations %>% dplyr::select(hhid, lga), by = "hhid")

# Group by lga - and summarise how many households belong to each lga: 
summary_lga <- nga_base_ai %>% 
  group_by(lga) %>% 
  summarise(n = n())

#-------------------------------------------------------------------------------

# FINALISE SCRIPT: 
nga_base_ai <- nga_base_ai %>% 
  dplyr::select(hhid, vb12_inadequate)

# Tidy environment: 
rm(list = setdiff(ls(), c("nga_base_ai")))

################################################################################
############################### END OF SCRIPT ##################################
################################################################################




