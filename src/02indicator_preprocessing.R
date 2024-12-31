################################################################################
############# SCRIPT FOR PRE-PROCESSING INDICATORS BEFORE ANALYSIS #############
################################################################################

# Author: Mo Osman
# Date created: 19-09-2024
# Last edited: 31-12-2024

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

# Read in functions required for compiling apparent intake estimates: 
source("src/01apparent_intake_functions.R")

#-------------------------------------------------------------------------------

# Get base case apparent intake for each MN: 
nga_base_ai <- apparent_intake("nga_lss1819")

# Tidy environment: 
rm(list = setdiff(ls(), c("nga_base_ai", "allen_ear")))

#-------------------------------------------------------------------------------

# Binarise Vitamin A, Folate and Vitamin B12 inadequacy: 
nga_base_ai$vita_rae_inadequate <- ifelse(nga_base_ai$vita_rae_mcg < allen_ear$ear_value[allen_ear$nutrient == "vita_rae_mcg"], 1, 0)
nga_base_ai$folate_inadequate <- ifelse(nga_base_ai$folate_mcg < allen_ear$ear_value[allen_ear$nutrient == "folate_mcg"], 1, 0)
nga_base_ai$vb12_inadequate <- ifelse(nga_base_ai$vitb12_mcg < allen_ear$ear_value[allen_ear$nutrient == "vitb12_mcg"], 1, 0)


# Vitamin A, Folate and Vitamin B12 inadequacy variables prepared.
# Zinc and Iron inadequacy to be done when Mo returns from BIS.

nga_base_ai <- nga_base_ai %>% 
  dplyr::select(hhid, vita_rae_inadequate, folate_inadequate, vb12_inadequate)

#-------------------------------------------------------------------------------

# Extract data on binarised household food group consumption. These variables can 
# be used for model validation: 

nga_food_groups <- read_csv("processed_data/nga_lss1819_food_consumption.csv")

# Animal-sourced food consumption (dairy, eggs, meat, poultry, fish):
nga_asf <- nga_food_groups |> 
  mutate(asf = ifelse(food_group %in% c("dairy", "eggs", "meat_poultry_fish"), 1, 0)) |> 
  group_by(hhid) |>
  summarise(asf = ifelse(sum(asf) > 0, 1, 0))

#-------------------------------------------------------------------------------

# Tidy environment: 
rm(list = setdiff(ls(), c("nga_base_ai", "nga_asf")))

################################################################################
############################### END OF SCRIPT ##################################
################################################################################




