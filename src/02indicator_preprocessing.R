################################################################################
############# SCRIPT FOR PRE-PROCESSING INDICATORS BEFORE ANALYSIS #############
################################################################################

# Author: Mo Osman
# Date created: 19-09-2024
# Last edited: 17-04-2024

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
source("src/01helper_functions.R")

#-------------------------------------------------------------------------------

# Get base case apparent intake for each MN: 
nga_base_ai <- apparent_intake("nga_lss1819")

# Tidy environment: 
rm(list = setdiff(ls(), c("nga_base_ai", "allen_ear")))

#-------------------------------------------------------------------------------

# Binarise Vitamin B12 and Folate inadequacy(Inadequate = 1, adequate = 0): 
nga_base_ai$folate_inadequate <- ifelse(nga_base_ai$folate_mcg < allen_ear$ear_value[allen_ear$nutrient == "folate_mcg"], 1, 0)
nga_base_ai$vb12_inadequate <- ifelse(nga_base_ai$vitb12_mcg < allen_ear$ear_value[allen_ear$nutrient == "vitb12_mcg"], 1, 0)

nga_base_ai <- nga_base_ai %>% 
  dplyr::select(hhid, folate_inadequate, vb12_inadequate)


# For purpose of this proof-of-concept study we will focus on these micronutrients, 
# however these methods can be replicated for any micronutrient.

# Tidy environment: 
rm(list = setdiff(ls(), c("nga_base_ai")))

################################################################################
############################### END OF SCRIPT ##################################
################################################################################




