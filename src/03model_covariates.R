################################################################################
###################### SCRIPT FOR PREPARING MODEL COVARIATES ###################
################################################################################

# Author: Mo Osman
# Date created: 11-11-2024
# Last edited: 12-11-2024

# In this script, I will prepare a dataframe of covariates for spatial modelling.

# INSTALL AND LOAD PACKAGES:

rq_packages <- c("readr", "tidyverse")

installed_packages <- rq_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(rq_packages[!installed_packages])
}

lapply(rq_packages, require, character.only = T)

rm(list= c("rq_packages", "installed_packages"))

#-------------------------------------------------------------------------------

# READ IN REQUIRED DATA

# household_info:
hh_info <- read_csv("processed_data/nga_lss1819_hh_info.csv")

# labour module:
labour <- read_csv("nlss_raw/sect4a1_labour.csv")

#-------------------------------------------------------------------------------

# Select required variables from the pre-processed hh_info dataframe: 
hh_info <- hh_info |>  
  dplyr::select(hhid, 
                sex_head, # Sex of head of household
                sep_quintile, # socio-economic quintile (1 = Poorest, 5 = Wealthiest)
                educ_head, # Highest education level of head of household
                res) # Household residence (urban or rural)

# educ_head variable requires feature engineering - transform into a factor 
# variable with levels:
# 1 = No school education, 2 = Primary education, 3 = Secondary education, 
# 4 = Higher education. 
# Requirement for the individual to have completed that education level for it 
# to be counted.
hh_info <- hh_info |>  
  mutate(educ_head = case_when(
    educ_head %in% c(0, 1, 2, 3, 11, 12, 13, 14, 15) ~ 1,
    educ_head %in% c(16, 21, 22, 23, 24, 25, 51, 52) ~ 2, 
    educ_head %in% c(26, 27, 28, 321) ~ 3,
    educ_head %in% c(31, 33, 34, 35, 41, 43, 61, 322, 411, 412, 421, 422, 423, 424) ~ 4,
    TRUE ~ NA_real_
  ))

# Note: a lot of missingness in the educ_head variable (n=4293)

#-------------------------------------------------------------------------------

# Create a variable that indicates if any member of the household has engaged in
# their own agricultural work:
# 0 = No, 1 = Yes

labour <- labour |> 
  group_by(hhid) |> 
  mutate(ag_work = ifelse(any(s04aq06 == 1, na.rm = TRUE), 1, 0)) |> 
  ungroup () |> 
  dplyr::select(hhid, ag_work) |> 
  distinct()

#-------------------------------------------------------------------------------

# Join dataframes: 
covariates <- hh_info |> 
  left_join(labour, by = "hhid")

rm(hh_info, labour)

# To generate this covariates dataframe in another script, run the following line 
# of code: 
# source("src/03model_covariates.R")

################################################################################
############################## END OF SCRIPT ###################################
################################################################################

