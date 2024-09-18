################################################################################
########### SCRIPT FOR GENERATING DIRECT AND SMOOTHED ADM2 ESTIMATES ###########
################################################################################

# Author: Mo Osman
# Date created: 17-09-2024
# Last edited: 

# In this script, I will generate both direct and smoothed ADM2 level estimates.

# INSTALL AND LOAD PACKAGES:

rq_packages <- c("readr", "tidyverse", "srvyr", "SUMMER", "sf", "spdep", "tmap",
                 "INLA")

installed_packages <- rq_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(rq_packages[!installed_packages])
}

lapply(rq_packages, require, character.only = T)

rm(list= c("rq_packages", "installed_packages"))

#-------------------------------------------------------------------------------

# READ AND MERGE DATA:

# Read in pre-processed indicators:
source("src/02indicator_preprocessing.R")

# Read in household locations that are matched to shapefiles: 
household_locations <- read_csv("shapefiles/household_locations.csv")

# Get further household info:
hh_info <- read_csv("processed_data/nga_lss1819_hh_info.csv")

# Merge data: 
nga_analysis_df <- nga_base_ai %>% 
  left_join(household_locations, by = "hhid") %>% 
  left_join(hh_info %>% dplyr::select(hhid, ea, survey_wgt), by = "hhid")

rm(nga_base_ai, household_locations, hh_info)

#-------------------------------------------------------------------------------

# READ SHAPEFILES: 

nigeria_0 <- st_read("shapefiles/nigeria_0")
nigeria_1 <- st_read("shapefiles/nigeria_1")
nigeria_2 <- st_read("shapefiles/nigeria_2")

#-------------------------------------------------------------------------------

# CREATE ADM2 SPATIAL ADJACENCY MATRIX: 

# Firstly extract vector of LGA names: 
row.names <- nigeria_2$lga

# Create neighbours list: 
nigeria_2_neigh <- poly2nb(nigeria_2, row.names = row.names, queen = TRUE)

# Use neightbours list to create adjacency matrix:
nigeria_2_adj <- nb2mat(nigeria_2_neigh, style = "B")

# Set column names of adjacency matrix to equal row names:
colnames(nigeria_2_adj) <- row.names

# Tidy environment:
rm(nigeria_2_neigh, row.names)

#-------------------------------------------------------------------------------

# CALCULATE HORVITZ-THOMPSON (DIRECT) ESTIMATES AT ADM2 LEVEL: 

# Firstly need to create a tbl_svy object to be used for analysis.
nga_analysis_df_svy <- nga_analysis_df %>% as_survey_design(ids = c("ea", "hhid"),
                                                            strata = "state",
                                                            weights = "survey_wgt",
                                                            nest = TRUE)

# Calculate vitamin B12 inadequacy stratified by LGA - include variance of estimates: 
vb12_inad_lga <- nga_analysis_df_svy %>% 
  group_by(lga) %>% 
  summarise(vb12_inad = survey_mean(vb12_inadequate, 
                                    na.rm = TRUE, 
                                    vartype = "var"))

#-------------------------------------------------------------------------------

# CALCULATE SMOOTHED ESTIMATES AT ADM2 LEVEL:
smoothed_vb12_inad <- fitGeneric(data = nga_analysis_df,
                                 geo = nigeria_2$geometry, 
                                 Amat = nigeria_2_adj,
                                 responseType = "binary",
                                 responseVar = "vb12_inadequate",
                                 regionVar = "lga",
                                 strataVar = "state",
                                 weightVar = "survey_wgt",
                                 clusterVar = "~ea+hhid",
                                 CI = 0.95)

# This function is currently not working due to the fact that some LGAs are not 
# included in the sample and therefore do not exist in the analysis df. Could 
# generate smoothed estimates ignoring these, but ideally would want to see if we
# can generate estimates for non sampled areas (TO REVISIT TOMORROW).
