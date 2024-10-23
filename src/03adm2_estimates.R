################################################################################
########### SCRIPT FOR GENERATING DIRECT AND SMOOTHED ADM2 ESTIMATES ###########
################################################################################

# Author: Mo Osman
# Date created: 19-09-2024
# Last edited: 22-10-2024

# In this script, I will generate both direct and smoothed ADM2 level estimates
# for inadequate intake of vitamin B12 in Nigeria. 

# This project uses pre-processed MIMI micro-nutrient intake data derived from:
# Nigeria Living Standards Survey (2018-19) & FAO/INFOODS West African food 
# composition table (2019)

# INSTALL AND LOAD PACKAGES:

rq_packages <- c("readr", "tidyverse", "srvyr", "sf", "spdep", "tmap", "INLA", 
                 "SUMMER", "wesanderson")

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
nigeria_1 <- st_read("shapefiles/nigeria_1")
nigeria_2 <- st_read("shapefiles/nigeria_2")

#-------------------------------------------------------------------------------

# CALCULATE NAIVE ESTIMATES AT THE ADM2 LEVEL: 
vb12_inad_lga <- nga_analysis_df %>% 
  group_by(lga) %>% 
  summarise(vb12_inad_naive = mean(vb12_inadequate, na.rm = TRUE),
            vb12_inad_naive_var = var(vb12_inadequate, na.rm = TRUE))

#-------------------------------------------------------------------------------

# CALCULATE DESIGN-BASED (DIRECT) ESTIMATES AT ADM2 LEVEL: 

# Firstly need to create a tbl_svy object to be used for analysis.
nga_analysis_df_svy <- nga_analysis_df %>% as_survey_design(ids = c("ea", "hhid"),
                                                            strata = "state",
                                                            weights = "survey_wgt",
                                                            nest = TRUE)

# Calculate vitamin B12 inadequacy stratified by LGA - include variance of estimates: 
vb12_inad_lga <- vb12_inad_lga %>% 
  left_join(nga_analysis_df_svy %>% 
              group_by(lga) %>% 
              summarise(vb12_inad_dir = survey_mean(vb12_inadequate, 
                                                    na.rm = TRUE, 
                                                    vartype = "var")), 
            by = "lga")
  

# Multiply vb12_inad by 100 to get percentage:
vb12_inad_lga$vb12_inad_naive <- vb12_inad_lga$vb12_inad_naive * 100
vb12_inad_lga$vb12_inad_dir <- vb12_inad_lga$vb12_inad_dir * 100

# Join spatial data: 
vb12_inad_lga <- nigeria_2 %>% 
  left_join(vb12_inad_lga, by = c("lga" = "lga"))

#-------------------------------------------------------------------------------

# MAP DIRECT ESTIMATES OF VITAMIN B12 INADEQUACY: 
plot_map <- function(data, col, title, metric, level) {
  
  # Create a map: 
  map <- tm_shape(data) + 
    tm_fill(col = col,
            title = metric, 
            style = "cont",
            breaks = seq(0, 100, by = 10),
            textNA = "Missing Data",
            legend.is.portrait = F,
            palette = wesanderson::wes_palette("Zissou1Continuous")) + 
    tm_layout(main.title = title, frame = F, main.title.size = 0.8, 
              main.title.position = "center", legend.outside.position = "bottom",
              legend.outside.size = 0.35) +
    tm_borders(lwd = 0) + 
    tm_legend(show = T) +
    tm_shape(nigeria_1) +
    tm_borders(col = "black", lwd = 0.8)
  
  return(map)
}

naive_map <- plot_map(data = vb12_inad_lga, 
                      col = "vb12_inad_naive", 
                      title = "Vitamin B12 (naive estimates)", 
                      metric = "Prevalence of inadequate intake (%)", 
                      level = "lga")

naive_map

# tmap_save(naive_map, "outputs/maps/vb12_inad_naive.png", width = 8, height = 8,
#           units = "in", dpi = 600)

direct_map <- plot_map(data = vb12_inad_lga, 
                       col = "vb12_inad_dir", 
                       title = "Vitamin B12 (design-based estimates)", 
                       metric = "Prevalence of inadequate intake (%)", 
                       level = "lga")

direct_map

# tmap_save(direct_map, "outputs/maps/vb12_inad_direct.png", width = 8, height = 8,
#           units = "in", dpi = 600)

rm(naive_map, direct_map)

#-------------------------------------------------------------------------------

# PREPARE INPUTS FOR BYM2 SMOOTHING MODEL: 

# Remove entries from nigeria_2 that are not in the sample - decision at present
# to only smooth estimates for areas that have observations. TO REVISIT: 
nigeria_2 <- nigeria_2[nigeria_2$lga %in% nga_analysis_df$lga,]

# Define spatial object:
geo <- as_Spatial(nigeria_2)

# Define spatial adjacency matrix: 
nb <- poly2nb(geo, queen = TRUE)
Amat <- nb2mat(nb, style = "B")
rownames(Amat) <- nigeria_2$lga
colnames(Amat) <- nigeria_2$lga

# Redfine nigeria_2 as a list containing the spatial object and the adjacency 
# matrix:
nigeria_2 <- list("geo" = geo, "Amat" = Amat)

# Tidy environment: 
rm(geo, Amat, nb)

# Convert from tbl_df to data.frame - This step is required for the fitGeneric()
# function to work: 
nga_analysis_df <- data.frame(nga_analysis_df)

#-------------------------------------------------------------------------------

# CALCULATE SMOOTHED ESTIMATES OF VITAMIN B12 INADEQUACY - BYM2 MODEL

# Fit generic smoothing model - naive to survey design:
smoothed_vb12_inad_naive <- fitGeneric(data = nga_analysis_df,
                                       Amat = nigeria_2$Amat,
                                       responseType = "binary",
                                       responseVar = "vb12_inadequate",
                                       regionVar = "lga",
                                       strataVar = NULL,
                                       weightVar = NULL,
                                       clusterVar = NULL,
                                       CI = 0.95)

head(smoothed_vb12_inad_naive$smooth)
head(smoothed_vb12_inad_naive$HT)
  
# # Fit generic smoothing model - acknowledging survey design:
# smoothed_vb12_inad <- fitGeneric(data = nga_analysis_df,
#                                  Amat = nigeria_2$Amat,
#                                  responseType = "binary",
#                                  responseVar = "vb12_inadequate",
#                                  regionVar = "lga",
#                                  strataVar = "state",
#                                  weightVar = "survey_wgt",
#                                  clusterVar = "~ea+hhid", 
#                                  nest = TRUE,
#                                  CI = 0.95)

# When specifying survey design, INLA program crashes: 
########### To revisit and identify issue ############

#-------------------------------------------------------------------------------

# MAP SMOOTHED ESTIMATES: 
vb12_inad_lga <- vb12_inad_lga %>% 
  left_join(smoothed_vb12_inad_naive$smooth %>% 
              dplyr::select(region, mean, var, lower, upper) %>% 
              mutate(vb12_inad_smoothed_naive = mean * 100,
                     var_smoothed_naive = var,
                     lower_smoothed_naive = lower * 100,
                     upper_smoothed_naive = upper * 100) %>% 
              rename(lga = region) %>% 
              dplyr::select(lga, vb12_inad_smoothed_naive, var_smoothed_naive, 
                            lower_smoothed_naive, upper_smoothed_naive), 
            by = "lga")

# Map smoothed (naive): 
vb12_smoothed_naive <- plot_map(data = vb12_inad_lga, 
                       col = "vb12_inad_smoothed_naive", 
                       title = "Vitamin B12 (naive smoothed estimates)", 
                       metric = "Prevalence of inadequate intake (%)", 
                       level = "lga")

vb12_smoothed_naive

# tmap_save(vb12_smoothed_naive, "outputs/maps/vb12_inad_smoothed_naive.png", width = 8, height = 8,
#           units = "in", dpi = 600)

#-------------------------------------------------------------------------------

# MAP LOWER AND UPPER BOUND OF CREDIBLE INTERVAL FOR SMOOTHED ESTIMATES: 

# Map lower bound of credible interval: 
vb12_smoothed_lower <- plot_map(data = vb12_inad_lga, 
                       col = "lower_smoothed_naive", 
                       title = "Vitamin B12 (smoothed estimates):\n lower bound of 95% credible interval", 
                       metric = "Prevalence of inadequate intake (%)", 
                       level = "lga")

vb12_smoothed_lower

# tmap_save(vb12_smoothed_lower, "outputs/maps/vb12_smoothed_lower.png", width = 8, height = 8,
#           units = "in", dpi = 600)

# Map upper bound of credible interval:
vb12_smoothed_upper <- plot_map(data = vb12_inad_lga, 
                       col = "upper_smoothed_naive", 
                       title = "Vitamin B12 (smoothed estimates):\n upper bound of 95% credible interval", 
                       metric = "Prevalence of inadequate intake (%)", 
                       level = "lga")

vb12_smoothed_upper

# tmap_save(vb12_smoothed_upper, "outputs/maps/vb12_smoothed_upper.png", width = 8, height = 8,
#           units = "in", dpi = 600)

#-------------------------------------------------------------------------------

# MAP VARIANCE ESTIMATES: 

# Firstly alter plot_map function to use a different colour palette (viridis): 

plot_map <- function(data, col, title, metric, level) {
  
  # Create a map: 
  map <- tm_shape(data) + 
    tm_fill(col = col,
            title = metric, 
            style = "cont",
            breaks = seq(0, 0.07, by = 0.01),
            textNA = "Missing Data",
            legend.is.portrait = F,
            palette = viridis::viridis(10)) + 
    tm_layout(main.title = title, frame = F, main.title.size = 0.8, 
              main.title.position = "center", legend.outside.position = "bottom",
              legend.outside.size = 0.35) +
    tm_borders(lwd = 0) + 
    tm_legend(show = T) +
    tm_shape(nigeria_1) +
    tm_borders(col = "black", lwd = 0.8)
  
  return(map)
}

# Variance of direct estimates: 
vb12_var_direct <- plot_map(data = vb12_inad_lga, 
                       col = "vb12_inad_dir_var", 
                       title = "Variance of direct HT estimates", 
                       metric = "Variance", 
                       level = "lga")

vb12_var_direct

# tmap_save(vb12_var_direct, "outputs/maps/var_direct.png", width = 8, height = 8,
#           units = "in", dpi = 600)

# Variance of smoothed estimates: 
vb12_var_smoothed <- plot_map(data = vb12_inad_lga, 
                       col = "var_smoothed_naive", 
                       title = "Variance of smoothed estimates", 
                       metric = "Variance", 
                       level = "lga")

vb12_var_smoothed

# tmap_save(vb12_var_smoothed, "outputs/maps/var_smoothed.png", width = 8, height = 8,
#           units = "in", dpi = 600)

#-------------------------------------------------------------------------------

# Clear environment: 
rm(list = ls())

################################################################################
############################# END OF SCRIPT ####################################
################################################################################
