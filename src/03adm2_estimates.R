################################################################################
########### SCRIPT FOR GENERATING DIRECT AND SMOOTHED ADM2 ESTIMATES ###########
################################################################################

# Author: Mo Osman
# Date created: 19-09-2024
# Last edited: 

# In this script, I will generate both direct and smoothed ADM2 level estimates.

# INSTALL AND LOAD PACKAGES:

rq_packages <- c("readr", "tidyverse", "srvyr", "sf", "spdep", "tmap", "INLA", 
                 "SUMMER")

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

# CALCULATE HORVITZ-THOMPSON (DIRECT) ESTIMATES AT ADM2 LEVEL: 

# Firstly need to create a tbl_svy object to be used for analysis.
nga_analysis_df_svy <- nga_analysis_df %>% as_survey_design(ids = c("ea", "hhid"),
                                                            strata = "state",
                                                            weights = "survey_wgt",
                                                            nest = TRUE)

# Calculate vitamin B12 inadequacy stratified by LGA - include variance of estimates: 
vb12_inad_lga <- nga_analysis_df_svy %>% 
  group_by(lga) %>% 
  summarise(vb12_inad_dir = survey_mean(vb12_inadequate, 
                                    na.rm = TRUE, 
                                    vartype = "var"))

# Multiply vb12_inad by 100 to get percentage:
vb12_inad_lga$vb12_inad_dir <- vb12_inad_lga$vb12_inad_dir * 100

# Join to spatial data: 
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

direct_map <- plot_map(data = vb12_inad_lga, 
                       col = "vb12_inad_dir", 
                       title = "Vitamin B12 (direct estimates)", 
                       metric = "Prevalence of inadequate intake (%)", 
                       level = "lga")

direct_map

# tmap_save(direct_map, "outputs/maps/vb12_inad_direct.png", width = 8, height = 8, 
#           units = "in", dpi = 600)

rm(direct_map)

#-------------------------------------------------------------------------------

# CALCULATE SMOOTHED ESTIMATES OF VITAMIN B12 INADEQUACY - BYM MODEL

# Remove entries from nigeria_2 that are not in the sample - decision at present
# to only smooth estimates for areas that have observations. Can revisit later: 
nigeria_2 <- nigeria_2[nigeria_2$lga %in% nga_analysis_df$lga,]

# Define spatial object:
geo <- as_Spatial(nigeria_2)

# Define spatial adjacency matrix: 
Amat <- poly2nb(geo, queen = TRUE)
Amat <- nb2mat(Amat, style = "B")
rownames(Amat) <- nigeria_2$lga
colnames(Amat) <- nigeria_2$lga

nigeria_2 <- list("geo" = geo, "Amat" = Amat)

# Tidy environment: 
rm(geo, Amat)

# Convert from tbl_df to data.frame - IMPORTANT, DO NOT REMOVE: 
nga_analysis_df <- data.frame(nga_analysis_df)

# Fit generic smoothing model:
smoothed_vb12_inad <- fitGeneric(data = nga_analysis_df,
                                 Amat = nigeria_2$Amat,
                                 responseType = "binary",
                                 responseVar = "vb12_inadequate",
                                 regionVar = "lga",
                                 strataVar = "state",
                                 weightVar = "survey_wgt",
                                 clusterVar = NULL,
                                 nest = TRUE,
                                 CI = 0.95, 
                                 verbose = TRUE)

head(smoothed_vb12_inad$smooth)

#-------------------------------------------------------------------------------

# MAP SMOOTHED ESTIMATES: 
vb12_smoothed <- smoothed_vb12_inad$smooth %>% 
  dplyr::select(region, mean, var, lower, upper) %>% 
  mutate(mean = mean * 100,
         lower = lower * 100, 
         upper = upper * 100) %>% 
  rename(lga = region,
         vb12_inad_smoothed = mean)


vb12_smoothed <- vb12_inad_lga %>%
  left_join(vb12_smoothed, by = "lga")



vitb12_lga_smooth <- plot_map(data = vb12_smoothed, 
                       col = "vb12_inad_smoothed", 
                       title = "Vitamin B12 (smoothed estimates)", 
                       metric = "Prevalence of inadequate intake (%)", 
                       level = "lga")

vitb12_lga_smooth

# tmap_save(vitb12_lga_smooth, "outputs/maps/vb12_inad_smoothed.png", width = 8, height = 8, 
#           units = "in", dpi = 600)

# Map lower bound of credible interval: 
vb12_smoothed_lower <- plot_map(data = vb12_smoothed, 
                       col = "lower", 
                       title = "Vitamin B12 (smoothed estimates):\n lower bound of 95% credible interval", 
                       metric = "Prevalence of inadequate intake (%)", 
                       level = "lga")

vb12_smoothed_lower

# tmap_save(vb12_smoothed_lower, "outputs/maps/vb12_smoothed_lower.png", width = 8, height = 8,
#           units = "in", dpi = 600)

# Map upper bound of credible interval:
vb12_smoothed_upper <- plot_map(data = vb12_smoothed, 
                       col = "upper", 
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
            breaks = seq(0, 0.062, by = 0.01),
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
vb12_var_smoothed <- plot_map(data = vb12_smoothed, 
                       col = "var", 
                       title = "Variance of smoothed estimates", 
                       metric = "Variance", 
                       level = "lga")

vb12_var_smoothed

# tmap_save(vb12_var_smoothed, "outputs/maps/var_smoothed.png", width = 8, height = 8,
#           units = "in", dpi = 600)
