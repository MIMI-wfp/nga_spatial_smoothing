################################################################################
########### SCRIPT FOR ESTIMATING ADM1 LEVEL FOOD GROUP CONSUMPTION ############
################################################################################

# Author: Mo Osman
# Date created: 31-12-2024
# Last edited: 

# In this script, I produce aggregated estimats of food group consumption at the
# ADM1 level (for food groups of interest). I will map these estimates, and these
# may be used for model validation. 

# INSTALL AND LOAD PACKAGES:

rq_packages <- c("readr", "tidyverse", "srvyr", "sf", "spdep", "tmap")

installed_packages <- rq_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(rq_packages[!installed_packages])
}

lapply(rq_packages, require, character.only = T)

rm(list= c("rq_packages", "installed_packages"))

#-------------------------------------------------------------------------------

# READ IN REQUIRED DATA: 

# Food group consumption data:
nga_food_groups <- read_csv("processed_data/nga_lss1819_food_consumption.csv")

# Nigeria ADM1 shapefile:
nigeria_1 <- st_read("shapefiles/nigeria_1")

# Household survey information data and state: 
nga_hh_info <- read_csv("processed_data/nga_lss1819_hh_info.csv") |> 
  left_join(read_csv("shapefiles/household_locations.csv") |> 
              dplyr::select(hhid, state),
            by = "hhid") |> 
  dplyr::select(hhid, state, ea, survey_wgt)

#-------------------------------------------------------------------------------

# BINARISE HOUSEHOLD FOOD GROUP CONSUMPTION:

# Animal-sourced food consumption (dairy, eggs, meat, poultry, fish):
nga_asf <- nga_food_groups |> 
  mutate(asf = ifelse(food_group %in% c("dairy", "eggs", "meat_poultry_fish"), 1, 0)) |> 
  group_by(hhid) |>
  summarise(asf = ifelse(sum(asf) > 0, 1, 0))

#-------------------------------------------------------------------------------

# ESTIMATE ADM1 LEVEL FOOD GROUP CONSUMPTION:: 

# Animal-sourced food: 
nga_asf <- nga_asf |> 
  left_join(nga_hh_info, by = "hhid")

# Create tbl_svy object for analysis: 
nga_asf <- nga_asf |> 
  as_survey_design(ids = c("ea", "hhid"),
                   strata = "state",
                   weights = "survey_wgt", 
                   nest = TRUE)

# Estimate prevalence of animal-sourced food consumption at the ADM1 level:
nga_asf <- nga_asf |> 
  group_by(state) |> 
  summarise(asf_consumption = survey_mean(asf, 
                                          na.rm = TRUE,
                                          vartype = "ci")) |> 
  # Multiply all values by 100: 
  mutate(across(where(is.numeric), ~ .x * 100))

# Join spatial data: 
nga_asf <- nigeria_1 |> 
  left_join(nga_asf, by = "state")

#-------------------------------------------------------------------------------

# MAP RESULTS:

plot_map <- function(data, col, title, metric) {
  
  # Create a map: 
  map <- tm_shape(data) + 
    tm_fill(col = col,
            title = metric, 
            style = "cont",
            breaks = seq(50, 100, by = 10),
            textNA = "Missing Data",
            legend.is.portrait = F,
            palette = "YlOrBr") + 
    tm_layout(main.title = title, frame = F, main.title.size = 0.8, 
              main.title.position = "center", legend.outside.position = "bottom",
              legend.outside.size = 0.35) +
    tm_borders(lwd = 0) + 
    tm_legend(show = T) +
    tm_shape(nigeria_1) +
    tm_borders(col = "black", lwd = 0.8)
  
  return(map)
}

nga_asf <- plot_map(data = nga_asf, 
                     col = "asf_consumption",
                     title = "Animal-sourced Food Consumption", 
                     metric = "% of households consumed")

# Save map:
tmap_save(nga_asf, filename = "outputs/maps/food_group_consumption/asf_adm1.png", 
          width = 10, height = 10, dpi = 300)

#-------------------------------------------------------------------------------

# Clear environment: 
rm(list = ls())

################################################################################
############################### END OF SCRIPT ##################################
################################################################################
