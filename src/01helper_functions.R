# FUNCTIONS TO COMPILE DATA FOR HOUSEHOLD APPARENT MICRONUTRIENT INTAKE:

#-------------------------------------------------------------------------------

# INSTALL AND LOAD PACKAGES:

rq_packages <- c("readr", "tidyverse", "here")

installed_packages <- rq_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(rq_packages[!installed_packages])
}

lapply(rq_packages, require, character.only = T)

rm(list= c("rq_packages", "installed_packages"))

#-------------------------------------------------------------------------------

path_to_file <- here::here("processed_data/")

allen_ear <- data.frame(
  # Data frame that specifies WHO E-AR values for each micronutrient:
  nutrient = c(
    "energy_kcal",
    "vita_rae_mcg",
    "thia_mg",
    "ribo_mg",
    "niac_mg",
    "vitb6_mg",
    "folate_mcg",
    "vitb12_mcg",
    "fe_mg",
    "ca_mg",
    "zn_mg"
  ),
  ear_value = c(
    2100,
    490, 
    0.9,
    1.3, 
    11, 
    1.3, 
    250, 
    2, 
    22.4, #low absorption
    860, 
    10.2# unrefined
  )
)

#-------------------------------------------------------------------------------

read_in_survey <- function(name_of_survey, path_to_file = here::here("processed_data/")){
  # given the name of the survey of country
  # the function reads in relevant data required for apparent intake calculations:
  
  hh_info <<-  read.csv(paste0(path_to_file, paste0(name_of_survey, "_hh_info.csv")))
  food_consumption<<- read.csv(paste0(path_to_file, paste0(name_of_survey, "_food_consumption.csv")))
  fc_table <<- read.csv(paste0(path_to_file, paste0(name_of_survey, "_fct.csv")))
}

#-------------------------------------------------------------------------------

apparent_intake <- function(name_of_survey, path_to_file = here::here("processed_data/")){
  # Estimates apparent intake of nutrients based on: 
  # 1 - quantity consumed of food items, 
  # 2 - nutrient composition of those food items,
  # 3 - adult female equivalent unit of the household
  read_in_survey(name_of_survey, path_to_file)
  
  
  # If NLSS survey then need to read in zone to do fct matches for milk by zone:
  if (name_of_survey == "nga_lss1819"){
    
    # Read in zone data - NB: only required for Nigeria
    cover <-  read.csv("nlss_raw/secta_cover.csv")
    
    # Left join zone to hh_info:
    hh_info <- hh_info %>% 
      dplyr::left_join(cover %>% dplyr::select(hhid, zone), by = "hhid")
    
    # Filter fc_table for milk:
    milk_fct <- fc_table %>% 
      filter(item_code == 110) %>% 
      filter(!is.na(zone))
    
    # Apparent nutrient intake from milk: 
    milk_ai <- food_consumption %>% 
      filter(item_code == 110) %>% # Filter to include only fresh milk
      left_join(hh_info %>% select(hhid, zone), 
                by = "hhid") %>%
      left_join(milk_fct, 
                by = c("item_code", "zone")) %>% # Join milk food composition table by item_code and zone
      mutate(across(-c(item_code, hhid, item_name, food_group, quantity_100g, 
                       quantity_g, zone),
                    ~.x*quantity_100g)) %>% # Multiply nutrient values by quantity consumed
      group_by(hhid) %>% # Aggregate by household id summing the values of consumption
      summarise(across(-c(item_code, item_name, quantity_100g, quantity_g, 
                          food_group, zone),
                       ~sum(., na.rm = T))) %>% 
      left_join(hh_info %>% select(hhid, afe), by = "hhid") %>% # Join afe
      mutate(across(-c(hhid, afe),~.x/afe)) %>% # Divide all nutrient values by afe
      ungroup()
    
    # Remove milk from the main food composition table: 
    fc_table <- fc_table %>% 
      filter(item_code != 110) %>% 
      dplyr::select(-zone)
    
    # Apparent nutrient intake from all other foods:
    base_ai <- food_consumption %>% 
      filter(item_code != 110) %>%
      left_join(fc_table, by = "item_code") %>% 
      mutate(
        across(
          -c(item_code, hhid,item_name ,food_group, quantity_100g, quantity_g),
          ~.x*quantity_100g
        )
      ) %>% 
      group_by(hhid) %>% 
      summarise(
        across(-c(item_code,item_name,quantity_100g,quantity_g, food_group),
               ~sum(.,na.rm = T))
      ) %>% 
      left_join(hh_info %>% select(hhid, afe), by = "hhid") %>% 
      mutate(
        across(
          -c(hhid,afe),
          ~.x/afe
        )
      ) %>% 
      ungroup()
    
    # Combine nutrient intake from the 2 apparent intake data-frames: 
    base_ai <- bind_rows(base_ai, milk_ai) %>% 
      group_by(hhid) %>%
      summarise(across(everything(), ~sum(., na.rm = T))) %>%
      ungroup() %>% 
      dplyr::select(-afe) 
    
    base_ai
  }
  
  else{
    x <- food_consumption %>% 
      left_join(fc_table, by = "item_code") %>% 
      mutate(
        across(
          -c(item_code, hhid,item_name , quantity_100g, quantity_g),
          ~.x*quantity_100g
        )
      ) %>% 
      group_by(hhid) %>% 
      summarise(
        across(-c(item_code,item_name,quantity_100g,quantity_g),
               ~sum(.,na.rm = T))
      ) %>% 
      left_join(hh_info %>% select(hhid, afe), by = "hhid") %>% 
      mutate(
        across(
          -c(hhid,afe),
          ~.x/afe
        )
      ) %>% 
      ungroup() %>% 
      select(-afe)
    x
  }
}

################################################################################
############################### END OF SCRIPT ##################################
################################################################################
