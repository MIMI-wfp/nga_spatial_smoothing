## FUNCTIONS TO COMPILE BASE MODELS

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
    2100,#who
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
  # the function reads in each part of the base model into general 
  # object names
  
  hh_info <<-  read.csv(paste0(path_to_file, paste0(name_of_survey, "_hh_info.csv")))
  food_consumption<<- read.csv(paste0(path_to_file, paste0(name_of_survey, "_food_consumption.csv")))
  fc_table <<- read.csv(paste0(path_to_file, paste0(name_of_survey, "_fct.csv")))
}

#-------------------------------------------------------------------------------

full_item_list <- function(name_of_survey, path_to_file = here::here("processed_data/")){
  # creates a data frame with a full list of food items for every
  # household. If food item is not consumed, quantity = 0
  # uesful for food group analyses
  
  hh_info <-   read.csv(paste0(path_to_file, paste0(name_of_survey, "_hh_info.csv")))
  food_consumption<- read.csv(paste0(path_to_file, paste0(name_of_survey, "_food_consumption.csv")))
  fc_table <- read.csv(paste0(path_to_file, paste0(name_of_survey, "_fct.csv")))
  
  x <- hh_info %>% 
    select(hhid,afe) %>% 
    cross_join(fc_table %>% 
                 select(item_code)) %>% 
    left_join(food_consumption %>% 
                group_by(hhid, item_code, food_group) %>% 
                summarise(across(
                  everything(),
                  ~sum(., na.rm = TRUE)
                )) %>% 
                ungroup(), 
              by = c("hhid", "item_code")) %>% 
    select(-food_group) %>% 
    mutate(
      across(
        c(quantity_100g, quantity_g),
        ~ifelse(is.na(.),0, .)
      )
    ) %>% 
    mutate(
      quantity_100g = quantity_100g/afe, 
      quantity_g = quantity_g/afe
    ) %>% 
    left_join(fc_table, by = "item_code") %>% 
    inner_join(food_consumption %>% 
                 select(item_code, food_group) %>% 
                 distinct(item_code, food_group),
               by = c('item_code')) %>% 
    dplyr::mutate(
      across(
        -c(hhid,afe, item_code,quantity_100g,quantity_g,item_name, food_group),
        ~.x*quantity_100g
      )
    )
  
  x
}

#-------------------------------------------------------------------------------

apparent_intake <- function(name_of_survey, path_to_file = here::here("processed_data/")){
  # Estimates apparent intake of nutrients based on consumed food items
  # and adult female equivalent unit of the household
  read_in_survey(name_of_survey, path_to_file)
  
  
  # If NLSS survey then need to read in zone to do fct matches for milk by zone:
  if (name_of_survey == "nga_lss1819"){
    
    # Read in zone data:
    # TODO :: fix this so we don't have to read in the data within the function
    cover <-  read.csv("processed_data/secta_cover.csv")
    
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
      ungroup() %>% 
      select(-afe)
    x
  }
}

#-------------------------------------------------------------------------------

household_data <- function(name_of_survey, path_to_file = here::here("processed_data/")){
  #reads in the household information data
  x <- read.csv(paste0(path_to_file, paste0(name_of_survey, "_hh_info.csv")))
  x
}

#-------------------------------------------------------------------------------

nutrient_density <- function(name_of_survey){
  # returns a data frame of nutrient density for each household
  # values are given in unit of mn per 1000kcal 
  #
  x <- food_consumption %>% 
    left_join(fc_table, by = "item_code") %>% 
    mutate(
      across(
        -c(item_code, hhid, item_name, food_group, quantity_100g, quantity_g),
        ~.x*quantity_100g
      )
    ) %>% 
    group_by(hhid) %>% 
    summarise(
      across(-c(item_code,item_name,quantity_100g,quantity_g, food_group),
             ~sum(.,na.rm = T))
    ) %>% 
    mutate(energy_1000kcal = energy_kcal/1000) %>% 
    mutate(
      across(
        -c(hhid),
        ~.x/energy_1000kcal,
        .names ="{.col}_1000kcal"
      )
    ) %>% 
    select(hhid, ends_with("1000kcal")) %>% 
    select(-energy_kcal_1000kcal)
  x
}

#-------------------------------------------------------------------------------

target_creation <- function(){
  eth_ess1819 <- apparent_intake("eth_ess1819")
  # eth_hices1516 <- apparent_intake("eth_hices1516")
  nga_lss1819 <- apparent_intake("nga_lss1819")
  ind_nss1112 <- apparent_intake("ind_nss1112")
  
  select_and_append <- function(survey, survey_id){
    survey <- survey %>% 
      select(hhid,vita_rae_mcg,folate_mcg,vitb12_mcg,
             fe_mg,zn_mg) %>% 
      rename(
        va_ai = vita_rae_mcg,
        fo_ai = folate_mcg,
        vb12_ai = vitb12_mcg,
        fe_ai = fe_mg,
        zn_ai = zn_mg
      ) %>% 
      mutate(
        hhid = as.character(hhid),
        va_ref = 
          allen_ear$ear_value[allen_ear$nutrient == "vita_rae_mcg"],
        fo_ref = allen_ear$ear_value[allen_ear$nutrient == "folate_mcg"],
        vb12_ref = allen_ear$ear_value[allen_ear$nutrient == "vitb12_mcg"],
        fe_ref =  allen_ear$ear_value[allen_ear$nutrient == "fe_mg"],
        zn_ref = allen_ear$ear_value[allen_ear$nutrient == "zn_mg"],
        va_nar = ifelse(va_ai<=va_ref, va_ai/va_ref,1),
        fo_nar = ifelse(fo_ai<=fo_ref, fo_ai/fo_ref,1),
        vb12_nar = ifelse(vb12_ai<=vb12_ref, vb12_ai/vb12_ref,1),
        fe_nar = ifelse(fe_ai<=fe_ref, fe_ai/fe_ref,1),
        zn_nar = ifelse(zn_ai<=zn_ref, zn_ai/zn_ref,1),
        mimi_simple = (va_nar+fo_nar+vb12_nar+fe_nar+zn_nar)/5,
        survey_id = survey_id
      ) %>% 
      select(-c(va_nar,fo_nar,vb12_nar,zn_nar,fe_nar))
    survey
  }
  x <-  select_and_append(eth_ess1819,"ETH_2018_ESS_v03_M") %>% 
    bind_rows(select_and_append(nga_lss1819,"NGA_2018_LSS_v01_M")) %>% 
    bind_rows(select_and_append(ind_nss1112,"DDI-IND-MOSPI-NSSO-68Rnd-Sch2.0-July2011-June2012"))
  x
}

