################################################################################
########### SCRIPT FOR GENERATING MODELLED (UNIT-LEVEL) ADM2 ESTIMATES #########
################################################################################

# Author: Sahoko Ishida
# Date created: 06-04-2025
# Last edited: 28-07-2025

# In this script, I will generate both smoothed ADM2 level estimates
# for inadequate intake of vitamin B12 in Nigeria. 

# This project uses pre-processed MIMI micro-nutrient intake data derived from:
# Nigeria Living Standards Survey (2018-19) & FAO/INFOODS West African food 
# composition table (2019)

# INSTALL AND LOAD PACKAGES:
rq_packages <- c("readr", "readxl","tidyverse","survey","srvyr", "sf", "spdep", "tmap", "INLA", 
                 "SUMMER", "wesanderson", "cmdstanr", "ggplot2", "ggrepel","ggpubr")

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
source("src/helper_functions.R")

# Read in household locations that are matched to shapefiles: 
household_locations <- read_csv("shapefiles/household_locations.csv")

# Get further household info:
hh_info <- read_csv("processed_data/nga_lss1819_hh_info.csv")

# Merge data: 
nga_analysis_df <- nga_base_ai %>% 
  left_join(household_locations, by = "hhid") %>% 
  left_join(hh_info %>% dplyr::select(hhid, ea, survey_wgt), by = "hhid")

rm(nga_base_ai, household_locations, hh_info)

# create a tbl_svy object to be used for analysis.
nga_analysis_df_svy <- nga_analysis_df %>% as_survey_design(ids = c("ea", "hhid"),
                                                            strata = "state",
                                                            weights = "survey_wgt",
                                                            nest = TRUE)

#-------------------------------------------------------------------------------
# READ SHAPEFILES: 
nigeria_1 <- st_read("shapefiles/nigeria_1")
nigeria_2 <- st_read("shapefiles/nigeria_2")

#-------------------------------------------------------------------------------
# main indicator 
nutrient <- "folate"
runname = "20250721-01"

sae_df <- function(df, nutrient, geo_level,geo_level_upper, add_phantom = FALSE){
  nutrient_colname <- paste0(nutrient,"_inadequate")
  df_aggregates <- df |> group_by(pick(c(geo_level_upper,geo_level))) |>
    summarise(
      naive_mean = mean(.data[[nutrient_colname]], na.rm = TRUE),
      naive_var = var(.data[[nutrient_colname]], na.rm = TRUE),
    )
  df_svy <- df |> 
    as_survey_design(ids = c("ea", "hhid"),
                     strata = all_of(geo_level_upper),
                     weights = "survey_wgt",
                     nest = TRUE)
  
  df_aggregates <- df_aggregates |>
    left_join(df_svy |>
                group_by(pick(geo_level)) |>
                summarise(
                  design_based = survey_mean(.data[[nutrient_colname]],
                                             na.rm = TRUE, vartype = "var")
                ) |>
                mutate(design_based_mean = design_based) |> select(-design_based)
    ) |>
    left_join(df |> group_by(pick(c(geo_level,"ea"))) |>
                summarise(n = n()) |> group_by(pick(geo_level)) |>
                summarise(degf = n()-1, n_obs = sum(n)) |> ungroup()
    )
  
  df_cluster_level <- df |> group_by(pick(c(geo_level_upper, geo_level, "ea"))) |> 
    summarise(cluster_mean = mean(.data[[nutrient_colname]], na.rm = TRUE),
              cluster_count = sum(.data[[nutrient_colname]], na.rm = TRUE),
              cluster_weight = sum(survey_wgt),
              #cluster_weight = mean(survey_wgt),
              n_hh = n(),
              n = sum(!is.na(.data[[nutrient_colname]]))
    ) |> ungroup() |>
    mutate(ea = as.character(ea))
  
  if (add_phantom==TRUE){
    df_geo_phantom <- df_cluster_level |>  group_by(pick(c(geo_level_upper, geo_level))) |>
      summarise(n_values = n_distinct(cluster_mean)) |>
      filter(n_values == 1)|>
      select(all_of(c(geo_level_upper,geo_level)))
    geo_upper_list <- unique(df_geo_phantom[[geo_level_upper]])
    i = 0
    for (each_upper in geo_upper_list){
      geo_phantom_list <- df_geo_phantom[[geo_level]][df_geo_phantom[[geo_level_upper]]==each_upper]
      df_tmp_to_add <- df_cluster_level |> 
        filter(.data[[geo_level_upper]]==each_upper) |> 
        group_by(pick(geo_level_upper)) |>
        summarise(cluster_mean = weighted.mean(cluster_mean,cluster_weight,na.rm = TRUE),
                  cluster_weight = sum(cluster_weight, na.rm = TRUE)/sum(n_hh),
                  n_hh = 1
                  #!!geo_level := each_geo
                  #ea = paste0("phantom",as.character(i))
        )
      for (each_geo in geo_phantom_list){
        i = i + 1
        df_cluster_level <- bind_rows(df_cluster_level,
                                      (df_tmp_to_add |> 
                                         mutate(!!geo_level := each_geo,
                                                ea = paste0("phantom",as.character(i)))
                                      )
        )
      }
    }
    df_cluster_level_svy <- df_cluster_level |> 
      as_survey_design(ids = c("ea"),
                       strata = all_of(geo_level_upper),
                       weights = "cluster_weight",
                       nest = TRUE)
    
    df_aggregates_tmp <- df_cluster_level_svy |>
      group_by(pick(geo_level)) |>
      summarise(design_based_ph = survey_mean(cluster_mean,na.rm = TRUE, vartype = "var"),
                degf_ph = n() - 1, n_obs_ph = sum(n_hh)) |>
      rename(design_based_ph_mean = design_based_ph)
    
    df_aggregates <- df_aggregates |> left_join(df_aggregates_tmp, by = geo_level)
    
  }
  return((list(df_aggregates, df_cluster_level)))
  #return(df_aggregates)
}

df_process <- sae_df(nga_analysis_df, nutrient=nutrient, geo_level="lga", geo_level_upper = "state", add_phantom = F)
inad_lga <- df_process[[1]]
inad_cluster <- df_process[[2]]


# Give adm2_index
nigeria_2$lga_id <- 1:nrow(nigeria_2)
# Join spatial data: 
inad_lga <- nigeria_2 %>% 
  left_join(inad_lga, by = c("lga" = "lga"))
inad_cluster <- nigeria_2 |>
  left_join(inad_cluster,by = c("lga" = "lga")) |> filter(lga%in%inad_cluster$lga)

# Prepare data for analysis
geo <- as_Spatial(nigeria_2)

# Define spatial adjacency matrix: 
nb <- poly2nb(geo, queen = TRUE)
Amat <- nb2mat(nb, style = "B")
rownames(Amat) <- nigeria_2$lga
colnames(Amat) <- nigeria_2$lga
prep <- prepare_bym2(Amat)
## quality check
qc = rep(0,length(prep$n1))
for (i in 1:length(prep$n1)){
  qc[i] = Amat[prep$n1[i],prep$n2[i]]
}
mean(qc) # should be 1

# inad_lga_complete <- inad_lga |> filter(!is.na(naive_mean),degf_ph!=0)
# inad_lga_complete$design_based_ph_var |> summary()
# # inad_lga[!inad_lga$lga_id%in%inad_lga_complete$lga_id,] |> view()

model_name = "cluster_level_BYM2_BetaBin" 
stanfile <- here('src',paste0(model_name,'.stan'))
mod <- cmdstan_model(stanfile)
data_list <- list(
  N = nrow(inad_lga),
  Nea = nrow(inad_cluster),
  adm2_index = inad_cluster$lga_id,
  n = inad_cluster$n,
  y = inad_cluster$cluster_count,
  N_edges = length(prep$n1),
  node1 = prep$n1,
  node2 = prep$n2,
  scaling_factor = prep$scaling_factor
)
# MCMC
folder_path = here("outputs","mcmc",nutrient, model_name,runname)
if (!dir.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}
fit <- mod$sample(
  data = data_list,
  seed = 4321,
  iter_warmup = 500,
  iter_sampling = 500,
  save_warmup = TRUE,
  chains = 2,
  parallel_chains = 2,
  refresh = 50 # print update every 50 iters
)

post_samples = fit$draws(format = "df",  inc_warmup = F)
post_samples$chain = as.character(post_samples$.chain)
colnames(post_samples) <- colnames(post_samples) %>%
  gsub("\\[", "_", .) %>%
  gsub("\\]", "", .) %>%
  gsub(",","_", .)
#params <- c('sigma_u',"rho","sigma_tau","gamma0","gamma1","gamma2")
params <-c('b0','sigma_u',"rho","d")
for (param in params){
  dat = post_samples 
  gg = ggplot(data=dat, aes_string(x=".iteration", y = param, color="chain"))+
    geom_line() + theme_minimal() 
  print(gg)
}
df_p = fit$draws(format = "df", variables=c('p'),  inc_warmup = F) |> mutate(chain=as.character(.chain))

inad_lga$post_mean <- df_p[,1:nrow(inad_lga)] |> colMeans() |> unname()
inad_lga$post_var <-  sapply(df_p[,1:nrow(inad_lga)],var) |> unname()
inad_lga$post_cv <- inad_lga$post_var / sqrt(inad_lga$post_mean)
inad_lga$design_based_cv <- inad_lga$design_based_var / sqrt(inad_lga$design_based_mean)
inad_lga$design_based_cv[!inad_lga$lga_id%in%inad_lga$lga_id] = NA

inad_lga$naive_mean100 <- inad_lga$naive_mean * 100
inad_lga$design_based_mean100 <-inad_lga$design_based_mean * 100
inad_lga$post_mean100 <- inad_lga$post_mean * 100

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

naive_map <- plot_map(data = inad_lga, 
                      col = "naive_mean100", 
                      title = paste0(str_to_title(nutrient)," (naive estimates)"), 
                      metric = "Prevalence of inadequate intake (%)", 
                      level = "lga")

naive_map

# naive_map <- plot_map(data = inad_cluster|>group_by(state,lga)|>summarise(prev100 = sum(cluster_count)/sum(n)*100), 
#                       col = "prev100", 
#                       title = paste0(str_to_title(nutrient)," (naive estimates)"), 
#                       metric = "Prevalence of inadequate intake (%)", 
#                       level = "lga")
# 
# naive_map
# 
# tmap_save(naive_map, "outputs/maps/vb12_inad_naive.png", width = 8, height = 8,
#           units = "in", dpi = 600)

direct_map <- plot_map(data = inad_lga , 
                       col = "design_based_mean100", 
                       title = paste0(str_to_title(nutrient)," (design-based estimates)"), 
                       metric = "Prevalence of inadequate intake (%)", 
                       level = "lga")

direct_map
# tmap_save(direct_map, "outputs/maps/vb12_inad_direct.png", width = 8, height = 8,
#           units = "in", dpi = 600)


FB_map <- plot_map(data = inad_lga , 
                   col = "post_mean100", 
                   title = paste0(str_to_title(nutrient)," (Full-Bayes Variance smoothing area level)"), 
                   metric = "Prevalence of inadequate intake (%)", 
                   level = "lga")

FB_map
# tmap_save(FB_map, here("outputs","maps",paste0("vb12_inad_FB_",model_name,".png")), width = 8, height = 8,
#        
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

library(ggplot2)
library(viridis)
library(sf)

plot_map_ggplot <- function(data, col, title, metric, 
                            limits = c(0, 0.085), 
                            breaks = seq(0, 0.085, by = 0.01)) {
  
  p <- ggplot() +
    geom_sf(data = data, aes_string(fill = col), color = NA) +
    geom_sf(data = nigeria_1, fill = NA, color = "black", lwd = 0.3) +
    scale_fill_viridis_c(
      name = metric,
      option = "viridis",
      limits = limits,
      breaks = breaks,
      na.value = "grey80",
      oob = scales::oob_squish,  # Keeps values within limits
      labels = scales::label_number(accuracy = 0.01),
      guide = guide_colorbar(
        title.position = "top",
        barwidth = unit(8, "cm"),
        barheight = unit(0.5, "cm"),
        label.position = "bottom"
      )
    ) +
    labs(title = title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom",
      legend.title.align = 0.5,
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  
  return(p)
}

# CV of direct estimates: 
vb12_cv_db <- plot_map_ggplot(data = inad_lga, 
                              col = "design_based_cv", 
                              title = "CV of direct HT estimates", 
                              metric = "Variance")

vb12_cv_db

# Variance of smoothed estimates: 
vb12_var_smoothed_FB <- plot_map_ggplot(data = inad_lga, 
                                        col = "post_cv", 
                                        title = "CV of smoothed estimates", 
                                        metric = "Variance")

vb12_var_smoothed_FB