################################################################################
################# SCRIPT FOR REPLICATING ANALYSES FOR VITAMIN A ################
################################################################################

# Author: Sahoko Ishida
# Adapted by: Mo Osman
# Date created: 05-02-2025
# Last edited: 

# In this script, I will generate both direct smoothed ADM2 level estimates
# for inadequate intake of vitamin A in Nigeria. 

# This project uses pre-processed MIMI micro-nutrient intake data derived from:
# Nigeria Living Standards Survey (2018-19) & FAO/INFOODS West African food 
# composition table (2019)

# If INLA pakackage has never been installed before, install and load: 


# INSTALL AND LOAD PACKAGES:

# INLA package is not available on CRAN, so it needs to be installed from a different repository:
if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dependencies = TRUE)
}

# CmdStanR package is not available on CRAN, so it needs to be installed from a different repository:
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
}

rq_packages <- c("readr", "readxl","tidyverse", "srvyr", "sf", "spdep", "tmap", "INLA", 
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
nga_analysis_df <- nga_base_ai |> 
  left_join(household_locations, by = "hhid") |> 
  left_join(hh_info |> dplyr::select(hhid, ea, survey_wgt), by = "hhid")

rm(nga_base_ai, household_locations, hh_info)

#-------------------------------------------------------------------------------

# READ SHAPEFILES: 
nigeria_1 <- st_read("shapefiles/nigeria_1")
nigeria_2 <- st_read("shapefiles/nigeria_2")

#-------------------------------------------------------------------------------

# CALCULATE NAIVE ESTIMATES AT THE ADM2 LEVEL: 
vita_inad_lga <- nga_analysis_df |> 
  group_by(lga) |> 
  summarise(vita_inad_naive = mean(vita_rae_inadequate, na.rm = TRUE),
            vita_inad_naive_var = var(vita_rae_inadequate, na.rm = TRUE))

#-------------------------------------------------------------------------------

# ADD A SYNTHETIC HOUSEHOLD TO LGAs WHERE ALL HOUSEHOLDS ARE IN ONE CATEGORY
# skip this step if no artificial data should be added
# 1 split data into two 
lga_01 <- (vita_inad_lga |> filter(vita_inad_naive%in%c(0,1)) |> select(lga))$lga
nga_analysis_01 <- nga_analysis_df |> filter(lga%in%lga_01)
nga_analysis_rest <- nga_analysis_df |> filter(!lga%in%lga_01)

set.seed(1234)

i = 0

for (eachlga in lga_01){
  i = i+1
  
  df_tmp <- nga_analysis_01 |> 
    filter(lga==eachlga) 
  
  df_synthetic <- df_tmp[sample(nrow(df_tmp),1),] |> 
    mutate(vita_rae_inadequate = abs(vita_rae_inadequate -1 ), hhid = i)
  
  nga_analysis_01 <- bind_rows(nga_analysis_01,df_synthetic)
}

nga_analysis_df <- bind_rows(nga_analysis_rest,nga_analysis_01)|> 
  arrange(ea, hhid)

#-------------------------------------------------------------------------------

# CALCULATE DESIGN-BASED (DIRECT) ESTIMATES AT ADM2 LEVEL

# Create tbl_svy object for analysis: 
nga_analysis_df_svy <- nga_analysis_df |> as_survey_design(ids = c("ea", "hhid"),
                                                            strata = "state",
                                                            weights = "survey_wgt",
                                                            nest = TRUE)

# Calculate vitamin A inadequacy stratified by LGA - include variance of estimates: 
vita_inad_lga <- vita_inad_lga |> 
  left_join(nga_analysis_df_svy |> 
              group_by(lga) |> 
              summarise(vita_inad_dir = survey_mean(vita_rae_inadequate, 
                                                    na.rm = TRUE, 
                                                    vartype = "var")),
            by = "lga")

# calculate degree of freedom
vita_inad_lga <- vita_inad_lga |> 
  left_join((nga_analysis_df |> 
               group_by(lga, ea) |> 
               summarise(n=n()) |> 
               group_by(lga) |> 
               summarise(d = n()-1, k = sum(n)) |> 
               ungroup()),
            by = "lga")

# Calculate normalised survey weight (ranging from 0 to 1): 
nga_analysis_df <- nga_analysis_df |> 
  left_join((nga_analysis_df|>
               group_by(lga) |>
               reframe(hhid = hhid,
                       survey_wgt=survey_wgt,
                       sum_wgt = sum(survey_wgt),
                       norm_survey_wgt =survey_wgt/sum_wgt) |>
               select(-lga, -survey_wgt,-sum_wgt)),
            by = 'hhid') 

# Calculate weighted proportion, (unbiased) variance of the estimate and degree of freedom:

vita_inad_lga <- vita_inad_lga |> 
  left_join( nga_analysis_df |> 
               group_by(lga) |>
               summarise(n_obs = n(),
                         n_eff = 1/sum(norm_survey_wgt^2), # n_eff denotes effective sample size
                         degf = (n_eff-1),
                         vita_inad_wdir = weighted.mean(vita_rae_inadequate, norm_survey_wgt),
                         vita_inad_wdir_var = vita_inad_wdir*(1-vita_inad_wdir)/degf),
             by = 'lga')

#-------------------------------------------------------------------------------

# PREPARE SPATIAL DATA

# Give adm2_index
nigeria_2$lga_id <- 1:nrow(nigeria_2)

# Join spatial data: 
vita_inad_lga <- nigeria_2 %>% 
  left_join(vita_inad_lga, by = c("lga" = "lga"))

# Define as spatial object:
geo <- as_Spatial(nigeria_2)

# Define spatial adjacency matrix: 
nb <- poly2nb(geo, queen = TRUE)
Amat <- nb2mat(nb, style = "B")
rownames(Amat) <- nigeria_2$lga
colnames(Amat) <- nigeria_2$lga
prep <- prepare_bym2(Amat)

# Perform a quality check: 
qc = rep(0,length(prep$n1))

for (i in 1:length(prep$n1)){
  qc[i] = Amat[prep$n1[i],prep$n2[i]]
}

mean(qc) # should be 1

# Filter to remove LGA's with only 1 EA, and those with zero variance:
vita_inad_lga_complete <- vita_inad_lga |> 
  filter(!is.na(vita_inad_wdir_var), degf!=0,vita_inad_wdir_var > 1e-10) 

#-------------------------------------------------------------------------------

# READ IN STANFILE AND FIT MODEL

stanfile <- here('src','areal_level_BYM2.stan')
mod <- cmdstan_model(stanfile)

# Specify data:
data_list <- list(
  N = nrow(vita_inad_lga),
  NS = nrow(vita_inad_lga_complete),
  adm2_index = vita_inad_lga_complete$lga_id,
  p_hat = vita_inad_lga_complete$vita_inad_wdir,
  v_hat = vita_inad_lga_complete$vita_inad_wdir_var,
  d = vita_inad_lga_complete$degf, 
  k = vita_inad_lga_complete$n_obs,
  N_edges = length(prep$n1),
  node1 = prep$n1,
  node2 = prep$n2,
  scaling_factor = prep$scaling_factor
)

# MCMC
fit <- mod$sample(
  data = data_list,
  seed = 123,
  iter_warmup = 500,
  iter_sampling = 500,
  save_warmup = TRUE,
  chains = 2,
  parallel_chains = 2,
  refresh = 50 # print update every 50 iters
)

fit$save_output_files(dir = here('outputs','mcmc'), basename = "stan_postsample_vitA_BYM2")
csv_files <- fit$output_files()
saveRDS(csv_files,here('outputs','mcmc','csv_files_postsample_vitA_BYM2.rds'))

################################################################################
# QUESTION - what are we looking for in these trace plots?

post_samples = fit$draws(format = "df",  inc_warmup = F)
post_samples$chain = as.character(post_samples$.chain)
colnames(post_samples) <- colnames(post_samples) %>%
  gsub("\\[", "_", .) %>%
  gsub("\\]", "", .) %>%
  gsub(",","_", .)

for (param in (c("sigma_u", "rho", "sigma_tau"))){
  dat = post_samples#[post_samples$chain==2,]
  gg = ggplot(data=dat, aes_string(x=".iteration", y = param, color="chain"))+
    geom_line() + theme_minimal() 
  print(gg)
}

df_p = fit$draws(format = "df", variables=c('p'),  inc_warmup = F) |> mutate(chain=as.character(.chain)) 
df_v = fit$draws(format = "df", variables=c('v'),  inc_warmup = F) |>mutate(chain=as.character(.chain))

colnames(df_p) <- colnames(df_p) %>%
  gsub("\\[", "", .) %>%
  gsub("\\]", "", .) %>%
  gsub(",","_", .)
colnames(df_v) <- colnames(df_v) %>%
  gsub("\\[", "", .) %>%
  gsub("\\]", "", .) %>%
  gsub(",","_", .)

for (param in colnames(df_p)[1:10]){
  dat = df_p
  gg = ggplot(data=dat, aes_string(x=".iteration", y = param, color="chain"))+
    geom_line() + theme_minimal() 
  print(gg)
}

vita_inad_lga$vita_inad_post_mean <- df_p[,1:nrow(vita_inad_lga)] |> colMeans() |> unname()
vita_inad_lga$vita_inad_post_var <-  sapply(df_p[,1:nrow(vita_inad_lga)],var) |> unname()

for (param in colnames(df_v)[1:10]){
  dat = df_v
  gg = ggplot(data=dat, aes_string(x=".iteration", y = param, color="chain"))+
    geom_line() + theme_minimal() 
  print(gg)
}

################################################################################

#-------------------------------------------------------------------------------

# CALCULATE PREVALENCE FROM DIFFERENT METHODS AND MAP:

vita_inad_lga$vita_inad_naive <- vita_inad_lga$vita_inad_naive * 100
vita_inad_lga$vita_inad_dir <- vita_inad_lga$vita_inad_dir * 100
vita_inad_lga$vita_inad_post_mean <- vita_inad_lga$vita_inad_post_mean * 100

# Define mapping funciton
plot_map <- function(data, col, title, metric, level) {
  
  map <- tm_shape(data, component.autoscale = FALSE) + 
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
    tm_shape(nigeria_1, component.autoscale = FALSE) +
    tm_borders(col = "black", lwd = 0.8)
  
  return(map)
}

# NOTE: tmap settings have changed in this new version which has altered the 
# appearance of the maps. Consider reverting to the previous version - REVISIT.

naive_map <- plot_map(data = vita_inad_lga, 
                      col = "vita_inad_naive", 
                      title = "Vitamin A (naive estimates)", 
                      metric = "Prevalence of inadequate intake (%)", 
                      level = "lga")

naive_map

# tmap_save(naive_map, "outputs/maps/vitamin_a/vita_inad_naive.png",
#           height = 8, width = 8, units = "in", dpi = 600)

direct_map <- plot_map(data = vita_inad_lga , 
                       col = "vita_inad_dir", 
                       title = "Vitamin A (design-based estimates)", 
                       metric = "Prevalence of inadequate intake (%)", 
                       level = "lga")

direct_map

# tmap_save(direct_map, "outputs/maps/vitamin_a/vita_inad_dir.png",
#           height = 8, width = 8, units = "in", dpi = 600)

FB_map <- plot_map(data = vita_inad_lga ,
                   col = "vita_inad_post_mean", 
                   title = "Vitamin A (Full-Bayes BYM2 smoothed estimates)", 
                   metric = "Prevalence of inadequate intake (%)", 
                   level = "lga")

FB_map

# tmap_save(FB_map, "outputs/maps/vitamin_a/vita_inad_BYM2.png",
#           height = 8, width = 8, units = "in", dpi = 600)

################################################################################
# QUESTION -- LGA's with missing data - have smoothed estimates with a very high
# prevalence that is not consistent with surrounding LGA's. these seems incorrect.
################################################################################

#-------------------------------------------------------------------------------

# VARIANCE MAPS: 

# redefine mapping function: 
plot_map <- function(data, col, title, metric, level) {
  
  # Create a map: 
  map <- tm_shape(data) + 
    tm_fill(col = col,
            title = metric, 
            style = "cont",
            breaks = seq(0, 0.2, by = 0.05),
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

# Direct: 
vita_var_direct <- plot_map(data = vita_inad_lga, 
                            col = "vita_inad_dir_var", 
                            title = "Variance of direct HT estimates (Vitamin A)", 
                            metric = "Variance", 
                            level = "lga")

vita_var_direct

# tmap_save(vita_var_direct, "outputs/maps/vitamin_a/vita_var_direct.png",
#           height = 8, width = 8, units = "in", dpi = 600)

# Smooth (BYM2): 
vita_var_smooth <- plot_map(data = vita_inad_lga, 
                            col = "vita_inad_post_var", 
                            title = "Variance of smoothed HT estimates (Vitamin A)", 
                            metric = "Variance", 
                            level = "lga")

vita_var_smooth

# tmap_save(vita_var_smooth, "outputs/maps/vitamin_a/vita_var_smooth.png",
#           height = 8, width = 8, units = "in", dpi = 600)

# MISSING LGAs have high variance smoothed estimates ###########################

#-------------------------------------------------------------------------------

# VALIDATION: 

# Clear environment and start validation steps
rm(list = ls())

source("src/02indicator_preprocessing.R")
source("src/helper_functions.R")

# --- Read data ----------------------
# Stan
stan_csv_filenames <- readRDS(here('outputs','mcmc', "csv_files_postsample_vitA_BYM2.rds"))

stanfit <- as_cmdstan_fit(stan_csv_filenames)

# Read in household locations that are matched to shapefiles: 
household_locations <- read_csv("shapefiles/household_locations.csv")

# Get further household info:
hh_info <- read_csv("processed_data/nga_lss1819_hh_info.csv")

# Merge data: 
nga_analysis_df <- nga_base_ai %>% 
  left_join(household_locations, by = "hhid") %>% 
  left_join(hh_info %>% dplyr::select(hhid, ea, survey_wgt), by = "hhid")

rm(nga_base_ai, household_locations, hh_info)

# READ SHAPEFILES: 
nigeria_1 <- st_read("shapefiles/nigeria_1")
nigeria_2 <- st_read("shapefiles/nigeria_2")

# Read adm2 population estimates
pop_estimates <- read.csv(here("population_estimates","adm2_pop_clean.csv"))
state_lga_matched <- read.csv(here("processed_data","state_lga_matched.csv"))
#### NEED THE STATE MATCHED LGA csv file that was used here ######

pop_estimates <- pop_estimates |> left_join(state_lga_matched,by='lga')