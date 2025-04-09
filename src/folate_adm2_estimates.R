################################################################################
########### SCRIPT FOR GENERATING DIRECT AND SMOOTHED ADM2 ESTIMATES ###########
################################################################################

# Author: Sahoko Ishida
# Date created: 11-03-2024
# Last edited: 20-03-2025

# In this script, I will generate both smoothed ADM2 level estimates
# for inadequate intake of vitamin B12 in Nigeria. 

# This project uses pre-processed MIMI micro-nutrient intake data derived from:
# Nigeria Living Standards Survey (2018-19) & FAO/INFOODS West African food 
# composition table (2019)

# INSTALL AND LOAD PACKAGES:

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
nga_analysis_df <- nga_base_ai %>% 
  left_join(household_locations, by = "hhid") %>% 
  left_join(hh_info %>% dplyr::select(hhid, ea, survey_wgt), by = "hhid")

nga_analysis_df |> 
  select(state, lga, ea) |> unique() |> 
  group_by(state,lga) |> 
  summarise(n_ea = n())

rm(nga_base_ai, household_locations, hh_info)

#-------------------------------------------------------------------------------

# READ SHAPEFILES: 
nigeria_1 <- st_read("shapefiles/nigeria_1")
nigeria_2 <- st_read("shapefiles/nigeria_2")

#-------------------------------------------------------------------------------

# CALCULATE NAIVE ESTIMATES AT THE ADM2 LEVEL: 
folate_inad_lga <- nga_analysis_df %>% 
  group_by(lga) %>% 
  summarise(folate_inad_naive = mean(folate_inadequate, na.rm = TRUE),
            folate_inad_naive_var = var(folate_inadequate, na.rm = TRUE))

#-------------------------------------------------------------------------------

# ADD A SYNTHETIC HOUSEHOLD TO LGAs WHERE ALL HOUSEHOLDS ARE IN ONE CATEGORY
# skip this step if no artificial data should be added
# 1 split data into two 
run = FALSE
if (run){
  lga_01 <- (folate_inad_lga |> filter(folate_inad_naive%in%c(0,1)) |> select(lga))$lga
  nga_analysis_01 <- nga_analysis_df |> filter(lga%in%lga_01)
  nga_analysis_rest <- nga_analysis_df |> filter(!lga%in%lga_01)
  set.seed(1234)
  i = 0
  for (eachlga in lga_01){
    i = i+1
    df_tmp <- nga_analysis_01 |> filter(lga==eachlga) 
    df_synthetic <- df_tmp[sample(nrow(df_tmp),1),] |> 
      mutate(folate_inadequate = abs(folate_inadequate -1 ),
             hhid = i)
    nga_analysis_01 <- bind_rows(nga_analysis_01,df_synthetic)
  }
  nga_analysis_df <- bind_rows(nga_analysis_rest,nga_analysis_01)|> arrange(ea, hhid)
}

#-------------------------------------------------------------------------------

# CALCULATE DESIGN-BASED (DIRECT) ESTIMATES AT ADM2 LEVEL: 
# Firstly need to create a tbl_svy object to be used for analysis.
nga_analysis_df_svy <- nga_analysis_df %>% as_survey_design(ids = c("ea", "hhid"),
                                                            strata = "state",
                                                            weights = "survey_wgt",
                                                            nest = TRUE)

# Calculate vitamin B12 inadequacy stratified by LGA - include variance of estimates: 
folate_inad_lga <- folate_inad_lga %>% 
  left_join(nga_analysis_df_svy %>% 
              group_by(lga) %>% 
              summarise(folate_inad_dir = survey_mean(folate_inadequate, 
                                                    na.rm = TRUE, 
                                                    vartype = "var")),
            by = "lga")

# calculate degree of freedom
folate_inad_lga <- folate_inad_lga |> 
  left_join(
    (nga_analysis_df |> group_by(lga, ea) |> 
       summarise(n=n()) |> group_by(lga) |> 
       summarise(d = n()-1, k = sum(n)) |> ungroup()),
    by = "lga"
  )

# Calculate variance & degree of freedom of weighted direct estimates (ignoring survey design)
## standardise weight
nga_analysis_df <- nga_analysis_df |> 
  left_join((nga_analysis_df|>
               group_by(lga) |>
               reframe(hhid = hhid,
                       survey_wgt=survey_wgt,
                       sum_wgt = sum(survey_wgt),
                       norm_survey_wgt =survey_wgt/sum_wgt) |>
               select(-lga, -survey_wgt,-sum_wgt)),
            by = 'hhid'
  ) 


## calculate weighted proportion, (unbiased) variance of the estimate and degree of freedom
folate_inad_lga <- folate_inad_lga |> 
  left_join( nga_analysis_df |> group_by(lga) |>
               summarise(n_obs = n(),
                         n_eff = 1/sum(norm_survey_wgt^2),
                         degf = (n_eff-1),
                         folate_inad_wdir = weighted.mean(folate_inadequate, norm_survey_wgt),
                         folate_inad_wdir_var = folate_inad_wdir*(1-folate_inad_wdir)/degf),
             by = 'lga'
  )

folate_inad_lga <- folate_inad_lga|> 
  left_join(nga_analysis_df |> 
              select(state, lga, ea) |> unique() |> 
              group_by(state,lga) |>  summarise(n_ea = n()),
            by = 'lga'
  )

# Give adm2_index
nigeria_2$lga_id <- 1:nrow(nigeria_2)

# Join spatial data: 
folate_inad_lga <- nigeria_2 %>% 
  left_join(folate_inad_lga, by = c("lga" = "lga"))

# Prepare data for 
# Define spatial object:
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
folate_inad_lga$folate_inad_naive <- folate_inad_lga$folate_inad_naive * 100
folate_inad_lga$folate_inad_dir <- folate_inad_lga$folate_inad_dir * 100


### removing all LGAs with only one EA 
#folate_inad_lga_complete<- folate_inad_lga |> filter(!is.na(folate_inad_dir_var),d!=0,folate_inad_dir_var>1e-10) # ALL LGA
folate_inad_lga_complete <- folate_inad_lga |> filter(!is.na(folate_inad_wdir)) 

stanfile <- here('src','areal_level_binom_BYM2.stan')
mod <- cmdstan_model(stanfile)
data_list <- list(
  N = nrow(folate_inad_lga),
  NS = nrow(folate_inad_lga_complete),
  adm2_index = folate_inad_lga_complete$lga_id,
  p_hat = folate_inad_lga_complete$folate_inad_wdir,
  n_star = round(folate_inad_lga_complete$n_eff),
  y_star = round(folate_inad_lga_complete$n_eff * folate_inad_lga_complete$folate_inad_wdir),
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
  iter_sampling = 2000,
  save_warmup = TRUE,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 50 iters
)

fit$save_output_files(dir = here('outputs','mcmc'), basename = "stan_postsample_folate_BYM2_arealbinom")
csv_files <- fit$output_files()
saveRDS(csv_files,here('outputs','mcmc','csv_files_postsample_folate_BYM2_arealbinom.rds'))

post_summary = fit$summary(variables=c("b0","rho","sigma_u"))
post_summary
post_samples = fit$draws(format = "df",  inc_warmup = F)
post_samples$chain = as.character(post_samples$.chain)
colnames(post_samples) <- colnames(post_samples) %>%
  gsub("\\[", "_", .) %>%
  gsub("\\]", "", .) %>%
  gsub(",","_", .)

for (param in c("b0","rho","sigma_u")){
  dat = post_samples#[post_samples$chain==2,]
  gg = ggplot(data=dat, aes_string(x=".iteration", y = param, color="chain"))+
    geom_line() + theme_minimal() 
  print(gg)
}

df_p = fit$draws(format = "df", variables=c('p'),  inc_warmup = F) |> mutate(chain=as.character(.chain)) 

colnames(df_p) <- colnames(df_p) %>%
  gsub("\\[", "", .) %>%
  gsub("\\]", "", .) %>%
  gsub(",","_", .)

for (param in colnames(df_p)[1:10]){
  dat = df_p
  gg = ggplot(data=dat, aes_string(x=".iteration", y = param, color="chain"))+
    geom_line() + theme_minimal() 
  print(gg)
}

folate_inad_lga$folate_inad_post_mean <- df_p[,1:nrow(folate_inad_lga)] |> colMeans() |> unname()
folate_inad_lga$folate_inad_post_var <-  sapply(df_p[,1:nrow(folate_inad_lga)],var) |> unname()

folate_inad_lga$folate_inad_post_mean <- folate_inad_lga$folate_inad_post_mean * 100

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

naive_map <- plot_map(data = folate_inad_lga, 
                      col = "folate_inad_naive", 
                      title = "Vitamin B12 (naive estimates)", 
                      metric = "Prevalence of inadequate intake (%)", 
                      level = "lga")

naive_map

# tmap_save(naive_map, "outputs/maps/folate/folate_inad_naive.png", width = 8, height = 8,
#           units = "in", dpi = 600)

direct_map <- plot_map(data = folate_inad_lga , 
                       col = "folate_inad_dir", 
                       title = "Vitamin B12 (design-based estimates)", 
                       metric = "Prevalence of inadequate intake (%)", 
                       level = "lga")

direct_map
# tmap_save(direct_map, "outputs/maps/folate_inad_direct.png", width = 8, height = 8,
#           units = "in", dpi = 600)

# folate_inad_lga_tmp <- folate_inad_lga
# folate_inad_lga_tmp$folate_inad_wdir[folate_inad_lga$degf==0|folate_inad_lga$folate_inad_wdir_var<1e-10] = NA
# folate_inad_lga_tmp$folate_inad_wdir = folate_inad_lga_tmp$folate_inad_wdir*100
# direct_map2 <- plot_map(data = folate_inad_lga_tmp , 
#                        col = "folate_inad_wdir", 
#                        title = "Vitamin B12 (design-based estimates)", 
#                        metric = "Prevalence of inadequate intake (%)", 
#                        level = "lga")
# 
# direct_map2
# tmap_save(direct_map2, "outputs/maps/folate_inad_direct_01rm.png", width = 8, height = 8,
#           units = "in", dpi = 600)


folate_inad_lga$folate_inad_post_mean_obs_only = folate_inad_lga$folate_inad_post_mean
folate_inad_lga$folate_inad_post_mean_obs_only[!folate_inad_lga$lga_id%in%folate_inad_lga_complete$lga_id] = NA
FB_map <- plot_map(data = folate_inad_lga,
                   col = "folate_inad_post_mean_obs_only", 
                   title = "Vitamin B12 (Full-Bayes BYM2 smoothed estimates)", 
                   metric = "Prevalence of inadequate intake (%)", 
                   level = "lga")

FB_map
tmap_save(FB_map, "outputs/maps/folate/folate_inad_FB_BYM2_arealbinom_obsonly.png", width = 8, height = 8,
          units = "in", dpi = 600)

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
vb12_var_direct <- plot_map(data = folate_inad_lga, 
                            col = "folate_inad_dir_var", 
                            title = "Variance of direct HT estimates", 
                            metric = "Variance", 
                            level = "lga")

vb12_var_direct

# tmap_save(vb12_var_direct, "outputs/maps/var_direct.png", width = 8, height = 8,
#           units = "in", dpi = 600)

# Variance of smoothed estimates: 
folate_inad_lga$folate_inad_post_var_obs_only = folate_inad_lga$folate_inad_post_var
folate_inad_lga$folate_inad_post_var_obs_only[!folate_inad_lga$lga_id%in%folate_inad_lga_complete$lga_id] = NA


vb12_var_smoothed_FB <- plot_map(data = folate_inad_lga, 
                                 col = "folate_inad_post_var", 
                                 title = "Variance of smoothed estimates", 
                                 metric = "Variance", 
                                 level = "lga")

vb12_var_smoothed_FB
# tmap_save(vb12_var_smoothed_FB, "outputs/maps/var_FB.png", width = 8, height = 8,
#           units = "in", dpi = 600)


#£ ---  End of Script---