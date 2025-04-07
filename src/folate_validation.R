################################################################################
########### SCRIPT FOR VALIDATING MODELLED ESTIMATES ###########################
################################################################################

# Author: Sahoko Ishida
# Date created: 07-11-2024
# Last edited: 31-01-2024

# This script is to check 

# This project uses pre-processed MIMI micro-nutrient intake data derived from:
# Nigeria Living Standards Survey (2018-19) & FAO/INFOODS West African food 
# composition table (2019)

# INSTALL AND LOAD PACKAGES:
rq_packages <- c("readr", "tidyverse", "srvyr", "sf", "spdep", "tmap", "INLA", 
                 "SUMMER", "wesanderson", "cmdstanr", "ggplot2")

installed_packages <- rq_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(rq_packages[!installed_packages])
}

lapply(rq_packages, require, character.only = T)

rm(list= c("rq_packages", "installed_packages"))

# ------------------- -----------
source("src/02indicator_preprocessing.R")
source("src/helper_functions.R")

# --- Read data ----------------------
# Stan
# stan_csv_filenames <- readRDS(here('outputs','mcmc', "csv_files_postsample_folate_BYM2.rds"))
stan_csv_filenames<- readRDS(here('outputs','mcmc','csv_files_postsample_folate_BYM2_arealbinom.rds'))
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
#pop_estimates <- read.csv(here("population_estimates","adm2_pop_clean.csv"))
state_lga_matched <- read.csv(here("processed_data","state_lga_matched.csv"))
#-------------------------------------------------------------------------------

# CALCULATE NAIVE ESTIMATES AT THE ADM2 LEVEL: 
fl_inad_lga <- nga_analysis_df |> 
  group_by(lga) |> 
  summarise(fl_inad_naive = mean(folate_inadequate, na.rm = TRUE),
            fl_inad_naive_var = var(folate_inadequate, na.rm = TRUE))

#-------------------------------------------------------------------------------

# ADD A SYNTHETIC HOUSEHOLD TO LGAs WHERE ALL HOUSEHOLDS ARE IN ONE CATEGORY
# skip this step if no artificial data should be added
# 1 split data into two 
run = FALSE
lga_01 <- (fl_inad_lga |> filter(fl_inad_naive%in%c(0,1)) |> select(lga))$lga
nga_analysis_01 <- nga_analysis_df |> filter(lga%in%lga_01)
nga_analysis_rest <- nga_analysis_df |> filter(!lga%in%lga_01)
set.seed(1234)
if (run){
  i = 0
  for (eachlga in lga_01){
    i = i+1
    df_tmp <- nga_analysis_01 |> 
      filter(lga==eachlga) 
    df_synthetic <- df_tmp[sample(nrow(df_tmp),1),] |> 
      mutate(folate_inadequate = abs(folate_inadequate -1), hhid = i)
    nga_analysis_01 <- bind_rows(nga_analysis_01,df_synthetic)
  }
  nga_analysis_df <- bind_rows(nga_analysis_rest,nga_analysis_01)|> 
    arrange(ea, hhid)
}


#-------------------------------------------------------------------------------

# CALCULATE DESIGN-BASED (DIRECT) ESTIMATES AT ADM2 LEVEL

# Create tbl_svy object for analysis: 
nga_analysis_df_svy <- nga_analysis_df |> as_survey_design(ids = c("ea", "hhid"),
                                                           strata = "state",
                                                           weights = "survey_wgt",
                                                           nest = TRUE)

# Calculate vitamin A inadequacy stratified by LGA - include variance of estimates: 
fl_inad_lga <- fl_inad_lga |> 
  left_join(nga_analysis_df_svy |> 
              group_by(lga) |> 
              summarise(fl_inad_dir = survey_mean(folate_inadequate, 
                                                  na.rm = TRUE, 
                                                  vartype = "var")),
            by = "lga")

# calculate degree of freedom
fl_inad_lga <- fl_inad_lga |> 
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
               select(-lga, -survey_wgt)),
            by = 'hhid') 

# Calculate weighted proportion, (unbiased) variance of the estimate and degree of freedom:

fl_inad_lga <- fl_inad_lga |> 
  left_join( nga_analysis_df |> 
               group_by(state,lga) |>
               summarise(n_obs = n(),
                         n_eff = 1/sum(norm_survey_wgt^2), # n_eff denotes effective sample size
                         degf = (n_eff-1),
                         sum_wgt = sum(survey_wgt),
                         fl_inad_wdir = weighted.mean(folate_inadequate, norm_survey_wgt),
                         fl_inad_wdir_var = fl_inad_wdir*(1-fl_inad_wdir)/degf),
             by = 'lga')

fl_inad_lga <- nigeria_2 %>% 
  left_join(fl_inad_lga, by = c("lga" = "lga"))
# ------------------------------------------------
# Model estimates 
df_p = stanfit$draws(format = "df", variables=c('p'),  inc_warmup = F) |> mutate(chain=as.character(.chain)) 
colnames(df_p) <- colnames(df_p) %>%
  gsub("\\[", "", .) %>%
  gsub("\\]", "", .) %>%
  gsub(",","_", .)
# --- compute mean and credible interval lga level -----
alpha = 0.1
fl_inad_lga$fl_inad_post_mean <- df_p[,1:nrow(fl_inad_lga)] |> colMeans() |> unname()
fl_inad_lga$fl_inad_post_var <-  sapply(df_p[,1:nrow(fl_inad_lga)],var) |> unname()
df_quantile =  apply(df_p[,1:nrow(fl_inad_lga)], 2 , quantile , probs = c(alpha/2,1-alpha/2) , na.rm = TRUE ) |> t()
fl_inad_lga$fl_inad_post_lower90 <- df_quantile[,1]
fl_inad_lga$fl_inad_post_upper90 <- df_quantile[,2]
fl_inad_lga$fl_inad_post_ci_length90<- fl_inad_lga$fl_inad_post_upper90 - fl_inad_lga$fl_inad_post_lower90 

# Coverage 
fl_inad_lga$fl_inad_wdir_lower90_wilson <- wilson_lower(fl_inad_lga$fl_inad_wdir, fl_inad_lga$n_eff,alpha = alpha)
fl_inad_lga$fl_inad_wdir_upper90_wilson <- wilson_upper(fl_inad_lga$fl_inad_wdir, fl_inad_lga$n_eff,alpha = alpha)
fl_inad_lga$fl_inad_wdir_ci_length90_wilson <- fl_inad_lga$fl_inad_wdir_upper90_wilson - fl_inad_lga$fl_inad_wdir_lower90_wilson 
mean(fl_inad_lga$fl_inad_wdir_ci_length90_wilson, na.rm = T)
fl_inad_lga$fl_inad_wdir_lower90_wald <- (fl_inad_lga$fl_inad_wdir - qnorm(1-alpha/2)*sqrt(fl_inad_lga$fl_inad_wdir_var))
fl_inad_lga$fl_inad_wdir_upper90_wald <- (fl_inad_lga$fl_inad_wdir + qnorm(1-alpha/2)*sqrt(fl_inad_lga$fl_inad_wdir_var))
fl_inad_lga$fl_inad_wdir_ci_length90_wald<- fl_inad_lga$fl_inad_wdir_upper90_wald - fl_inad_lga$fl_inad_wdir_lower90_wald 
# int-lgength
mean(fl_inad_lga$fl_inad_wdir_ci_length90_wald, na.rm = T)
mean(fl_inad_lga$fl_inad_post_ci_length90, na.rm = T)
# coverage
fl_inad_lga$int_overlap <-check_overlap(fl_inad_lga$fl_inad_wdir_lower90_wilson,fl_inad_lga$fl_inad_wdir_upper90_wilson,
                                          fl_inad_lga$fl_inad_post_lower90,fl_inad_lga$fl_inad_post_upper90)
mean(fl_inad_lga$int_overlap[!fl_inad_lga$fl_inad_wdir%in%c(0.0,1.0)], na.rm=T)
mean(fl_inad_lga$int_overlap, na.rm=T)
# --- State level Validation -----------------
# compute mean and credible interval at state level
# lga_pop <- pop_estimates |> group_by(state) |> 
#   reframe(lga = lga,pop_lga = population,pop_state = sum(pop_lga))  |>
#   mutate(geo_weight = pop_lga/pop_state) 

nga_analysis_df <- nga_analysis_df |> 
  left_join((nga_analysis_df|>
               group_by(state) |>
               reframe(hhid = hhid,
                       survey_wgt=survey_wgt,
                       sum_wgt_state = sum(survey_wgt),
                       norm_survey_wgt_state = survey_wgt/sum_wgt_state) |>
               select(-state, -survey_wgt)),
            by = 'hhid'
  ) 

lga_in_survey <- nga_analysis_df$lga |> unique()
fl_inad_state <- nga_analysis_df |> group_by(state) |>
  summarise(n_obs = n(),
            n_eff = 1/sum(norm_survey_wgt_state^2),
            degf = (n_eff-1),
            fl_inad_wdir =  weighted.mean(folate_inadequate, norm_survey_wgt_state),
            fl_inad_wdir_var = fl_inad_wdir*(1-fl_inad_wdir)/degf
  ) 

fl_inad_state$fl_inad_wdir_lower90_wilson <- wilson_lower(fl_inad_state$fl_inad_wdir, fl_inad_state$n_eff,alpha = alpha)
fl_inad_state$fl_inad_wdir_upper90_wilson <- wilson_upper(fl_inad_state$fl_inad_wdir, fl_inad_state$n_eff,alpha = alpha)
fl_inad_state$fl_inad_wdir_ci_length90_wilson <- fl_inad_state$fl_inad_wdir_upper90_wilson - fl_inad_state$fl_inad_wdir_lower90_wilson 
mean(fl_inad_state$fl_inad_wdir_ci_length90_wilson, na.rm = T)
fl_inad_state$fl_inad_wdir_lower90_wald <- (fl_inad_state$fl_inad_wdir - qnorm(1-alpha/2)*sqrt(fl_inad_state$fl_inad_wdir_var))
fl_inad_state$fl_inad_wdir_upper90_wald <- (fl_inad_state$fl_inad_wdir + qnorm(1-alpha/2)*sqrt(fl_inad_state$fl_inad_wdir_var))
fl_inad_state$fl_inad_wdir_ci_length90_wald<- fl_inad_state$fl_inad_wdir_upper90_wald - fl_inad_state$fl_inad_wdir_lower90_wald 
mean(fl_inad_state$fl_inad_wdir_ci_length90_wald, na.rm = T)

df_post <- cbind(fl_inad_lga |> select(state,lga,sum_wgt),
                 t(df_p[,1:nrow(fl_inad_lga)])) |> 
  filter(!is.na(sum_wgt))

df_post_state <- df_post |> as.data.frame() |>
  select(-geometry) |>
  group_by(state) |>
  summarise(across(all_of(colnames(df_post)[grepl('X',colnames(df_post))]), ~ weighted.mean(., w = sum_wgt)))


# df_post <- cbind(fl_inad_lga |> select(lga),t(df_p[,1:nrow(fl_inad_lga)])) |>
#   left_join(lga_pop  |> select(state, lga, geo_weight) ,
#             by = 'lga') |> filter(lga!="bakassi")
# 
# df_post_state <- df_post |> as.data.frame() |>
#   select(-geometry) |> #filter(lga%in%lga_in_survey) |>  un-comment to exclude LGAs not included in the original survey
#   group_by(state) |>
#   summarise(across(all_of(colnames(df_post)[grepl('X',colnames(df_post))]), ~ weighted.mean(., w = geo_weight)))

df_post_state$fl_inad_post_mean <- df_post_state[,2:(nrow(df_p)+1)] |> rowMeans()
df_quantile =  apply(df_post_state[,2:nrow(df_p)], 1 , quantile , probs = c(alpha/2,1-alpha/2) , na.rm = TRUE ) |> t()
df_post_state$fl_inad_post_lower90 <- df_quantile[,1]
df_post_state$fl_inad_post_upper90 <- df_quantile[,2]
df_post_state$fl_inad_posr_ci_length90 <- df_post_state$fl_inad_post_upper90 - df_post_state$fl_inad_post_lower90 
mean(df_post_state$fl_inad_posr_ci_length90)
fl_inad_state  <- fl_inad_state |> left_join((df_post_state|> select(-colnames(df_post)[grepl('X',colnames(df_post))])), by = 'state')

mean(abs(fl_inad_state$fl_inad_wdir - fl_inad_state$fl_inad_post_mean))
cor(fl_inad_state$fl_inad_wdir, fl_inad_state$fl_inad_post_mean, method = "pearson")
fl_inad_state$int_overlap <-check_overlap(fl_inad_state$fl_inad_wdir_lower90_wald,fl_inad_state$fl_inad_wdir_upper90_wald,
                                            fl_inad_state$fl_inad_post_lower90,fl_inad_state$fl_inad_post_upper90)
fl_inad_state$error = fl_inad_state$fl_inad_post_mean - fl_inad_state$fl_inad_wdir 

mean(fl_inad_state$int_overlap)
mean(abs(fl_inad_state$error))
#write_csv(fl_inad_state, here("processed_data","state_level_estimates_withSynthetic.csv"))

#-------------------------------------------------------------------------------
# Plots
#-------------------------------------------------------------------------------
fl_inad_state <- fl_inad_state |> left_join(nigeria_1|> select(state, geometry), by='state')
fl_inad_state <- st_as_sf(fl_inad_state)
p = ggplot(fl_inad_state |> st_as_sf()
           , aes(fill = fl_inad_post_mean)) +
  geom_sf() +
  #facet_wrap(~Age) +
  scale_fill_gradientn(values = c("1","0.8","0.6","0.4","0.2","0"), 
                       colors = c("red4","red3","#FC4E07","#E7B800","green4"),
                       limits=c(0,1))+
  labs(title='Full Bayse estimate',fill = "VB12 deficiency")+
  theme_minimal(base_size = 14) 
p
#ggsave(here("outputs","maps","folate",paste0("fl_inad_state_FB.png")), p, width = 10, height = 10,bg = "white")

p = ggplot(fl_inad_state |> st_as_sf()
           , aes(fill = fl_inad_wdir)) +
  geom_sf() +
  #facet_wrap(~Age) +
  scale_fill_gradientn(values = c("1","0.8","0.6","0.4","0.2","0"), 
                       colors = c("red4","red3","#FC4E07","#E7B800","green4"),
                       limits=c(0,1))+
  labs(title='Direct estimate',fill = "VB12 deficiency")+
  theme_minimal(base_size = 14) 
p
#ggsave(here("outputs","maps","folate",paste0("fl_inad_state_direct.png")), p, width = 10, height = 10,bg = "white")

p = ggplot(fl_inad_state |> st_as_sf()
           , aes(fill = int_overlap,label = str_to_title(state))) +
  geom_sf() +
  #geom_sf_text(fun.geometry = sf::st_centroid) + # make Manhattan behave itself
  labs(title='Coverage')+
  theme_minimal(base_size = 14) 
p
#ggsave(here("outputs","maps",paste0("fl_inad_state_coverage.png")), p, width = 10, height = 10,bg = "white")


p = ggplot(fl_inad_state |> st_as_sf()
           , aes(fill = error)) +
  geom_sf() +
  #facet_wrap(~Age) +
  scale_fill_gradient2(  high = "tomato",
                         mid = "white",
                         low = "royalblue") +
  # scale_fill_gradientn(values = c("-0.3","-0.15","0","0.15","0.3"), 
  #                      colors = c("red4","red3","#FC4E07","#E7B800","green4"))+
  labs(title='Error')+
  theme_minimal(base_size = 14) 
p
#ggsave(here("outputs","maps",paste0("fl_inad_state_error.png")), p, width = 10, height = 10,bg = "white")

p <- ggplot(fl_inad_state,aes(x=fl_inad_post_mean,y=fl_inad_wdir,label = str_to_title(state))) +
  geom_point(aes(color = int_overlap))+
  geom_abline(slope = 1,intercept = 0, linetype=2)+
  ggrepel::geom_text_repel()+
  xlim(0,1) + ylim(0,1)+
  labs(x='FB modelled estimates',y='Direct Estimates')+
  theme_minimal(base_size = 14)
p
ggsave(here("outputs","maps","folate",paste0("fl_inad_state_scatter.png")), p, width = 10, height = 10,bg = "white")

df_plot <- data.frame(
  y_est = c(fl_inad_state$fl_inad_wdir,fl_inad_state$fl_inad_post_mean),
  lower = c(fl_inad_state$fl_inad_wdir_lower90_wald, fl_inad_state$fl_inad_post_lower90),
  upper = c(fl_inad_state$fl_inad_wdir_upper90_wald, fl_inad_state$fl_inad_post_upper90),
  type = rep(c("Direct","Modelled"), each = nrow(fl_inad_state)),
  state = rep(fl_inad_state$state, times=2)
)
for (i in 1:37){
  p = ggplot(df_plot |> filter(state==unique(df_plot$state[i])), aes(y = y_est,x= type)) +
    geom_linerange(aes(ymin=upper,ymax=lower))+
    geom_point()+
    ggtitle(str_to_title(df_plot$state[i]))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme_minimal(base_size = 14) 
  print(p)
}

