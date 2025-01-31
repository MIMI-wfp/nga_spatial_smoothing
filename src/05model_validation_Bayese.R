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
stan_csv_filenames <- readRDS(here('outputs','mcmc','csv_files_postsample_vb12_BYM2.rds'))
# stan_csv_filenames <- c(here('outputs','mcmc','stan_postsample_BYM2_synthetic-202501311100-1-959a29.csv'),
#                         here('outputs','mcmc','stan_postsample_BYM2_synthetic-202501311100-2-959a29.csv'))

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
pop_estimates <- pop_estimates |> left_join(state_lga_matched,by='lga')
# --- Further prosessing ----
vb12_inad_lga <- nga_analysis_df %>% 
  group_by(lga) %>% 
  summarise(vb12_inad_naive = mean(vb12_inadequate, na.rm = TRUE),
            vb12_inad_naive_var = var(vb12_inadequate, na.rm = TRUE))
#-------------------------------------------------------------------------------
# 1 split data into two 
lga_01 <- (vb12_inad_lga |> filter(vb12_inad_naive%in%c(0,1)) |> select(lga))$lga

nga_analysis_01 <- nga_analysis_df |> filter(lga%in%lga_01)
nga_analysis_rest <- nga_analysis_df |> filter(!lga%in%lga_01)
i = 0
for (eachlga in lga_01){
  i = i+1
  df_tmp <- nga_analysis_01 |> filter(lga==eachlga) 
  df_synthetic <- df_tmp[sample(nrow(df_tmp),1),] |> 
    mutate(vb12_inadequate = abs(vb12_inadequate -1 ),
           hhid = i)
  nga_analysis_01 <- bind_rows(nga_analysis_01,df_synthetic)
}
nga_analysis_01 <- nga_analysis_01 |> arrange(ea, hhid, vb12_inadequate)
nga_analysis_df <- bind_rows(nga_analysis_rest,nga_analysis_01)|> arrange(ea, hhid, vb12_inadequate)
#-------------------------------------------------------------------------------

# CALCULATE DESIGN-BASED (DIRECT) ESTIMATES AT ADM2 LEVEL: 
# Firstly need to create a tbl_svy object to be used for analysis.
# Survey data
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
vb12_inad_lga <- vb12_inad_lga |> 
  left_join( nga_analysis_df |> group_by(lga) |>
               summarise(n_obs = n(),
                         n_eff = 1/sum(norm_survey_wgt^2),
                         degf = (n_eff-1),
                         vb12_inad_wdir = weighted.mean(vb12_inadequate, norm_survey_wgt),
                         vb12_inad_wdir_var = vb12_inad_wdir*(1-vb12_inad_wdir)/degf),
             by = 'lga'
  )
vb12_inad_lga <- nigeria_2 %>% 
  left_join(vb12_inad_lga, by = c("lga" = "lga"))
# ------------------------------------------------
# Model estimates 
df_p = stanfit$draws(format = "df", variables=c('p'),  inc_warmup = F) |> mutate(chain=as.character(.chain)) 
colnames(df_p) <- colnames(df_p) %>%
  gsub("\\[", "", .) %>%
  gsub("\\]", "", .) %>%
  gsub(",","_", .)
# --- compute mean and credible interval lga level -----
alpha = 0.1
vb12_inad_lga$vb12_inad_post_mean <- df_p[,1:nrow(vb12_inad_lga)] |> colMeans() |> unname()
vb12_inad_lga$vb12_inad_post_var <-  sapply(df_p[,1:nrow(vb12_inad_lga)],var) |> unname()
df_quantile =  apply(df_p[,1:nrow(vb12_inad_lga)], 2 , quantile , probs = c(alpha/2,1-alpha/2) , na.rm = TRUE ) |> t()
vb12_inad_lga$vb12_inad_post_lower90 <- df_quantile[,1]
vb12_inad_lga$vb12_inad_post_upper90 <- df_quantile[,2]
vb12_inad_lga$vb12_inad_post_ci_length90<- vb12_inad_lga$vb12_inad_post_upper90 - vb12_inad_lga$vb12_inad_post_lower90 

# Coverage 
vb12_inad_lga$vb12_inad_wdir_lower90_wilson <- wilson_lower(vb12_inad_lga$vb12_inad_wdir, vb12_inad_lga$n_eff,alpha = alpha)
vb12_inad_lga$vb12_inad_wdir_upper90_wilson <- wilson_upper(vb12_inad_lga$vb12_inad_wdir, vb12_inad_lga$n_eff,alpha = alpha)
vb12_inad_lga$vb12_inad_wdir_ci_length90_wilson <- vb12_inad_lga$vb12_inad_wdir_upper90_wilson - vb12_inad_lga$vb12_inad_wdir_lower90_wilson 
mean(vb12_inad_lga$vb12_inad_wdir_ci_length90_wilson, na.rm = T)
vb12_inad_lga$vb12_inad_wdir_lower90_wald <- (vb12_inad_lga$vb12_inad_wdir - qnorm(1-alpha/2)*sqrt(vb12_inad_lga$vb12_inad_wdir_var))
vb12_inad_lga$vb12_inad_wdir_upper90_wald <- (vb12_inad_lga$vb12_inad_wdir + qnorm(1-alpha/2)*sqrt(vb12_inad_lga$vb12_inad_wdir_var))
vb12_inad_lga$vb12_inad_wdir_ci_length90_wald<- vb12_inad_lga$vb12_inad_wdir_upper90_wald - vb12_inad_lga$vb12_inad_wdir_lower90_wald 
# int-lgength
mean(vb12_inad_lga$vb12_inad_wdir_ci_length90_wald, na.rm = T)
mean(vb12_inad_lga$vb12_inad_post_ci_length90, na.rm = T)
# coverage
vb12_inad_lga$int_overlap <-check_overlap(vb12_inad_lga$vb12_inad_wdir_lower90_wilson,vb12_inad_lga$vb12_inad_wdir_upper90_wilson,
                                          vb12_inad_lga$vb12_inad_post_lower90,vb12_inad_lga$vb12_inad_post_upper90)
mean(vb12_inad_lga$int_overlap[!vb12_inad_lga$vb12_inad_wdir%in%c(0.0,1.0)], na.rm=T)
mean(vb12_inad_lga$int_overlap, na.rm=T)
# --- State level Validation -----------------
# compute mean and credible interval at state level
lga_pop <- pop_estimates |> group_by(state) |> 
  reframe(lga = lga,pop_lga = population,pop_state = sum(pop_lga))  |>
  mutate(geo_weight = pop_lga/pop_state) 
nga_analysis_df <- nga_analysis_df |> 
  left_join((nga_analysis_df|>
               group_by(state) |>
               reframe(hhid = hhid,
                       survey_wgt=survey_wgt,
                       sum_wgt = sum(survey_wgt),
                       norm_survey_wgt_state = survey_wgt/sum_wgt) |>
               select(-state, -survey_wgt,-sum_wgt)),
            by = 'hhid'
  ) 

lga_in_survey <- nga_analysis_df$lga |> unique()
vb12_inad_state <- nga_analysis_df |> filter(lga!="bakassi") |> group_by(state) |>
  summarise(n_obs = n(),
            n_eff = 1/sum(norm_survey_wgt_state^2),
            degf = (n_eff-1),
            vb12_inad_wdir =  weighted.mean(vb12_inadequate, norm_survey_wgt_state),
            vb12_inad_wdir_var = vb12_inad_wdir*(1-vb12_inad_wdir)/degf
  ) |>
  left_join((nga_analysis_df |> filter(lga!="bakassi") |> 
              select(state, lga) |> unique() |> 
              group_by(state) |> summarise(n_lga = n())),
            by = 'state'
            )

vb12_inad_state$vb12_inad_wdir_lower90_wilson <- wilson_lower(vb12_inad_state$vb12_inad_wdir, vb12_inad_state$n_eff,alpha = alpha)
vb12_inad_state$vb12_inad_wdir_upper90_wilson <- wilson_upper(vb12_inad_state$vb12_inad_wdir, vb12_inad_state$n_eff,alpha = alpha)
vb12_inad_state$vb12_inad_wdir_ci_length90_wilson <- vb12_inad_state$vb12_inad_wdir_upper90_wilson - vb12_inad_state$vb12_inad_wdir_lower90_wilson 
mean(vb12_inad_state$vb12_inad_wdir_ci_length90_wilson, na.rm = T)
vb12_inad_state$vb12_inad_wdir_lower90_wald <- (vb12_inad_state$vb12_inad_wdir - qnorm(1-alpha/2)*sqrt(vb12_inad_state$vb12_inad_wdir_var))
vb12_inad_state$vb12_inad_wdir_upper90_wald <- (vb12_inad_state$vb12_inad_wdir + qnorm(1-alpha/2)*sqrt(vb12_inad_state$vb12_inad_wdir_var))
vb12_inad_state$vb12_inad_wdir_ci_length90_wald<- vb12_inad_state$vb12_inad_wdir_upper90_wald - vb12_inad_state$vb12_inad_wdir_lower90_wald 
mean(vb12_inad_state$vb12_inad_wdir_ci_length90_wald, na.rm = T)

df_post <- cbind(vb12_inad_lga |> select(lga),t(df_p[,1:nrow(vb12_inad_lga)])) |> 
  left_join(lga_pop  |> select(state, lga, geo_weight) , 
            by = 'lga') |> filter(lga!="bakassi") 

df_post_state <- df_post |> as.data.frame() |> 
  select(-geometry) |> #filter(lga%in%lga_in_survey) |>  un-comment to exclude LGAs not included in the original survey
  group_by(state) |>
  summarise(across(all_of(colnames(df_post)[grepl('X',colnames(df_post))]), ~ weighted.mean(., w = geo_weight))) 
df_post_state$vb12_inad_post_mean <- df_post_state[,2:(nrow(df_p)+1)] |> rowMeans()
df_quantile =  apply(df_post_state[,2:nrow(df_p)], 1 , quantile , probs = c(alpha/2,1-alpha/2) , na.rm = TRUE ) |> t()
df_post_state$vb12_inad_post_lower90 <- df_quantile[,1]
df_post_state$vb12_inad_post_upper90 <- df_quantile[,2]
df_post_state$vb12_inad_posr_ci_length90 <- df_post_state$vb12_inad_post_upper90 - df_post_state$vb12_inad_post_lower90 
mean(df_post_state$vb12_inad_posr_ci_length90)
vb12_inad_state  <- vb12_inad_state |> left_join((df_post_state|> select(-colnames(df_post)[grepl('X',colnames(df_post))])), by = 'state')

mean(abs(vb12_inad_state$vb12_inad_wdir - vb12_inad_state$vb12_inad_post_mean))
vb12_inad_state$int_overlap <-check_overlap(vb12_inad_state$vb12_inad_wdir_lower90_wald,vb12_inad_state$vb12_inad_wdir_upper90_wald,
                                            vb12_inad_state$vb12_inad_post_lower90,vb12_inad_state$vb12_inad_post_upper90)
vb12_inad_state$error = vb12_inad_state$vb12_inad_post_mean - vb12_inad_state$vb12_inad_wdir 

mean(vb12_inad_state$int_overlap)
mean(abs(vb12_inad_state$error))
write_csv(vb12_inad_state, here("processed_data","state_level_estimates_withSynthetic.csv"))

#-------------------------------------------------------------------------------
# Plots
#-------------------------------------------------------------------------------
vb12_inad_state <- vb12_inad_state |> left_join(nigeria_1|> select(state, geometry), by='state')
vb12_inad_state <- st_as_sf(vb12_inad_state)
p = ggplot(vb12_inad_state |> st_as_sf()
           , aes(fill = vb12_inad_post_mean)) +
  geom_sf() +
  #facet_wrap(~Age) +
  scale_fill_gradientn(values = c("1","0.8","0.6","0.4","0.2","0"), 
                       colors = c("red4","red3","#FC4E07","#E7B800","green4"),
                       limits=c(0,1))+
  labs(title='Full Bayse estimate',fill = "VB12 deficiency")+
  theme_minimal(base_size = 14) 
p
ggsave(here("outputs","maps",paste0("vb12_inad_state_FB.png")), p, width = 10, height = 10,bg = "white")

p = ggplot(vb12_inad_state |> st_as_sf()
           , aes(fill = vb12_inad_wdir)) +
  geom_sf() +
  #facet_wrap(~Age) +
  scale_fill_gradientn(values = c("1","0.8","0.6","0.4","0.2","0"), 
                       colors = c("red4","red3","#FC4E07","#E7B800","green4"),
                       limits=c(0,1))+
  labs(title='Direct estimate',fill = "VB12 deficiency")+
  theme_minimal(base_size = 14) 
p
ggsave(here("outputs","maps",paste0("vb12_inad_state_direct.png")), p, width = 10, height = 10,bg = "white")

p = ggplot(vb12_inad_state |> st_as_sf()
           , aes(fill = int_overlap,label = str_to_title(state))) +
  geom_sf() +
  #geom_sf_text(fun.geometry = sf::st_centroid) + # make Manhattan behave itself
  labs(title='Coverage')+
  theme_minimal(base_size = 14) 
p
ggsave(here("outputs","maps",paste0("vb12_inad_state_coverage.png")), p, width = 10, height = 10,bg = "white")

p = ggplot(vb12_inad_state |> st_as_sf()
           , aes(fill = n_lga)) +
  geom_sf() +
  labs(title='Number of LGAs')+
  theme_minimal(base_size = 14) 
p

p = ggplot(vb12_inad_state |> st_as_sf()
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
ggsave(here("outputs","maps",paste0("vb12_inad_state_error.png")), p, width = 10, height = 10,bg = "white")

p <- ggplot(vb12_inad_state,aes(x=vb12_inad_post_mean,y=vb12_inad_wdir,label = str_to_title(state))) +
  geom_point(aes(color = int_overlap))+
  geom_abline(slope = 1,intercept = 0, linetype=2)+
  ggrepel::geom_text_repel()+
  xlim(0,1) + ylim(0,1)+
  labs(x='FB modelled estimates',y='Direct Estimates')+
  theme_minimal(base_size = 14)
p
ggsave(here("outputs","maps",paste0("vb12_inad_state_scatter.png")), p, width = 10, height = 10,bg = "white")

df_plot <- data.frame(
  y_est = c(vb12_inad_state$vb12_inad_wdir,vb12_inad_state$vb12_inad_post_mean),
  lower = c(vb12_inad_state$vb12_inad_wdir_lower90_wald, vb12_inad_state$vb12_inad_post_lower90),
  upper = c(vb12_inad_state$vb12_inad_wdir_upper90_wald, vb12_inad_state$vb12_inad_post_upper90),
  type = rep(c("Direct","Modelled"), each = nrow(vb12_inad_state)),
  state = rep(vb12_inad_state$state, times=2)
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

