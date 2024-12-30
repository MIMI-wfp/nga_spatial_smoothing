################################################################################
#################### SCRIPT FOR PREPARING POPULATION ESTIMATES #################
################################################################################

# Author: Mo Osman
# Date created: 27-12-2024
# Last edited: 30-12-2024

# In this script, I will extract 2020 population estimates at the ADM2 level in 
# Nigeria. Raw data provided by WFP Nigeria country office (Benjamin UMEH)

# INSTALL AND LOAD PACKAGES:

rq_packages <- c("readr", "tidyverse", "readxl", "sf")

installed_packages <- rq_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(rq_packages[!installed_packages])
}

lapply(rq_packages, require, character.only = T)

rm(list= c("rq_packages", "installed_packages"))

#-------------------------------------------------------------------------------

# READ IN DATA: 

# Get lga names as they are listed in the shapefiles:
lga <- st_read("shapefiles/nigeria_2") |> 
  dplyr::select(lga)

adm2_pop <- read_excel("population_estimates/nigeria_Admin2 Population.xlsx", 
                       sheet = "2000 - 2020")

# Remove first 2 rows from data-frame, and use 3rd row as column names: 
adm2_pop <- adm2_pop[-c(1:2),]
colnames(adm2_pop) <- adm2_pop[1,]
adm2_pop <- adm2_pop[-1,]

# Sense-check 2020 population estimates: 
adm2_pop <- adm2_pop |> 
  dplyr::select(COUNTRY, ADM1_NAME, ADM2_NAME, BTOTL_2020)
# Can see here that country population estimate for 2020 is 204,909,220
# This is approximately consistent with estimates from the WorldBank. 

#-------------------------------------------------------------------------------

# Select required columns and clean data: 
adm2_pop <- adm2_pop |>
  dplyr::select(ADM2_NAME, BTOTL_2020) |> 
  rename(lga = ADM2_NAME, population = BTOTL_2020) |>
  filter(!is.na(lga))

# Change entries in lga column to lowercase: 
adm2_pop$lga <- tolower(adm2_pop$lga)

# Check to see if any LGA's do not match:
lga |> 
  anti_join(adm2_pop, by = c("lga" = "lga"))

# There are 114 LGA's which did not match the data-set provided, will need to 
# look into each of these individually as there may be different reasons for not 
# matching (e.g. spelling errors, name changes, alternative names)

# Rename LGA's to match shapefile names:
adm2_pop$lga[adm2_pop$lga == "abaji"] <- "abaji area council"
adm2_pop$lga[adm2_pop$lga == "abua/odual"] <- "abua odua"
adm2_pop$lga[adm2_pop$lga == "ado-ekiti"] <- "ado ekiti"
adm2_pop$lga[adm2_pop$lga == "ado-odo/ota"] <- "ado odo/ota"
adm2_pop$lga[adm2_pop$lga == "ado"] <- "ador"
adm2_pop$lga[adm2_pop$lga == "ajeromi-ifelodun"] <- "ajeromi/ifelodun"
adm2_pop$lga[adm2_pop$lga == "akamkpa"] <- "akamkpa buyo"
adm2_pop$lga[adm2_pop$lga == "akoko-edo"] <- "akoko edo"
adm2_pop$lga[adm2_pop$lga == "akuku-toru"] <- "akuku toru"
adm2_pop$lga[adm2_pop$lga == "aliero"] <- "aleiro"
adm2_pop$lga[adm2_pop$lga == "amuwo-odofin"] <- "amuwo odofin"
adm2_pop$lga[adm2_pop$lga == "aninri"] <- "aniniri"
adm2_pop$lga[adm2_pop$lga == "arewa-dandi"] <- "arewa"
adm2_pop$lga[adm2_pop$lga == "asari-toru"] <- "asari toru"
adm2_pop$lga[adm2_pop$lga == "askira/uba"] <- "asikira/uba"
adm2_pop$lga[adm2_pop$lga == "atakunmosa east"] <- "atakumosa east"
adm2_pop$lga[adm2_pop$lga == "atakunmosa west"] <- "atakumosa west"
adm2_pop$lga[adm2_pop$lga == "atisbo"] <- "atigbo"
adm2_pop$lga[adm2_pop$lga == "ayedaade"] <- "ayedade"
adm2_pop$lga[adm2_pop$lga == "bassa" & adm2_pop$population == 192142] <- "bassa (kogi)"
adm2_pop$lga[adm2_pop$lga == "bassa"] <- "bassa (plateau)"
adm2_pop$lga[adm2_pop$lga == "birnin magaji/kiyaw"] <- "birnin magaji"
adm2_pop$lga[adm2_pop$lga == "boluwaduro"] <- "bolowaduro"
adm2_pop$lga[adm2_pop$lga == "bwari"] <- "bwari area council"
adm2_pop$lga[adm2_pop$lga == "calabar municipality"] <- "calabar municipal"
adm2_pop$lga[adm2_pop$lga == "damboa"] <- "damaboa"
adm2_pop$lga[adm2_pop$lga == "dambam"] <- "damban"
adm2_pop$lga[adm2_pop$lga == "dambatta"] <- "danbatta"
adm2_pop$lga[adm2_pop$lga == "dange-shuni"] <- "dange shuni"
adm2_pop$lga[adm2_pop$lga == "wasagu/danko"] <- "danko wasagu"
adm2_pop$lga[adm2_pop$lga == "dutsin-ma"] <- "dutsin ma"
adm2_pop$lga[adm2_pop$lga == "edati"] <- "edati idati"
adm2_pop$lga[adm2_pop$lga == "yewa north"] <- "egbado north/yewa"
adm2_pop$lga[adm2_pop$lga == "yewa south"] <- "egbado south/"
adm2_pop$lga[adm2_pop$lga == "ekiti south west"] <- "ekiti south"
adm2_pop$lga[adm2_pop$lga == "emuoha"] <- "emohu"
adm2_pop$lga[adm2_pop$lga == "ese-odo"] <- "ese odo"
adm2_pop$lga[adm2_pop$lga == "eti-osa"] <- "eti osa"
adm2_pop$lga[adm2_pop$lga == "ezinihitte"] <- "ezinihitte mbaise"
adm2_pop$lga[adm2_pop$lga == "fufore"] <- "fufore/gurin"
adm2_pop$lga[adm2_pop$lga == "garun malam"] <- "garum mallam"
adm2_pop$lga[adm2_pop$lga == "gwagwalada"] <- "gwagwalada area council"
adm2_pop$lga[adm2_pop$lga == "ibeju-lekki"] <- "ibeju lekki"
adm2_pop$lga[adm2_pop$lga == "ido-osi"] <- "ido/osi"
adm2_pop$lga[adm2_pop$lga == "ifako-ijaye"] <- "ifako ijaye"
adm2_pop$lga[adm2_pop$lga == "ifelodun" & adm2_pop$population == 235490] <- "ifelodun (kwara)"
adm2_pop$lga[adm2_pop$lga == "igalamela/odolu"] <- "igalamela odolu"
adm2_pop$lga[adm2_pop$lga == "igueben"] <- "igugben"
adm2_pop$lga[adm2_pop$lga == "ihitte/uboma"] <- "ihitte uboma"
adm2_pop$lga[adm2_pop$lga == "ijebu-ode"] <- "ijebu ode"
adm2_pop$lga[adm2_pop$lga == "ikorodu"] <- "ikorordu"
adm2_pop$lga[adm2_pop$lga == "ikpoba-okha"] <- "ikpooba okha"
adm2_pop$lga[adm2_pop$lga == "ile-oluji/oke-igbo"] <- "ileoluji/okeigbo"
adm2_pop$lga[adm2_pop$lga == "imeko-afon"] <- "imeko/afon"
adm2_pop$lga[adm2_pop$lga == "irepodun" & adm2_pop$population == 176401] <- "irepodun (kwara)"
adm2_pop$lga[adm2_pop$lga == "irepodun"] <- "irepodun (osun)"
adm2_pop$lga[adm2_pop$lga == "isuikwuato"] <- "isuikwato"
adm2_pop$lga[adm2_pop$lga == "itas/gadau"] <- "itas gadau"
adm2_pop$lga[adm2_pop$lga == "jama’are"] <- "jama'are"
adm2_pop$lga[adm2_pop$lga == "kabba bunu"] <- "kabba/bunu"
adm2_pop$lga[adm2_pop$lga == "kibiya"] <- "kabiya"
adm2_pop$lga[adm2_pop$lga == "kala/balge"] <- "kala balge"
adm2_pop$lga[adm2_pop$lga == "katsina-ala"] <- "katsina ala"
adm2_pop$lga[adm2_pop$lga == "kaura-namoda"] <- "kaura namoda"
adm2_pop$lga[adm2_pop$lga == "kiri kasamma"] <- "kirika kasamma"
adm2_pop$lga[adm2_pop$lga == "kogi"] <- "kogi(k.k)"
adm2_pop$lga[adm2_pop$lga == "koko/besse"] <- "koko besse"
adm2_pop$lga[adm2_pop$lga == "kuje"] <- "kuje area council"
adm2_pop$lga[adm2_pop$lga == "kwali"] <- "kwali area council"
adm2_pop$lga[adm2_pop$lga == "langtang north"] <- "lantang north"
adm2_pop$lga[adm2_pop$lga == "langtang south"] <- "lantang south"
adm2_pop$lga[adm2_pop$lga == "maiduguri"] <- "maiduguri metropolitan"
adm2_pop$lga[adm2_pop$lga == "maigatari"] <- "maigatar"
adm2_pop$lga[adm2_pop$lga == "malumfashi"] <- "malunfashi"
adm2_pop$lga[adm2_pop$lga == "mayo-belwa"] <- "mayo belwa"
adm2_pop$lga[adm2_pop$lga == "mopa-muro"] <- "mopamuro"
adm2_pop$lga[adm2_pop$lga == "abuja municipal"] <- "municipal area council" # Appears to be Abuja Municipal area council
adm2_pop$lga[adm2_pop$lga == "nassarawa"] <- "nasarawa (kano)"
adm2_pop$lga[adm2_pop$lga == "nasarawa"] <- "nasarawa (nasarawa)"
adm2_pop$lga[adm2_pop$lga == "nasarawa egon"] <- "nasarawa eggon"
adm2_pop$lga[adm2_pop$lga == "ngor-okpala"] <- "ngor okpala"
adm2_pop$lga[adm2_pop$lga == "nkwerre"] <- "nkwere"
adm2_pop$lga[adm2_pop$lga == "esit eket"] <- "nsit eket"
adm2_pop$lga[adm2_pop$lga == "obafemi-owode"] <- "obafemi owode"
adm2_pop$lga[adm2_pop$lga == "obi" & adm2_pop$population == 134000] <- "obi (benue)"
adm2_pop$lga[adm2_pop$lga == "obi"] <- "obi (nasarawa)"
adm2_pop$lga[adm2_pop$lga == "obio/akpor"] <- "obio akpor"
adm2_pop$lga[adm2_pop$lga == "ogu/bolo"] <- "ogu bolo"
adm2_pop$lga[adm2_pop$lga == "ohaji/egbema"] <- "ohaji egbema"
adm2_pop$lga[adm2_pop$lga == "oke-ero"] <- "oke ero"
adm2_pop$lga[adm2_pop$lga == "okrika"] <- "okirika"
adm2_pop$lga[adm2_pop$lga == "omuma"] <- "omumma"
adm2_pop$lga[adm2_pop$lga == "onuimo"] <- "ono imo"
adm2_pop$lga[adm2_pop$lga == "oorelope"] <- "orelope"
adm2_pop$lga[adm2_pop$lga == "oriire"] <- "ori ire"
adm2_pop$lga[adm2_pop$lga == "oshodi-isolo"] <- "oshodi/isolo"
adm2_pop$lga[adm2_pop$lga == "osisioma ngwa"] <- "osisioma north"
adm2_pop$lga[adm2_pop$lga == "sabon birni"] <- "sabon birnin"
adm2_pop$lga[adm2_pop$lga == "shagamu"] <- "sagamu"
adm2_pop$lga[adm2_pop$lga == "somolu"] <- "shomolu"
adm2_pop$lga[adm2_pop$lga == "sule tankarkar"] <- "sule tankar kar"
adm2_pop$lga[adm2_pop$lga == "surulere" & adm2_pop$population == 1047767] <- "surulere (lagos)"
adm2_pop$lga[adm2_pop$lga == "surulere"] <- "surulere (oyo)"
adm2_pop$lga[adm2_pop$lga == "tambuwal"] <- "tambuwai"
adm2_pop$lga[adm2_pop$lga == "abi"] <- "ugep south abi"
adm2_pop$lga[adm2_pop$lga == "uhunmwonde"] <- "uhunmuonde"
adm2_pop$lga[adm2_pop$lga == "ukwuani"] <- "ukwani"
adm2_pop$lga[adm2_pop$lga == "wamako"] <- "wamakko"
adm2_pop$lga[adm2_pop$lga == "yabo"] <- "yabo bodingo"
adm2_pop$lga[adm2_pop$lga == "yakurr"] <- "yakurr ugep north"
adm2_pop$lga[adm2_pop$lga == "yala"] <- "yalla"
adm2_pop$lga[adm2_pop$lga == "yamaltu/deba"] <- "yamaltu deba"
adm2_pop$lga[adm2_pop$lga == "yankwashi"] <- "yan kwashi"



#-------------------------------------------------------------------------------

# Left join adm2_pop to lga to see if there are any LGAs that fail to join: 
lga <- left_join(lga, adm2_pop, by = "lga")

# How many NA's in the population column: 
lga <- lga |> filter(is.na(population))

# Only LGA which has not joined is Bakassi, this is as expected due to the fact 
# that it is not included in the population data because sovreignty of this LGA
# has been transferred to Cameroon as of 2007.

sum(is.na(adm2_pop$population))

rm(lga)

#-------------------------------------------------------------------------------

# WRITE DATA: 
write_csv(adm2_pop, "population_estimates/adm2_pop_clean.csv")

# Clear environment: 
rm(list = ls())

################################################################################
################################ END OF SCRIPT #################################
################################################################################











