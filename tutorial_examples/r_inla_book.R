# install.packages("SpatialEpi")
library(SpatialEpi)

map <- pennLC$spatial.polygon
plot(map)

# Get neightbours list: 
library(spdep)
nb <- poly2nb(map, queen = TRUE)
head(nb)

# Aasasign county names
d <- data.frame(county = names(map), neigh = rep(0, length(map)))
rownames(d) <- names(map)
map <- SpatialPolygonsDataFrame(map, d, match.ID = TRUE)

# Add longitude and lattitude coordinates:
coord <- coordinates(map)
map$long <- coord[, 1]
map$lat <- coord[, 2]
map$ID <- 1:dim(map@data)[1]

# Convert the map to a simple feature object: 
library(sf)
mapsf <- st_as_sf(map)

# Summarise counts of incidence of lung cancer in Pensylvania by county:
library(dplyr)
d <- group_by(pennLC$data, county) %>% summarize(Y = sum(cases))
head(d)

# Specify strata: 
pennLC$data <- pennLC$data[order(
  pennLC$data$county,
  pennLC$data$race,
  pennLC$data$gender,
  pennLC$data$age
), ]

# Obtain expected counts in each county:
E <- expected(
  population = pennLC$data$population,
  cases = pennLC$data$cases, n.strata = 16
)

# List expected and observed counts for each county, making sure that values 
# correspond to the correct county: 
d$E <- E[match(d$county, unique(pennLC$data$county))]
head(d)

# Compute SIR (standardised incidence ratio) values: 
d$SIR <- d$Y / d$E
head(d)

# Merge SIR values with the map:
map <- merge(map, d)
mapsf <- st_as_sf(map)

# Produce map: 
ggplot(mapsf) + geom_sf(aes(fill = SIR)) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red"
  ) +
  theme_bw()

# Specify BYM model forumla: 
formula <- Y ~
  f(idareau, model = "besag", graph = g, scale.model = TRUE) +
  f(idareav, model = "iid")

map$idareau <- 1:nrow(map@data)
map$idareav <- 1:nrow(map@data)

# Compute neighbourhood matrix: 
library(spdep)
library(INLA)
nb <- poly2nb(map)
head(nb)

nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")

res <- inla(formula,
            family = "poisson", data = map@data,
            E = E, control.predictor = list(compute = TRUE)
)

summary(res)
