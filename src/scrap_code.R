#-------------------------------------------------------------------------------

# CREATE ADM2 SPATIAL ADJACENCY MATRIX: 

# Amat <- getAmat(nigeria_2, nigeria_2$lga)

# # Firstly extract vector of LGA names: 
# row.names <- nigeria_2$lga
# 
# # Create neighbours list: 
# nigeria_2_neigh <- poly2nb(nigeria_2, 
#                            row.names = row.names,
#                            queen = TRUE)
# 
# # Use neightbours list to create adjacency matrix:
# nigeria_2_adj <- nb2mat(nigeria_2_neigh, style = "B")
# nigeria_2_adj <- as.matrix(nigeria_2_adj)
# 
# # Set column names of adjacency matrix to equal row names:
# colnames(nigeria_2_adj) <- row.names
# 
# # Tidy environment:
# rm(nigeria_2_neigh, row.names)

#-------------------------------------------------------------------------------

# CALCULATE SMOOTHED ESTIMATES AT ADM2 LEVEL - BYM2 MODEL (R-INLA):

# Filter vb12_inad_lga to only include LGAs with estimates: 
vb12_inad_lga <- vb12_inad_lga %>% 
  filter(!is.na(vb12_inad))

# Set vector of indices of spatially structured component:
vb12_inad_lga$idareau <- 1:nrow(vb12_inad_lga)

# Set a vector of indices of unstructured component:
vb12_inad_lga$idareav <- 1:nrow(vb12_inad_lga)


# # Assign penalised complexity (PC) priors:
# prior <- list(prec = list(prior = "pc.prec",
#                           param = c(0.5 / 0.31, 0.01)),
#               phi = list(prior = "pc",
#                          param = c(0.5, 2 / 3)))

# Can revisit these priors later - refer to Simpson et al. (2017) for guidance.

# Compute a neighbourhood matrix - in the format required by R-INLA:
vb12_inad_lga <- poly2nb(vb12_inad_lga, queen = TRUE)
nb2INLA("outputs/nigeria_2.adj", vb12_inad_lga)
g <- inla.read.graph(filename = "outputs/nigeria_2.adj")

# Specify the formula for the BYM2 model:
formula <- vb12_inadequate ~ f(idareau, model = "besag", graph = g, scale.model = T) +
  f(idareav, model = "iid")

# Fit
smoothed <- inla(formula = formula,
                 family = "poisson",
                 data = vb12_inad_lga,
                 control.predictor = list(compute = TRUE))

# ?? vb12_inadequate needs to be a vector of survey weighted prevalence. 
# Currently, the model takes issue with the response variable not being the same
# length as idarea. Need to ensure that the indexing is correct.

# # Fit the BYM2 model:
# fit <- inla(formula = formula,
#             family = "poisson",
#             data = nga_analysis_df,
#             control.predictor = list(compute = TRUE))

#-------------------------------------------------------------------------------