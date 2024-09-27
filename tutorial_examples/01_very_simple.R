#' Fitting a simple glm model with inlabru
#'=========================================================
#'
#+results="hide",warning=FALSE,message=FALSE

library(inlabru)
library(lme4)

#' Set default INLA values for tutorial session. (We will explain these later - 
#' just take them as given for now.)

init.tutorial()

#' The aim of the practical
#'-----------------------------------
#' In this practical we are going to fit a simple Poisson regression model to count data.

#' We have a simple data set of the number of awards earned by students at a high school in a year, math is
#' continuous predictor variable and represents students’ scores on their math final exam, and prog is a
#' categorical predictor variable with three levels indicating the type of program in which the students were
#' enrolled. It is coded as 1 = “General”, 2 = “Academic” and 3 = “Vocational”. L
#'
#' Get the data
#'-----------------------------------

load("data/awards.RData")

#' We need to speficy a model in order to fit it.
#' 
#' Our response variable in the data frame `awards` is called `num_awards` so the model specification
#' needs to have that on the left of the `~`.
#' The right has to have the covariate(s) (here `math`) and we also explicitly state `+ Intercept`
#' (most, but not all, of the models we use have explicit intercept parameters),

#' Specify the model:
#' 
#+warning=FALSE,message=FALSE

cmp1 <- num_awards ~ math + Intercept
# cmp1 <- num_awards ~ math + 1 # This works too

#' Now fit the model with `bru()`.
#+warning=FALSE,message=FALSE

fit.glm.bru <- bru(cmp1, family = "poisson", data = awards)
summary(fit.glm.bru)

#' Now fit the same model with `glm()` for comparison.
#+warning=FALSE,message=FALSE

cmp2 <- num_awards ~ math
fit2.glm <- glm(cmp2, family="poisson", data = awards)
summary(fit2.glm)

#' There may be a slight issue with overdispersion here, so we include
#' an individual level random effect here.
#+warning=FALSE,message=FALSE

cmp3 <- num_awards ~ math + Intercept + rand.eff(map = 1:200, model = "iid", n = 200)

#' Now fit the model with `bru()`.
#+warning=FALSE,message=FALSE
fit3.glmm.bru <- bru(cmp3, family = "poisson", data = awards )

#' Use
#+results=FALSE

summary(fit3.glmm.bru) ## Experimental inlabru summary output
INLA:::summary.inla(fit3.glmm.bru) ## Standard inla summary output

#' to look at the results.  Note that for the random effect, only some summary
#' information is shown.  Later, we will use the `predict()` function to analyse
#' estimates in more detail.
#'
#' The random effects `inlabru` estimate is similar to the following:
#+warning=FALSE,message=FALSE,results=FALSE

cmp4<- num_awards ~ math + (1|id)
fit4.glmm<-glmer(cmp4, family = poisson, data = awards)
summary(fit4.glmm)

#' You'll notice that the random effect precision estimate from INLA is very
#' different from the `glmer` estimate. However, if instead of the default
#' prior precision model $Gamma(1,10^{-5})$ we use a PC-prior, the results
#' are comparable. This shows that it's important to consider priors, and
#' not always unquestioningly use default software parameters; the priors
#' are part of the model, and general purpose software authors cannot know
#' what is appropriate for your particular problem, so the user of the
#' software is ultimately responsible for modelling choices.
#+warning=FALSE,message=FALSE,results=FALSE

cmp5 <- num_awards ~ math + Intercept + rand.eff(map = 1:200, model = "iid", n = 200,
                                                 hyper=list(prec=list(param=c(10,0.1),
                                                                      prior="pc.prec")))
fit5.glmm.bru <- bru(cmp5, family = "poisson", data = awards )
summary(fit5.glmm.bru)
INLA:::summary.inla(fit5.glmm.bru)


#' Additional questions:
#' Estimate 95% credible interval for (a) the linear predictor, (b) the Poisson parameter, and 
#' also predict the number of awards of someone with a math score of 55
# No random effect:
lp <- predict(fit.glm.bru, data.frame(math=55), ~ math + Intercept, n=500) # linear predictor
lp
lambda <- predict(fit.glm.bru, data.frame(math=55), ~ exp(math + Intercept)) # Poisson parameter
lambda
numawards = 0:10
post = data.frame(numawards = numawards, pdf = dpois(numawards, lambda=lambda$mean))
ggplot(post) +
  geom_line(aes(x = numawards, y = pdf), lty=2) +
  geom_point(aes(x = numawards, y = pdf))
# Here's some code to give you all posterior samples (see "help(generate.bru)"):
vals = generate(fit.glm.bru, data.frame(math=55), ~ exp(math + Intercept),n=500)
vals = unlist(vals)
hist(vals,nclass=30)
quantile(vals,probs=c(0.025,0.5,0.975))



#' with PC prior random effect (predicting with random effect: 
#' Refit model fit5.gll.bru, but using column "id" as the random effect variable:
cmp5a <- num_awards ~ math + Intercept + rand.eff(map = id, model = "iid", n = 200,
                                                  hyper=list(prec=list(param=c(10,0.1),
                                                                       prior="pc.prec")))
fit5a.glmm.bru <- bru(cmp5a, family = "poisson", data = awards )
# Check that you get the same model (you do!)
summary(fit5a.glmm.bru)
INLA:::summary.inla(fit5a.glmm.bru)
# Now predict for all n values of the id variable:
n=dim(awards)[1]
lambda5a <- predict(fit5a.glmm.bru, data.frame(math=rep(55,n), id=1:n), ~ exp(math + Intercept + rand.eff))
lambda5a

