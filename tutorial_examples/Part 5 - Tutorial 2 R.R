
## ------------------------------------------------------------------------
library(SUMMER)
library(sp)
library(ggplot2)
library(gridExtra)
data(DemoData2)
head(DemoData2)


## ------------------------------------------------------------------------
naive <- aggregate(tobacco.use ~ region, data = DemoData2, FUN = mean)


## ------------------------------------------------------------------------
data(DemoMap2)
smoothed <- fitGeneric(data = DemoData2, geo = DemoMap2$geo, Amat = DemoMap2$Amat, responseType = "binary", responseVar = "tobacco.use", regionVar = "region", strataVar = NULL, weightVar = NULL, clusterVar = NULL, CI = 0.95)


## ------------------------------------------------------------------------
head(smoothed$smooth)
head(smoothed$HT)


## ------------------------------------------------------------------------
library(survey)
design <- svydesign(ids = ~clustid+id, weights = ~weights, strata = ~strata, data = DemoData2) 
direct <- svyby(~tobacco.use, ~region, design, svymean)
head(direct)


## ------------------------------------------------------------------------
smoothweighted <- fitGeneric(data = DemoData2, geo = DemoMap2$geo, Amat = DemoMap2$Amat, responseType = "binary", responseVar = "tobacco.use", regionVar = "region", strataVar = "strata", weightVar = "weights", clusterVar = "~clustid+id", CI = 0.95)


## ------------------------------------------------------------------------
head(smoothweighted$HT)
head(smoothweighted$smooth)


## ---- -------------------------------------------------------
prev <- NULL
prev <- rbind(prev, data.frame(region = smoothed$HT$region, 
	                           mean = smoothed$HT$HT.est.original,
	                           var = smoothed$HT$HT.variance.original,
	                           type = "Naive"))
prev <- rbind(prev, data.frame(region = smoothed$smooth$region, 
							   mean = smoothed$smooth$mean.original, 
							   var = smoothed$smooth$variance.original,
							   type = "Smoothed"))
prev <- rbind(prev, data.frame(region = smoothweighted$HT$region, 
							   mean = smoothweighted$HT$HT.est.original,
							   var = smoothweighted$HT$HT.variance.original,
							   type = "Weighted"))
prev <- rbind(prev, data.frame(region = smoothweighted$smooth$region, 
							   mean = smoothweighted$smooth$mean.original, 
							   var = smoothweighted$smooth$variance.original,
							   type = "Smooth Weighted"))


## ------------------------------------------------------------------------
g1 <- mapPlot(prev, geo = DemoMap2$geo, by.data = "region", by.geo = "REGNAME", variables = "type", values = "mean", is.long=TRUE, legend.label = "Estimates", ncol = 4)
g2 <- mapPlot(prev, geo = DemoMap2$geo, by.data = "region", by.geo = "REGNAME", variables = "type", values = "var", is.long=TRUE, legend.label = "Variance", ncol = 4)
grid.arrange(g1, g2, ncol = 1)


## ------------------------------------------------------------------------
data(DemoData)
data <- DemoData[[1]]
head(data)


## ------------------------------------------------------------------------
data <- subset(data, age == 0)


## ------------------------------------------------------------------------
data$time.id <- match(data$time, levels(data$time))
fit <- fitGeneric(data = data, geo = DemoMap$geo, Amat = DemoMap$Amat, responseType = "binary", responseVar = "died", regionVar = "region", strataVar = "strata", weightVar = "weights", clusterVar = "~clustid+id", CI = 0.95, timeVar = "time.id", time.model = "rw1",type.st = 4)


## ---- -------------------------------------------------------
nmr <- NULL
nmr <- rbind(nmr, data.frame(region = fit$HT$region, 
							   time = fit$HT$time,
							   mean = fit$HT$HT.est.original, 
							   type = "Direct"))
nmr <- rbind(nmr, data.frame(region = fit$smooth$region, 
							   time = fit$smooth$time, 
							   mean = fit$smooth$mean.original, 
							   type = "Smoothed"))
nmr$time <- levels(data$time)[nmr$time]
nmr$time <- factor(nmr$time, levels = levels(data$time))


## ---- fig.width = 10, fig.height = 4-------------------------------------
ggplot(nmr, aes(x = time, y = mean, color = type, group = type)) + geom_point() + geom_line() + facet_wrap(~region, ncol = 4) + theme(legend.position = "bottom")


## ---- warning=FALSE------------------------------------------------------
years <- levels(DemoData[[1]]$time)
data_multi <- getDirectList(births = DemoData, years = years,regionVar = "region",  timeVar = "time", clusterVar = "~clustid+id", ageVar = "age", weightsVar = "weights", geo.recode = NULL)


## ------------------------------------------------------------------------
data <- aggregateSurvey(data_multi)
dim(data)


## --------------------------------------------------------
years.all <- c(years, "15-19")
fit1 <- fitINLA(data = data, geo = NULL, Amat = NULL, year_label = years.all, rw = 2, is.yearly=FALSE)


## --------------------------------------------------------
fit2 <- fitINLA(data = data, geo = NULL, Amat = NULL, year_label = years.all, year_range = c(1985, 2019), rw = 2, is.yearly=TRUE, m = 5)


## ------------------------------------------------------------------------
out1 <- getSmoothed(fit1)
out2 <- getSmoothed(fit2)


## -----------------------------------------------------------
g <- NULL
ylim <- range(c(out2$lower, out2$upper))
g[[1]] <- plot(out1, is.subnational=FALSE) + ggtitle("National period model") + ylim(ylim)
g[[2]] <- plot(out2, is.subnational=FALSE) + ggtitle("National yearly model") + ylim(ylim)
grid.arrange(grobs=g, ncol = 2)


## --------------------------------------------------------
fit3 <- fitINLA(data = data, geo = DemoMap$geo, Amat = DemoMap$Amat, year_label = years.all, rw = 2, is.yearly=FALSE, type.st = 4)
out3 <- getSmoothed(fit3, Amat = DemoMap$Amat)


## -------------------------------------------------------
fit4 <- fitINLA(data = data, geo = DemoMap$geo, Amat = DemoMap$Amat, year_label = years.all, year_range = c(1985, 2019), rw = 2, is.yearly=TRUE, m = 5, type.st = 4)
out4 <- getSmoothed(fit4, Amat = DemoMap$Amat)


## -----------------------------------------------------------
g2 <- NULL
ylim <- range(c(out3$median, out4$median))
g2[[1]] <- plot(out3, is.yearly=FALSE, is.subnational=TRUE) + 
			ggtitle("Subnational period model") + ylim(ylim)
g2[[2]] <- plot(out4, is.yearly=TRUE, is.subnational=TRUE) + 
			ggtitle("Subnational yearly model")+ ylim(ylim)
grid.arrange(grobs=g2, ncol = 2)


## ------------------------------------------------------------------------
plot(out4, is.yearly=TRUE, is.subnational=FALSE,  data.add = data_multi, option.add = list(point = "mean", by = "surveyYears"), color.add = "blue", plot.CI = TRUE, alpha.CI = 0.2) + facet_wrap(~region) 



## ----------------------------------------------
mapPlot(data = subset(out4, is.yearly==F), geo = DemoMap$geo,
        variables=c("years"), values = c("median"), by.data = "region", by.geo = "NAME_final", is.long=TRUE, ncol = 4, per1000 = TRUE, legend.label = "U5MR")

