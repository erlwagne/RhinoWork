## RHAU CODE

## Set graphics device
if(.Platform$OS.type == "windows") options(device = windows) else options(device = quartz)

library(rstanarm)
library(bayesplot)
library(Hmisc)
library(sm)
library(loo)
library(denstrip)
library(yarrr)
library(corrplot)

#--------------------------
# DATA  
#--------------------------

data.path <- 

# Nest success data
nest_data <- read.csv(file.path("~","R","RhinoWork","data","RhAu2.csv"), fileEncoding="UTF-8-BOM") ## Nest data at site level
rhau_agg <- read.csv(file.path("~","R","RhinoWork","data","RhAu4.csv"), header = TRUE)

# Merge environmental covariates into site-level nest data
env_data <- subset(rhau_agg, Island == "DI", 
                   select = c(Year, sst_spring, sst_summer, mei_avg, cu_spring, 
                              cu_summer, st_onset, st_length, pdo_index))
names(env_data) <- gsub("sst_", "sst_DI_", names(env_data))
env_data <- data.frame(env_data, 
                       sst_PI_spring = rhau_agg$sst_spring[rhau_agg$Island=="PI"],
                       sst_PI_summer = rhau_agg$sst_summer[rhau_agg$Island=="PI"])

#---------------------------------
# PCA of Oceanographic Predictors
#---------------------------------

pca_env <- princomp(~ sst_DI_spring + sst_PI_spring + mei_avg + cu_spring + st_onset + pdo_index, 
                    data = env_data, cor = TRUE)
dev.new()
par(mfrow = c(3,1))
plot(pca_env)  # scree plot
biplot(pca_env) # biplot of PCAs and oceanographic indicators
barplot(pca_env$loadings[,1], names.arg = dimnames(pca_env$loadings)[[1]], main = "PC1 loadings")
pca_env$loadings 
pca_env$sdev / sum(pca_env$sdev)

# Add PC1 and PC2 environmental data
env_data <- data.frame(env_data, PC1 = scale(pca_env$scores[,1]), PC2 = scale(pca_env$scores[,2]))
rhau <- merge(nest_data, env_data, all.x = TRUE)
summary(rhau)

#---------------------------------
# GLMMs
#---------------------------------

## Hierarchical Bayesian binomial models with rstanarm ##
## Reproductive Success ##

## Random effects of site and year

# Intercept-only model 
mod0 <- stan_glmer(cbind(LastCheck,Egg-LastCheck) ~ (1 | Site) + (1 | Year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(mod0, prob = c(0.025,0.5,0.975))
launch_shinystan(mod0)

# Inter-island differences, constant across years
mod1 <- stan_glmer(cbind(LastCheck,Egg-LastCheck) ~ Island + (1 | Site) + (1 | Year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(mod1, prob = c(0.025,0.5,0.975))
launch_shinystan(mod1)

# Inter-island differences, varying among years
mod2 <- stan_glmer(cbind(LastCheck,Egg-LastCheck) ~ Island + (1 | Site) + (Island | Year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(mod2, prob = c(0.025,0.5,0.975))
launch_shinystan(mod2)

# Inter-island differences, varying among years, plus first PC
mod3 <- stan_glmer(cbind(LastCheck,Egg-LastCheck) ~ Island + PC1 + (1 | Site) + (Island | Year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(mod3, prob = c(0.025,0.5,0.975))
launch_shinystan(mod3)

mod4 <- stan_glmer(cbind(LastCheck,Egg-LastCheck) ~ Island + PC2 + (1 | Site) + (Island | Year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(mod4, prob = c(0.025,0.5,0.975))
launch_shinystan(mod4)

mod5 <- stan_glmer(cbind(LastCheck,Egg-LastCheck) ~ Island + PC1 + (1 | Site),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(mod5, prob = c(0.025,0.5,0.975))
launch_shinystan(mod5)

mod6 <- stan_glmer(cbind(LastCheck,Egg-LastCheck) ~ Island + PC2 + (1 | Site),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(mod6, prob = c(0.025,0.5,0.975))
launch_shinystan(mod6)

## Model selection by approximate leave-one-out cross-validation

loos <-  list(mod0 = loo(mod0), mod1 = loo(mod1), mod2 = loo(mod2), mod3 = loo(mod3), mod4 = loo(mod4),
              mod5 = loo(mod5), mod6 = loo(mod6))
mod_sel <- loo::compare(loos$mod0, loos$mod1, loos$mod2)
mod_sel_full <- loo::compare(loos$mod0, loos$mod1, loos$mod2, loos$mod3, loos$mod4, loos$mod5, loos$mod6)
mod2vs1 <- loo::compare(loos$mod2, loos$mod1)
mod1vs0 <- loo::compare(loos$mod1, loos$mod0)
mod2vs0 <- loo::compare(loos$mod2, loos$mod0)
mod2vs3 <- loo::compare(loos$mod2, loos$mod3)
mod2vs4 <- loo::compare(loos$mod2, loos$mod4)
mod3vs4 <- loo::compare(loos$mod3, loos$mod4)
mod3vs5 <- loo::compare(loos$mod3, loos$mod5)
mod4vs6 <- loo::compare(loos$mod4, loos$mod6)

mod_sel
mod2vs1
mod1vs0
mod2vs0
mod2vs3
mod2vs4
mod3vs4
mod3vs5
mod4vs6

#-------------------
# FIGURES
#-------------------

## Correlation plot of covariates
dev.new()
<<<<<<< HEAD:analysis/RhAu_BayesGLMM_script.R
dat <- subset(rhau, select = c(sst_spring, sst_summer, mei_avg, cu_spring, 
                                 cu_summer, st_onset, st_length, pdo_index))
=======
dat <- subset(rhau, select = c(sst.DI.spring, sst.PI.spring, mei.avg, cu.spring, 
                              st.onset, pdo.index))
colnames(dat) <- c("SST spring (DI)", "SST spring (PI)", "MEI", "Coastal Upwelling (Spring)", 
                   "Spring Transition (Onset)", "PDO Index")
>>>>>>> master:analysis/rhau_BayesGLMM2.R
corrplot(cor(dat), method = "ellipse")

## Model validation by graphical posterior predictive checking
# default plot: marginal distribution of data and posterior predictive distribution
pp_check(mod2)

## Posterior distribution of average between-island difference in full model,
## with median and 95% credible interval
dev.new()
par(mar = c(5.1,4.3,4.1,1))

x <- as.matrix(mod2, pars = "IslandPI")
px <- sm.density(x, display = "none")
px <- sm.density(x, eval.points = sort(c(px$eval.points, quantile(x, c(0.025,0.5,0.975)))), display = "none")
ci <- which(names(px$eval.points)=="2.5%"):which(names(px$eval.points)=="97.5%")
plot(px$eval.points, px$estimate, col = "darkgray", type = "l", lwd = 3,
     las = 1, cex.lab = 1.5, cex.axis = 1.2, bty = "n", 
     ylim = c(0,max(px$estimate)), yaxs = "i",
     xlab = "Fledging log-odds ratio: Protection vs. Destruction", ylab = "Probability density",
     main = "Posterior median and 95% credible interval")
polygon(c(px$eval.points[ci], rev(px$eval.points[ci])), c(px$estimate[ci], rep(0,length(ci))),
        col = transparent("darkgray",0.7), border = NA)
segments(x0 = px$eval.points[names(px$eval.points)=="50%"], y0 = 0, 
         y1 = px$estimate[names(px$eval.points)=="50%"], col = "darkgray", lwd = 2)


## Time series of fitted and observed fledging success at DI and PI
## (This will take some time to render. For alternate version, 
## comment out both calls to densregion() and uncomment polygon()).
dev.new(width = 10, height = 7)
par(mar = c(5.1,4.3,4.1,1))

YY <- sort(unique(rhau$Year))
newdata <- data.frame(Year = rep(YY, each = 2), Site = "newsite",
                      Island = rep(c("DI","PI"), length(YY)))
pfit <- posterior_linpred(mod2, transform = TRUE, re.form = ~ (1 | Year), newdata = newdata)
pobs <- aggregate(cbind(Egg, LastCheck) ~ Year + Island, data = rhau, sum)
pobs_ci <- binconf(pobs$LastCheck, n = pobs$Egg)
eval.points <- range(pobs_ci, apply(pfit,2,quantile,c(0.025,0.975)))
eval.points <- seq(min(eval.points), max(eval.points), length = 300)
xyDI <- pfit[,newdata$Island=="DI"]
pxyDI <- t(apply(xyDI, 2, function(x) 
  sm.density(x, eval.points = eval.points, display = "none")$estimate))
xyPI <- pfit[,newdata$Island=="PI"]
pxPI <- sm.density(as.vector(xyPI), display = "none")
pxyPI <- t(apply(xyPI, 2, function(x) 
  sm.density(x, eval.points = eval.points, display = "none")$estimate))

plot(newdata$Year[newdata$Island=="PI"], apply(pfit[,newdata$Island=="PI"],2,median), 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, type = "l", lwd = 3, col = "cornflowerblue",
     xlab = "Year", ylab = "Probability of fledging", 
     ylim = range(pobs_ci, apply(pfit,2,quantile,c(0.025,0.975))), 
     xaxs = "i", xaxt = "n", yaxs = "i")
lines(newdata$Year[newdata$Island=="DI"], apply(pfit[,newdata$Island=="DI"],2,median),
      lwd = 3, col = "darkgray")
axis(1, at = min(rhau$Year):max(rhau$Year))

#densregion(newdata$Year[newdata$Island=="PI"], eval.points, pxyPI, 
#          colmax = transparent("darkgray",0.1), colmin = "transparent")
polygon(c(newdata$Year[newdata$Island=="PI"], rev(newdata$Year[newdata$Island=="PI"])),
        c(apply(pfit[,newdata$Island=="PI"],2,quantile,0.025), 
          rev(apply(pfit[,newdata$Island=="PI"],2,quantile,0.975))), 
        col = transparent("cornflowerblue",0.7), border = NA)
#densregion(newdata$Year[newdata$Island=="DI"], eval.points, pxyDI, 
#           colmax = transparent("black",0.3), colmin = "transparent")
polygon(c(newdata$Year[newdata$Island=="DI"], rev(newdata$Year[newdata$Island=="DI"])),
        c(apply(pfit[,newdata$Island=="DI"],2,quantile,0.025), 
          rev(apply(pfit[,newdata$Island=="DI"],2,quantile,0.975))), 
        col = transparent("darkgray",0.7), border = NA)

YYPI <- pobs$Year[pobs$Island=="PI"]
YYPI <- YYPI + c(0.05, rep(-0.05, length(YYPI) - 2), -0.1)
points(YYPI, pobs.ci[pobs$Island=="PI","PointEst"], 
       col = "cornflowerblue", pch = 15, cex = 1.5)
segments(x0 = YYPI, y0 = pobs.ci[pobs$Island=="PI","Lower"],
         y1 = pobs.ci[pobs$Island=="PI","Upper"], col = "cornflowerblue")
YYDI <- pobs$Year[pobs$Island=="DI"]
YYDI <- YYDI + c(0.1, rep(0.05, length(YYDI) - 2), -0.05)
points(YYDI, pobs.ci[pobs$Island=="DI","PointEst"], 
       col = "darkgray", pch = 16, cex = 1.5)
segments(x0 = YYDI, y0 = pobs_ci[pobs$Island=="DI","Lower"],
         y1 = pobs_ci[pobs$Island=="DI","Upper"], col = "darkgray")

legend("topleft", c("Protection","Destruction"), lwd = 3, pch = c(15,16), 
       col = c("cornflowerblue","darkgray"))