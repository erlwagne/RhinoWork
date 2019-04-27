## RHAU CODE

## Set graphics device
if(.Platform$OS.type == "windows") options(device = windows) else options(device = quartz)

library(here)
library(rstanarm)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(sm)
library(loo)
library(denstrip)
library(yarrr)
library(corrplot)
library(lubridate)
library(tidyr)
library(dplyr)
source(here::here("analysis","loo_compair.R"))

#--------------------------
# DATA  
#--------------------------

# Nest success data
nest_data <- read.csv(here::here("data","RhAu2.csv"), fileEncoding="UTF-8-BOM")
names(nest_data) <- gsub("lastcheck", "last_check", tolower(names(nest_data)))

# PDO from ERSST V3b https://www.esrl.noaa.gov/psd/pdo/ Using EOF from 1920 to 2014 for N Pacific
# Monthly 1854-2019
pdo <- read.csv(here::here("data","pdo.timeseries.ersstv3b.csv"), na.strings = "-9999.000")
colnames(pdo)[2] <- "pdo"
pdo$Date <- ymd(pdo$Date)
pdo <- mutate(pdo, year = ifelse(month(Date) > 8, year(Date) + 1, year(Date)),
              month = month(Date, label = TRUE), Date = NULL)
pdo <- spread(pdo, month, pdo)
pdo <- select(pdo, c(year, Sep:Dec, Jan:Aug))
pdo <- data.frame(pdo, pdo_index = rowMeans(select(pdo, Nov:Mar)))

# MEI v.2 1979-2018
# Wide format: rows are years, columns are months
mei <- read.csv(here::here("data","MEIv2.csv"))
colnames(mei) <- tolower(colnames(mei))
mei <- gather(mei, months, mei, decjan:novdec)
mei <- mutate(mei, year = ifelse(months %in% c("sepoct","octnov","novdec"), year + 1, year))
mei <- spread(mei, months, mei)
mei <- select(mei, c(year,sepoct,octnov,novdec,decjan,
                     janfeb,febmar,marapr,aprmay,mayjun,junjul,julaug,augsep))
mei <- data.frame(mei, mei_avg = rowMeans(select(mei, sepoct:augsep)))

# Average monthly SST from DFO stations
# Wide format: year x month
sst_amph <- read.csv(here::here("data",grep("Amphitrite", list.files(here::here("data")), 
                                            value = TRUE)),
                     skip = 1, na.strings = "99.99")
colnames(sst_amph) <- tolower(colnames(sst_amph))
sst_amph <- mutate(sst_amph, sst_amph_spring = rowMeans(select(sst_amph, apr:jun)))

sst_race <- read.csv(here::here("data",grep("Race_Rocks", list.files(here::here("data")), 
                                            value = TRUE)),
                     skip = 1, na.strings = "99.99")
colnames(sst_race) <- tolower(colnames(sst_race))
sst_race <- mutate(sst_race, sst_race_spring = rowMeans(select(sst_race, apr:jun)))

# Coastal Upwelling Index 48N 125W 1946-2018
# Wide format: rows are years, columns are months
cui <- read.csv(here::here("data","CoastalUpwellingIndex.csv"))[,-1]
colnames(cui)[1] <- "year"
cui <- mutate(cui, cui_spring = rowMeans(select(cui, Apr:Jun)))

# Biological spring transition from NWFSC Ocean Ecosystem Indicators
# Long format, 1970-2018
biol_trans <- read.csv(here::here("data","biological_spring_transition_NWFSC.csv"),
                       skip = 10, na.strings = "Never ", stringsAsFactors = FALSE)
colnames(biol_trans) <- c("year","st_onset","st_end","st_duration")
biol_trans <- mutate(biol_trans, st_onset = yday(ydm(paste(year, st_onset, sep = "-"))))
biol_trans <- mutate(biol_trans, st_onset = replace_na(st_onset, 365), 
                     st_end = replace_na(st_end, 365))

# Area-Averaged of Sea Surface Temperature at 11 microns (Day) monthly 4 km [MODIS-Aqua ()
# at Protection and Destruction Island
# Monthly 2002-2018
sst_DI <- read.csv(here::here("data","SST_DI_2002_2019.csv"), skip = 8, 
                   na.strings = "-32767")[,1:2]
colnames(sst_DI) <- c("date","sst")
sst_DI$date <- mdy_hm(sst_DI$date)
sst_DI <- mutate(sst_DI, year = year(date), month = month(date, label = TRUE), date = NULL)
sst_DI <- spread(sst_DI, month, sst)
sst_DI <- mutate(sst_DI, sst_DI_spring = rowMeans(select(sst_DI, Apr:Jun)))

sst_PI <- read.csv(here::here("data","SST_PI_2002_2019.csv"), skip = 8, 
                   na.strings = "-32767")[,1:2]
colnames(sst_PI) <- c("date","sst")
sst_PI$date <- ymd_hms(sst_PI$date)
sst_PI <- mutate(sst_PI, year = year(date), month = month(date, label = TRUE), date = NULL)
sst_PI <- spread(sst_PI, month, sst)
sst_PI <- mutate(sst_PI, sst_PI_spring = rowMeans(select(sst_PI, Apr:Jun)))

#  Area-Averaged of Chlorophyll a concentration monthly 4 km [MODIS-Aqua MODISA_L3m_CHL v2018]
# at Protection and Destruction Island
# Monthly 2002-2019
chla_DI <- read.csv(here::here("data","Chla_DI_2002_2019.csv"), skip = 8, na.strings = "-32767")[,1:2]
colnames(chla_DI) <- c("date","chla")
chla_DI$date <- mdy_hm(chla_DI$date)
chla_DI <- mutate(chla_DI, year = year(date), month = month(date, label = TRUE), date = NULL)
chla_DI <- spread(chla_DI, month, chla)
chla_DI <- mutate(chla_DI, chla_DI_spring = rowMeans(select(chla_DI, Apr:Jun)))

chla_PI <- read.csv(here::here("data","Chla_PI_2002_2019.csv"), skip = 8, na.strings = "-32767")[,1:2]
colnames(chla_PI) <- c("date","chla")
chla_PI$date <- mdy_hm(chla_PI$date)
chla_PI <- mutate(chla_PI, year = year(date), month = month(date, label = TRUE), date = NULL)
chla_PI <- spread(chla_PI, month, chla)
chla_PI <- mutate(chla_PI, chla_PI_spring = rowMeans(select(chla_PI, Apr:Jun)))

# Merge covariate data
env_data <- Reduce(inner_join, list(select(pdo, c(year, pdo_index)), 
                                    select(mei, c(year, mei_avg)), 
                                    select(sst_DI, c(year, sst_DI_spring)), 
                                    select(sst_PI, c(year, sst_PI_spring)),
                                    select(sst_amph, c(year, sst_amph_spring)),
                                    select(sst_race, c(year, sst_race_spring)),
                                    select(cui, c(year, cui_spring)),
                                    select(chla_DI, c(year, chla_DI_spring)),
                                    select(chla_PI, c(year, chla_PI_spring)),
                                    select(biol_trans, c(year, st_onset, st_duration))))

#---------------------------------
# PCA of Oceanographic Indicators
#---------------------------------

## Correlation plot of indicators
dat <- select(env_data, c(pdo_index, mei_avg, sst_DI_spring, sst_PI_spring, 
                          cui_spring, chla_DI_spring, chla_PI_spring, 
                          st_onset, st_duration))
# colnames(dat) <- c("SST spring (DI)", "SST spring (PI)", "MEI", "Coastal Upwelling (Spring)", 
#                    "Spring Transition (Onset)", "PDO Index")
dev.new()
corrplot(cor(dat, use = "pairwise"), method = "ellipse", diag = FALSE)

## Pairs plot of indicators
dev.new(width = 12, height = 12)
pairs(select(env_data, c(pdo_index, mei_avg, sst_DI_spring, sst_PI_spring, 
                         cui_spring, chla_DI_spring, chla_PI_spring, st_onset, st_duration)), 
      gap = 0.2, pch = 16, cex = 1.2, col = transparent("slategray",0.6))

# PCA
pca_env <- prcomp(~ pdo_index + mei_avg + sst_DI_spring + sst_PI_spring + 
                    cui_spring + chla_DI_spring + chla_PI_spring + st_onset, 
                  data = env_data, scale = TRUE)

pca_env            # rotation matrix gives the loadings
summary(pca_env)   # proportion of variance associated with each PC

## PCA plots
dev.new(width = 12, height = 12)
par(mfcol = c(2,2))
# scree plot
imp <- summary(pca_env)$importance["Proportion of Variance",]
barplot(imp, xlab = "", ylab = "Proportion of variance", names.arg = names(imp))   
# biplot of PCAs and oceanographic indicators
biplot(pca_env) 
# PC1 loadings
xloc <- barplot(pca_env$rotation[,"PC1"], xaxt = "n", main = "PC1 loadings")
text(xloc, par("usr")[3], labels = dimnames(pca_env$rotation)[[1]], adj = c(1,1), srt = 45, xpd = TRUE)
# PC2 loadings
xloc <- barplot(pca_env$rotation[,"PC2"], xaxt = "n", main = "PC2 loadings")
text(xloc, par("usr")[3], labels = dimnames(pca_env$rotation)[[1]], adj = c(1,1), srt = 45, xpd = TRUE)


## Add PC1 and PC2 to covariate data
scores <- predict(pca_env, newdata = env_data)
env_data <- data.frame(env_data, PC1 = scale(scores[,"PC1"]), PC2 = scale(scores[,"PC2"]))

## Merge PC1 and PC2 into nest data
rhau <- left_join(nest_data, select(env_data, c(year, PC1, PC2)))

## Time-series plots of indicators and PCs
## ggplot2
dev.new(width = 7, height = 10)
pdo.gg <- ggplot(data = env_data, aes(x = year, y = pdo_index)) + geom_line() +
  labs(x = "", y = "PDO", size = 5) + theme_gray()
mei.gg <- ggplot(data = env_data, aes(x = year, y = mei_avg)) + geom_line()+
  labs(x = "", y = "MEI") + theme_gray()
sstdi.gg <- ggplot(data = env_data, aes(x = year, y = sst_DI_spring)) + geom_line() +
  labs(x = "", y = "SST (DI)") + theme_gray()
sstpi.gg <- ggplot(data = env_data, aes(x = year, y = sst_PI_spring)) + geom_line() +
  labs(x = "", y = "SST (PI)") + theme_gray()
cui.gg <- ggplot(data = env_data, aes(x = year, y = cui_spring)) + geom_line() +
  labs(x = "", y = "Coastal Upwelling") + theme_gray()
sto.gg <- ggplot(data = env_data, aes(x = year, y = st_onset)) + geom_line() +
  labs(x = "", y = "Spring Transition (d)") + theme_gray()
chldi.gg <- ggplot(data = env_data, aes(x = year, y = chla_DI_spring)) + geom_line() +
  labs(x = "", y = "Chl a (DI)") + theme_gray()
chlpi.gg <- ggplot(data = env_data, aes(x = year, y = chla_PI_spring)) + geom_line() +
  labs(x = "", y = "Chl a (PI)") + theme_gray()
pc1.gg <- ggplot(data = env_data, aes(x = year, y = PC1)) + geom_line() +
  labs(x = "Year", y = "PC1") + theme_gray()
pc2.gg <- ggplot(data = env_data, aes(x = year, y = PC2)) + geom_line() +
  labs(x = "Year", y = "PC2") + theme_gray()
grid.arrange(pdo.gg, mei.gg, sstdi.gg, sstpi.gg, cui.gg, sto.gg, chldi.gg, chlpi.gg, pc1.gg, pc2.gg,
             nrow = 5, ncol = 2)

# base graphics
dev.new(width = 10, height = 10)
par(mfrow = c(5,2), mar = c(2,5,1,1), oma = c(3,0,0,0))
vars <- c("pdo_index","mei_avg","sst_DI_spring","sst_PI_spring","cui_spring","st_onset",
  "chla_DI_spring","chla_PI_spring","PC1","PC2")
ylabs <- c("PDO","MEI","SST (DI)", "SST (PI)", "Coastal upwelling", "Spring transition (d)",
           bquote("Chl " * italic(a) * " (DI)"), bquote("Chl " * italic(a) * " (PI)"), "PC1", "PC2")
for(i in 1:length(vars)) {
  plot(env_data$year, env_data[,vars[i]], pch = "", las = 1, cex.lab = 1.5, cex.axis = 1,
       xlab = "", ylab = ylabs[i], xpd = NA)
  mtext(ifelse(par("mfg")[1] == 5, "Year", ""), side = 1, line = 3, cex = 1.5*par("cex"))
  # grid(nx = NA, ny = NULL, col = "lightgray", lty = 1)
  abline(v = env_data$year, col = "lightgray")
  rug(env_data$year[env_data$year %% 5 != 0], ticksize = -0.04)
  lines(env_data$year, env_data[,vars[i]], type = "l", lwd = 2)
}

# # much nicer
# dev.new(width = 7, height = 10)
# par(mfrow = c(7,1), mar = c(1,5,1,1), oma = c(3,0,0,0))
# vars <- c("pdo_index","mei_avg","sst_DI_spring","sst_PI_spring","cui_spring","st_onset",
#           "chla_DI_spring","chla_PI_spring","PC1","PC2")
# cols <- c("purple","orangered","darkred","salmon","darkblue","green","darkgreen","lightgreen",
#           "black","darkgray")
# ylabs <- c("PDO","MEI","SST", "SST", "CUI", "Transition", 
#            bquote("Chl " * italic(a)), bquote("Chl " * italic(a)), "PC", "PC")
# for(i in unique(ylabs)) {
#   plot(env_data$year, env_data[,vars[grep(i,ylabs)[1]]], type = "l", lwd = 2,
#        col = cols[grep(i,ylabs)[1]], las = 1, cex.axis = 1.2, cex.lab = 1.5, xlab = "", ylab = i,
#        ylim = range(env_data[,vars[grep(i,ylabs)]], na.rm = TRUE), xaxt = "n", yaxt = "n", bty = "n")
#   axis(side = 1, at = env_data$year, cex.axis = 1.2)
#   axis(side = 2, at = axisTicks(par("usr")[3:4], log = FALSE, nint = 3), cex.axis = 1.2)
#   rug(env_data$year[env_data$year %% 5 != 0], ticksize = -0.04)
#   mtext(ifelse(par("mfg")[1] == 5, "Year", ""), side = 1, line = 3, cex = 1.5*par("cex"))
#   if(sum(ylabs == i) == 2) {
#     lines(env_data$year, env_data[,vars[grep(i,ylabs)[2]]], lwd = 2, col = cols[grep(i,ylabs)[2]])
#     legend("topright", c("DI","PI"), lwd = 2, col = cols[grep(i,ylabs)])
#   }
# }




#---------------------------------
# GLMMs
#---------------------------------

## Hierarchical Bayesian binomial models with rstanarm ##
## Burrow Occupancy ##

## Random effects of site and year

# Intercept-only model 
occ0 <- stan_glmer(cbind(egg, viable - egg) ~ (1 | site) + (1 | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(occ0, prob = c(0.025,0.5,0.975), digits = 2)
launch_shinystan(occ0)

# Inter-island differences, constant across years
occ1 <- stan_glmer(cbind(egg, viable - egg) ~ island + (1 | site) + (1 | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(occ1, prob = c(0.025,0.5,0.975), digits = 2)
launch_shinystan(occ1)

# Inter-island differences, varying among years
occ2 <- stan_glmer(cbind(egg, viable - egg) ~ island + (1 | site) + (island | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(occ2, prob = c(0.025,0.5,0.975), digits = 2)
launch_shinystan(occ2)

# Inter-island differences, varying among years, plus PC1 + PC2
occ3 <- stan_glmer(cbind(egg, viable - egg) ~ island + PC1 + PC2 + (1 | site) + (island | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(occ3, prob = c(0.025,0.5,0.975), digits = 2)
launch_shinystan(occ3)

# Inter-island differences plus PC1 + PC2, no random time-variation
occ4 <- stan_glmer(cbind(egg, viable - egg) ~ island + PC1 + PC2 + (1 | site),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(occ4, prob = c(0.025,0.5,0.975), digits = 2)
launch_shinystan(occ4)

## Model selection by approximate leave-one-out cross-validation
occ_loos <-  lapply(list(occ0 = occ0, occ1 = occ1, occ2 = occ2, occ3 = occ3, occ4 = occ4), loo)
occ_loos

# unpaired comparisons
occ_compare <- loo_compare(occ_loos)
print(occ_compare, simplify = FALSE)

# pairwise comparisons
occ_compair <- loo_compair(occ_loos)
occ_compair

## Cache stanfits and loo objects
save(list = grep("occ", ls(), value = TRUE), file = here::here("analysis","cache","occ_stanfit_loo.RData"))


#-------------------------
# BURROW OCCUPANCY FIGURES
#-------------------------

## Model validation by graphical posterior predictive checking
# default plot: marginal distribution of data and posterior predictive distribution
pp_check(occ3)

## Posterior distribution of average between-island difference in full model,
## with median and 95% credible interval
dev.new()
par(mar = c(5.1,4.3,4.1,1))

x <- as.matrix(occ2, pars = "islandPI")
px <- sm.density(x, display = "none")
px <- sm.density(x, eval.points = sort(c(px$eval.points, quantile(x, c(0.025,0.5,0.975)))), display = "none")
ci <- which(names(px$eval.points)=="2.5%"):which(names(px$eval.points)=="97.5%")
plot(px$eval.points, px$estimate, col = "darkgray", type = "l", lwd = 3,
     las = 1, cex.lab = 1.5, cex.axis = 1.2, bty = "n", 
     ylim = c(0,max(px$estimate)), yaxs = "i",
     xlab = "Occupancy log-odds ratio: Protection vs. Destruction", 
     ylab = "Probability density",
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

YY <- sort(unique(rhau$year))
newdata <- data.frame(year = rep(YY, each = 2), site = "newsite",
                      island = rep(c("DI","PI"), length(YY)))
pfit <- posterior_linpred(occ2, transform = TRUE, re.form = ~ (1 | year), newdata = newdata)
pobs <- aggregate(cbind(viable, egg) ~ year + island, data = rhau, sum)
pobs_ci <- binconf(pobs$egg, n = pobs$viable)
eval.points <- range(pobs_ci, apply(pfit,2,quantile,c(0.025,0.975)))
eval.points <- seq(min(eval.points), max(eval.points), length = 300)
xyDI <- pfit[,newdata$island=="DI"]
pxyDI <- t(apply(xyDI, 2, function(x) 
  sm.density(x, eval.points = eval.points, display = "none")$estimate))
xyPI <- pfit[,newdata$island=="PI"]
pxPI <- sm.density(as.vector(xyPI), display = "none")
pxyPI <- t(apply(xyPI, 2, function(x) 
  sm.density(x, eval.points = eval.points, display = "none")$estimate))

plot(newdata$year[newdata$island=="PI"], apply(pfit[,newdata$island=="PI"],2,median), 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, type = "l", lwd = 3, col = "cornflowerblue",
     xlab = "Year", ylab = "Probability of Burrow Occupancy", 
     ylim = range(pobs_ci, apply(pfit,2,quantile,c(0.025,0.975))), 
     xaxs = "i", xaxt = "n", yaxs = "i")
lines(newdata$year[newdata$island=="DI"], apply(pfit[,newdata$island=="DI"],2,median),
      lwd = 3, col = "darkgray")
axis(1, at = min(rhau$year):max(rhau$year))

#densregion(newdata$year[newdata$ysland=="PI"], eval.points, pxyPI, 
#          colmax = transparent("darkgray",0.1), colmin = "transparent")
polygon(c(newdata$year[newdata$island=="PI"], rev(newdata$year[newdata$island=="PI"])),
        c(apply(pfit[,newdata$island=="PI"],2,quantile,0.025), 
          rev(apply(pfit[,newdata$island=="PI"],2,quantile,0.975))), 
        col = transparent("cornflowerblue",0.7), border = NA)
#densregion(newdata$Year[newdata$Island=="DI"], eval.points, pxyDI, 
#           colmax = transparent("black",0.3), colmin = "transparent")
polygon(c(newdata$year[newdata$island=="DI"], rev(newdata$year[newdata$island=="DI"])),
        c(apply(pfit[,newdata$island=="DI"],2,quantile,0.025), 
          rev(apply(pfit[,newdata$island=="DI"],2,quantile,0.975))), 
        col = transparent("darkgray",0.7), border = NA)

YYPI <- pobs$year[pobs$island=="PI"]
YYPI <- YYPI + c(0.05, rep(-0.05, length(YYPI) - 2), -0.1)
points(YYPI, pobs_ci[pobs$island=="PI","PointEst"], 
       col = "cornflowerblue", pch = 15, cex = 1.5)
segments(x0 = YYPI, y0 = pobs_ci[pobs$island=="PI","Lower"],
         y1 = pobs_ci[pobs$island=="PI","Upper"], col = "cornflowerblue")
YYDI <- pobs$year[pobs$island=="DI"]
YYDI <- YYDI + c(0.1, rep(0.05, length(YYDI) - 2), -0.05)
points(YYDI, pobs_ci[pobs$island=="DI","PointEst"], 
       col = "darkgray", pch = 16, cex = 1.5)
segments(x0 = YYDI, y0 = pobs_ci[pobs$island=="DI","Lower"],
         y1 = pobs_ci[pobs$island=="DI","Upper"], col = "darkgray")

legend("topleft", c("Protection","Destruction"), lwd = 3, pch = c(15,16), 
       col = c("cornflowerblue","darkgray"))

#----------------------------
## REPRODUCTIVE SUCCESS ##
#----------------------------

## Random effects of site and year

# Intercept-only model 
suc0 <- stan_glmer(cbind(last_check, egg - last_check) ~ (1 | site) + (1 | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(suc0, prob = c(0.025,0.5,0.975), digits = 2)
launch_shinystan(suc0)

# Inter-island differences, constant across years
suc1 <- stan_glmer(cbind(last_check, egg - last_check) ~ island + (1 | site) + (1 | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(suc1, prob = c(0.025,0.5,0.975), digits = 2)
launch_shinystan(suc1)

# Inter-island differences, varying among years
suc2 <- stan_glmer(cbind(last_check, egg - last_check) ~ island + (1 | site) + (island | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(suc2, prob = c(0.025,0.5,0.975), digits = 2)
launch_shinystan(suc2)

# Inter-island differences, varying among years, plus PC1 + PC2
suc3 <- stan_glmer(cbind(last_check, egg - last_check) ~ island + PC1 + PC2 + (1 | site) + (island | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(suc3, prob = c(0.025,0.5,0.975), digits = 2)
launch_shinystan(suc3)

# Inter-island differences plus PC1 + PC2, no random time-variation
suc4 <- stan_glmer(cbind(last_check, egg - last_check) ~ island + PC1 + PC2 + (1 | site),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)

summary(suc4, prob = c(0.025,0.5,0.975), digits = 2)
launch_shinystan(suc4)

## Model selection by approximate leave-one-out cross-validation
suc_loos <- lapply(list(suc0 = suc0, suc1 = suc1, suc2 = suc2, suc3 = suc3, suc4 = suc4), 
                    loo, k_threshold = 0.7)
suc_loos

# unpaired comparisons
suc_compare <- loo_compare(suc_loos)
print(suc_compare, simplify = FALSE)

# pairwise comparisons
suc_compair <- loo_compair(suc_loos)
suc_compair

## Cache stanfits and loo objects
save(list = grep("suc", ls(), value = TRUE), file = here::here("analysis","cache","suc_stanfit_loo.RData"))

#-----------------------------
# REPRODUCTIVE SUCCESS FIGURES
#-----------------------------

## Model validation by graphical posterior predictive checking
# default plot: marginal distribution of data and posterior predictive distribution
pp_check(suc2)

## Posterior distribution of average between-island difference in full model,
## with median and 95% credible interval
dev.new()
par(mar = c(5.1,4.3,4.1,1))

x <- as.matrix(suc2, pars = "islandPI")
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

YY <- sort(unique(rhau$year))
newdata <- data.frame(year = rep(YY, each = 2), site = "newsite",
                      island = rep(c("DI","PI"), length(YY)))
pfit <- posterior_linpred(suc2, transform = TRUE, re.form = ~ (1 | year), newdata = newdata)
pobs <- aggregate(cbind(egg, last_check) ~ year + island, data = rhau, sum)
pobs_ci <- binconf(pobs$last_check, n = pobs$egg)
eval.points <- range(pobs_ci, apply(pfit,2,quantile,c(0.025,0.975)))
eval.points <- seq(min(eval.points), max(eval.points), length = 300)
xyDI <- pfit[,newdata$island=="DI"]
pxyDI <- t(apply(xyDI, 2, function(x) 
  sm.density(x, eval.points = eval.points, display = "none")$estimate))
xyPI <- pfit[,newdata$island=="PI"]
pxPI <- sm.density(as.vector(xyPI), display = "none")
pxyPI <- t(apply(xyPI, 2, function(x) 
  sm.density(x, eval.points = eval.points, display = "none")$estimate))

plot(newdata$year[newdata$island=="PI"], apply(pfit[,newdata$island=="PI"],2,median), 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, type = "l", lwd = 3, col = "cornflowerblue",
     xlab = "Year", ylab = "Probability of fledging", 
     ylim = range(pobs_ci, apply(pfit,2,quantile,c(0.025,0.975))), 
     xaxs = "i", xaxt = "n", yaxs = "i")
lines(newdata$year[newdata$island=="DI"], apply(pfit[,newdata$island=="DI"],2,median),
      lwd = 3, col = "darkgray")
axis(1, at = min(rhau$year):max(rhau$year))

#densregion(newdata$year[newdata$island=="PI"], eval.points, pxyPI, 
#          colmax = transparent("darkgray",0.1), colmin = "transparent")
polygon(c(newdata$year[newdata$island=="PI"], rev(newdata$year[newdata$island=="PI"])),
        c(apply(pfit[,newdata$island=="PI"],2,quantile,0.025), 
          rev(apply(pfit[,newdata$island=="PI"],2,quantile,0.975))), 
        col = transparent("cornflowerblue",0.7), border = NA)
#densregion(newdata$year[newdata$island=="DI"], eval.points, pxyDI, 
#           colmax = transparent("black",0.3), colmin = "transparent")
polygon(c(newdata$year[newdata$island=="DI"], rev(newdata$year[newdata$island=="DI"])),
        c(apply(pfit[,newdata$island=="DI"],2,quantile,0.025), 
          rev(apply(pfit[,newdata$island=="DI"],2,quantile,0.975))), 
        col = transparent("darkgray",0.7), border = NA)

YYPI <- pobs$year[pobs$island=="PI"]
YYPI <- YYPI + c(0.05, rep(-0.05, length(YYPI) - 2), -0.1)
points(YYPI, pobs_ci[pobs$island=="PI","PointEst"], 
       col = "cornflowerblue", pch = 15, cex = 1.5)
segments(x0 = YYPI, y0 = pobs_ci[pobs$island=="PI","Lower"],
         y1 = pobs_ci[pobs$island=="PI","Upper"], col = "cornflowerblue")
YYDI <- pobs$year[pobs$island=="DI"]
YYDI <- YYDI + c(0.1, rep(0.05, length(YYDI) - 2), -0.05)
points(YYDI, pobs_ci[pobs$island=="DI","PointEst"], 
       col = "darkgray", pch = 16, cex = 1.5)
segments(x0 = YYDI, y0 = pobs_ci[pobs$island=="DI","Lower"],
         y1 = pobs_ci[pobs$island=="DI","Upper"], col = "darkgray")

legend("topleft", c("Protection","Destruction"), lwd = 3, pch = c(15,16), 
       col = c("cornflowerblue","darkgray"))

