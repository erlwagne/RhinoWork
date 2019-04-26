Nest Success of Two Rhinoceros Auklet Colonies in the Salish Sea and
California Current
================
Eric Wagner, Eric Buhle, …
April 21, 2019

  - [Overview](#overview)
  - [Setup and data](#setup-and-data)
  - [Principal Components Analysis of Oceanographic
    Indicators](#principal-components-analysis-of-oceanographic-indicators)
  - [GLMMs of Burrow Occupancy](#glmms-of-burrow-occupancy)

# Overview

This is a Rhinoceros Auklet from the Protection Island
colony:

![](https://www.eopugetsound.org/sites/default/files/styles/magazinewidth_592px/public/topical_articles/images/16343412270_5cfaa5c480_o.jpg?itok=WjzI1_2K)

# Setup and data

First we load the libraries and functions we’ll use.

``` r
if(!require(here)) install.packages("here::here")
library(here)
if(!require(rstanarm)) install.package("rstanarm")
library(rstanarm)
if(!require(bayesplot)) install.package("bayesplot")
library(bayesplot)
if(!require(Hmisc)) install.package("Hmisc")
library(Hmisc)
if(!require(sm)) install.package("sm")
library(sm)
if(!require(loo)) install.package("loo")
library(loo)
if(!require(denstrip)) install.package("denstrip")
library(denstrip)
if(!require(yarrr)) install.package("yarrr")
library(yarrr)
if(!require(corrplot)) install.package("corrplot")
library(corrplot)
if(!require(lubridate)) install.package("lubridate")
library(lubridate)
if(!require(tidyr)) install.package("tidyr")
library(tidyr)
if(!require(dplyr)) install.package("dplyr")
library(dplyr)
if(!require(gridExtra)) install.package("gridExtra")
library(gridExtra)
source(here::here("analysis","loo_compair.R"))
```

Next we’ll read in the datasets from various sources and manipulate them
into a usable format.

``` r
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
# Monthly 2009-2018
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
```

# Principal Components Analysis of Oceanographic Indicators

Let’s explore the patterns of (a)synchrony among the oceanographic
indicators to see how severe the multicollinearity might be if they were
used as raw regression
inputs.

``` r
dat <- select(env_data, c(pdo_index, mei_avg, sst_DI_spring, sst_PI_spring, 
                          cui_spring, chla_DI_spring, chla_PI_spring, 
                          st_onset, st_duration))
corrplot(cor(dat, use = "pairwise"), method = "ellipse", diag = FALSE)
```

![](RhAu_BayesGLMM_files/figure-gfm/env_corrplot-1.png)<!-- -->

Now we perform a principal components analysis to extract the major
trends in the suite of oceanographic indicators for use as regression
inputs.

``` r
pca_env <- prcomp(~ pdo_index + mei_avg + sst_DI_spring + sst_PI_spring + 
                    cui_spring + chla_DI_spring + chla_PI_spring + st_onset, 
                  data = env_data, scale = TRUE)

pca_env            # rotation matrix gives the loadings
```

    Standard deviations (1, .., p=8):
    [1] 2.0728549 1.3624867 0.9924571 0.6192212 0.4579104 0.3927643 0.2679374 0.2067856
    
    Rotation (n x k) = (8 x 8):
                          PC1        PC2         PC3         PC4         PC5        PC6         PC7         PC8
    pdo_index      -0.4032537 -0.1131162  0.37581649 -0.31536032 -0.60799403  0.3409443 -0.15072067  0.27438037
    mei_avg        -0.4465664 -0.1275237  0.07999949 -0.22583164  0.28890643 -0.5977609  0.27744056  0.45733801
    sst_DI_spring  -0.2374015 -0.4532955 -0.56717794  0.22677323 -0.24954683 -0.2470813 -0.49144976 -0.01418353
    sst_PI_spring  -0.4277657 -0.2711304 -0.16615976  0.07713533 -0.08299182  0.1986218  0.64087140 -0.50288314
    cui_spring     -0.2704591  0.4816178  0.29890905  0.56249483 -0.31357726 -0.3900527 -0.03809154 -0.19294126
    chla_DI_spring -0.2859470  0.4816069 -0.25504438 -0.60163238  0.07664705 -0.1002367 -0.26588956 -0.41552439
    chla_PI_spring -0.2834287  0.4465721 -0.47089350  0.25744514  0.10822038  0.4240045  0.10530644  0.47920966
    st_onset       -0.4041983 -0.1691934  0.35626848  0.21801839  0.60149643  0.2897772 -0.40599671 -0.15144634

``` r
summary(pca_env)   # proportion of variance associated with each PC
```

    Importance of components:
                              PC1    PC2    PC3     PC4     PC5     PC6     PC7     PC8
    Standard deviation     2.0729 1.3625 0.9925 0.61922 0.45791 0.39276 0.26794 0.20679
    Proportion of Variance 0.5371 0.2321 0.1231 0.04793 0.02621 0.01928 0.00897 0.00535
    Cumulative Proportion  0.5371 0.7691 0.8923 0.94019 0.96640 0.98568 0.99465 1.00000

``` r
## Add PC1 and PC2 to covariate data
scores <- predict(pca_env, newdata = env_data)
env_data <- data.frame(env_data, PC1 = scale(scores[,"PC1"]), PC2 = scale(scores[,"PC2"]))

## Merge PC1 and PC2 into nest data
rhau <- left_join(nest_data, select(env_data, c(year, PC1, PC2)))
```

    Joining, by = "year"

``` r
par(mfcol = c(2,2))
# scree plot
imp <- summary(pca_env)$importance["Proportion of Variance",]
barplot(imp, xlab = "", ylab = "Proportion of variance", names.arg = names(imp))   
# biplot of PCAs and oceanographic indicators
biplot(pca_env) 
# PC1 loadings
xloc <- barplot(pca_env$rotation[,"PC1"], xaxt = "n", main = "PC1 loadings")
text(xloc, par("usr")[3], labels = dimnames(pca_env$rotation)[[1]], adj = c(1,1), srt = 45, 
     xpd = TRUE)
# PC2 loadings
xloc <- barplot(pca_env$rotation[,"PC2"], xaxt = "n", main = "PC2 loadings")
text(xloc, par("usr")[3], labels = dimnames(pca_env$rotation)[[1]], adj = c(1,1), srt = 45, 
     xpd = TRUE)
```

![](RhAu_BayesGLMM_files/figure-gfm/PCA_plots-1.png)<!-- -->

``` r
## Time-series plots of indicators and PCs
pdo.gg <- ggplot(data = env_data, aes(x = year, y = pdo_index)) + geom_line() +
  labs(x = "Year", y = "PDO", size = 5) + theme_gray()
mei.gg <- ggplot(data = env_data, aes(x = year, y = mei_avg)) + geom_line()+
  labs(x = "Year", y = "MEI") + theme_gray()
sstdi.gg <- ggplot(data = env_data, aes(x = year, y = sst_DI_spring)) + geom_line() +
  labs(x = "Year", y = "SST (DI)") + theme_gray()
sstpi.gg <- ggplot(data = env_data, aes(x = year, y = sst_PI_spring)) + geom_line() +
  labs(x = "Year", y = "SST (PI)") + theme_gray()
cui.gg <- ggplot(data = env_data, aes(x = year, y = cui_spring)) + geom_line() +
  labs(x = "Year", y = "Coastal Upwelling") + theme_gray()
sto.gg <- ggplot(data = env_data, aes(x = year, y = st_onset)) + geom_line() +
  labs(x = "Year", y = "Spring Transition (d)") + theme_gray()
chldi.gg <- ggplot(data = env_data, aes(x = year, y = chla_DI_spring)) + geom_line() +
  labs(x = "Year", y = "Chl a (DI)") + theme_gray()
chlpi.gg <- ggplot(data = env_data, aes(x = year, y = chla_PI_spring)) + geom_line() +
  labs(x = "Year", y = "Chl a (PI)") + theme_gray()
pc1.gg <- ggplot(data = env_data, aes(x = year, y = PC1)) + geom_line() +
  labs(x = "Year", y = "PC1") + theme_gray()
pc2.gg <- ggplot(data = env_data, aes(x = year, y = PC2)) + geom_line() +
  labs(x = "Year", y = "PC2") + theme_gray()
grid.arrange(pdo.gg, mei.gg, sstdi.gg, sstpi.gg, cui.gg, sto.gg, chldi.gg, chlpi.gg, pc1.gg, pc2.gg,
             nrow = 5, ncol = 2)
```

![](RhAu_BayesGLMM_files/figure-gfm/env_timeseries_plots-1.png)<!-- -->

# GLMMs of Burrow Occupancy

Now we are ready to start fitting some GLMMs to the burrow occupancy
data…\[more on modeling, explanations of each model\]

``` r
## Random effects of site and year

# Intercept-only model 
occ0 <- stan_glmer(cbind(egg, viable - egg) ~ (1 | site) + (1 | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)
```

``` r
summary(occ0, prob = c(0.025,0.5,0.975), digits = 2)
```

``` 

Model Info:

 function:     stan_glmer
 family:       binomial [logit]
 formula:      cbind(egg, viable - egg) ~ (1 | site) + (1 | year)
 algorithm:    sampling
 priors:       see help('prior_summary')
 sample:       3000 (posterior sample size)
 observations: 139
 groups:       site (17), year (9)

Estimates:
                                      mean    sd      2.5%    50%     97.5%
(Intercept)                            0.54    0.11    0.33    0.54    0.75
b[(Intercept) site:Boathouse]         -0.26    0.16   -0.58   -0.26    0.04
b[(Intercept) site:Catwalk.Salmon]     0.05    0.18   -0.31    0.04    0.39
b[(Intercept) site:Catwalk.Willow]    -0.30    0.18   -0.68   -0.29    0.04
b[(Intercept) site:Grassy.Knoll]      -0.10    0.20   -0.51   -0.09    0.28
b[(Intercept) site:s741]              -0.18    0.18   -0.55   -0.18    0.16
b[(Intercept) site:s742]               0.20    0.18   -0.12    0.20    0.58
b[(Intercept) site:s743]               0.22    0.19   -0.13    0.22    0.63
b[(Intercept) site:s744]              -0.02    0.17   -0.38   -0.02    0.31
b[(Intercept) site:s746]              -0.20    0.19   -0.58   -0.19    0.15
b[(Intercept) site:s747]               0.23    0.19   -0.12    0.22    0.63
b[(Intercept) site:s748]              -0.16    0.21   -0.58   -0.15    0.25
b[(Intercept) site:s749]               0.25    0.18   -0.09    0.25    0.62
b[(Intercept) site:s750]               0.32    0.19   -0.03    0.31    0.74
b[(Intercept) site:s751]               0.12    0.21   -0.27    0.11    0.55
b[(Intercept) site:Salmon.Sequel]     -0.16    0.19   -0.52   -0.16    0.20
b[(Intercept) site:SW.Point]          -0.24    0.16   -0.58   -0.24    0.06
b[(Intercept) site:Willow.Gully]       0.18    0.17   -0.13    0.18    0.55
b[(Intercept) year:2010]               0.10    0.12   -0.10    0.08    0.38
b[(Intercept) year:2011]               0.04    0.12   -0.18    0.03    0.32
b[(Intercept) year:2012]              -0.02    0.11   -0.25   -0.01    0.19
b[(Intercept) year:2013]              -0.04    0.12   -0.30   -0.02    0.19
b[(Intercept) year:2014]              -0.09    0.12   -0.37   -0.06    0.13
b[(Intercept) year:2015]              -0.03    0.11   -0.28   -0.02    0.18
b[(Intercept) year:2016]               0.16    0.14   -0.04    0.14    0.47
b[(Intercept) year:2017]              -0.09    0.12   -0.36   -0.07    0.10
b[(Intercept) year:2018]              -0.05    0.11   -0.30   -0.03    0.15
Sigma[site:(Intercept),(Intercept)]    0.09    0.05    0.02    0.08    0.22
Sigma[year:(Intercept),(Intercept)]    0.03    0.04    0.00    0.02    0.13
mean_PPD                               7.39    0.19    7.00    7.40    7.76
log-posterior                       -312.31    5.47 -323.92 -311.90 -302.88

Diagnostics:
                                    mcse Rhat n_eff
(Intercept)                         0.00 1.00 1860 
b[(Intercept) site:Boathouse]       0.00 1.00 3105 
b[(Intercept) site:Catwalk.Salmon]  0.00 1.00 4383 
b[(Intercept) site:Catwalk.Willow]  0.00 1.00 3144 
b[(Intercept) site:Grassy.Knoll]    0.00 1.00 4893 
b[(Intercept) site:s741]            0.00 1.00 4308 
b[(Intercept) site:s742]            0.00 1.00 3873 
b[(Intercept) site:s743]            0.00 1.00 4759 
b[(Intercept) site:s744]            0.00 1.00 4092 
b[(Intercept) site:s746]            0.00 1.00 3859 
b[(Intercept) site:s747]            0.00 1.00 4244 
b[(Intercept) site:s748]            0.00 1.00 4871 
b[(Intercept) site:s749]            0.00 1.00 4322 
b[(Intercept) site:s750]            0.00 1.00 3287 
b[(Intercept) site:s751]            0.00 1.00 5190 
b[(Intercept) site:Salmon.Sequel]   0.00 1.00 4436 
b[(Intercept) site:SW.Point]        0.00 1.00 3267 
b[(Intercept) site:Willow.Gully]    0.00 1.00 4089 
b[(Intercept) year:2010]            0.00 1.00 2019 
b[(Intercept) year:2011]            0.00 1.00 3397 
b[(Intercept) year:2012]            0.00 1.00 4451 
b[(Intercept) year:2013]            0.00 1.00 3837 
b[(Intercept) year:2014]            0.00 1.00 2747 
b[(Intercept) year:2015]            0.00 1.00 3806 
b[(Intercept) year:2016]            0.00 1.00 1598 
b[(Intercept) year:2017]            0.00 1.00 2620 
b[(Intercept) year:2018]            0.00 1.00 2764 
Sigma[site:(Intercept),(Intercept)] 0.00 1.00 1459 
Sigma[year:(Intercept),(Intercept)] 0.00 1.00 1300 
mean_PPD                            0.00 1.00 2935 
log-posterior                       0.20 1.00  777 

For each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).
```

``` r
# Inter-island differences, constant across years
occ1 <- stan_glmer(cbind(egg, viable - egg) ~ island + (1 | site) + (1 | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)
```

``` r
summary(occ1, prob = c(0.025,0.5,0.975), digits = 2)
```

``` 

Model Info:

 function:     stan_glmer
 family:       binomial [logit]
 formula:      cbind(egg, viable - egg) ~ island + (1 | site) + (1 | year)
 algorithm:    sampling
 priors:       see help('prior_summary')
 sample:       3000 (posterior sample size)
 observations: 139
 groups:       site (17), year (9)

Estimates:
                                      mean    sd      2.5%    50%     97.5%
(Intercept)                            0.35    0.13    0.08    0.35    0.62
islandPI                               0.33    0.16    0.01    0.33    0.63
b[(Intercept) site:Boathouse]         -0.10    0.16   -0.44   -0.10    0.19
b[(Intercept) site:Catwalk.Salmon]     0.13    0.18   -0.21    0.11    0.53
b[(Intercept) site:Catwalk.Willow]    -0.15    0.17   -0.53   -0.13    0.16
b[(Intercept) site:Grassy.Knoll]       0.00    0.18   -0.35    0.00    0.36
b[(Intercept) site:s741]              -0.22    0.18   -0.60   -0.21    0.08
b[(Intercept) site:s742]               0.10    0.17   -0.20    0.08    0.46
b[(Intercept) site:s743]               0.12    0.17   -0.21    0.11    0.48
b[(Intercept) site:s744]              -0.09    0.16   -0.42   -0.08    0.22
b[(Intercept) site:s746]              -0.23    0.19   -0.64   -0.21    0.09
b[(Intercept) site:s747]               0.12    0.18   -0.20    0.11    0.50
b[(Intercept) site:s748]              -0.17    0.20   -0.61   -0.15    0.17
b[(Intercept) site:s749]               0.14    0.17   -0.16    0.13    0.51
b[(Intercept) site:s750]               0.20    0.19   -0.12    0.18    0.63
b[(Intercept) site:s751]               0.05    0.19   -0.32    0.04    0.44
b[(Intercept) site:Salmon.Sequel]     -0.04    0.17   -0.38   -0.03    0.30
b[(Intercept) site:SW.Point]          -0.09    0.16   -0.44   -0.09    0.22
b[(Intercept) site:Willow.Gully]       0.25    0.19   -0.06    0.23    0.68
b[(Intercept) year:2010]               0.10    0.12   -0.10    0.08    0.39
b[(Intercept) year:2011]               0.04    0.12   -0.20    0.03    0.31
b[(Intercept) year:2012]              -0.02    0.12   -0.29   -0.01    0.20
b[(Intercept) year:2013]              -0.04    0.12   -0.31   -0.03    0.17
b[(Intercept) year:2014]              -0.09    0.12   -0.37   -0.07    0.11
b[(Intercept) year:2015]              -0.03    0.11   -0.28   -0.02    0.19
b[(Intercept) year:2016]               0.17    0.15   -0.04    0.15    0.50
b[(Intercept) year:2017]              -0.09    0.12   -0.36   -0.07    0.12
b[(Intercept) year:2018]              -0.05    0.11   -0.30   -0.04    0.16
Sigma[site:(Intercept),(Intercept)]    0.06    0.05    0.00    0.05    0.19
Sigma[year:(Intercept),(Intercept)]    0.03    0.04    0.00    0.02    0.13
mean_PPD                               7.38    0.20    6.99    7.39    7.76
log-posterior                       -313.63    5.79 -326.20 -313.16 -303.44

Diagnostics:
                                    mcse Rhat n_eff
(Intercept)                         0.00 1.00 1582 
islandPI                            0.00 1.00 1793 
b[(Intercept) site:Boathouse]       0.00 1.00 2059 
b[(Intercept) site:Catwalk.Salmon]  0.00 1.00 3061 
b[(Intercept) site:Catwalk.Willow]  0.00 1.00 2391 
b[(Intercept) site:Grassy.Knoll]    0.00 1.00 3063 
b[(Intercept) site:s741]            0.00 1.00 1881 
b[(Intercept) site:s742]            0.00 1.00 3179 
b[(Intercept) site:s743]            0.00 1.00 2738 
b[(Intercept) site:s744]            0.00 1.00 2606 
b[(Intercept) site:s746]            0.00 1.00 1931 
b[(Intercept) site:s747]            0.00 1.00 2441 
b[(Intercept) site:s748]            0.00 1.00 2934 
b[(Intercept) site:s749]            0.00 1.00 3044 
b[(Intercept) site:s750]            0.00 1.00 2248 
b[(Intercept) site:s751]            0.00 1.00 3163 
b[(Intercept) site:Salmon.Sequel]   0.00 1.00 3054 
b[(Intercept) site:SW.Point]        0.00 1.00 2277 
b[(Intercept) site:Willow.Gully]    0.00 1.00 2017 
b[(Intercept) year:2010]            0.00 1.00 2230 
b[(Intercept) year:2011]            0.00 1.00 3426 
b[(Intercept) year:2012]            0.00 1.00 2905 
b[(Intercept) year:2013]            0.00 1.00 2193 
b[(Intercept) year:2014]            0.00 1.00 2299 
b[(Intercept) year:2015]            0.00 1.00 2709 
b[(Intercept) year:2016]            0.00 1.00 1642 
b[(Intercept) year:2017]            0.00 1.00 2058 
b[(Intercept) year:2018]            0.00 1.00 2567 
Sigma[site:(Intercept),(Intercept)] 0.00 1.00 1125 
Sigma[year:(Intercept),(Intercept)] 0.00 1.00 1088 
mean_PPD                            0.00 1.00 3294 
log-posterior                       0.23 1.00  621 

For each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).
```

``` r
# Inter-island differences, varying among years
occ2 <- stan_glmer(cbind(egg, viable - egg) ~ island + (1 | site) + (island | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)
```

``` r
summary(occ2, prob = c(0.025,0.5,0.975), digits = 2)
```

``` 

Model Info:

 function:     stan_glmer
 family:       binomial [logit]
 formula:      cbind(egg, viable - egg) ~ island + (1 | site) + (island | year)
 algorithm:    sampling
 priors:       see help('prior_summary')
 sample:       3000 (posterior sample size)
 observations: 139
 groups:       site (17), year (9)

Estimates:
                                      mean    sd      2.5%    50%     97.5%
(Intercept)                            0.34    0.14    0.07    0.34    0.62
islandPI                               0.34    0.18   -0.02    0.34    0.68
b[(Intercept) site:Boathouse]         -0.10    0.15   -0.43   -0.09    0.19
b[(Intercept) site:Catwalk.Salmon]     0.12    0.17   -0.18    0.11    0.49
b[(Intercept) site:Catwalk.Willow]    -0.14    0.17   -0.52   -0.12    0.18
b[(Intercept) site:Grassy.Knoll]       0.00    0.18   -0.35    0.00    0.36
b[(Intercept) site:s741]              -0.22    0.18   -0.61   -0.21    0.09
b[(Intercept) site:s742]               0.10    0.17   -0.21    0.08    0.47
b[(Intercept) site:s743]               0.11    0.17   -0.22    0.10    0.49
b[(Intercept) site:s744]              -0.08    0.16   -0.41   -0.07    0.23
b[(Intercept) site:s746]              -0.22    0.19   -0.62   -0.20    0.08
b[(Intercept) site:s747]               0.12    0.18   -0.20    0.10    0.51
b[(Intercept) site:s748]              -0.16    0.19   -0.60   -0.14    0.18
b[(Intercept) site:s749]               0.13    0.17   -0.18    0.11    0.49
b[(Intercept) site:s750]               0.19    0.18   -0.11    0.17    0.60
b[(Intercept) site:s751]               0.04    0.19   -0.32    0.03    0.43
b[(Intercept) site:Salmon.Sequel]     -0.04    0.17   -0.42   -0.03    0.29
b[(Intercept) site:SW.Point]          -0.09    0.16   -0.44   -0.08    0.19
b[(Intercept) site:Willow.Gully]       0.25    0.18   -0.06    0.23    0.64
b[(Intercept) year:2010]               0.07    0.13   -0.15    0.05    0.38
b[islandPI year:2010]                  0.06    0.18   -0.31    0.04    0.45
b[(Intercept) year:2011]               0.04    0.13   -0.21    0.02    0.33
b[islandPI year:2011]                  0.01    0.18   -0.37    0.00    0.40
b[(Intercept) year:2012]              -0.03    0.12   -0.30   -0.01    0.20
b[islandPI year:2012]                  0.03    0.17   -0.31    0.01    0.43
b[(Intercept) year:2013]               0.01    0.13   -0.22    0.00    0.32
b[islandPI year:2013]                 -0.14    0.20   -0.64   -0.09    0.15
b[(Intercept) year:2014]              -0.10    0.14   -0.45   -0.07    0.10
b[islandPI year:2014]                  0.09    0.19   -0.23    0.05    0.55
b[(Intercept) year:2015]              -0.06    0.13   -0.37   -0.04    0.15
b[islandPI year:2015]                  0.13    0.20   -0.16    0.08    0.62
b[(Intercept) year:2016]               0.12    0.14   -0.10    0.09    0.45
b[islandPI year:2016]                  0.13    0.20   -0.21    0.09    0.59
b[(Intercept) year:2017]              -0.02    0.13   -0.29   -0.02    0.25
b[islandPI year:2017]                 -0.17    0.20   -0.67   -0.13    0.11
b[(Intercept) year:2018]               0.00    0.12   -0.24    0.00    0.28
b[islandPI year:2018]                 -0.14    0.19   -0.61   -0.10    0.13
Sigma[site:(Intercept),(Intercept)]    0.06    0.05    0.00    0.05    0.18
Sigma[year:(Intercept),(Intercept)]    0.03    0.04    0.00    0.02    0.14
Sigma[year:islandPI,(Intercept)]      -0.01    0.03   -0.10    0.00    0.03
Sigma[year:islandPI,islandPI]          0.06    0.08    0.00    0.03    0.28
mean_PPD                               7.37    0.20    6.96    7.37    7.77
log-posterior                       -330.96    6.55 -344.63 -330.64 -318.95

Diagnostics:
                                    mcse Rhat n_eff
(Intercept)                         0.00 1.00 1844 
islandPI                            0.00 1.00 2148 
b[(Intercept) site:Boathouse]       0.00 1.00 2661 
b[(Intercept) site:Catwalk.Salmon]  0.00 1.00 2870 
b[(Intercept) site:Catwalk.Willow]  0.00 1.00 2520 
b[(Intercept) site:Grassy.Knoll]    0.00 1.00 4430 
b[(Intercept) site:s741]            0.00 1.00 2648 
b[(Intercept) site:s742]            0.00 1.00 4153 
b[(Intercept) site:s743]            0.00 1.00 3711 
b[(Intercept) site:s744]            0.00 1.00 4191 
b[(Intercept) site:s746]            0.00 1.00 2342 
b[(Intercept) site:s747]            0.00 1.00 3296 
b[(Intercept) site:s748]            0.00 1.00 3291 
b[(Intercept) site:s749]            0.00 1.00 3899 
b[(Intercept) site:s750]            0.00 1.00 2770 
b[(Intercept) site:s751]            0.00 1.00 5149 
b[(Intercept) site:Salmon.Sequel]   0.00 1.00 3534 
b[(Intercept) site:SW.Point]        0.00 1.00 2445 
b[(Intercept) site:Willow.Gully]    0.00 1.00 1749 
b[(Intercept) year:2010]            0.00 1.00 2537 
b[islandPI year:2010]               0.00 1.00 3098 
b[(Intercept) year:2011]            0.00 1.00 2850 
b[islandPI year:2011]               0.00 1.00 2938 
b[(Intercept) year:2012]            0.00 1.00 3767 
b[islandPI year:2012]               0.00 1.00 3694 
b[(Intercept) year:2013]            0.00 1.00 2419 
b[islandPI year:2013]               0.00 1.00 1792 
b[(Intercept) year:2014]            0.00 1.00 2251 
b[islandPI year:2014]               0.00 1.00 2355 
b[(Intercept) year:2015]            0.00 1.00 2377 
b[islandPI year:2015]               0.00 1.00 2038 
b[(Intercept) year:2016]            0.00 1.00 2049 
b[islandPI year:2016]               0.00 1.00 2240 
b[(Intercept) year:2017]            0.00 1.00 2436 
b[islandPI year:2017]               0.00 1.00 1774 
b[(Intercept) year:2018]            0.00 1.00 2112 
b[islandPI year:2018]               0.00 1.00 1755 
Sigma[site:(Intercept),(Intercept)] 0.00 1.00 1168 
Sigma[year:(Intercept),(Intercept)] 0.00 1.00 1100 
Sigma[year:islandPI,(Intercept)]    0.00 1.00 1756 
Sigma[year:islandPI,islandPI]       0.00 1.00 1309 
mean_PPD                            0.00 1.00 3287 
log-posterior                       0.25 1.00  678 

For each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).
```

``` r
# Inter-island differences, varying among years, plus PC1 + PC2
occ3 <- stan_glmer(cbind(egg, viable - egg) ~ island + PC1 + PC2 + (1 | site) + (island | year),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)
```

``` r
summary(occ3, prob = c(0.025,0.5,0.975), digits = 2)
```

``` 

Model Info:

 function:     stan_glmer
 family:       binomial [logit]
 formula:      cbind(egg, viable - egg) ~ island + PC1 + PC2 + (1 | site) + 
       (island | year)
 algorithm:    sampling
 priors:       see help('prior_summary')
 sample:       3000 (posterior sample size)
 observations: 139
 groups:       site (17), year (9)

Estimates:
                                      mean    sd      2.5%    50%     97.5%
(Intercept)                            0.34    0.14    0.07    0.34    0.62
islandPI                               0.34    0.18   -0.03    0.34    0.70
PC1                                   -0.04    0.09   -0.21   -0.04    0.14
PC2                                   -0.09    0.15   -0.40   -0.09    0.20
b[(Intercept) site:Boathouse]         -0.10    0.16   -0.44   -0.09    0.20
b[(Intercept) site:Catwalk.Salmon]     0.13    0.17   -0.19    0.12    0.50
b[(Intercept) site:Catwalk.Willow]    -0.15    0.17   -0.52   -0.13    0.16
b[(Intercept) site:Grassy.Knoll]       0.00    0.18   -0.38    0.00    0.36
b[(Intercept) site:s741]              -0.23    0.18   -0.63   -0.21    0.08
b[(Intercept) site:s742]               0.09    0.17   -0.23    0.08    0.48
b[(Intercept) site:s743]               0.11    0.18   -0.21    0.10    0.50
b[(Intercept) site:s744]              -0.09    0.16   -0.43   -0.08    0.22
b[(Intercept) site:s746]              -0.23    0.19   -0.63   -0.22    0.09
b[(Intercept) site:s747]               0.12    0.18   -0.22    0.10    0.50
b[(Intercept) site:s748]              -0.17    0.20   -0.64   -0.15    0.18
b[(Intercept) site:s749]               0.14    0.17   -0.17    0.13    0.52
b[(Intercept) site:s750]               0.20    0.19   -0.14    0.18    0.59
b[(Intercept) site:s751]               0.04    0.19   -0.35    0.03    0.42
b[(Intercept) site:Salmon.Sequel]     -0.04    0.17   -0.40   -0.03    0.31
b[(Intercept) site:SW.Point]          -0.09    0.16   -0.41   -0.08    0.20
b[(Intercept) site:Willow.Gully]       0.26    0.19   -0.07    0.24    0.66
b[(Intercept) year:2010]               0.07    0.15   -0.17    0.05    0.42
b[islandPI year:2010]                  0.05    0.18   -0.32    0.03    0.45
b[(Intercept) year:2011]               0.06    0.17   -0.23    0.03    0.45
b[islandPI year:2011]                  0.01    0.19   -0.39    0.01    0.40
b[(Intercept) year:2012]              -0.01    0.15   -0.34    0.00    0.31
b[islandPI year:2012]                  0.05    0.18   -0.31    0.03    0.44
b[(Intercept) year:2013]               0.03    0.15   -0.24    0.01    0.39
b[islandPI year:2013]                 -0.14    0.20   -0.65   -0.10    0.17
b[(Intercept) year:2014]              -0.12    0.16   -0.53   -0.09    0.09
b[islandPI year:2014]                  0.10    0.20   -0.26    0.06    0.57
b[(Intercept) year:2015]              -0.06    0.20   -0.56   -0.03    0.28
b[islandPI year:2015]                  0.17    0.22   -0.15    0.13    0.68
b[(Intercept) year:2016]               0.10    0.18   -0.20    0.08    0.51
b[islandPI year:2016]                  0.14    0.21   -0.23    0.11    0.60
b[(Intercept) year:2017]              -0.08    0.18   -0.47   -0.05    0.27
b[islandPI year:2017]                 -0.22    0.23   -0.73   -0.17    0.12
b[(Intercept) year:2018]               0.00    0.13   -0.26   -0.01    0.27
b[islandPI year:2018]                 -0.17    0.20   -0.64   -0.13    0.13
Sigma[site:(Intercept),(Intercept)]    0.06    0.05    0.00    0.05    0.20
Sigma[year:(Intercept),(Intercept)]    0.04    0.07    0.00    0.02    0.21
Sigma[year:islandPI,(Intercept)]      -0.01    0.05   -0.14    0.00    0.05
Sigma[year:islandPI,islandPI]          0.08    0.11    0.00    0.05    0.34
mean_PPD                               7.38    0.20    7.00    7.39    7.76
log-posterior                       -331.87    6.89 -346.71 -331.45 -319.35

Diagnostics:
                                    mcse Rhat n_eff
(Intercept)                         0.00 1.00 1692 
islandPI                            0.00 1.00 1410 
PC1                                 0.00 1.00 1348 
PC2                                 0.00 1.00 1600 
b[(Intercept) site:Boathouse]       0.00 1.00 2317 
b[(Intercept) site:Catwalk.Salmon]  0.00 1.00 2543 
b[(Intercept) site:Catwalk.Willow]  0.00 1.00 2567 
b[(Intercept) site:Grassy.Knoll]    0.00 1.00 3300 
b[(Intercept) site:s741]            0.00 1.00 2234 
b[(Intercept) site:s742]            0.00 1.00 2828 
b[(Intercept) site:s743]            0.00 1.00 2781 
b[(Intercept) site:s744]            0.00 1.00 2540 
b[(Intercept) site:s746]            0.00 1.00 2277 
b[(Intercept) site:s747]            0.00 1.00 2721 
b[(Intercept) site:s748]            0.00 1.00 2554 
b[(Intercept) site:s749]            0.00 1.00 2163 
b[(Intercept) site:s750]            0.00 1.00 2113 
b[(Intercept) site:s751]            0.00 1.00 3914 
b[(Intercept) site:Salmon.Sequel]   0.00 1.00 3052 
b[(Intercept) site:SW.Point]        0.00 1.00 2566 
b[(Intercept) site:Willow.Gully]    0.00 1.00 1705 
b[(Intercept) year:2010]            0.00 1.00 2350 
b[islandPI year:2010]               0.00 1.00 2476 
b[(Intercept) year:2011]            0.00 1.00 1809 
b[islandPI year:2011]               0.00 1.00 3186 
b[(Intercept) year:2012]            0.00 1.00 2099 
b[islandPI year:2012]               0.00 1.00 2469 
b[(Intercept) year:2013]            0.00 1.00 2401 
b[islandPI year:2013]               0.00 1.00 1797 
b[(Intercept) year:2014]            0.00 1.00 1220 
b[islandPI year:2014]               0.00 1.00 1827 
b[(Intercept) year:2015]            0.01 1.00 1384 
b[islandPI year:2015]               0.01 1.00 1613 
b[(Intercept) year:2016]            0.00 1.00 1428 
b[islandPI year:2016]               0.00 1.00 2004 
b[(Intercept) year:2017]            0.00 1.00 1386 
b[islandPI year:2017]               0.01 1.00 1488 
b[(Intercept) year:2018]            0.00 1.00 2498 
b[islandPI year:2018]               0.00 1.00 1757 
Sigma[site:(Intercept),(Intercept)] 0.00 1.00 1021 
Sigma[year:(Intercept),(Intercept)] 0.00 1.00  768 
Sigma[year:islandPI,(Intercept)]    0.00 1.00 1223 
Sigma[year:islandPI,islandPI]       0.00 1.00 1087 
mean_PPD                            0.00 1.00 3240 
log-posterior                       0.30 1.00  536 

For each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).
```

``` r
# Inter-island differences plus PC1 + PC2, no random time-variation
occ4 <- stan_glmer(cbind(egg, viable - egg) ~ island + PC1 + PC2 + (1 | site),
                   data =  rhau,
                   family = binomial(link = logit),
                   prior = normal(0,5),
                   prior_intercept = normal(0,5),
                   prior_covariance = decov(),
                   chains = 3, iter = 2000, warmup = 1000, cores = 3)
```

``` r
summary(occ4, prob = c(0.025,0.5,0.975), digits = 2)
```

``` 

Model Info:

 function:     stan_glmer
 family:       binomial [logit]
 formula:      cbind(egg, viable - egg) ~ island + PC1 + PC2 + (1 | site)
 algorithm:    sampling
 priors:       see help('prior_summary')
 sample:       3000 (posterior sample size)
 observations: 139
 groups:       site (17)

Estimates:
                                      mean    sd      2.5%    50%     97.5%
(Intercept)                            0.33    0.12    0.09    0.33    0.58
islandPI                               0.33    0.16    0.02    0.33    0.64
PC1                                   -0.05    0.05   -0.15   -0.05    0.04
PC2                                   -0.06    0.08   -0.22   -0.06    0.10
b[(Intercept) site:Boathouse]         -0.09    0.16   -0.43   -0.08    0.21
b[(Intercept) site:Catwalk.Salmon]     0.12    0.18   -0.20    0.11    0.53
b[(Intercept) site:Catwalk.Willow]    -0.13    0.17   -0.48   -0.11    0.18
b[(Intercept) site:Grassy.Knoll]      -0.01    0.18   -0.39   -0.01    0.38
b[(Intercept) site:s741]              -0.21    0.18   -0.60   -0.19    0.08
b[(Intercept) site:s742]               0.09    0.17   -0.21    0.07    0.44
b[(Intercept) site:s743]               0.11    0.17   -0.20    0.09    0.47
b[(Intercept) site:s744]              -0.08    0.16   -0.43   -0.07    0.21
b[(Intercept) site:s746]              -0.22    0.19   -0.63   -0.20    0.09
b[(Intercept) site:s747]               0.11    0.17   -0.20    0.09    0.48
b[(Intercept) site:s748]              -0.14    0.19   -0.59   -0.11    0.19
b[(Intercept) site:s749]               0.12    0.17   -0.17    0.11    0.49
b[(Intercept) site:s750]               0.18    0.18   -0.12    0.16    0.57
b[(Intercept) site:s751]               0.03    0.17   -0.33    0.02    0.38
b[(Intercept) site:Salmon.Sequel]     -0.03    0.17   -0.38   -0.02    0.30
b[(Intercept) site:SW.Point]          -0.08    0.16   -0.41   -0.07    0.21
b[(Intercept) site:Willow.Gully]       0.24    0.19   -0.05    0.22    0.65
Sigma[site:(Intercept),(Intercept)]    0.06    0.05    0.00    0.05    0.18
mean_PPD                               7.39    0.20    6.99    7.39    7.78
log-posterior                       -303.11    4.82 -313.65 -302.71 -294.75

Diagnostics:
                                    mcse Rhat n_eff
(Intercept)                         0.00 1.00 1505 
islandPI                            0.00 1.00 1505 
PC1                                 0.00 1.00 3724 
PC2                                 0.00 1.00 4066 
b[(Intercept) site:Boathouse]       0.00 1.00 1948 
b[(Intercept) site:Catwalk.Salmon]  0.00 1.00 2207 
b[(Intercept) site:Catwalk.Willow]  0.00 1.00 1757 
b[(Intercept) site:Grassy.Knoll]    0.00 1.00 3539 
b[(Intercept) site:s741]            0.00 1.00 1590 
b[(Intercept) site:s742]            0.00 1.00 2944 
b[(Intercept) site:s743]            0.00 1.00 2840 
b[(Intercept) site:s744]            0.00 1.00 2486 
b[(Intercept) site:s746]            0.00 1.00 1475 
b[(Intercept) site:s747]            0.00 1.00 3196 
b[(Intercept) site:s748]            0.00 1.00 2609 
b[(Intercept) site:s749]            0.00 1.00 2597 
b[(Intercept) site:s750]            0.00 1.00 2197 
b[(Intercept) site:s751]            0.00 1.00 3782 
b[(Intercept) site:Salmon.Sequel]   0.00 1.00 2652 
b[(Intercept) site:SW.Point]        0.00 1.00 2010 
b[(Intercept) site:Willow.Gully]    0.00 1.00 1430 
Sigma[site:(Intercept),(Intercept)] 0.00 1.00  867 
mean_PPD                            0.00 1.00 3212 
log-posterior                       0.21 1.00  552 

For each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).
```

We compare the expected out-of-sample predictive performance of the
candidate models using the Bayesian Pareto-smoothed approxiamte
leave-one-out cross-validation score…

``` r
## Model selection by approximate leave-one-out cross-validation
occ_loos <-  lapply(list(occ0 = occ0, occ1 = occ1, occ2 = occ2, occ3 = occ3, occ4 = occ4), loo)

# unpaired comparisons
occ_compare <- loo_compare(occ_loos)

# pairwise comparisons
occ_compair <- loo_compair(occ_loos)
```

``` r
occ_loos
```

    $occ0
    
    Computed from 3000 by 139 log-likelihood matrix
    
             Estimate   SE
    elpd_loo   -275.7  8.0
    p_loo        15.9  1.7
    looic       551.4 16.0
    ------
    Monte Carlo SE of elpd_loo is 0.1.
    
    Pareto k diagnostic values:
                             Count Pct.    Min. n_eff
    (-Inf, 0.5]   (good)     137   98.6%   1073      
     (0.5, 0.7]   (ok)         2    1.4%   1284      
       (0.7, 1]   (bad)        0    0.0%   <NA>      
       (1, Inf)   (very bad)   0    0.0%   <NA>      
    
    All Pareto k estimates are ok (k < 0.7).
    See help('pareto-k-diagnostic') for details.
    
    $occ1
    
    Computed from 3000 by 139 log-likelihood matrix
    
             Estimate   SE
    elpd_loo   -275.7  7.9
    p_loo        15.3  1.6
    looic       551.4 15.9
    ------
    Monte Carlo SE of elpd_loo is 0.1.
    
    All Pareto k estimates are good (k < 0.5).
    See help('pareto-k-diagnostic') for details.
    
    $occ2
    
    Computed from 3000 by 139 log-likelihood matrix
    
             Estimate   SE
    elpd_loo   -275.1  7.8
    p_loo        17.4  1.8
    looic       550.2 15.6
    ------
    Monte Carlo SE of elpd_loo is 0.1.
    
    All Pareto k estimates are good (k < 0.5).
    See help('pareto-k-diagnostic') for details.
    
    $occ3
    
    Computed from 3000 by 139 log-likelihood matrix
    
             Estimate   SE
    elpd_loo   -275.4  7.8
    p_loo        18.8  1.9
    looic       550.8 15.6
    ------
    Monte Carlo SE of elpd_loo is 0.1.
    
    Pareto k diagnostic values:
                             Count Pct.    Min. n_eff
    (-Inf, 0.5]   (good)     137   98.6%   461       
     (0.5, 0.7]   (ok)         2    1.4%   1683      
       (0.7, 1]   (bad)        0    0.0%   <NA>      
       (1, Inf)   (very bad)   0    0.0%   <NA>      
    
    All Pareto k estimates are ok (k < 0.7).
    See help('pareto-k-diagnostic') for details.
    
    $occ4
    
    Computed from 3000 by 139 log-likelihood matrix
    
             Estimate   SE
    elpd_loo   -277.7  8.1
    p_loo        13.1  1.4
    looic       555.5 16.1
    ------
    Monte Carlo SE of elpd_loo is 0.1.
    
    All Pareto k estimates are good (k < 0.5).
    See help('pareto-k-diagnostic') for details.

``` r
print(occ_compare, simplify = FALSE)
```

``` 
     elpd_diff se_diff elpd_loo se_elpd_loo p_loo  se_p_loo looic  se_looic
occ2    0.0       0.0  -275.1      7.8        17.4    1.8    550.2   15.6  
occ3   -0.3       0.8  -275.4      7.8        18.8    1.9    550.8   15.6  
occ0   -0.6       2.2  -275.7      8.0        15.9    1.7    551.4   16.0  
occ1   -0.6       1.3  -275.7      7.9        15.3    1.6    551.4   15.9  
occ4   -2.7       2.2  -277.7      8.1        13.1    1.4    555.5   16.1  
```

``` r
occ_compair
```

    $`occ3 vs. occ2`
    looic_diff         se 
           0.6        1.5 
    
    $`occ0 vs. occ3`
    looic_diff         se 
           0.6        4.7 
    
    $`occ1 vs. occ0`
    looic_diff         se 
           0.0        3.1 
    
    $`occ4 vs. occ1`
    looic_diff         se 
           4.0        3.7 

If we have refit the models, we cache the results to save time…

``` r
## Cache stanfits and loo objects
save(list = grep("occ", ls(), value = TRUE), file = here::here("analysis","cache","occ_stanfit_loo.RData"))
```
