Reproductive Success of Two Rhinoceros Auklet Colonies in the Salish Sea
and California Current
================
Eric Wagner, Eric Buhle, …
2019-04-28

  - [Overview](#overview)
  - [Setup and data](#setup-and-data)
  - [Principal Components Analysis of Oceanographic
    Indicators](#principal-components-analysis-of-oceanographic-indicators)
  - [GLMMs of Burrow Occupancy](#glmms-of-burrow-occupancy)
      - [Model Fitting and Assessment](#model-fitting-and-assessment)
      - [Model Selection](#model-selection)
      - [Plots](#plots)
  - [GLMMs of Apparent Fledging
    Success](#glmms-of-apparent-fledging-success)
      - [Model Fitting and Assessment](#model-fitting-and-assessment-1)
      - [Model Selection](#model-selection-1)
      - [Plots](#plots-1)

# Overview

This is a Rhinoceros Auklet from the Protection Island
colony:

<img src="https://www.eopugetsound.org/sites/default/files/styles/magazinewidth_592px/public/topical_articles/images/16343412270_5cfaa5c480_o.jpg?itok=WjzI1_2K" width="50%" style="display: block; margin: auto auto auto 0;" />

\[Cap’n: Add elegant prose about the study system, questions, and brief
methods here; add subsections if desired\]

# Setup and data

First we load the libraries and functions we’ll use.

Next we’ll read in the datasets from various sources and manipulate them
into a usable format.

Let’s take a look at the nest
data…

| year | island | site           | viable | occupied | egg | hatch | last\_check |
| ---: | :----- | :------------- | -----: | -------: | --: | ----: | ----------: |
| 2010 | DI     | Boathouse      |     22 |       15 |  15 |    13 |          10 |
| 2010 | DI     | Catwalk.Salmon |     13 |        7 |   5 |     5 |           4 |
| 2010 | DI     | Catwalk.Willow |     13 |        8 |   8 |     8 |           8 |
| 2010 | DI     | Salmon.Sequel  |     11 |        7 |   7 |     5 |           3 |
| 2010 | DI     | SW.Point       |     14 |       10 |   8 |     8 |           7 |
| 2010 | DI     | Willow.Gully   |     11 |       10 |  10 |    10 |          10 |

And here are the environmental covariates (describe oceanographic
indicators)…

| year | pdo\_index | mei\_avg | sst\_DI\_spring | sst\_PI\_spring | sst\_amph\_spring | sst\_race\_spring | cui\_spring | chla\_DI\_spring | chla\_PI\_spring | st\_onset | st\_duration |
| ---: | ---------: | -------: | --------------: | --------------: | ----------------: | ----------------: | ----------: | ---------------: | ---------------: | --------: | -----------: |
| 2002 |     \-0.90 |     0.11 |              NA |              NA |             10.36 |              9.60 |       12.33 |               NA |               NA |       108 |          197 |
| 2003 |       1.06 |     0.32 |           12.11 |           10.35 |             10.90 |             10.08 |       10.33 |             8.76 |             4.56 |       156 |          113 |
| 2004 |     \-0.28 |     0.08 |           13.38 |           10.59 |             11.45 |             10.22 |       11.00 |             5.49 |             5.35 |       146 |          142 |
| 2005 |     \-0.20 |     0.27 |           13.32 |           10.64 |             11.74 |             10.47 |      \-0.33 |             2.70 |             3.43 |       238 |           33 |
| 2006 |     \-0.16 |   \-0.35 |           11.81 |           10.22 |             11.03 |              9.90 |       16.67 |             7.01 |             4.82 |       158 |          146 |
| 2007 |     \-0.42 |   \-0.12 |           11.83 |            9.46 |             10.26 |              9.72 |       10.67 |             8.67 |             4.78 |        81 |          263 |
| 2008 |     \-1.05 |   \-1.11 |           10.88 |            8.71 |             10.34 |              9.07 |       17.33 |             6.53 |             4.84 |        64 |          237 |
| 2009 |     \-1.40 |   \-0.51 |           11.84 |            9.96 |             10.92 |              9.43 |       15.67 |             7.37 |             6.04 |        65 |          270 |
| 2010 |     \-0.08 |   \-0.16 |           11.90 |            9.90 |             10.93 |              9.59 |       10.00 |             4.52 |             3.60 |       169 |          159 |
| 2011 |     \-1.30 |   \-1.52 |           11.69 |            9.12 |             10.77 |              9.40 |        7.67 |             3.29 |             3.82 |        82 |          190 |
| 2012 |     \-1.70 |   \-0.60 |           11.40 |            9.26 |             10.27 |              9.44 |      \-0.67 |             5.93 |             4.52 |       125 |          174 |
| 2013 |     \-1.14 |   \-0.38 |           12.63 |            9.92 |             10.76 |              9.01 |        8.00 |            10.08 |             5.47 |        91 |          178 |
| 2014 |     \-0.44 |   \-0.14 |           12.92 |           10.55 |             11.41 |              9.61 |        4.00 |             9.19 |             6.62 |       126 |          142 |
| 2015 |       1.44 |     0.85 |           12.00 |           10.97 |             11.64 |             10.15 |       33.33 |            10.84 |             7.15 |       365 |            0 |
| 2016 |       0.76 |     1.07 |           12.75 |           10.96 |             11.81 |             10.36 |       20.00 |             8.98 |             5.65 |       365 |            0 |
| 2017 |       0.16 |   \-0.45 |           12.18 |           10.26 |             11.19 |             10.18 |      \-5.33 |             4.96 |             3.68 |       179 |           98 |
| 2018 |     \-0.29 |   \-0.52 |           12.86 |           10.14 |             11.20 |              9.92 |        8.00 |             6.66 |             5.84 |       149 |          182 |

# Principal Components Analysis of Oceanographic Indicators

Let’s explore the patterns of (a)synchrony among the oceanographic
indicators to see how severe the multicollinearity might be if they were
used as raw regression
inputs.

<img src="RhAu_BayesGLMM_files/figure-gfm/env_corrplot-1.png" width="50%" style="display: block; margin: auto;" />

Now we perform a principal components analysis to extract the major
trends in the suite of oceanographic indicators for use as regression
inputs.

    Importance of components:
                              PC1    PC2    PC3     PC4     PC5     PC6     PC7     PC8
    Standard deviation     2.0729 1.3625 0.9925 0.61922 0.45791 0.39276 0.26794 0.20679
    Proportion of Variance 0.5371 0.2321 0.1231 0.04793 0.02621 0.01928 0.00897 0.00535
    Cumulative Proportion  0.5371 0.7691 0.8923 0.94019 0.96640 0.98568 0.99465 1.00000

    Joining, by = "year"

<img src="RhAu_BayesGLMM_files/figure-gfm/PCA_plots-1.png" width="75%" style="display: block; margin: auto;" />

To help interpret the principal components, let’s see how they track the
time series of oceanographic
indicators…

<img src="RhAu_BayesGLMM_files/figure-gfm/env_timeseries_plots-1.png" width="50%" style="display: block; margin: auto;" />

# GLMMs of Burrow Occupancy

## Model Fitting and Assessment

Now we are ready to start fitting some GLMMs to the burrow occupancy
data…\[more on modeling, explanations of each model (starting with model
2 as a baseline description of all known variance components in the
study design)\]

Model 1

    stan_glmer
     family:       binomial [logit]
     formula:      cbind(egg, viable - egg) ~ island + (1 | site) + (island | year)
     observations: 139
    ------
                Median MAD_SD
    (Intercept) 0.3    0.1   
    islandPI    0.3    0.2   
    
    Error terms:
     Groups Name        Std.Dev. Corr 
     site   (Intercept) 0.25          
     year   (Intercept) 0.17          
            islandPI    0.24     -0.24
    Num. levels: site 17, year 9 
    
    Sample avg. posterior predictive distribution of y:
             Median MAD_SD
    mean_PPD 7.4    0.2   
    
    ------
    * For help interpreting the printed output see ?print.stanreg
    * For info on the priors used see ?prior_summary.stanreg

Model 2

    stan_glmer
     family:       binomial [logit]
     formula:      cbind(egg, viable - egg) ~ island + PC1 + PC2 + (1 | site) + 
           (island | year)
     observations: 139
    ------
                Median MAD_SD
    (Intercept)  0.3    0.1  
    islandPI     0.3    0.2  
    PC1          0.0    0.1  
    PC2         -0.1    0.1  
    
    Error terms:
     Groups Name        Std.Dev. Corr 
     site   (Intercept) 0.25          
     year   (Intercept) 0.20          
            islandPI    0.27     -0.11
    Num. levels: site 17, year 9 
    
    Sample avg. posterior predictive distribution of y:
             Median MAD_SD
    mean_PPD 7.4    0.2   
    
    ------
    * For help interpreting the printed output see ?print.stanreg
    * For info on the priors used see ?prior_summary.stanreg

Model 3

    stan_glmer
     family:       binomial [logit]
     formula:      cbind(egg, viable - egg) ~ island + PC1 + PC2 + (1 | site)
     observations: 139
    ------
                Median MAD_SD
    (Intercept)  0.3    0.1  
    islandPI     0.3    0.2  
    PC1         -0.1    0.0  
    PC2         -0.1    0.1  
    
    Error terms:
     Groups Name        Std.Dev.
     site   (Intercept) 0.23    
    Num. levels: site 17 
    
    Sample avg. posterior predictive distribution of y:
             Median MAD_SD
    mean_PPD 7.4    0.2   
    
    ------
    * For help interpreting the printed output see ?print.stanreg
    * For info on the priors used see ?prior_summary.stanreg

We can assess the adequacy of the “full” model `occ3` by generating
replicate pseudo-data from its posterior predictive distribution and
comparing them to the real
data…

<img src="RhAu_BayesGLMM_files/figure-gfm/plot_ppd_occ-1.png" width="50%" style="display: block; margin: auto;" />

## Model Selection

We compare the expected out-of-sample predictive performance of the
candidate models using the Bayesian Pareto-smoothed approxiamte
leave-one-out cross-validation
score…

``` 
     elpd_diff se_diff elpd_loo se_elpd_loo p_loo  se_p_loo looic  se_looic
occ2    0.0       0.0  -274.8      7.8        17.0    1.8    549.5   15.5  
occ0   -0.9       2.1  -275.7      8.0        15.9    1.7    551.4   16.0  
occ1   -1.1       1.3  -275.9      7.9        15.0    1.6    551.8   15.9  
occ3   -1.2       0.8  -275.9      7.9        19.0    2.0    551.8   15.7  
occ4   -3.1       2.2  -277.8      8.1        13.0    1.4    555.6   16.2  
```

    $`occ0 vs. occ2`
    looic_diff         se 
           1.9        4.3 
    
    $`occ1 vs. occ0`
    looic_diff         se 
           0.4        3.3 
    
    $`occ3 vs. occ1`
    looic_diff         se 
           0.1        3.3 
    
    $`occ4 vs. occ3`
    looic_diff         se 
           3.8        4.9 

## Plots

We can explore the predictions of the selected model by plotting the
posterior distribution of fitted values as time series for each island
and comparing them to the sample proportion of occupied burrows (with
binomial 95% confidence
intervals).

<img src="RhAu_BayesGLMM_files/figure-gfm/occ_fitted_timeseries-1.png" width="60%" style="display: block; margin: auto;" />

# GLMMs of Apparent Fledging Success

## Model Fitting and Assessment

Model 1

    stan_glmer
     family:       binomial [logit]
     formula:      cbind(last_check, egg - last_check) ~ island + (1 | site) + (island | 
           year)
     observations: 139
    ------
                Median MAD_SD
    (Intercept)  1.9    0.3  
    islandPI    -0.6    0.3  
    
    Error terms:
     Groups Name        Std.Dev. Corr 
     site   (Intercept) 0.38          
     year   (Intercept) 0.61          
            islandPI    0.59     -0.23
    Num. levels: site 17, year 9 
    
    Sample avg. posterior predictive distribution of y:
             Median MAD_SD
    mean_PPD 5.9    0.1   
    
    ------
    * For help interpreting the printed output see ?print.stanreg
    * For info on the priors used see ?prior_summary.stanreg

Model 2

    stan_glmer
     family:       binomial [logit]
     formula:      cbind(last_check, egg - last_check) ~ island + PC1 + PC2 + (1 | 
           site) + (island | year)
     observations: 139
    ------
                Median MAD_SD
    (Intercept)  1.9    0.3  
    islandPI    -0.6    0.3  
    PC1          0.1    0.2  
    PC2          0.0    0.4  
    
    Error terms:
     Groups Name        Std.Dev. Corr 
     site   (Intercept) 0.38          
     year   (Intercept) 0.79          
            islandPI    0.62     -0.31
    Num. levels: site 17, year 9 
    
    Sample avg. posterior predictive distribution of y:
             Median MAD_SD
    mean_PPD 5.9    0.1   
    
    ------
    * For help interpreting the printed output see ?print.stanreg
    * For info on the priors used see ?prior_summary.stanreg

Model 3

    stan_glmer
     family:       binomial [logit]
     formula:      cbind(last_check, egg - last_check) ~ island + PC1 + PC2 + (1 | 
           site)
     observations: 139
    ------
                Median MAD_SD
    (Intercept)  1.9    0.2  
    islandPI    -0.6    0.3  
    PC1          0.2    0.1  
    PC2          0.1    0.1  
    
    Error terms:
     Groups Name        Std.Dev.
     site   (Intercept) 0.38    
    Num. levels: site 17 
    
    Sample avg. posterior predictive distribution of y:
             Median MAD_SD
    mean_PPD 5.9    0.1   
    
    ------
    * For help interpreting the printed output see ?print.stanreg
    * For info on the priors used see ?prior_summary.stanreg

Marginal posterior predictive density of largest
model

<img src="RhAu_BayesGLMM_files/figure-gfm/plot_ppd_suc-1.png" width="50%" style="display: block; margin: auto;" />

## Model Selection

``` 
     elpd_diff se_diff elpd_loo se_elpd_loo p_loo  se_p_loo looic  se_looic
suc2    0.0       0.0  -200.1      9.8        21.4    2.6    400.1   19.5  
suc3   -0.4       0.6  -200.5      9.9        22.7    2.7    401.0   19.8  
suc1   -2.1       2.2  -202.1     10.3        18.6    2.3    404.3   20.7  
suc0   -3.5       3.2  -203.6     10.3        20.9    2.8    407.1   20.7  
suc4  -14.5       5.7  -214.6     11.2        16.5    2.2    429.2   22.4  
```

    $`suc3 vs. suc2`
    looic_diff         se 
           0.8        1.2 
    
    $`suc1 vs. suc3`
    looic_diff         se 
           3.3        4.5 
    
    $`suc0 vs. suc1`
    looic_diff         se 
           2.8        4.5 
    
    $`suc4 vs. suc0`
    looic_diff         se 
          22.1       11.1 

## Plots

We can explore the predictions of the selected model by plotting the
posterior distribution of fitted values as time series for each island
and comparing them to the sample proportion of occupied burrows (with
binomial 95% confidence
intervals).

<img src="RhAu_BayesGLMM_files/figure-gfm/suc_fitted_timeseries-1.png" width="60%" style="display: block; margin: auto;" />
