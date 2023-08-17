# Longitudinal MR Simulation Study

Code and settings for simulation of longitudinal exposure and MVMR of its mean and variability

## Idea

Using longitudinal data of an exposure *X*, one can reduce its dimensions by using the mean and standard deviation (SD) of *X* per individual as proxies. Both might have an causal effect on an outcome *Y*

Here, I want to store the code to run the simulation, and all relevant setting files.

## Global settings

I want to test 4 global settings

1) shared SNP set for the exposure (likely) and no missing data points (unlikely)
2) shared SNP set for the exposure (likely) and missing data points (likely)
2) distinct SNP sets for the exposure (unlikely) and no missing data points (unlikely)
2) distinct SNP sets for the exposure (unlikely) and missing data points (likely)

## GX association models

1) Mean and SD per individual over time-series (2 exposures, using SNP effect on mean and SD as instruments)
2) Weights of the eigenfunctions per individual (1-3 exposures, depending on fPCA, using SNP effect on weights as instruments)
3) Linear mixed model over time-series (2 exposures, using SNP effect on intercept and slope (time-interaction) as instruments)
4) Generalized additive models for location, scale and shape (gamlss, 2 exposures, using SNP effect on $\mu$ and $\sigma$ as instruments)
