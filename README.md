
Code of the Paper "Unlocking the Power of Time-Since-Infection Models: Data Augmentation for Improved Instantaneous Reproduction Number Estimation"
==============================================
  
  
## Outline
1. Description
2. Sudo Algorithm
3. Package requirements
4. Code for Simulation 
5. Code for Real Data
6. How to run the Algorithm on your data?
  
## Description
This README is prepared for journal peer review of the "Unlocking the Power of Time-Since-Infection Models: Data Augmentation for Improved Instantaneous Reproduction Number Estimation" paper. 

The proposed refined TSI model integrates hospitalization data and breaks the reliance on incidence data only. The correspondingly proposed MCEM algorithm builds on a composite likelihood approach and is proposed for simultaneously estimating the instantaneous reproduction number and the generation time (infectiousness profile and hospitalization propensity) of a transmission disease that meets the basic assumptions of the time-since-infection model with daily incident cases, covariates data (risk factors), and hospitalization data. It contributes to the field by unlocking the potential of TSI models through data augmentation, which refines real-time estimation of Rt, enables the estimation of hospitalization-related parameters previously accessible primarily through contact tracing, and facilitates analysis of associations between risk factors and disease transmission dynamics. 

## Data structure 
![](Figure1_3row.png)

## Package Requirements
- A database with clear and consistent variable names
- R version: R (>= 4.3.1)
- On Mac: download and install [EpiEstim](https://CRAN.R-project.org/package=EpiEstim), [dlnm](https://CRAN.R-project.org/package=dlnm), [tsModel](https://CRAN.R-project.org/package=tsModel), [ggplot2](https://CRAN.R-project.org/package=ggplot2), [nlme](https://CRAN.R-project.org/package=nlme), [sjstats](https://CRAN.R-project.org/package=sjstats), [corrplot](https://CRAN.R-project.org/package=corrplot), [pracma](https://CRAN.R-project.org/package=pracma), [dplyr](https://CRAN.R-project.org/package=dplyr)
## Run QSOEID example with code

##  Code for Simulation 
Code for simulation can be found in the file "R Code" and "sampled data"

##  Code for Real Data
Code for simulation can be found in the file "Covid_code_real_data_organize", where the precleaned data (only part of the real data is shared due to privacy) contains data from four metropolitans, 
  #Miami-Dade 12086
  #Wayne MI 26163
  #NYC NY 36061
  #COOK IL 17031

## How to run this MCEM algorithm on your data?

* Get the data ready, which requires no missing values and clear variable names. 
* Set the tunning parameters `NoCov, T, R_0, I_0` so the model is specified. The input data includes covariates Z, being a (T \times NoCov) matrix, incident cases I, and hospitalization H, both being a (1 \times T) matrix. 
