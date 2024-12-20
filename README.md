
Code of the Paper "Unlocking the Power of Time-Since-Infection Models: Data Augmentation for Improved Instantaneous Reproduction Number Estimation"
==============================================
  
  
## Outline
1. Description
2. Data Structure and Results
3. Package Requirements
4. Code for Simulation 
5. Code for Real Data
6. How to run this MCEM algorithm on your data?
  
## Description
This README is prepared for journal peer review of the "Unlocking the Power of Time-Since-Infection Models: Data Augmentation for Improved Instantaneous Reproduction Number Estimation" paper. 

![Proposed Algorithm](hospitalization_RNN.jpg)

The proposed refined TSI model integrates hospitalization data and breaks the reliance on incidence data only. The correspondingly proposed MCEM algorithm builds on a composite likelihood approach and is proposed for simultaneously estimating the instantaneous reproduction number and the generation time (infectiousness profile and hospitalization propensity) of a transmission disease that meets the basic assumptions of the time-since-infection model with daily incident cases, covariates data (risk factors), and hospitalization data. It contributes to the field by unlocking the potential of TSI models through data augmentation, which refines real-time estimation of Rt, enables the estimation of hospitalization-related parameters previously accessible primarily through contact tracing, and facilitates analysis of associations between risk factors and disease transmission dynamics. 

## Data Structure and Results
![Data Structure](Figure1_3row.png)

![Results](Case_1_revision_wdw7.jpeg)

![Results](Case4EstR.jpeg)

## Package Requirements
- A database with clear and consistent variable names
- R version: R (>= 4.3.1)
- On Mac: download and install [EpiEstim](https://CRAN.R-project.org/package=EpiEstim), [dlnm](https://CRAN.R-project.org/package=dlnm), [tsModel](https://CRAN.R-project.org/package=tsModel), [ggplot2](https://CRAN.R-project.org/package=ggplot2), [nlme](https://CRAN.R-project.org/package=nlme), [sjstats](https://CRAN.R-project.org/package=sjstats), [corrplot](https://CRAN.R-project.org/package=corrplot), [pracma](https://CRAN.R-project.org/package=pracma), [dplyr](https://CRAN.R-project.org/package=dplyr), [mvtnorm](https://CRAN.R-project.org/package=mvtnorm), [mc2d](https://CRAN.R-project.org/package=mc2d), [TTR](https://CRAN.R-project.org/package=TTR), [numDeriv](https://CRAN.R-project.org/package=numDeriv)
## Run QSOEID example with code

##  Code for Simulation 
Code for simulation can be found in the file "R Code" and "sampled data", one illustrational code is "hospitalization_Case1_main_function.R".

##  Code for Real Data
Code for simulation can be found in the file "Covid_code_real_data_organize", where the precleaned data (only part of the real data is shared due to privacy) contains data from four metropolitans, 
  #Miami-Dade 12086
  #Wayne MI 26163
  #NYC NY 36061
  #COOK IL 17031

## How to run this MCEM algorithm on your data?

* Get the data ready, which requires no missing values and clear variable names, one of the exampled data can be found in the file "sampled data". 
* Set the tunning parameters `NoCov, T, R_0, I_0` so the model is specified. The input data includes covariates Z, being a (T \times NoCov) matrix, incident cases I, and hospitalization H, both being a (1 \times T) matrix. 
