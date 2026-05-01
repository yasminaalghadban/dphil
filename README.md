# DPhil Code Repository
An R-based workflow for gestational diabetes cohort construction, predictor derivation, model external validation, and survival modelling using linked primary care and hospital datasets.

This repository has been cleaned for GitHub use by removing private local file paths, hardcoded usernames, passwords, and server details; and using project-relative paths where possible. The analytical logic has been preserved as closely as possible, but the code is now easier to review, rerun, and maintain.

### CPRD Dataset Creation: Build the dataset used for Chapters 3 & 5
The dataset creation pipeline:
- identifies women with first-ever GDM from primary care and hospital data
- restricts to first GDM diagnosis from 2010 onward
- links GDM diagnosis to a likely index pregnancy
- derives covariates including maternal age, deprivation, substance use, gestational age, parity, SMM, prior pregnancy complications, hypertensive disorders, BMI, ethnicity, smoking, alcohol use, treatment method, and delivery method
- determines post-delivery diabetes and complication outcomes
- generates the exclusion dataset for women with evidence of diabetes before the index pregnancy
- exports the final analysis dataset

### Chapter 3 analysis
The analysis pipeline:
- applies study exclusions
- standardises variables and categories
- builds time-to-event datasets
- creates spline terms for maternal age and gestational age
- generates included-versus-excluded summaries
- estimates incidence rates
- fits Cox proportional hazards models
- calculates discrimination and recalibration summaries
- builds risk-group tables

### Chapter 5 analysis
The analysis pipeline:
- reconstructs primary care versus secondary care GDM ascertainment groups
- compares baseline characteristics by ascertainment source
- examines overlap between PCR and HES cohorts
- evaluates patterns of missing data
- performs multiple imputation using `mice`
- fits pooled logistic and Cox models across imputed datasets
- generates adjusted survival curves and proportional hazards diagnostics
