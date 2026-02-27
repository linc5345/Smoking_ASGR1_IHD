# Proteome-wide genetic analyses identify ASGR1 as a potential mediator of smoking on ischaemic heart disease
Alexander C Tinworth, Pang Yao, Andri Iona, Alfred Pozarickij, Adam Von Ende, Iona Millwood, Robin G Walters, Robert Clarke, Fiona Bragg, Zhengming Chen

This repository contains code to perform the primary MR analyses and ABF colocalization in the above study.

## Stage 1: Smoking and plasma levels of 2922 proteins from the Olink panel
Requirements: 
  1. GSCAN and/or UKB summary statistics for Smoking intitiation, cigarettes per day, and the lifetime smoking index
  2. Summary statistics for circulating plasma protein levels (UKB-PPP used in our analyses).
  3. A file containing gene coordinate information for proteins (to define cis regions for exclusion)
  4. Perform MR function (https://github.com/adamkvonende/mr_pipeline)

Run: Smok_to_prot.R

## Stage 1: Reverse _cis_-MR of protein levels to smoking variables
Requirements:
  1. Summary statistics for circulating plasma protein levels (UKB-PPP used in our analyses)
  2. A file containing gene coordinate information for proteins (to define cis regions for inclusion)
  3. GSCAN and/or UKB summary statistics for Smoking intitiation, cigarettes per day, and the lifetime smoking index
  4. Perform MR function (https://github.com/adamkvonende/mr_pipeline)

Run: Prot_to_smok.R

## Stage 2: _cis_-MR of protein levels to disease risk
Requirements:
  1. Summary statistics for circulating plasma protein levels (UKB-PPP used in our analyses)
  2. A file containing gene coordinate information for proteins (to define cis regions for inclusion)
  3. Disease outcome summary statistics (see manuscript for details)
  4. Perform MR function (https://github.com/adamkvonende/mr_pipeline)

Run: Prot_to_disease.R

## Stage 2: Colocalization ABF
Requirements:
   1. Summary statistics for circulating plasma protein levels (UKB-PPP used in our analyses)
   2. Disease outcome summary statistics (see manuscript for details)

Run: Coloc_abf.R

## Stage 3: MVMR of smoking to protein levels adjusting for BMI and alcohol consumption
Requirments:
  1. GSCAN and/or UKB summary statistics for Smoking intitiation or the lifetime smoking index
  2. Summary statistics for circulating plasma protein levels (UKB-PPP used in our analyses)
  3. Summary statistics for BMI and/or alcohol consumption (see manuscript for details)

Run: MVMR.R

## Stage 4: Example calculation of proportion mediated
Requirements:
  1. File containing beta and se for exposure -> mediator, mediator -> outcome, and exposure -> outcome associations
Run: Mediation_calc.R
