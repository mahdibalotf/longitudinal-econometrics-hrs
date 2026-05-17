# Longitudinal Analysis of Cognitive Trajectories (HRS)

## Overview
This repository implements an econometric and data engineering pipeline using 11 waves of the RAND HRS dataset to analyze cognitive aging. It addresses look-ahead bias, temporal violations, and informative attrition.

## Methodological Highlights
* **Exposure Engineering:** Built a time-updated absorbing state for surgical menopause and a mode/median resolution algorithm for recall errors.
* **Longitudinal Modeling:** Fit Linear Mixed-Effects Models (LMM) with random intercepts and slopes using `lme4`.
* **Bias Correction:** Applied Inverse Probability of Censoring Weighting (IPCW) to adjust for attrition and mortality.
* **Visualizations:** Generated standardized predictive margins (`ggeffects`) and Kaplan-Meier retention curves.

## Repository Structure
* `phase1_phase2_integrity_cohort.R`: Data cleaning, temporal checks, and cohort selection.
* `phase3_phase4_descriptives_lmm.R`: Baseline SMD diagnostics and primary LMM growth models.
* `phase5_phase6_ipcw_sensitivity.R`: IPCW weight construction, weighted LMMs, and quadratic sensitivity.
* `phase7_final_predictions_visualizations.R`: Marginal predictions, clinical contrasts, and survival curves.

## Tech Stack
* R (`tidyverse`, `lme4`, `lmerTest`, `ggeffects`, `survival`)

*Note: The scripts are designed to process the RAND HRS longitudinal files. Raw datasets are not included in this repository to comply with HRS data use agreements.*
