#  Clinical Trial Biostatistics & Epidemiology Engine
**Advanced R Pipeline for Longitudinal Data, Survival Analysis, and Missing Data Imputation**

##  Project Overview
This repository contains a programmatic biostatistics pipeline designed to simulate and analyze Phase III clinical trial data. The engine models a retrospective cohort of 5,000 Type 2 Diabetes Mellitus (T2DM) patients comparing the cardiovascular and glycemic outcomes of **SGLT2 inhibitors vs. DPP4 inhibitors**.

The pipeline demonstrates advanced handling of real-world, noisy clinical data using statistical methodologies frequently required in doctoral epidemiological research and FDA/EMA trial submissions.

---

##  Advanced Methodologies Implemented

### 1. Missing Data Handling: Multiple Imputation (MICE)
Real-world clinical data is rarely complete. To comply with ICH E9(R1) guidelines, this engine simulates "Missing At Random" (MAR) mechanics (e.g., older patients missing BMI records). It utilizes **Multivariate Imputation by Chained Equations (MICE)** via Predictive Mean Matching (PMM) to preserve statistical power and reduce attrition bias without relying on flawed listwise deletion.

### 2. Propensity Score Matching (PSM)
To mimic the covariate balance of a randomized controlled trial (RCT) within observational data, the engine applies Nearest-Neighbor PSM. It matches patients 1:1 based on Age, Sex, Baseline eGFR, HbA1c, and duration of diabetes using a strict caliper to eliminate selection bias.

### 3. Survival Analysis (Cox Proportional Hazards)
Analyzes the time-to-event for Major Adverse Cardiovascular Events (MACE). The engine computes hazard ratios adjusted for baseline confounders and generates Kaplan-Meier survival curves to visualize freedom-from-MACE over a 36-month follow-up period.

### 4. Linear Mixed-Effects Models (LMM)
Evaluates longitudinal efficacy. Rather than simple pre/post t-tests, the pipeline fits an LMM to track HbA1c trajectories over 24 months (Visits at Month 0, 6, 12, 18, 24). It incorporates random intercepts for individual patients to account for intra-patient correlation over time.

---

##  Technical Stack
* **Language:** R
* **Missing Data:** `mice`, `VIM`
* **Matching:** `MatchIt`
* **Survival Analysis:** `survival`, `survminer`
* **Longitudinal Modeling:** `lme4`, `broom.mixed`
* **Data Wrangling:** `tidyverse`

##  Execution
To run the full simulation and analysis pipeline:
```R
source("survival_psm_lmm_analysis.R")
