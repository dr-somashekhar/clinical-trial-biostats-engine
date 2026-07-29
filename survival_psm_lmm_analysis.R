# =========================================================================================
# 📊 ADVANCED CLINICAL BIOSTATISTICS & EPIDEMIOLOGY ENGINE
# 🧬 Focus: Survival Analysis, Propensity Score Matching (PSM), and Linear Mixed Models
# 
# Description: This script processes longitudinal clinical trial data comparing 
# SGLT2 inhibitors vs. DPP4 inhibitors in Type 2 Diabetes Mellitus (T2DM). 
# It utilizes advanced biostatistical methods required for PhD-level research.
# =========================================================================================

# -----------------------------------------------------------------------------------------
# SECTION 1: LIBRARY INITIALIZATION
# -----------------------------------------------------------------------------------------
# Required packages for advanced statistical modeling
library(tidyverse)    # Data manipulation and visualization (ggplot2, dplyr)
library(survival)     # Core survival analysis functions
library(survminer)    # Elegant Kaplan-Meier plotting
library(MatchIt)      # Propensity score matching to reduce selection bias
library(lme4)         # Linear mixed-effects models for longitudinal data
library(broom.mixed)  # Tidy outputs for mixed models
library(gtsummary)    # Publication-ready summary tables

# -----------------------------------------------------------------------------------------
# SECTION 2: LARGE-SCALE COHORT SIMULATION (N = 5000)
# -----------------------------------------------------------------------------------------
set.seed(2026) # Setting seed for reproducibility

cat("\n[1] INITIALIZING CLINICAL COHORT SIMULATION...\n")
n_patients <- 5000

# Generating baseline covariates with intentional confounding
clinical_cohort <- data.frame(
  Patient_ID = 1:n_patients,
  Age = rnorm(n_patients, mean = 62, sd = 8),
  Sex = sample(c("Male", "Female"), n_patients, replace = TRUE, prob = c(0.55, 0.45)),
  Duration_DM_Years = rpois(n_patients, lambda = 8),
  Baseline_HbA1c = rnorm(n_patients, mean = 8.2, sd = 1.1),
  Baseline_eGFR = rnorm(n_patients, mean = 75, sd = 15)
)

# Confounding assignment: Older, sicker patients are more likely to get DPP4i
propensity_true <- plogis(-2 + (clinical_cohort$Age * 0.05) - (clinical_cohort$Baseline_eGFR * 0.02))
clinical_cohort$Treatment <- ifelse(runif(n_patients) < propensity_true, "DPP4i", "SGLT2i")
clinical_cohort$Treatment <- as.factor(clinical_cohort$Treatment)

# Simulating Survival Outcomes: Time to Major Adverse Cardiovascular Event (MACE)
# SGLT2i reduces hazard; Age increases hazard
clinical_cohort$Hazard <- exp(0.04 * clinical_cohort$Age - 0.02 * clinical_cohort$Baseline_eGFR + 
                              ifelse(clinical_cohort$Treatment == "SGLT2i", -0.4, 0.1))
clinical_cohort$Time_To_MACE <- rexp(n_patients, rate = 0.01 * clinical_cohort$Hazard)
clinical_cohort$MACE_Event <- ifelse(clinical_cohort$Time_To_MACE < 36, 1, 0) # 36 month follow up
clinical_cohort$Time_To_MACE <- pmin(clinical_cohort$Time_To_MACE, 36) # Censoring at 36 months

# -----------------------------------------------------------------------------------------
# SECTION 3: PROPENSITY SCORE MATCHING (PSM)
# -----------------------------------------------------------------------------------------
cat("\n[2] PERFORMING NEAREST-NEIGHBOR PROPENSITY SCORE MATCHING...\n")

# Objective: Balance baseline covariates between Treatment groups to mimic randomization
psm_model <- matchit(Treatment ~ Age + Sex + Duration_DM_Years + Baseline_HbA1c + Baseline_eGFR,
                     data = clinical_cohort, 
                     method = "nearest", 
                     distance = "glm",
                     ratio = 1,          # 1:1 matching
                     caliper = 0.2)      # Strict caliper to ensure close matches

# Extract the perfectly balanced matched dataset
matched_cohort <- match.data(psm_model)

cat("Matching Complete. Original Cohort: ", nrow(clinical_cohort), 
    " | Matched Cohort: ", nrow(matched_cohort), "\n")

# -----------------------------------------------------------------------------------------
# SECTION 4: SURVIVAL ANALYSIS (KAPLAN-MEIER & COX PROPORTIONAL HAZARDS)
# -----------------------------------------------------------------------------------------
cat("\n[3] EXECUTING TIME-TO-EVENT SURVIVAL ANALYSIS (MACE)...\n")

# 4.1 Kaplan-Meier Survival Curve
km_fit <- survfit(Surv(Time_To_MACE, MACE_Event) ~ Treatment, data = matched_cohort)

# Code to plot the Kaplan-Meier curve (output suppressed for console)
km_plot <- ggsurvplot(km_fit, 
                      data = matched_cohort, 
                      pval = TRUE, 
                      risk.table = TRUE, 
                      conf.int = TRUE,
                      palette = c("#E7B800", "#2E9FDF"),
                      title = "Kaplan-Meier Curve: Freedom from MACE (SGLT2i vs. DPP4i)",
                      xlab = "Time in Months",
                      ylab = "MACE-Free Survival Probability")

# 4.2 Multivariable Cox Proportional Hazards Model
cox_model <- coxph(Surv(Time_To_MACE, MACE_Event) ~ Treatment + Age + Sex + Baseline_eGFR, 
                   data = matched_cohort)

cat("\n--- Cox Proportional Hazards Model Summary ---\n")
print(summary(cox_model))

# -----------------------------------------------------------------------------------------
# SECTION 5: LINEAR MIXED-EFFECTS MODELING (LONGITUDINAL DATA)
# -----------------------------------------------------------------------------------------
cat("\n[4] EXECUTING LONGITUDINAL LINEAR MIXED-EFFECTS MODEL (LMM)...\n")
# Objective: Analyze repeated HbA1c measures over 24 months accounting for intra-patient correlation

# Simulating longitudinal data (Visits at Month 0, 6, 12, 18, 24)
timepoints <- c(0, 6, 12, 18, 24)
longitudinal_data <- expand_grid(Patient_ID = matched_cohort$Patient_ID, Time_Month = timepoints) %>%
  left_join(matched_cohort %>% select(Patient_ID, Treatment, Baseline_HbA1c), by = "Patient_ID")

# SGLT2i shows a steeper initial drop in HbA1c that sustains over time
longitudinal_data <- longitudinal_data %>%
  mutate(
    HbA1c_Drop = ifelse(Treatment == "SGLT2i", 
                        -0.08 * Time_Month,  # SGLT2i trajectory
                        -0.03 * Time_Month), # DPP4i trajectory
    Random_Effect = rnorm(n(), mean = 0, sd = 0.3),
    Current_HbA1c = Baseline_HbA1c + HbA1c_Drop + Random_Effect
  )

# 5.1 Fit the Linear Mixed-Effects Model
# Formula: Current_HbA1c ~ Treatment * Time + (1 | Patient_ID)
# The (1 | Patient_ID) adds a random intercept for each patient
lmm_fit <- lmer(Current_HbA1c ~ Treatment * Time_Month + (1 | Patient_ID), data = longitudinal_data)

cat("\n--- Linear Mixed-Effects Model (LMM) Summary ---\n")
print(summary(lmm_fit))

cat("\n[5] BIOSTATISTICAL PIPELINE EXECUTION COMPLETE.\n")
