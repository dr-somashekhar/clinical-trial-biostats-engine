# =========================================================================================
# ADVANCED CLINICAL BIOSTATISTICS & EPIDEMIOLOGY ENGINE
# Focus: Multiple Imputation, Propensity Score Matching (PSM), Survival Analysis, and
#        Linear Mixed Models (LMM)
#
# Description: Simulates a retrospective cohort comparing SGLT2 inhibitors vs. DPP4
# inhibitors in Type 2 Diabetes Mellitus (T2DM) and analyses it end to end.
# All data are synthetic; no real patient information is used.
#
# Usage:  Rscript survival_psm_lmm_analysis.R      (from the repository folder)
#   or:   source("survival_psm_lmm_analysis.R")    (from an R / RStudio session)
# Outputs are written to the ./output folder.
# =========================================================================================

# -----------------------------------------------------------------------------------------
# SECTION 1: LIBRARY INITIALIZATION
# -----------------------------------------------------------------------------------------
required_packages <- c("tidyverse", "survival", "survminer", "MatchIt", "mice",
                       "lme4", "lmerTest", "broom.mixed", "gtsummary", "gt", "smd")

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)    # Data manipulation and visualization (ggplot2, dplyr)
  library(survival)     # Core survival analysis functions
  library(survminer)    # Kaplan-Meier plotting
  library(MatchIt)      # Propensity score matching to reduce selection bias
  library(lmerTest)     # lme4 mixed models with p-values (masks lme4::lmer)
  library(broom.mixed)  # Tidy outputs for mixed models
  library(gtsummary)    # Publication-ready summary tables
})
# mice is called as mice:: to avoid masking dplyr/stats functions.

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------------------
# SECTION 2: LARGE-SCALE COHORT SIMULATION (N = 5000)
# -----------------------------------------------------------------------------------------
set.seed(2026) # Seed for reproducibility

cat("\n[1] INITIALIZING CLINICAL COHORT SIMULATION...\n")
n_patients <- 5000
follow_up_months <- 36

clinical_cohort <- data.frame(
  Patient_ID = 1:n_patients,
  Age = rnorm(n_patients, mean = 62, sd = 8),
  Sex = sample(c("Male", "Female"), n_patients, replace = TRUE, prob = c(0.55, 0.45)),
  Duration_DM_Years = rpois(n_patients, lambda = 8),
  Baseline_HbA1c = rnorm(n_patients, mean = 8.2, sd = 1.1),
  Baseline_eGFR = rnorm(n_patients, mean = 75, sd = 15),
  Baseline_BMI = rnorm(n_patients, mean = 30, sd = 5)
)

# Confounding by indication: older patients and those with lower eGFR are more likely
# to receive DPP4i.
propensity_true <- plogis(-2 + 0.05 * clinical_cohort$Age - 0.02 * clinical_cohort$Baseline_eGFR)
clinical_cohort$Treatment <- factor(
  ifelse(runif(n_patients) < propensity_true, "DPP4i", "SGLT2i"),
  levels = c("DPP4i", "SGLT2i")   # DPP4i = reference (comparator), SGLT2i = treated
)

# Time to Major Adverse Cardiovascular Event (MACE): SGLT2i lowers the hazard, age raises it.
clinical_cohort$Hazard <- exp(0.04 * clinical_cohort$Age - 0.02 * clinical_cohort$Baseline_eGFR +
                              ifelse(clinical_cohort$Treatment == "SGLT2i", -0.4, 0.1))
latent_time <- rexp(n_patients, rate = 0.01 * clinical_cohort$Hazard)
clinical_cohort$MACE_Event <- as.integer(latent_time < follow_up_months)
clinical_cohort$Time_To_MACE <- pmin(latent_time, follow_up_months) # Administrative censoring
clinical_cohort$Hazard <- NULL  # Latent simulation parameter, not an observable variable

# Missing At Random (MAR) mechanism: older patients are more likely to lack a BMI record.
p_missing_bmi <- plogis(-2.5 + 0.08 * (clinical_cohort$Age - 62))
bmi_missing <- runif(n_patients) < p_missing_bmi
clinical_cohort$Baseline_BMI[bmi_missing] <- NA
cat(sprintf("Baseline_BMI missing for %d of %d patients (%.1f%%).\n",
            sum(bmi_missing), n_patients, 100 * mean(bmi_missing)))

# -----------------------------------------------------------------------------------------
# SECTION 3: MULTIPLE IMPUTATION (MICE, PREDICTIVE MEAN MATCHING)
# -----------------------------------------------------------------------------------------
cat("\n[2] IMPUTING MISSING DATA (MICE / PMM, m = 5)...\n")
n_imputations <- 5

# Patient_ID is excluded as a predictor. The outcome (event indicator and Nelson-Aalen
# cumulative hazard) is included in the imputation model, as recommended for survival analyses.
imputation_data <- clinical_cohort
imputation_data$Cum_Hazard <- mice::nelsonaalen(imputation_data, Time_To_MACE, MACE_Event)
imputation_data$Sex <- factor(imputation_data$Sex)

predictor_matrix <- mice::make.predictorMatrix(imputation_data)
predictor_matrix[, "Patient_ID"] <- 0
predictor_matrix[, "Time_To_MACE"] <- 0   # Cum_Hazard carries the time information

imputed <- mice::mice(imputation_data, m = n_imputations, method = "pmm",
                      predictorMatrix = predictor_matrix, seed = 2026, printFlag = FALSE)

# -----------------------------------------------------------------------------------------
# SECTION 4: PROPENSITY SCORE MATCHING (PSM) - WITHIN EACH IMPUTED DATASET
# -----------------------------------------------------------------------------------------
cat("\n[3] PERFORMING NEAREST-NEIGHBOR PROPENSITY SCORE MATCHING...\n")

# Objective: balance baseline covariates between treatment groups to mimic randomization.
ps_formula <- Treatment ~ Age + Sex + Duration_DM_Years + Baseline_HbA1c +
  Baseline_eGFR + Baseline_BMI

match_one <- function(data) {
  matchit(ps_formula, data = data, method = "nearest", distance = "glm",
          ratio = 1,          # 1:1 matching
          caliper = 0.1)      # Caliper = 0.1 SD of the propensity score (logit scale); tight enough for |SMD| < 0.1
}

psm_models <- lapply(seq_len(n_imputations), function(i) match_one(mice::complete(imputed, i)))
matched_sets <- lapply(psm_models, match.data)

# The first imputed dataset is used for descriptive output (balance table, KM curve, LMM).
psm_model <- psm_models[[1]]
matched_cohort <- matched_sets[[1]]

cat("Matching complete (imputation 1). Original cohort:", n_patients,
    "| Matched cohort:", nrow(matched_cohort), "\n")

cat("\n--- Covariate balance (standardized mean differences; |SMD| < 0.1 is good) ---\n")
balance <- summary(psm_model)
print(round(balance$sum.matched[, c("Means Treated", "Means Control", "Std. Mean Diff.")], 3))

png(file.path(output_dir, "psm_balance_plot.png"), width = 900, height = 700)
plot(balance, var.order = "unmatched")
dev.off()

baseline_table <- matched_cohort %>%
  select(Treatment, Age, Sex, Duration_DM_Years, Baseline_HbA1c, Baseline_eGFR, Baseline_BMI) %>%
  tbl_summary(by = Treatment) %>%
  add_difference()
tryCatch(
  gt::gtsave(as_gt(baseline_table), file.path(output_dir, "baseline_table_matched.html")),
  error = function(e) message("Could not save baseline table: ", conditionMessage(e))
)

# -----------------------------------------------------------------------------------------
# SECTION 5: SURVIVAL ANALYSIS (KAPLAN-MEIER & COX PROPORTIONAL HAZARDS)
# -----------------------------------------------------------------------------------------
cat("\n[4] EXECUTING TIME-TO-EVENT SURVIVAL ANALYSIS (MACE)...\n")

# 5.1 Kaplan-Meier survival curve (imputation 1 matched cohort)
km_fit <- survfit(Surv(Time_To_MACE, MACE_Event) ~ Treatment, data = matched_cohort)

km_plot <- ggsurvplot(km_fit,
                      data = matched_cohort,
                      pval = TRUE,
                      risk.table = TRUE,
                      conf.int = TRUE,
                      palette = c("#E7B800", "#2E9FDF"),   # DPP4i, SGLT2i
                      legend.title = "Treatment", legend.labs = c("DPP4i", "SGLT2i"),
                      title = "Kaplan-Meier Curve: Freedom from MACE (SGLT2i vs. DPP4i)",
                      xlab = "Time in Months",
                      ylab = "MACE-Free Survival Probability")
png(file.path(output_dir, "kaplan_meier_mace.png"), width = 1000, height = 800)
print(km_plot)
dev.off()

# 5.2 Multivariable Cox model, fitted in each matched imputation. Robust standard errors
# clustered on matched pair (subclass) account for the matched design.
fit_cox <- function(data) {
  coxph(Surv(Time_To_MACE, MACE_Event) ~ Treatment + Age + Sex + Baseline_eGFR,
        data = data, cluster = subclass)
}
cox_models <- lapply(matched_sets, fit_cox)

cat("\n--- Cox Proportional Hazards Model Summary (imputation 1) ---\n")
print(summary(cox_models[[1]]))

cat("\n--- Proportional hazards assumption (Schoenfeld residuals; p < 0.05 = violation) ---\n")
print(cox.zph(cox_models[[1]]))

# 5.3 Pool the treatment effect across imputations with Rubin's rules
treatment_term <- "TreatmentSGLT2i"
q_hat <- sapply(cox_models, function(m) coef(m)[treatment_term])
u_hat <- sapply(cox_models, function(m) vcov(m)[treatment_term, treatment_term])
pooled <- mice::pool.scalar(q_hat, u_hat, n = nrow(matched_cohort))
pooled_se <- sqrt(pooled$t)
cat("\n--- Pooled SGLT2i vs DPP4i effect on MACE (Rubin's rules, m =", n_imputations, ") ---\n")
cat(sprintf("Hazard ratio: %.3f (95%% CI %.3f - %.3f)\n",
            exp(pooled$qbar),
            exp(pooled$qbar - qnorm(0.975) * pooled_se),
            exp(pooled$qbar + qnorm(0.975) * pooled_se)))

# -----------------------------------------------------------------------------------------
# SECTION 6: LINEAR MIXED-EFFECTS MODELING (LONGITUDINAL DATA)
# -----------------------------------------------------------------------------------------
cat("\n[5] EXECUTING LONGITUDINAL LINEAR MIXED-EFFECTS MODEL (LMM)...\n")
# Objective: analyze repeated HbA1c measures over 24 months accounting for intra-patient
# correlation.

timepoints <- c(0, 6, 12, 18, 24)  # Visits at months 0, 6, 12, 18, 24
patient_baseline <- matched_cohort %>%
  select(Patient_ID, Treatment, Baseline_HbA1c) %>%
  mutate(Patient_Intercept = rnorm(n(), mean = 0, sd = 0.5))  # one true random intercept per patient

longitudinal_data <- expand_grid(Patient_ID = patient_baseline$Patient_ID, Time_Month = timepoints) %>%
  left_join(patient_baseline, by = "Patient_ID") %>%
  mutate(
    # SGLT2i shows a steeper HbA1c decline that is sustained over time
    HbA1c_Drop = ifelse(Treatment == "SGLT2i", -0.08 * Time_Month, -0.03 * Time_Month),
    Residual_Error = rnorm(n(), mean = 0, sd = 0.3),  # visit-level measurement noise
    Current_HbA1c = Baseline_HbA1c + Patient_Intercept + HbA1c_Drop + Residual_Error
  )

# Formula: Current_HbA1c ~ Treatment * Time + (1 | Patient_ID)
# (1 | Patient_ID) adds a random intercept for each patient.
lmm_fit <- lmer(Current_HbA1c ~ Treatment * Time_Month + (1 | Patient_ID),
                data = longitudinal_data)

cat("\n--- Linear Mixed-Effects Model (LMM) Summary ---\n")
print(summary(lmm_fit))

lmm_tidy <- broom.mixed::tidy(lmm_fit, effects = "fixed", conf.int = TRUE)
write.csv(lmm_tidy, file.path(output_dir, "lmm_fixed_effects.csv"), row.names = FALSE)

cat("\n[6] BIOSTATISTICAL PIPELINE EXECUTION COMPLETE. Results saved in ./", output_dir, "\n", sep = "")
