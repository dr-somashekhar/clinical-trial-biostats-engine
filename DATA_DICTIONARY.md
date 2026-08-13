#  Data Dictionary & Clinical Codebook

This document describes the simulated variables and biological features utilized within the modeling pipelines in this repository.

## 1. Patient Demographics & Anthropometrics
| Variable Name | Data Type | Description | Unit / Coding |
| :--- | :--- | :--- | :--- |
| `Patient_ID` | Integer | Unique cryptographic identifier for each patient | `1` to `N` |
| `Age` | Numeric | Age of the patient at baseline enrollment | Years |
| `Sex` | Categorical | Biological sex assigned at birth | `Male`, `Female` |
| `Baseline_BMI` | Numeric | Body Mass Index at baseline | kg/m^2 |

## 2. Clinical Biomarkers & Laboratory Values
| Variable Name | Data Type | Description | Unit / Coding |
| :--- | :--- | :--- | :--- |
| `Baseline_HbA1c` | Numeric | Glycated hemoglobin level at enrollment | Percentage (%) |
| `Baseline_eGFR` | Numeric | Estimated Glomerular Filtration Rate | mL/min/1.73m^2 |
| `ALT` | Numeric | Alanine Aminotransferase | U/L |
| `AST` | Numeric | Aspartate Aminotransferase | U/L |

## 3. Pharmacotherapy & Interventions
| Variable Name | Data Type | Description | Unit / Coding |
| :--- | :--- | :--- | :--- |
| `Treatment` | Categorical | Assigned antidiabetic pharmacotherapy | `SGLT2i`, `DPP4i`, `Vehicle` |
| `Duration_DM_Years`| Numeric | Years since initial Type 2 Diabetes diagnosis | Years |

## 4. Outcomes & Toxicological Metrics
| Variable Name | Data Type | Description | Unit / Coding |
| :--- | :--- | :--- | :--- |
| `Time_To_MACE` | Numeric | Time from enrollment to Major Adverse Cardiovascular Event | Months |
| `MACE_Event` | Binary | Occurrence of MACE during the follow-up period | `0` = Censored, `1` = Event |
| `Toxicity_Score` | Numeric | Continuous algorithmic severity score for Drug-Induced Liver Injury | Scale `0` - `100` |

*Note: All datasets are synthetically generated for pipeline validation and contain no real patient Protected Health Information (PHI).*
