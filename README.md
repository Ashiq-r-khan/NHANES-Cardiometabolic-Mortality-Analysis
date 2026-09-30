# Cardiometabolic Risk and Mortality in US Adults

Survey-weighted survival analysis of NHANES 1999-2018 linked to National Death Index mortality, with competing risks, causal inference, a 10-year risk model and a Power BI dashboard.

The cohort is 50,819 US adults who had a full NHANES exam between 1999 and 2018, followed for death up to the end of 2019. 8,279 of them died. The project asks five questions:

1. How common are diabetes, prediabetes, hypertension, obesity and chronic kidney disease (CKD), and did that change over 20 years?
2. How much higher is the risk of death with each condition, after adjusting for age and everything else?
3. What do people with each condition die of?
4. How many deaths can be attributed to each risk factor?
5. Can 10-year risk of death be predicted from one exam?

Python and DuckDB SQL do the download, cleaning and descriptive work. All statistics are in R with the full survey design (weights, strata, PSUs) and 20 multiple imputations, so the estimates describe the 213.7 million US adults the sample represents. The full write-up with code, figures and tables is in [`NHANES Cardiometabolic Mortality Analysis Report.pdf`](NHANES%20Cardiometabolic%20Mortality%20Analysis%20Report.pdf).

## Key results

| Question | Result |
|---|---|
| Cohort | 50,819 adults, 8,279 deaths, 479,971 person-years, median follow-up 9 years (max 20.8) |
| Diabetes prevalence | 11.8% age-standardized (95% CI 11.4 to 12.2), up from 9.0% in 1999-2000 to 14.4% in 2017-2018 |
| Undiagnosed diabetes | 23.4% of all diabetes, no significant change over 20 years |
| Prediabetes | No excess risk of death once age is adjusted for: HR 0.99 (0.92 to 1.06) |
| Diabetes | HR 1.44 (1.34 to 1.55) fully adjusted. Total effect with IPW: HR 1.83 (1.58 to 2.13), 6.3 extra deaths per 100 people over 10 years |
| CKD | HR 1.59 (1.50 to 1.68). Largest share of deaths of any risk factor: 14.7% (13.2 to 16.1) |
| Hypertension | HR 1.22 (1.14 to 1.31). 11.8% of all deaths and 22.9% of CVD deaths attributable |
| Cause of death | Diabetes raised CVD death (HR 1.50) and other causes (1.62), not cancer (1.03). Diabetes was on the death certificate of only 34.5% of people who died with diagnosed diabetes |
| Competing risks | 1 minus Kaplan-Meier overstated cause-specific risk by up to 28.5% compared with Aalen-Johansen |
| 10-year risk model | C-index 0.872 on 2007-2010 (trained on 1999-2006), observed 10.4% vs predicted 10.2%, calibration slope 0.98 |

HR = hazard ratio from survey-weighted Cox models, pooled over 20 imputations, 95% CI in brackets.

## Dashboard

Eight report pages, three drill-through pages and tooltip pages, built on 42 tables written by the last R script. The dashboard only looks up and formats the R results, so its numbers match the report exactly. Open [`NHANES Analysis Dashboard.pbix`](NHANES%20Analysis%20Dashboard.pbix) in Power BI Desktop, or view [`NHANES Analysis Dashboard.pdf`](NHANES%20Analysis%20Dashboard.pdf).

| | |
|---|---|
| ![Executive summary](images/dashboard_01_executive_summary.png) | ![Prevalence and trends](images/dashboard_02_prevalence_trends.png) |
| Executive summary | Prevalence and trends |
| ![Who dies](images/dashboard_03_who_dies.png) | ![Survival and hazard ratios](images/dashboard_04_survival_hazard_ratios.png) |
| Who dies: crude vs age-standardized | Survival and hazard ratios |
| ![Competing risks](images/dashboard_05_competing_risks.png) | ![Causal and burden](images/dashboard_06_causal_burden.png) |
| Competing risks and causes of death | Causal effects and attributable burden |
| ![Model performance](images/dashboard_07_model_performance.png) | ![Risk and what-if](images/dashboard_08_risk_whatif.png) |
| Model performance | Risk tiers and what-if |

## Report

50 pages: executive summary, then one section per notebook or script with the code, figures and result tables, then discussion, limitations and references. Every number in the text comes from the output CSVs.

![Report preview](images/report_preview.png)

## Workflow

| Step | Tool | File | What it does |
|---|---|---|---|
| 01 | Python | [`01_download_nhanes_mortality.ipynb`](notebooks/01_download_nhanes_mortality.ipynb) | Downloads 14 NHANES components for 10 cycles and the NCHS Linked Mortality Files |
| 02 | Python | [`02_data_audit.ipynb`](notebooks/02_data_audit.ipynb) | Harmonizes variable names across cycles, derives conditions, eGFR (CKD-EPI 2021), cohort rules, combined weights |
| 03 | Python, SQL | [`03_duckdb_sql.ipynb`](notebooks/03_duckdb_sql.ipynb) | 11 descriptive questions in DuckDB SQL, including a KDIGO grid and a life table built with window functions |
| 04 | Python | [`04_eda.ipynb`](notebooks/04_eda.ipynb) | Skew, correlation, missing data vs outcome, risk shapes, crude vs age-standardized rates |
| 05 | R | [`05_survey_design_imputation.R`](r/05_survey_design_imputation.R) | Survey design, multiple imputation (mice, m = 20), weighted prevalence, trend tests, Table 1 |
| 06 | R | [`06_survival_cox.R`](r/06_survival_cox.R) | Weighted death rates, Kaplan-Meier, Cox models M1 to M3, splines, proportional hazards, 6 sensitivity analyses |
| 07 | R | [`07_competing_risks.R`](r/07_competing_risks.R) | Cause-of-death mix, Aalen-Johansen cumulative incidence, cause-specific and Fine-Gray models |
| 08 | R | [`08_causal_attributable.R`](r/08_causal_attributable.R) | IPW marginal hazard ratios, covariate balance, E-values, population attributable fractions |
| 09 | R | [`09_risk_prediction_powerbi.R`](r/09_risk_prediction_powerbi.R) | Three 10-year risk models with temporal validation, calibration, decision curves, risk tiers, Power BI tables |

## Methods

**Cohort.** Adults 20+, examined at the mobile exam centre, eligible for mortality linkage, not pregnant. Follow-up starts at the exam. Diabetes uses ADA lab criteria plus self-report, hypertension uses 140/90 or medication, CKD is eGFR below 60 or ACR 30+ mg/g. The 10 cycles are combined with the CDC weight rule (4-year weight x 2/10 for 1999-2002, 2-year weight x 1/10 after).

**SQL.** DuckDB directly on Parquet. CTEs, window functions, `RANK()`, `CONCAT_WS` profiles, and an actuarial life table with running sums. All crude and unweighted, used to check the data before modelling.

**EDA.** Age standardization showed that the crude prediabetes excess was age, and that crude rates hid the risk of current smoking. People missing BMI, ACR or waist died at 3 to 4 times the rate of people with values, so the R work uses multiple imputation instead of complete cases. U and J shaped risk curves led to splines.

**Survey design and imputation.** `survey` package with PSU, strata and weights (`nest = TRUE`). All 301 PSUs are in the cohort, so the design gives the same standard errors as a full-sample domain analysis. `mice` with 20 imputations, the death indicator and Nelson-Aalen hazard in the imputation model, pooled with Rubin's rules. Prevalence is age-standardized to the 2000 US standard population.

**Survival.** Survey-weighted Cox models (`svycoxph`) with age as a natural spline, in three adjustment steps. Dose-response splines for HbA1c, BMI, SBP, eGFR, HDL and ACR. Schoenfeld test plus time-split models, and sensitivity analyses: complete cases, excluding early deaths, excluding baseline CVD or cancer, without survey cycle.

**Competing risks.** Weighted Aalen-Johansen cumulative incidence, cause-specific Cox models, and Fine-Gray subdistribution models for CVD and cancer death.

**Causal.** Separate confounder sets per exposure, leaving mediators out. Stabilized IPW times survey weight, balance checked with standardized mean differences, truncation sensitivity, E-values, and population attributable fractions with Miettinen's formula.

**Prediction.** Temporal split: develop on 1999-2006, validate on 2007-2010. Three weighted Cox models (basic, conditions, full with lab splines), pooled over imputations. C-index, IPCW Brier score, O/E ratio, calibration slope, decision curves.

## Project structure

```
NHANES Cardiometabolic Mortality Project/
├── data/
│   ├── raw/                download_manifest.csv (raw .xpt and .dat files not in the repo)
│   ├── parquet/            one Parquet file per NHANES component, plus mortality
│   └── processed/          analytic_cohort.parquet (50,819 rows x 58 columns)
├── notebooks/              01 to 04, Python (Google Colab)
├── r/                      05 to 09, R scripts (RStudio)
├── outputs/
│   ├── sql_results/        q01 to q11
│   ├── tables/             02_* and 04_* tables
│   ├── figures/            04_* figures
│   ├── R/                  one folder per R script with its CSVs, figures and session info
│   └── powerbi/            42 tables, _relationships.csv, _table_dictionary.csv
├── images/                 README screenshots
├── NHANES Cardiometabolic Mortality Analysis Report.pdf
├── NHANES Analysis Dashboard.pbix
├── NHANES Analysis Dashboard.pdf
└── requirements.txt
```

## How to run

**What is not in the repo.** The raw CDC files (`data/raw/nhanes`, `data/raw/mortality`, about 260 MB) and the R `.rds` files (imputed data and model caches, about 65 MB) are left out. The Parquet files and the analytic cohort are included, so notebooks 02 to 04 and all R scripts run without downloading anything.

**Python (notebooks 01 to 04).**

1. In Google Colab: copy the folder to `MyDrive/NHANES Cardiometabolic Mortality Project` and run the notebooks in order. The first cell mounts Drive.
2. Locally: `pip install -r requirements.txt`, then run the notebooks from inside `notebooks/`. The first cell falls back to the parent folder when it is not in Colab.
3. Notebook 01 needs internet access to `wwwn.cdc.gov` and `ftp.cdc.gov`. Skip it unless you want to rebuild `data/raw` and `data/parquet`.

**R (scripts 05 to 09).**

1. R 4.5 with `tidyverse`, `arrow`, `survey`, `mitools`, `mice` and `survival` installed.
2. Change `root` at the top of each script to your project path.
3. Run 05 to 09 in order. Each script reads the outputs of the ones before it and saves into its own folder under `outputs/R/`.
4. The first run takes time: the imputation in 05 (10 to 30 minutes), about 440 Cox fits in 06 (10 to 20 minutes) and the Fine-Gray models in 07. Slow steps are cached in `.rds` files, so a rerun picks up where it stopped. The imputation seed is 2026.

**Power BI.** The data is embedded in the `.pbix`. To refresh from the CSVs, change the `DataFolder` parameter to your `outputs/powerbi` path.

## Data

- Centers for Disease Control and Prevention, National Center for Health Statistics. [National Health and Nutrition Examination Survey](https://wwwn.cdc.gov/nchs/nhanes/), 1999-2018.
- National Center for Health Statistics. [Public-use Linked Mortality Files](https://www.cdc.gov/nchs/data-linkage/mortality-public.htm), follow-up to 31 December 2019.

Both are public-use files, freely available from CDC/NCHS.

## Limitations

- Every exposure is measured once, at the exam. Changes during follow-up are not seen.
- Diabetes rests on a single HbA1c or fasting glucose value, not a repeat test.
- The HbA1c lab method changed in 2007-2008, which explains part of the rise in prediabetes. Survey cycle is in every model.
- Prior CVD and cancer are self-reported.
- NCHS perturbs the date and cause of death for some records in the public-use mortality files.
- Physical activity, diet and access to care are not measured. E-values show how strong a hidden confounder would need to be.
- Hazard ratios are associations from observational data. The IPW and attributable fractions assume no unmeasured confounding.
- The risk model is validated on 2007-2010 participants, not on a separate cohort.

## Author

**Md. Ashiqur Rahman Khan**, B.Sc. in Statistics, Mawlana Bhashani Science and Technology University

[LinkedIn](https://www.linkedin.com/in/md-ashiqur-rahman-khan-b475b1316/) · ashiqurrahmankhan04@gmail.com
