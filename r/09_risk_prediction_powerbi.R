# ============================================================================
# Script 09: 10-year mortality risk prediction and Power BI export
# Part A: develop on cycles 1999-2006, validate on 2007-2010 (temporal
#         validation), compare three models, calibration, decision curves,
#         final model refit on all cycles, risk tiers and what-if profiles
# Part B: write the Power BI data model (dim_, fact_, stat_, model_ CSVs + relationships)
# Input : outputs/R/05-08 result folders
# Output: outputs/R/09_risk_prediction/ and outputs/powerbi/
# ============================================================================

library(tidyverse)
library(survival)
library(splines)

options(dplyr.summarise.inform = FALSE, scipen = 999)

# Set paths
root <- "E:/NHANES Cardiometabolic Mortality Project"
rdir <- file.path(root, "outputs", "R")
out <- file.path(rdir, "09_risk_prediction")
pbi <- file.path(root, "outputs", "powerbi")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(pbi, recursive = TRUE, showWarnings = FALSE)

show_save <- function(x, name) {
  write_csv(x, file.path(out, paste0(name, ".csv")), na = "")
  print(x, n = Inf, width = Inf)
  invisible(x)
}
plot_save <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "#fcfcfb")
  print(p)
}
header <- function(x) cat("\n============================================================\n", x, "\n============================================================\n", sep = "")
wmean <- function(x, w) sum(w * x) / sum(w)

pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
                  plot.title = element_text(face = "bold", colour = "#0b0b0b"), plot.subtitle = element_text(colour = "#52514e"),
                  axis.text = element_text(colour = "#52514e"), strip.text = element_text(face = "bold", hjust = 0, colour = "#0b0b0b"),
                  legend.position = "top", legend.justification = "left", plot.background = element_rect(fill = "#fcfcfb", colour = NA)))

# ============================================================================
# PART A: RISK PREDICTION
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Data and temporal split
# Development: 1999-2006 (4 cycles, 13-21 years of follow-up)
# Validation : 2007-2010 (2 cycles, 9-13 years of follow-up)
# Survey cycle is never a predictor: a model used today cannot know it.
# Survey weights are used so performance refers to the US adult population.
# ----------------------------------------------------------------------------
header("1. DATA AND TEMPORAL SPLIT")

horizon <- 10
long <- readRDS(file.path(rdir, "05_survey_design_imputation", "05_imputed_long.rds")) %>%
  mutate(time_y = time_months / 12, w = wt_mec / mean(wt_mec), htn_med0 = coalesce(htn_med, 0),
         glycemic4 = factor(case_when(glycemic == "diabetes" & dm_dx == 1 ~ "diagnosed diabetes", glycemic == "diabetes" ~ "undiagnosed diabetes", TRUE ~ as.character(glycemic)),
                            c("normal", "prediabetes", "undiagnosed diabetes", "diagnosed diabetes")),
         bmi_cat = factor(bmi_cat, c("normal", "underweight", "overweight", "obese")),
         sample = if_else(cycle_num <= 4, "development", if_else(cycle_num <= 6, "validation", "later cycles")))
imps <- split(long, long$.imp)
m <- length(imps)

split_summary <- imps[[1]] %>%
  group_by(sample) %>%
  summarise(cycles = paste(range(as.character(cycle)), collapse = " to "), n = n(), deaths = sum(death),
            deaths_within_10y = sum(death == 1 & time_y <= horizon), min_followup_y = round(min(time_y[death == 0]), 1))
show_save(split_summary, "09_split_summary")

# ----------------------------------------------------------------------------
# 2. Three candidate models (Cox, survey-weighted, fitted in every imputation;
# coefficients and baseline hazard averaged across imputations)
# A Basic     : age, sex, race/ethnicity, smoking
# B Conditions: A + glycemic status (with diagnosis), hypertension, CKD, BMI group, prior CVD, cancer
# C Full      : A + education, income, continuous biomarkers as splines, treatment flags, prior CVD, cancer
# ----------------------------------------------------------------------------
header("2. MODEL DEVELOPMENT")

dev1 <- filter(imps[[1]], sample == "development")
spline_term <- function(v, d) {
  x <- d[[v]]
  sprintf("ns(%s, knots = c(%s), Boundary.knots = c(%s))", v, paste(round(quantile(x, c(0.35, 0.65)), 3), collapse = ", "),
          paste(round(quantile(x, c(0.05, 0.95)), 3), collapse = ", "))
}
age_term <- "ns(age, knots = c(35, 50, 65), Boundary.knots = c(20, 85))"
models <- list(
  `A: Basic` = c(age_term, "sex", "race_eth", "smoking"),
  `B: Conditions` = c(age_term, "sex", "race_eth", "smoking", "glycemic4", "hypertension", "ckd", "bmi_cat", "cvd_history", "cancer"),
  `C: Full` = c(age_term, "sex", "race_eth", "smoking", "education", "pir", map_chr(c("bmi", "sbp", "hba1c", "egfr", "log_acr", "hdl"), spline_term, d = dev1),
                "tc", "htn_med0", "dm_dx", "cvd_history", "cancer"))

fit_pooled <- function(terms, data_list) {
  f <- reformulate(terms, "Surv(time_y, death)")
  fits <- map(data_list, function(d) {
    fit <- coxph(f, data = d, weights = w, model = TRUE)
    bh <- basehaz(fit, centered = FALSE)
    list(beta = coef(fit), H0 = max(c(0, bh$hazard[bh$time <= horizon])), tt = delete.response(terms(fit)), xlevels = fit$xlevels)
  })
  list(terms = terms, tt = fits[[1]]$tt, xlevels = fits[[1]]$xlevels,
       beta = rowMeans(sapply(fits, function(x) x$beta)), H0 = mean(map_dbl(fits, "H0")))
}
predict_risk <- function(model, d) {
  X <- model.matrix(model$tt, model.frame(model$tt, d, xlev = model$xlevels))[, -1, drop = FALSE]
  lp <- drop(X[, names(model$beta)] %*% model$beta)
  tibble(lp = lp, risk = 1 - exp(-model$H0 * exp(lp)))
}

dev_models <- map(models, function(terms) fit_pooled(terms, map(imps, function(d) filter(d, sample == "development"))))
walk2(names(dev_models), dev_models, function(n, mod) cat(sprintf("%-14s %2d coefficients, baseline 10-year cumulative hazard %.4f\n", n, length(mod$beta), mod$H0)))

# ----------------------------------------------------------------------------
# 3. Temporal validation (2007-2010), metrics averaged over imputations
# C-index      : ranking; restricted to the first 10 years
# Brier        : 10-year prediction error with censoring weights (IPCW)
# Scaled Brier : improvement over predicting the average risk for everyone
# O/E          : observed (weighted KM) / mean predicted 10-year risk
# Cal. slope   : 1 = perfect spread of predictions; < 1 = predictions too extreme
# ----------------------------------------------------------------------------
header("3. TEMPORAL VALIDATION (2007-2010)")

km_risk <- function(time, status, w, t = horizon) {
  s <- summary(survfit(Surv(time, status) ~ 1, weights = w), times = t, extend = TRUE)
  1 - s$surv
}
brier_ipcw <- function(d, risk, t = horizon) {
  cens <- survfit(Surv(d$time_y, 1 - d$death) ~ 1, weights = d$w)
  G <- stepfun(cens$time, c(1, cens$surv))
  event <- d$time_y <= t & d$death == 1
  wt <- case_when(event ~ 1 / G(pmax(d$time_y - 1e-6, 0)), d$time_y > t ~ 1 / G(t), TRUE ~ 0)
  wmean(wt * ((as.numeric(event) - risk)^2), d$w)
}
metrics <- function(d, pr) {
  obs <- km_risk(d$time_y, d$death, d$w)
  bs <- brier_ipcw(d, pr$risk)
  bs0 <- brier_ipcw(d, rep(obs, nrow(d)))
  tibble(c_index = concordance(Surv(time_y, death) ~ pr$risk, data = d, weights = w, reverse = TRUE, ymax = horizon)$concordance,
         brier = bs, scaled_brier = 1 - bs / bs0, observed_risk = obs, mean_predicted_risk = wmean(pr$risk, d$w),
         oe_ratio = obs / wmean(pr$risk, d$w),
         calibration_slope = unname(coef(coxph(Surv(pmin(time_y, horizon), death * (time_y <= horizon)) ~ pr$lp, data = d, weights = w))))
}

val <- map(imps, function(d) filter(d, sample == "validation"))
performance <- imap_dfr(dev_models, function(mod, name) map_dfr(val, function(d) metrics(d, predict_risk(mod, d))) %>%
                          summarise(across(everything(), list(mean = mean, min = min, max = max))) %>% mutate(model = name)) %>%
  select(model, everything())
performance_tbl <- performance %>%
  transmute(model, c_index = round(c_index_mean, 3), c_index_range = sprintf("%.3f-%.3f", c_index_min, c_index_max),
            brier = round(brier_mean, 4), scaled_brier = round(scaled_brier_mean, 3), observed_risk_pct = round(100 * observed_risk_mean, 1),
            mean_predicted_risk_pct = round(100 * mean_predicted_risk_mean, 1), oe_ratio = round(oe_ratio_mean, 3), calibration_slope = round(calibration_slope_mean, 3))
show_save(performance_tbl, "09_validation_performance")

calibration <- imap_dfr(dev_models, function(mod, name) map_dfr(val, function(d) {
  pr <- predict_risk(mod, d)
  d %>% mutate(risk = pr$risk, decile = ntile(risk, 10)) %>%
    group_by(decile) %>%
    summarise(predicted = wmean(risk, w), observed = km_risk(time_y, death, w), n = n())
}, .id = "imp") %>% group_by(decile) %>% summarise(across(c(predicted, observed, n), mean)) %>% mutate(model = name))
write_csv(calibration, file.path(out, "09_calibration_deciles.csv"))

p_cal <- ggplot(calibration, aes(100 * predicted, 100 * observed, colour = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "#52514e", linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = pal[1:3]) +
  coord_equal() +
  labs(title = "Calibration in the 2007-2010 validation cohort", subtitle = "Predicted vs observed (weighted Kaplan-Meier) risk by decile of predicted risk",
       x = "Predicted 10-year risk of death (%)", y = "Observed 10-year risk of death (%)", colour = NULL)
plot_save(p_cal, "09_calibration", w = 8, h = 7.5)

# Decision curve: net benefit = TP/n - FP/n x pt/(1 - pt), with TP and FP from
# weighted Kaplan-Meier among people above the threshold (handles censoring)
thresholds <- seq(0.02, 0.5, by = 0.01)
net_benefit <- function(d, risk, pt) {
  hi <- risk >= pt
  p_hi <- wmean(hi, d$w)
  if (sum(hi) < 20) return(NA_real_)
  r_hi <- km_risk(d$time_y[hi], d$death[hi], d$w[hi])
  r_hi * p_hi - (1 - r_hi) * p_hi * pt / (1 - pt)
}
dca <- map_dfr(val, function(d) {
  all_risk <- km_risk(d$time_y, d$death, d$w)
  bind_rows(imap_dfr(dev_models, function(mod, name) { r <- predict_risk(mod, d)$risk; tibble(model = name, threshold = thresholds, net_benefit = map_dbl(thresholds, net_benefit, d = d, risk = r)) }),
            tibble(model = "Treat all", threshold = thresholds, net_benefit = all_risk - (1 - all_risk) * thresholds / (1 - thresholds)),
            tibble(model = "Treat none", threshold = thresholds, net_benefit = 0))
}, .id = "imp") %>% group_by(model, threshold) %>% summarise(net_benefit = mean(net_benefit, na.rm = TRUE)) %>% ungroup()
write_csv(dca, file.path(out, "09_decision_curve.csv"))

p_dca <- dca %>%
  mutate(model = factor(model, c(names(models), "Treat all", "Treat none"))) %>%
  ggplot(aes(100 * threshold, net_benefit, colour = model, linetype = model)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(pal[1:3], "#52514e", "#9a9994")) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "dashed", "dotted")) +
  coord_cartesian(ylim = c(-0.01, max(dca$net_benefit, na.rm = TRUE) * 1.05)) +
  labs(title = "Decision curve analysis (2007-2010 validation cohort)", subtitle = "Net benefit of using each model to flag people above a 10-year risk threshold",
       x = "Risk threshold (%)", y = "Net benefit", colour = NULL, linetype = NULL)
plot_save(p_dca, "09_decision_curve", h = 5.5)

# ----------------------------------------------------------------------------
# 4. Final models: refit on all 10 cycles (standard after validation)
# ----------------------------------------------------------------------------
header("4. FINAL MODELS (ALL CYCLES) AND HAZARD RATIOS")

final_B <- fit_pooled(models$`B: Conditions`, imps)
final_C <- fit_pooled(models$`C: Full`, imps)
coefficients <- bind_rows(tibble(model = "B: Conditions", term = names(final_B$beta), beta = final_B$beta),
                          tibble(model = "C: Full", term = names(final_C$beta), beta = final_C$beta)) %>%
  mutate(hr = round(exp(beta), 3), beta = round(beta, 5)) %>%
  bind_rows(tibble(model = c("B: Conditions", "C: Full"), term = "baseline_cumhaz_10y", beta = round(c(final_B$H0, final_C$H0), 6)))
write_csv(coefficients, file.path(out, "09_final_model_coefficients.csv"))
coefficients %>% filter(model == "B: Conditions", !startsWith(term, "ns(")) %>% print(n = Inf, width = Inf)

# ----------------------------------------------------------------------------
# 5. Risk tiers in the US adult population (final model C, imputation-averaged)
# ----------------------------------------------------------------------------
header("5. RISK TIERS (FINAL MODEL C)")

tier_breaks <- c(0, 0.05, 0.10, 0.20, 0.30, 1)
tier_labels <- c("Low (<5%)", "Borderline (5-10%)", "Intermediate (10-20%)", "High (20-30%)", "Very high (30%+)")
risk_all <- map_dfr(imps, function(d) d %>% transmute(SEQN, .imp, risk = predict_risk(final_C, d)$risk)) %>%
  group_by(SEQN) %>% summarise(risk_10y = mean(risk)) %>% ungroup()
base1 <- imps[[1]] %>% left_join(risk_all, by = "SEQN") %>%
  mutate(risk_tier = cut(risk_10y, tier_breaks, tier_labels, right = FALSE))

risk_tiers <- base1 %>%
  group_by(risk_tier) %>%
  summarise(n = n(), population_pct = 100 * sum(wt_mec), mean_age = wmean(age, wt_mec), mean_predicted_pct = 100 * wmean(risk_10y, wt_mec),
            observed_10y_pct = 100 * km_risk(time_y, death, w), diabetes_pct = 100 * wmean(diabetes, wt_mec), hypertension_pct = 100 * wmean(hypertension, wt_mec),
            ckd_pct = 100 * wmean(ckd, wt_mec), current_smoker_pct = 100 * wmean(smoking == "current", wt_mec), deaths = sum(death)) %>%
  mutate(population_pct = population_pct / sum(population_pct) * 100) %>%
  mutate(across(where(is.double), function(x) round(x, 1)))
show_save(risk_tiers, "09_risk_tiers")

p_tiers <- risk_tiers %>%
  select(risk_tier, Predicted = mean_predicted_pct, Observed = observed_10y_pct) %>%
  pivot_longer(-risk_tier) %>%
  mutate(name = factor(name, c("Predicted", "Observed"))) %>%
  ggplot(aes(risk_tier, value, fill = name)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", value)), position = position_dodge(width = 0.75), vjust = -0.4, size = 3.3, colour = "#52514e") +
  scale_fill_manual(values = c(Predicted = pal[1], Observed = pal[2])) +
  labs(title = "Predicted vs observed 10-year risk of death by risk tier", subtitle = "Final model C applied to all US adults 20+, NHANES 1999-2018, survey-weighted",
       x = NULL, y = "10-year risk of death (%)", fill = NULL)
plot_save(p_tiers, "09_risk_tiers", h = 5.5)

# ----------------------------------------------------------------------------
# 6. What-if risk profiles for the Power BI calculator (final model B)
# Every combination of the slicer values, other factors at: no prior CVD,
# no cancer. Power BI only looks up the precomputed risk.
# ----------------------------------------------------------------------------
header("6. WHAT-IF RISK PROFILES")

whatif <- expand_grid(age = c(30, 40, 50, 60, 70, 80), sex = levels(long$sex), race_eth = levels(long$race_eth), smoking = levels(long$smoking),
                      glycemic4 = levels(long$glycemic4), hypertension = 0:1, ckd = 0:1, bmi_cat = levels(long$bmi_cat), cvd_history = 0, cancer = 0) %>%
  mutate(across(c(sex, race_eth, smoking, glycemic4, bmi_cat), function(x) factor(x, levels(long[[cur_column()]])))) %>%
  mutate(risk_10y_pct = round(100 * predict_risk(final_B, pick(everything()))$risk, 2))
cat(sprintf("What-if grid: %s profiles\n", format(nrow(whatif), big.mark = ",")))
whatif %>% filter(race_eth == "Non-Hispanic White", smoking == "never", bmi_cat == "normal", sex == "male") %>%
  select(age, glycemic4, hypertension, ckd, risk_10y_pct) %>%
  pivot_wider(names_from = age, values_from = risk_10y_pct, names_prefix = "age_") %>% print(n = Inf, width = Inf)

# ============================================================================
# PART B: POWER BI DATA MODEL (written to outputs/powerbi)
# Star-style model: 7 dimension tables with one key column each, shared by the
# fact table and the result tables. Key columns have the same name in every
# table, so Power BI can auto-detect most relationships. _relationships.csv
# lists every relationship and is checked here for unmatched keys.
# ============================================================================
header("7. POWER BI DATA MODEL")

rd <- function(folder, file) read_csv(file.path(rdir, folder, file), show_col_types = FALSE)
parse_hr <- function(x) { m <- str_match(x, "([0-9.]+) \\(([0-9.]+)-([0-9.]+)\\)"); tibble(hr = as.numeric(m[, 2]), lcl = as.numeric(m[, 3]), ucl = as.numeric(m[, 4])) }
f05 <- "05_survey_design_imputation"; f06 <- "06_survival_cox"; f07 <- "07_competing_risks"; f08 <- "08_causal_attributable"

# ---- Dimensions ----
cycles <- levels(long$cycle)
group_lookup <- c(`Glycemic status` = "glycemic", Hypertension = "hypertension", `Chronic kidney disease` = "ckd", `BMI category` = "bmi",
                  Smoking = "smoking", age3 = "age", sex = "sex", race_eth = "race", education = "education")
dim_group_level <- tribble(
  ~group_key, ~group_label, ~level, ~level_label,
  "glycemic", "Glycemic status", "normal", "Normal glucose",
  "glycemic", "Glycemic status", "prediabetes", "Prediabetes",
  "glycemic", "Glycemic status", "undiagnosed diabetes", "Undiagnosed diabetes",
  "glycemic", "Glycemic status", "diagnosed diabetes", "Diagnosed diabetes",
  "hypertension", "Hypertension", "no hypertension", "No hypertension",
  "hypertension", "Hypertension", "hypertension", "Hypertension",
  "ckd", "Chronic kidney disease", "no CKD", "No CKD",
  "ckd", "Chronic kidney disease", "CKD", "CKD",
  "bmi", "BMI category", "underweight", "Underweight (<18.5)",
  "bmi", "BMI category", "normal", "Normal (18.5-24.9)",
  "bmi", "BMI category", "overweight", "Overweight (25-29.9)",
  "bmi", "BMI category", "obese", "Obese (30+)",
  "smoking", "Smoking", "never", "Never smoker",
  "smoking", "Smoking", "former", "Former smoker",
  "smoking", "Smoking", "current", "Current smoker",
  "age", "Age group", "20-39", "20-39 years",
  "age", "Age group", "40-59", "40-59 years",
  "age", "Age group", "60+", "60+ years",
  "sex", "Sex", "male", "Male",
  "sex", "Sex", "female", "Female",
  "race", "Race/ethnicity", "Non-Hispanic White", "Non-Hispanic White",
  "race", "Race/ethnicity", "Non-Hispanic Black", "Non-Hispanic Black",
  "race", "Race/ethnicity", "Mexican American", "Mexican American",
  "race", "Race/ethnicity", "Other Hispanic", "Other Hispanic",
  "race", "Race/ethnicity", "Other/Multiracial", "Other/Multiracial",
  "education", "Education", "<HS", "Less than high school",
  "education", "Education", "HS/GED", "High school/GED",
  "education", "Education", ">HS", "More than high school") %>%
  mutate(level_key = paste(group_key, level, sep = "|"), group_sort = match(group_key, unique(group_key)), level_sort = row_number()) %>%
  select(level_key, group_key, group_label, level, level_label, group_sort, level_sort)
add_level_key <- function(x) mutate(x, group_key = unname(group_lookup[group]), level_key = paste(group_key, level, sep = "|"))

prev_overall <- rd(f05, "05_prevalence_overall.csv")
dims <- list(
  dim_cycle = tibble(cycle = cycles, cycle_num = seq_along(cycles), start_year = as.integer(substr(cycles, 1, 4)), end_year = as.integer(substr(cycles, 6, 9))),
  dim_condition = prev_overall %>% distinct(condition, condition_label = label) %>%
    mutate(domain = case_when(condition %in% c("diabetes", "dm_dx", "dm_undiagnosed", "prediabetes") ~ "Glycemic",
                              condition %in% c("hypertension", "hypertension_acc") ~ "Blood pressure",
                              condition %in% c("ckd", "albuminuria") ~ "Kidney", TRUE ~ "Other risk factors"), condition_sort = row_number()),
  dim_group_level = dim_group_level,
  dim_cause = tibble(cause = c("All causes", "CVD", "Cancer", "Other causes"),
                     cause_description = c("Any cause of death", "Heart disease and stroke (NCHS codes 1, 5)", "Malignant neoplasms (code 2)", "All other underlying causes"),
                     cause_sort = 0:3),
  dim_risk_tier = tibble(risk_tier = tier_labels, lower_pct = 100 * head(tier_breaks, -1), upper_pct = 100 * tail(tier_breaks, -1), tier_sort = seq_along(tier_labels)),
  dim_prediction_model = tibble(model = c(names(models), "Treat all", "Treat none"),
                                model_description = c("Age, sex, race/ethnicity, smoking", "Basic + glycemic status, hypertension, CKD, BMI group, prior CVD, cancer",
                                                      "Basic + education, income, biomarker splines (BMI, SBP, HbA1c, eGFR, ACR, HDL), cholesterol, treatment, prior CVD, cancer",
                                                      "Reference: flag everyone", "Reference: flag no one"),
                                model_type = c("Model", "Model", "Model", "Reference", "Reference"), model_sort = 1:5),
  dim_exposure = tibble(exposure = c("Glycemic status", "Glycemic status (diagnosis)", "Diabetes", "Hypertension", "Chronic kidney disease", "BMI category"), exposure_sort = 1:6))

# ---- Fact and result tables ----
facts <- list(
  fact_participants = base1 %>% transmute(SEQN, cycle = as.character(cycle), weight = wt_mec, age, age_group = as.character(age3), sex = as.character(sex),
                                          race_eth = as.character(race_eth), education = as.character(education), smoking = as.character(smoking),
                                          level_key = paste("glycemic", glycemic4, sep = "|"), hypertension, ckd, bmi_cat = as.character(bmi_cat),
                                          cvd_history, cancer, followup_years = round(time_y, 2), died = death,
                                          cause_of_death = if_else(death == 1, as.character(cause), "alive"), risk_10y_pct = round(100 * risk_10y, 2), risk_tier = as.character(risk_tier)),
  stat_kpis = tibble(kpi_sort = 1:10,
                     metric = c("Participants", "Deaths", "Person-years", "Median follow-up (years)", "US adults represented (millions)",
                                "Diabetes prevalence, age-standardized (%)", "Undiagnosed share of diabetes (%)", "Diabetes HR, fully adjusted",
                                "Deaths attributable to CKD (%)", "Model C C-index (temporal validation)"),
                     value = c(nrow(base1), sum(base1$death), round(sum(base1$time_y)), round(median(base1$time_y), 1), round(sum(base1$wt_mec) / 1e6, 1),
                               prev_overall %>% filter(condition == "diabetes", estimate == "Age-standardized (MI)") %>% pull(prevalence_pct),
                               prev_overall %>% filter(estimate == "Crude (MI)") %>% { round(100 * .$prevalence_pct[.$condition == "dm_undiagnosed"] / .$prevalence_pct[.$condition == "diabetes"], 1) },
                               rd(f06, "06_hazard_ratios.csv") %>% filter(exposure == "Glycemic status", level == "diabetes", startsWith(model, "M3")) %>% pull(hr),
                               rd(f08, "08_population_attributable_fractions.csv") %>% filter(outcome == "All-cause", exposure == "Chronic kidney disease") %>% pull(paf_pct),
                               performance_tbl %>% filter(model == "C: Full") %>% pull(c_index))),
  stat_design = rd(f05, "05_design_summary.csv"),
  stat_prevalence_overall = prev_overall %>% select(-label),
  stat_prevalence_by_cycle = rd(f05, "05_prevalence_by_cycle.csv") %>% select(-label),
  stat_prevalence_by_group = rd(f05, "05_prevalence_by_group.csv") %>% select(-label) %>% add_level_key(),
  stat_prevalence_trend_tests = rd(f05, "05_trend_tests.csv") %>% select(-label),
  stat_table1 = rd(f05, "05_table1_weighted.csv") %>% mutate(row_sort = row_number()),
  stat_death_rates = rd(f06, "06_death_rates.csv") %>% add_level_key(),
  stat_km_curves = rd(f06, "06_km_curves.csv") %>% add_level_key(),
  stat_survival_at_years = rd(f06, "06_survival_at_years.csv") %>% add_level_key(),
  stat_hazard_ratios = rd(f06, "06_hazard_ratios.csv") %>% rename(adjustment = model) %>% mutate(adjustment_short = word(adjustment, 1, sep = ":")),
  stat_hr_full_model = rd(f06, "06_m3_full_model.csv"),
  stat_spline_curves = rd(f06, "06_spline_curves.csv"),
  stat_hr_sensitivity = rd(f06, "06_sensitivity.csv"),
  stat_cause_mix = rd(f07, "07_cause_mix.csv") %>% add_level_key(),
  stat_cause_detail = rd(f07, "07_cause_detail.csv") %>% rename(underlying_cause = cause) %>% mutate(level_key = paste("glycemic", glycemic_status, sep = "|")),
  stat_cause_specific_rates = rd(f07, "07_cause_specific_rates.csv") %>% add_level_key(),
  stat_cumulative_incidence = rd(f07, "07_cumulative_incidence.csv") %>% add_level_key(),
  stat_cumulative_incidence_15y = rd(f07, "07_cumulative_incidence_15y.csv") %>% add_level_key(),
  stat_cause_specific_hr = rd(f07, "07_cause_specific_hr.csv"),
  stat_finegray = rd(f07, "07_cause_specific_vs_finegray.csv") %>%
    { bind_rows(bind_cols(select(., cause, term), parse_hr(.[["Cause-specific HR"]])) %>% mutate(hr_type = "Cause-specific"),
                bind_cols(select(., cause, term), parse_hr(.[["Subdistribution HR (Fine-Gray)"]])) %>% mutate(hr_type = "Fine-Gray")) },
  stat_dm_death_certificate = rd(f07, "07_diabetes_on_death_certificate.csv") %>% mutate(level_key = paste("glycemic", glycemic_status, sep = "|")),
  stat_paf = rd(f08, "08_population_attributable_fractions.csv") %>% mutate(cause = recode(outcome, `All-cause` = "All causes")) %>% select(-outcome) %>% rename(risk_factor = exposure),
  stat_ipw_hr = rd(f08, "08_ipw_hazard_ratios.csv") %>%
    { bind_rows(bind_cols(select(., exposure), parse_hr(.[["IPW marginal HR"]])) %>% mutate(hr_type = "IPW marginal"),
                bind_cols(select(., exposure), parse_hr(.[["Conditional HR (06, M3)"]])) %>% mutate(hr_type = "Conditional (M3)")) },
  stat_ipw_survival = rd(f08, "08_ipw_adjusted_survival.csv"),
  stat_ipw_10y_risk = rd(f08, "08_ipw_10y_risk.csv"),
  stat_e_values = rd(f08, "08_e_values.csv"),
  stat_covariate_balance = rd(f08, "08_covariate_balance.csv"),
  model_performance = performance_tbl,
  model_calibration = calibration %>% mutate(across(c(predicted, observed), function(x) round(100 * x, 2))) %>% rename(predicted_pct = predicted, observed_pct = observed),
  model_decision_curve = dca %>% mutate(threshold_pct = 100 * threshold) %>% select(model, threshold_pct, net_benefit),
  model_risk_tiers = risk_tiers %>% mutate(risk_tier = as.character(risk_tier)),
  model_whatif_profiles = whatif %>% mutate(across(where(is.factor), as.character)) %>% rename(glycemic_status = glycemic4),
  model_coefficients = coefficients)
tables <- c(dims, facts)

# ---- Relationships: many (fact/result table) to one (dimension), single direction ----
relationships <- tribble(
  ~from_table, ~from_column, ~to_table, ~to_column,
  "fact_participants", "cycle", "dim_cycle", "cycle",
  "fact_participants", "level_key", "dim_group_level", "level_key",
  "fact_participants", "risk_tier", "dim_risk_tier", "risk_tier",
  "stat_prevalence_by_cycle", "cycle", "dim_cycle", "cycle",
  "stat_prevalence_overall", "condition", "dim_condition", "condition",
  "stat_prevalence_by_cycle", "condition", "dim_condition", "condition",
  "stat_prevalence_by_group", "condition", "dim_condition", "condition",
  "stat_prevalence_trend_tests", "condition", "dim_condition", "condition",
  "stat_prevalence_by_group", "level_key", "dim_group_level", "level_key",
  "stat_death_rates", "level_key", "dim_group_level", "level_key",
  "stat_km_curves", "level_key", "dim_group_level", "level_key",
  "stat_survival_at_years", "level_key", "dim_group_level", "level_key",
  "stat_cause_mix", "level_key", "dim_group_level", "level_key",
  "stat_cause_detail", "level_key", "dim_group_level", "level_key",
  "stat_cause_specific_rates", "level_key", "dim_group_level", "level_key",
  "stat_cumulative_incidence", "level_key", "dim_group_level", "level_key",
  "stat_cumulative_incidence_15y", "level_key", "dim_group_level", "level_key",
  "stat_dm_death_certificate", "level_key", "dim_group_level", "level_key",
  "stat_cause_mix", "cause", "dim_cause", "cause",
  "stat_cause_specific_rates", "cause", "dim_cause", "cause",
  "stat_cumulative_incidence", "cause", "dim_cause", "cause",
  "stat_cumulative_incidence_15y", "cause", "dim_cause", "cause",
  "stat_cause_specific_hr", "cause", "dim_cause", "cause",
  "stat_finegray", "cause", "dim_cause", "cause",
  "stat_paf", "cause", "dim_cause", "cause",
  "stat_hazard_ratios", "exposure", "dim_exposure", "exposure",
  "stat_ipw_hr", "exposure", "dim_exposure", "exposure",
  "stat_ipw_survival", "exposure", "dim_exposure", "exposure",
  "stat_ipw_10y_risk", "exposure", "dim_exposure", "exposure",
  "stat_covariate_balance", "exposure", "dim_exposure", "exposure",
  "model_risk_tiers", "risk_tier", "dim_risk_tier", "risk_tier",
  "model_performance", "model", "dim_prediction_model", "model",
  "model_calibration", "model", "dim_prediction_model", "model",
  "model_decision_curve", "model", "dim_prediction_model", "model",
  "model_coefficients", "model", "dim_prediction_model", "model") %>%
  mutate(cardinality = "Many to one", cross_filter = "Single",
         unmatched_keys = pmap_int(list(from_table, from_column, to_table, to_column), function(ft, fc, tt, tc) length(setdiff(unique(tables[[ft]][[fc]]), tables[[tt]][[tc]]))))

cat("Referential integrity (every key in a result table must exist in its dimension):\n")
print(relationships %>% select(from_table, from_column, to_table, unmatched_keys), n = Inf, width = Inf)
stopifnot("Some keys have no match in their dimension table" = all(relationships$unmatched_keys == 0),
          "A dimension key is not unique" = all(map_lgl(unique(relationships$to_table), function(t) !anyDuplicated(tables[[t]][[relationships$to_column[relationships$to_table == t][1]]]))))

# ---- Write ----
unlink(list.files(pbi, pattern = "\\.csv$", full.names = TRUE))
iwalk(tables, function(x, name) write_csv(x, file.path(pbi, paste0(name, ".csv")), na = ""))
write_csv(select(relationships, -unmatched_keys), file.path(pbi, "_relationships.csv"))

linked <- unique(c(relationships$from_table, relationships$to_table))
dictionary <- tibble(table = names(tables), rows = map_int(tables, nrow), columns = map_chr(tables, function(x) paste(names(x), collapse = ", "))) %>%
  mutate(role = case_when(startsWith(table, "dim_") ~ "Dimension", startsWith(table, "fact_") ~ "Fact", startsWith(table, "model_") ~ "Model result", TRUE ~ "Statistical result"),
         in_relationships = table %in% linked,
         source = case_when(startsWith(table, "dim_") ~ "Script 09", table %in% c("fact_participants", "stat_kpis") ~ "Scripts 05-09",
                            grepl("prevalence|table1|design", table) ~ "Script 05", grepl("death_rates|km|survival_at|hazard|hr_full|spline|hr_sens", table) ~ "Script 06",
                            grepl("cause|incidence|finegray|certificate", table) ~ "Script 07", startsWith(table, "model_") ~ "Script 09", TRUE ~ "Script 08"))
write_csv(dictionary, file.path(pbi, "_table_dictionary.csv"))

cat(sprintf("\nPower BI model written to %s\n", pbi))
cat(sprintf("%d tables: %d dimensions, %d fact, %d statistical results, %d model results | %d relationships (%d tables linked, %d standalone)\n",
            nrow(dictionary), sum(dictionary$role == "Dimension"), sum(dictionary$role == "Fact"), sum(dictionary$role == "Statistical result"),
            sum(dictionary$role == "Model result"), nrow(relationships), sum(dictionary$in_relationships), sum(!dictionary$in_relationships)))
print(dictionary %>% select(table, role, rows, in_relationships, source), n = Inf, width = Inf)

writeLines(capture.output(sessionInfo()), file.path(out, "09_session_info.txt"))
cat(sprintf("\n✓ Part A outputs saved in: %s\n", out))
print(list.files(out))
