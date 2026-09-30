# ============================================================================
# Script 05: Survey design, multiple imputation and weighted prevalence
# NHANES 1999-2018 linked to NDI mortality through 2019
# Input : data/processed/analytic_cohort.parquet (notebook 02)
# Output: outputs/R/05_survey_design_imputation/
# ============================================================================

library(tidyverse)
library(arrow)
library(survey)
library(mitools)
library(mice)
options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE, scipen = 999)

# Set paths
root <- "E:/NHANES Cardiometabolic Mortality Project"
out <- file.path(root, "outputs", "R", "05_survey_design_imputation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

show_save <- function(x, name) {
  write_csv(x, file.path(out, paste0(name, ".csv")), na = "")
  print(x, n = Inf, width = Inf)
  invisible(x)
}
plot_save <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "#fcfcfb")
  print(p)
}
mi_tidy <- function(results) {
  r <- MIcombine(results)
  tibble(term = names(coef(r)), est = unname(coef(r)), se = sqrt(unname(diag(vcov(r))))) %>%
    mutate(lcl = est - qnorm(0.975) * se, ucl = est + qnorm(0.975) * se, p = 2 * pnorm(-abs(est / se)))
}
header <- function(x) cat("\n============================================================\n", x, "\n============================================================\n", sep = "")

# 2000 US standard population (NCHS Statistical Notes No. 20), thousands: weights 0.3966, 0.3718, 0.2316
std_pop <- tibble(age3 = c("20-39", "40-59", "60+"), pop = c(77670, 72816, 45364))

pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
                  plot.title = element_text(face = "bold", colour = "#0b0b0b"), plot.subtitle = element_text(colour = "#52514e"),
                  axis.text = element_text(colour = "#52514e"), strip.text = element_text(face = "bold", hjust = 0, colour = "#0b0b0b"),
                  legend.position = "top", legend.justification = "left", plot.background = element_rect(fill = "#fcfcfb", colour = NA)))

# ----------------------------------------------------------------------------
# 1. Load the cohort
# derive() rebuilds status variables with the notebook 02 rules. Observed
# statuses are kept; only missing ones are filled (used again after imputation).
# ----------------------------------------------------------------------------
header("1. LOAD COHORT")

derive <- function(x) mutate(x,
  acr = exp(log_acr),
  bmi_cat = cut(bmi, c(0, 18.5, 25, 30, Inf), right = FALSE, labels = c("underweight", "normal", "overweight", "obese")),
  obesity = as.numeric(bmi >= 30),
  dm_dx = as.numeric(dm_told == 1 | dm_insulin %in% 1 | dm_pills %in% 1),
  diabetes = coalesce(diabetes, as.numeric(dm_dx == 1 | hba1c >= 6.5 | coalesce(fpg >= 126, FALSE))),
  prediabetes = if_else(diabetes %in% 1, 0, coalesce(prediabetes, as.numeric((hba1c >= 5.7 & hba1c < 6.5) | coalesce(fpg >= 100 & fpg < 126, FALSE)))),
  glycemic = factor(case_when(diabetes == 1 ~ "diabetes", prediabetes == 1 ~ "prediabetes", diabetes == 0 & prediabetes == 0 ~ "normal"), c("normal", "prediabetes", "diabetes")),
  dm_undiagnosed = as.numeric(diabetes == 1 & dm_dx == 0),
  hypertension = coalesce(hypertension, as.numeric(sbp >= 140 | dbp >= 90 | htn_med %in% 1)),
  hypertension_acc = coalesce(hypertension_acc, as.numeric(sbp >= 130 | dbp >= 80 | htn_med %in% 1)),
  albuminuria = coalesce(albuminuria, as.numeric(acr >= 30)),
  ckd = coalesce(ckd, as.numeric(egfr < 60 | acr >= 30)),
  current_smoker = as.numeric(smoking == "current"))

cohort <- read_parquet(file.path(root, "data", "processed", "analytic_cohort.parquet")) %>%
  transmute(SEQN, cycle = factor(cycle), cycle_num = as.integer(cycle), psu, strata, wt_mec,
            age, age3 = cut(age, c(20, 40, 60, Inf), right = FALSE, labels = std_pop$age3), age_group = factor(age_group),
            sex = factor(sex_lbl, c("male", "female")),
            race_eth = factor(race_eth_lbl, c("Non-Hispanic White", "Non-Hispanic Black", "Mexican American", "Other Hispanic", "Other/Multiracial")),
            education = factor(education_lbl, c("<HS", "HS/GED", ">HS")), pir, smoking = factor(smoking, c("never", "former", "current")),
            bmi, waist, sbp, dbp, hba1c, fpg, tc, hdl, non_hdl, tg, ldl, log_acr = log(acr), egfr,
            dm_told, dm_insulin, dm_pills, htn_told, htn_med, chol_told, chol_med,
            diabetes, prediabetes, hypertension, hypertension_acc, albuminuria, ckd,
            chf, chd, angina, mi, stroke, cvd_history, cancer,
            time_months, death, cause = factor(if_else(cause == "unknown", "other", cause), c("alive", "cvd", "cancer", "other")),
            ucod_leading, mcod_diabetes, mcod_hypertension) %>%
  derive()

cat(sprintf("Participants: %s | Deaths: %s | Person-years: %s | Cycles: %d\n", format(nrow(cohort), big.mark = ","),
            format(sum(cohort$death), big.mark = ","), format(round(sum(cohort$time_months) / 12), big.mark = ","), nlevels(cohort$cycle)))
print(table(cohort$cause))

# ----------------------------------------------------------------------------
# 2. Survey design
# Strata are unique across cycles, so 10 cycles stack with nest = TRUE.
# wt_mec is the combined 10-cycle MEC weight from notebook 02. If the cohort
# keeps every PSU of the full sample, a design on the cohort gives the same
# standard errors as a full-sample domain analysis.
# ----------------------------------------------------------------------------
header("2. SURVEY DESIGN")

psu_full <- read_parquet(file.path(root, "data", "parquet", "demographics.parquet"), col_select = c("SDMVSTRA", "SDMVPSU")) %>% distinct() %>% nrow()
des_obs <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt_mec, nest = TRUE, data = cohort)

design_summary <- tibble(participants = nrow(cohort), strata = n_distinct(cohort$strata), psu_in_cohort = nrow(distinct(cohort, strata, psu)),
                         psu_in_full_nhanes = psu_full, min_participants_per_psu = min(count(cohort, strata, psu)$n),
                         design_df = degf(des_obs), weighted_population_millions = round(sum(cohort$wt_mec) / 1e6, 1))
show_save(design_summary, "05_design_summary")

# ----------------------------------------------------------------------------
# 3. Missing data
# fpg, tg and ldl are fasting-subsample only (~55% missing by design) and are
# not imputed; medication questions are skip patterns and are not imputed.
# ----------------------------------------------------------------------------
header("3. MISSING DATA")

imp_vars <- c("pir", "education", "smoking", "bmi", "waist", "sbp", "dbp", "hba1c", "tc", "hdl", "log_acr", "egfr", "dm_told", "cvd_history", "cancer")

missing_summary <- tibble(variable = imp_vars, n_missing = colSums(is.na(cohort[imp_vars]))) %>%
  mutate(pct_missing = round(100 * n_missing / nrow(cohort), 2)) %>% arrange(desc(n_missing))
show_save(missing_summary, "05_missing_summary")

complete_vs_incomplete <- cohort %>%
  mutate(record = if_else(complete.cases(pick(all_of(imp_vars))), "complete", "incomplete")) %>%
  group_by(record) %>%
  summarise(n = n(), pct = round(100 * n / nrow(cohort), 1), mean_age = round(mean(age), 1), deaths = sum(death),
            deaths_per_1000_py = round(1000 * sum(death) / sum(time_months / 12), 1))
show_save(complete_vs_incomplete, "05_complete_vs_incomplete")
cat("Incomplete records are older and die at a higher rate, so complete-case analysis would bias results. Using multiple imputation.\n")

# ----------------------------------------------------------------------------
# 4. Multiple imputation (mice, m = 20, 10 iterations, parallel seed 2026)
# pmm for continuous/binary, polyreg for education and smoking. Death indicator
# and Nelson-Aalen cumulative hazard included (White & Royston 2009).
# Cached: delete 05_mice.rds to impute again. First run ~10-30 minutes.
# ----------------------------------------------------------------------------
header("4. MULTIPLE IMPUTATION")

imp_data <- cohort %>%
  mutate(na_hazard = nelsonaalen(cohort, time_months, death), log_wt = log(wt_mec)) %>%
  select(all_of(imp_vars), age, sex, race_eth, cycle, death, na_hazard, log_wt)

imp_file <- file.path(out, "05_mice.rds")
if (file.exists(imp_file)) {
  cat("Loading cached imputation...\n")
  imp <- readRDS(imp_file)
} else {
  cat("Running imputation in parallel (this takes a while)...\n")
  imp <- futuremice(imp_data, m = 20, maxit = 10, parallelseed = 2026)
  saveRDS(imp, imp_file)
}
print(imp$method[imp$method != ""])
if (is.null(imp$loggedEvents)) cat("No logged events: no predictor dropped for collinearity.\n") else print(imp$loggedEvents)

cont_vars <- c("pir", "bmi", "waist", "sbp", "dbp", "hba1c", "tc", "hdl", "log_acr", "egfr")

p_conv <- as.data.frame.table(imp$chainMean, responseName = "mean") %>%
  setNames(c("variable", "iteration", "chain", "mean")) %>%
  filter(variable %in% cont_vars) %>%
  mutate(iteration = as.integer(iteration)) %>%
  ggplot(aes(iteration, mean, group = chain)) +
  geom_line(colour = pal[1], alpha = 0.35, linewidth = 0.5) +
  facet_wrap(~variable, scales = "free_y", ncol = 5) +
  labs(title = "Imputation chains mix without trend", subtitle = paste0("Mean of imputed values per iteration, one line per imputation (m = ", imp$m, ")"), x = "Iteration", y = NULL)
plot_save(p_conv, "05_imputation_convergence")

obs_imp <- map_dfr(cont_vars, function(v) bind_rows(
  tibble(variable = v, source = "Observed", value = na.omit(cohort[[v]])),
  tibble(variable = v, source = "Imputed", value = unlist(imp$imp[[v]], use.names = FALSE)))) %>%
  mutate(source = factor(source, c("Observed", "Imputed")))

p_dens <- ggplot(obs_imp, aes(value, colour = source)) +
  geom_density(linewidth = 0.7) +
  facet_wrap(~variable, scales = "free", ncol = 5) +
  scale_colour_manual(values = c(Observed = pal[1], Imputed = pal[2])) +
  labs(title = "Imputed values stay within the observed range", subtitle = paste0("Observed values vs values imputed across all ", imp$m, " imputations"), x = NULL, y = NULL, colour = NULL) +
  theme(axis.text.y = element_blank())
plot_save(p_dens, "05_imputation_density")

imputation_check <- obs_imp %>%
  group_by(variable, source) %>% summarise(mean = mean(value)) %>% ungroup() %>%
  pivot_wider(names_from = source, values_from = mean) %>%
  mutate(across(c(Observed, Imputed), function(x) round(x, 2)), difference = Imputed - Observed)
show_save(imputation_check, "05_imputation_check")

# ----------------------------------------------------------------------------
# 5. Completed data: 20 stacked datasets with statuses re-derived.
# Saved as 05_imputed_long.rds, the input for scripts 06-09.
# ----------------------------------------------------------------------------
header("5. COMPLETED DATA")

long <- complete(imp, "long") %>%
  select(.imp, .id, all_of(imp_vars)) %>%
  mutate(.id = as.integer(.id)) %>%
  left_join(cohort %>% select(-all_of(imp_vars)) %>% mutate(.id = row_number()), by = ".id") %>%
  derive()
saveRDS(long, file.path(out, "05_imputed_long.rds"))

cat(sprintf("Rows: %s | Imputations: %d | Participants: %s | Missing glycemic/hypertension/CKD: %d/%d/%d\n",
            format(nrow(long), big.mark = ","), n_distinct(long$.imp), format(n_distinct(long$SEQN), big.mark = ","),
            sum(is.na(long$glycemic)), sum(is.na(long$hypertension)), sum(is.na(long$ckd))))

# ----------------------------------------------------------------------------
# 6. Weighted prevalence, pooled over imputations with Rubin's rules.
# Age-standardized: direct method, 2000 US standard population (20-39, 40-59, 60+).
# ----------------------------------------------------------------------------
header("6. WEIGHTED PREVALENCE")

des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt_mec, nest = TRUE, data = imputationList(split(long, long$.imp)))
conds <- c(diabetes = "Diabetes (total)", dm_dx = "Diagnosed diabetes", dm_undiagnosed = "Undiagnosed diabetes", prediabetes = "Prediabetes",
           hypertension = "Hypertension (140/90 or meds)", hypertension_acc = "Hypertension (130/80 or meds)", obesity = "Obesity (BMI 30+)",
           ckd = "Chronic kidney disease", albuminuria = "Albuminuria (ACR 30+)", current_smoker = "Current smoking", cvd_history = "Self-reported CVD")
f_conds <- reformulate(names(conds))
pct <- function(x) mutate(x, across(c(est, lcl, ucl), function(v) round(100 * v, 1)), label = unname(conds[condition]))

cat("\nOverall, 1999-2018 (crude and age-standardized with MI, available cases for comparison):\n")
prev_overall <- bind_rows(
  with(des, fun = function(d) svymean(f_conds, d)) %>% mi_tidy() %>% mutate(condition = term, estimate = "Crude (MI)"),
  with(des, fun = function(d) svymean(f_conds, svystandardize(d, ~age3, ~1, std_pop$pop))) %>% mi_tidy() %>% mutate(condition = term, estimate = "Age-standardized (MI)"),
  map_dfr(names(conds), function(v) {
    s <- svymean(reformulate(v), des_obs, na.rm = TRUE)
    tibble(condition = v, est = unname(coef(s)), lcl = confint(s)[1], ucl = confint(s)[2], n_available = sum(!is.na(cohort[[v]])))
  }) %>% mutate(estimate = "Crude (available cases)")) %>%
  pct() %>% select(label, condition, estimate, prevalence_pct = est, lcl, ucl, n_available)
write_csv(prev_overall, file.path(out, "05_prevalence_overall.csv"), na = "")
prev_overall %>%
  mutate(value = sprintf("%.1f (%.1f-%.1f)", prevalence_pct, lcl, ucl)) %>%
  select(label, estimate, value) %>%
  pivot_wider(names_from = estimate, values_from = value) %>%
  print(n = Inf, width = Inf)

cat("\nBy survey cycle:\n")
by_cycle <- bind_rows(
  with(des, fun = function(d) svyby(f_conds, ~cycle, d, svymean)) %>% mi_tidy() %>% mutate(type = "Crude"),
  with(des, fun = function(d) svyby(f_conds, ~cycle, svystandardize(d, ~age3, ~cycle, std_pop$pop), svymean)) %>% mi_tidy() %>% mutate(type = "Age-standardized")) %>%
  separate(term, c("cycle", "condition"), sep = ":") %>%
  pct() %>% select(label, condition, cycle, type, prevalence_pct = est, lcl, ucl)
write_csv(by_cycle, file.path(out, "05_prevalence_by_cycle.csv"), na = "")
by_cycle %>%
  filter(type == "Age-standardized") %>%
  select(label, cycle, prevalence_pct) %>%
  pivot_wider(names_from = cycle, values_from = prevalence_pct) %>%
  print(n = Inf, width = Inf)

cat("\nTrend tests (survey logistic regression on cycle 1-10, adjusted for age group, sex, race/ethnicity; OR per 2-year cycle):\n")
trend_tests <- map_dfr(names(conds), function(v)
  with(des, fun = function(d) svyglm(reformulate(c("cycle_num", "age3", "sex", "race_eth"), v), d, family = quasibinomial())) %>%
    mi_tidy() %>% filter(term == "cycle_num") %>% mutate(condition = v)) %>%
  transmute(label = unname(conds[condition]), condition, or_per_cycle = round(exp(est), 3), lcl = round(exp(lcl), 3), ucl = round(exp(ucl), 3),
            p_trend = format.pval(p, digits = 2, eps = 0.001))
show_save(trend_tests, "05_trend_tests")

key <- c("diabetes", "dm_undiagnosed", "prediabetes", "hypertension", "obesity", "ckd")

p_trend <- by_cycle %>%
  filter(condition %in% key) %>%
  mutate(label = factor(label, conds[key]), cycle = factor(cycle), type = factor(type, c("Age-standardized", "Crude"))) %>%
  ggplot(aes(cycle, prevalence_pct, colour = type, fill = type, group = type)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~label, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = pal[1:2]) +
  scale_fill_manual(values = pal[1:2]) +
  scale_x_discrete(labels = function(x) substr(x, 1, 4)) +
  labs(title = "Cardiometabolic conditions in US adults, 1999-2018", subtitle = paste0("Weighted prevalence (%) with 95% CI, pooled over ", imp$m, " imputations"),
       x = "Survey cycle (start year)", y = "Prevalence (%)", colour = NULL, fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_save(p_trend, "05_prevalence_trends", h = 7)

cat("\nBy population group (age-specific by age group; age-standardized by sex, race/ethnicity, education):\n")
prev_by <- function(by, std = TRUE) with(des, fun = function(d) svyby(f_conds, by, if (std) svystandardize(d, ~age3, by, std_pop$pop) else d, svymean)) %>%
  mi_tidy() %>% separate(term, c("level", "condition"), sep = ":") %>% mutate(group = all.vars(by)[1], type = if (std) "Age-standardized" else "Age-specific")

by_group <- bind_rows(prev_by(~age3, std = FALSE), prev_by(~sex), prev_by(~race_eth), prev_by(~education)) %>%
  pct() %>% select(label, condition, group, level, type, prevalence_pct = est, lcl, ucl)
write_csv(by_group, file.path(out, "05_prevalence_by_group.csv"), na = "")
by_group %>%
  filter(condition %in% key) %>%
  select(group, level, label, prevalence_pct) %>%
  pivot_wider(names_from = label, values_from = prevalence_pct) %>%
  print(n = Inf, width = Inf)

p_race <- by_group %>%
  filter(group == "race_eth", condition %in% key) %>%
  mutate(label = factor(label, conds[key]), level = factor(level, rev(levels(cohort$race_eth)))) %>%
  ggplot(aes(prevalence_pct, level)) +
  geom_linerange(aes(xmin = lcl, xmax = ucl), colour = pal[1], linewidth = 0.6) +
  geom_point(colour = pal[1], size = 2.5) +
  facet_wrap(~label, scales = "free_x", ncol = 3) +
  labs(title = "Age-standardized prevalence by race and ethnicity", subtitle = "US adults 20+, NHANES 1999-2018, weighted, 95% CI", x = "Prevalence (%)", y = NULL) +
  theme(panel.spacing.x = unit(1.5, "lines"))
plot_save(p_race, "05_prevalence_by_race", h = 7)

# ----------------------------------------------------------------------------
# 7. Weighted baseline characteristics (Table 1) by glycemic status
# ACR as geometric mean; n and deaths are unweighted, averaged over imputations.
# ----------------------------------------------------------------------------
header("7. TABLE 1 (WEIGHTED)")

t1_vars <- ~age + sex + race_eth + education + pir + smoking + bmi + waist + sbp + dbp + hba1c + tc + hdl + egfr + log_acr +
  hypertension + ckd + cvd_history + cancer
cont <- c("age", "pir", "bmi", "waist", "sbp", "dbp", "hba1c", "tc", "hdl", "egfr", "log_acr")
t1_labels <- c(n = "Participants, n (unweighted)", deaths = "Deaths, n (unweighted)", age = "Age, years", pir = "Income-to-poverty ratio",
               bmi = "BMI, kg/m2", waist = "Waist circumference, cm", sbp = "Systolic BP, mmHg", dbp = "Diastolic BP, mmHg", hba1c = "HbA1c, %",
               tc = "Total cholesterol, mg/dL", hdl = "HDL cholesterol, mg/dL", egfr = "eGFR, mL/min/1.73m2", log_acr = "Urine ACR, mg/g (geometric mean)",
               hypertension = "Hypertension", ckd = "Chronic kidney disease", cvd_history = "Self-reported CVD", cancer = "Self-reported cancer")
prefixes <- c(sex = "Sex", race_eth = "Race/ethnicity", education = "Education", smoking = "Smoking")
label_term <- function(t) {
  p <- names(prefixes)[startsWith(t, names(prefixes))]
  if (!is.na(t1_labels[t])) t1_labels[[t]] else paste0(prefixes[[p]], ": ", substring(t, nchar(p) + 1))
}

t1 <- bind_rows(
  with(des, fun = function(d) svyby(t1_vars, ~glycemic, d, svymean)) %>% mi_tidy() %>% separate(term, c("group", "term"), sep = ":"),
  with(des, fun = function(d) svymean(t1_vars, d)) %>% mi_tidy() %>% mutate(group = "overall")) %>%
  mutate(across(c(est, lcl, ucl), function(v) case_when(term == "log_acr" ~ exp(v), term %in% cont ~ v, TRUE ~ 100 * v)),
         value = if_else(term %in% cont, sprintf("%.1f (%.1f-%.1f)", est, lcl, ucl), sprintf("%.1f%%", est)))

counts <- long %>%
  group_by(.imp, group = as.character(glycemic)) %>% summarise(n = n(), deaths = sum(death)) %>%
  bind_rows(long %>% group_by(.imp) %>% summarise(n = n(), deaths = sum(death)) %>% mutate(group = "overall")) %>%
  group_by(group) %>% summarise(across(c(n, deaths), function(x) as.character(round(mean(x))))) %>%
  pivot_longer(c(n, deaths), names_to = "term", values_to = "value")

table1 <- bind_rows(counts, t1 %>% select(group, term, value)) %>%
  filter(term != "sexmale") %>%
  pivot_wider(names_from = group, values_from = value) %>%
  mutate(characteristic = map_chr(term, label_term)) %>%
  select(characteristic, normal, prediabetes, diabetes, overall)
show_save(table1, "05_table1_weighted")

# ----------------------------------------------------------------------------
# 8. Summary
# ----------------------------------------------------------------------------
header("8. SUMMARY: AGE-STANDARDIZED PREVALENCE AND TREND")

prev_overall %>%
  filter(estimate == "Age-standardized (MI)") %>%
  select(label, prevalence_pct, lcl, ucl) %>%
  left_join(trend_tests %>% select(label, or_per_cycle, p_trend), by = "label") %>%
  print(n = Inf, width = Inf)

writeLines(capture.output(sessionInfo()), file.path(out, "05_session_info.txt"))
cat(sprintf("\n✓ All outputs saved in: %s\n", out))
print(list.files(out))
