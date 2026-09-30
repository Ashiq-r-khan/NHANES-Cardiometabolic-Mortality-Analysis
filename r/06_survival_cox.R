# ============================================================================
# Script 06: Survival analysis, weighted Kaplan-Meier and Cox models
# NHANES 1999-2018 linked to NDI mortality through 2019
# Input : outputs/R/05_survey_design_imputation/05_imputed_long.rds, 05_mice.rds
# Output: outputs/R/06_survival_cox/
# ============================================================================

library(tidyverse)
library(survey)
library(mitools)
library(splines)

options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE, scipen = 999)

# Set paths
root <- "E:/NHANES Cardiometabolic Mortality Project"
inp <- file.path(root, "outputs", "R", "05_survey_design_imputation")
out <- file.path(root, "outputs", "R", "06_survival_cox")
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
fmt_p <- function(p) if_else(p < 0.001, "<0.001", sprintf("%.3f", p))
fmt_hr <- function(est, lcl, ucl) sprintf("%.2f (%.2f-%.2f)", exp(est), exp(lcl), exp(ucl))

std_prop <- c(`20-39` = 0.3966, `40-59` = 0.3718, `60+` = 0.2316)

pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
                  plot.title = element_text(face = "bold", colour = "#0b0b0b"), plot.subtitle = element_text(colour = "#52514e"),
                  axis.text = element_text(colour = "#52514e"), strip.text = element_text(face = "bold", hjust = 0, colour = "#0b0b0b"),
                  legend.position = "top", legend.justification = "left", plot.background = element_rect(fill = "#fcfcfb", colour = NA)))

# Cox fits are slow (~440 survey Cox models over 20 imputations, ~10-20 min on first run).
# Each result is cached in 06_cox_cache.rds as soon as it finishes, so a rerun or an
# interrupted run picks up where it stopped. Delete the file to refit everything.
cache_file <- file.path(out, "06_cox_cache.rds")
cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()
cached <- function(key, expr) {
  if (is.null(cache[[key]])) {
    cat(sprintf("  fitting: %s\n", key))
    cache[[key]] <<- expr
    saveRDS(cache, cache_file)
  }
  cache[[key]]
}

# ----------------------------------------------------------------------------
# 1. Load imputed data and build the survey design
# ----------------------------------------------------------------------------
header("1. LOAD DATA")

long <- readRDS(file.path(inp, "05_imputed_long.rds")) %>%
  mutate(time_y = time_months / 12,
         glycemic4 = factor(case_when(glycemic == "diabetes" & dm_dx == 1 ~ "diagnosed diabetes", glycemic == "diabetes" ~ "undiagnosed diabetes", TRUE ~ as.character(glycemic)),
                            c("normal", "prediabetes", "undiagnosed diabetes", "diagnosed diabetes")),
         bmi_cat = factor(bmi_cat, c("normal", "underweight", "overweight", "obese")),
         hypertension_f = factor(hypertension, 0:1, c("no hypertension", "hypertension")),
         ckd_f = factor(ckd, 0:1, c("no CKD", "CKD")),
         t10 = pmin(time_y, 10), d10 = as.numeric(death == 1 & time_y <= 10), t_late = time_y - 10) %>%
  group_by(SEQN) %>% mutate(any_cvd_cancer = max(cvd_history == 1 | cancer == 1)) %>% ungroup()

des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt_mec, nest = TRUE, data = imputationList(split(long, long$.imp)))
m <- n_distinct(long$.imp)
cat(sprintf("Imputations: %d | Participants: %s | Deaths: %s | Person-years: %s | Max follow-up: %.1f years\n", m,
            format(n_distinct(long$SEQN), big.mark = ","), format(sum(long$death[long$.imp == 1]), big.mark = ","),
            format(round(sum(long$time_y[long$.imp == 1])), big.mark = ","), max(long$time_y)))

# ----------------------------------------------------------------------------
# 2. Weighted death rates per 1,000 person-years
# Crude: survey ratio of deaths to person-years. Age-standardized: age-specific
# rates (20-39, 40-59, 60+) weighted to the 2000 US standard population.
# ----------------------------------------------------------------------------
header("2. WEIGHTED DEATH RATES PER 1,000 PERSON-YEARS")

rate_groups <- c(glycemic4 = "Glycemic status", hypertension_f = "Hypertension", ckd_f = "Chronic kidney disease", bmi_cat = "BMI category", smoking = "Smoking")

death_rates <- map_dfr(names(rate_groups), function(v) {
  lv <- levels(long[[v]])
  crude <- with(des, fun = function(d) svyby(~death, reformulate(v), denominator = ~time_y, design = d, FUN = svyratio)) %>% mi_tidy() %>% mutate(level = term, type = "Crude")
  std <- with(des, fun = function(d) {
    r <- svyby(~death, reformulate(c(v, "age3")), denominator = ~time_y, design = d, FUN = svyratio, covmat = TRUE)
    svycontrast(r, setNames(lapply(lv, function(l) setNames(std_prop, paste(l, names(std_prop), sep = "."))), lv))
  }) %>% mi_tidy() %>% mutate(level = term, type = "Age-standardized")
  bind_rows(crude, std) %>% mutate(group = rate_groups[[v]], level = factor(level, lv))
}) %>%
  transmute(group, level, type, rate_per_1000_py = round(1000 * est, 1), lcl = round(1000 * lcl, 1), ucl = round(1000 * ucl, 1))
show_save(death_rates, "06_death_rates")

p_rates <- death_rates %>%
  mutate(type = factor(type, c("Age-standardized", "Crude")), group = factor(group, rate_groups)) %>%
  ggplot(aes(rate_per_1000_py, fct_rev(level), colour = type)) +
  geom_linerange(aes(xmin = lcl, xmax = ucl), linewidth = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.3, position = position_dodge(width = 0.6)) +
  facet_wrap(~group, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = pal[1:2]) +
  labs(title = "All-cause death rates by cardiometabolic status", subtitle = "Survey-weighted deaths per 1,000 person-years with 95% CI, pooled over imputations",
       x = "Deaths per 1,000 person-years", y = NULL, colour = NULL)
plot_save(p_rates, "06_death_rates", h = 8)

# ----------------------------------------------------------------------------
# 3. Weighted Kaplan-Meier survival, averaged over imputations
# Curves stop at 18 years: only the 1999-2000 cycle is followed longer, so the
# tail rests on few people and jumps.
# ----------------------------------------------------------------------------
header("3. WEIGHTED KAPLAN-MEIER SURVIVAL")

grid <- seq(0, 18, by = 0.1)
km_groups <- c(glycemic4 = "Glycemic status", hypertension_f = "Hypertension", ckd_f = "Chronic kidney disease", bmi_cat = "BMI category")

km <- map_dfr(names(km_groups), function(v) map_dfr(des$designs, function(d) {
  k <- svykm(reformulate(v, "Surv(time_y, death)"), d, se = FALSE)
  map_dfr(names(k), function(g) tibble(level = g, time = grid, surv = approx(k[[g]]$time, k[[g]]$surv, grid, method = "constant", f = 0, rule = 2, ties = last)$y))
}) %>%
  group_by(level, time) %>% summarise(surv = mean(surv)) %>% ungroup() %>%
  mutate(group = km_groups[[v]], level = factor(level, levels(long[[v]]))))
write_csv(km, file.path(out, "06_km_curves.csv"))

survival_at <- km %>%
  filter(time %in% c(5, 10, 15)) %>%
  mutate(surv = round(100 * surv, 1), time = paste0("survival_", time, "y_pct")) %>%
  pivot_wider(names_from = time, values_from = surv) %>%
  select(group, level, everything()) %>% arrange(factor(group, km_groups), level)
show_save(survival_at, "06_survival_at_years")

walk2(names(km_groups), km_groups, function(v, g) {
  p <- km %>%
    filter(group == g) %>%
    ggplot(aes(time, 100 * surv, colour = level)) +
    geom_step(linewidth = 0.8) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(breaks = seq(0, 18, 3)) +
    labs(title = paste("Survey-weighted survival by", tolower(g)), subtitle = "Kaplan-Meier estimates, US adults 20+, NHANES 1999-2018, mortality follow-up to 2019",
         x = "Years of follow-up", y = "Survival (%)", colour = NULL)
  plot_save(p, paste0("06_km_", sub("4$|_f$|_cat$", "", v)), w = 9, h = 5.5)
})

# ----------------------------------------------------------------------------
# 4. Survey-weighted Cox models, pooled with Rubin's rules
# M1: age (natural spline) + sex
# M2: M1 + race/ethnicity, education, income-to-poverty ratio, smoking, survey cycle
# M3: M2 + all conditions together (glycemic, hypertension, CKD, BMI, CVD, cancer)
# ----------------------------------------------------------------------------
header("4. COX PROPORTIONAL HAZARDS MODELS")

age_term <- "ns(age, knots = c(35, 50, 65), Boundary.knots = c(20, 85))"
m1 <- c(age_term, "sex")
m2 <- c(m1, "race_eth", "education", "pir", "smoking", "cycle")
m3 <- c(m2, "glycemic", "hypertension", "ckd", "bmi_cat", "cvd_history", "cancer")
cox_f <- function(terms, time = "time_y", event = "death") reformulate(unique(terms), sprintf("Surv(%s, %s)", time, event))
fit_mi <- function(f, design = des) with(design, fun = function(d) svycoxph(f, design = d)) %>% mi_tidy()

exposures <- c(glycemic = "Glycemic status", glycemic4 = "Glycemic status (diagnosis)", hypertension = "Hypertension", ckd = "Chronic kidney disease", bmi_cat = "BMI category")

hr_models <- map_dfr(names(exposures), function(e) bind_rows(
  cached(paste0("M1_", e), fit_mi(cox_f(c(e, m1)))) %>% mutate(model = "M1: age, sex"),
  cached(paste0("M2_", e), fit_mi(cox_f(c(e, m2)))) %>% mutate(model = "M2: + sociodemographics, smoking, cycle"),
  if (e != "glycemic4") cached("M3", fit_mi(cox_f(m3))) %>% mutate(model = "M3: + all conditions")) %>%
  filter(startsWith(term, e), !(e == "glycemic" & startsWith(term, "glycemic4"))) %>%
  mutate(exposure = exposures[[e]], level = if_else(term == e, "yes", sub(paste0("^", e), "", term))))

hazard_ratios <- hr_models %>%
  transmute(exposure, level, model, hr = round(exp(est), 2), lcl = round(exp(lcl), 2), ucl = round(exp(ucl), 2), p = fmt_p(p))
show_save(hazard_ratios, "06_hazard_ratios")

cat("\nHazard ratios, wide view (reference: normal glucose, no hypertension, no CKD, normal BMI):\n")
hr_models %>%
  mutate(value = fmt_hr(est, lcl, ucl)) %>%
  select(exposure, level, model, value) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  print(n = Inf, width = Inf)

cat("\nFull M3 model (all non-spline terms):\n")
m3_full <- cache[["M3"]] %>%
  filter(!startsWith(term, "ns(")) %>%
  transmute(term, hr = round(exp(est), 2), lcl = round(exp(lcl), 2), ucl = round(exp(ucl), 2), p = fmt_p(p))
show_save(m3_full, "06_m3_full_model")

p_forest <- hr_models %>%
  mutate(label = paste0(exposure, ": ", level), label = factor(label, rev(unique(label))),
         model = factor(model, unique(hr_models$model))) %>%
  ggplot(aes(exp(est), label, colour = model)) +
  geom_vline(xintercept = 1, colour = "#52514e", linewidth = 0.4) +
  geom_linerange(aes(xmin = exp(lcl), xmax = exp(ucl)), linewidth = 0.6, position = position_dodge(width = 0.7)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.7)) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4)) +
  scale_colour_manual(values = pal[1:3]) +
  labs(title = "All-cause mortality hazard ratios", subtitle = "Survey-weighted Cox models, 95% CI, pooled over imputations (log scale)",
       x = "Hazard ratio", y = NULL, colour = NULL) +
  guides(colour = guide_legend(nrow = 3))
plot_save(p_forest, "06_hazard_ratios_forest", h = 7)

# ----------------------------------------------------------------------------
# 5. Dose-response with natural cubic splines (M2 adjustment)
# The EDA showed U/J shapes, so each measure enters as a spline with knots at
# fixed percentiles; HR is relative to a clinical reference value.
# ----------------------------------------------------------------------------
header("5. DOSE-RESPONSE SPLINES")

spl <- tibble(var = c("hba1c", "bmi", "sbp", "egfr", "hdl", "log_acr"),
              label = c("HbA1c (%)", "BMI (kg/m2)", "Systolic BP (mmHg)", "eGFR (mL/min/1.73m2)", "HDL cholesterol (mg/dL)", "Urine ACR (log10 mg/g)"),
              ref = c(5.4, 25, 120, 95, 55, log(10)))
d1 <- des$designs[[1]]$variables

spline_curves <- pmap_dfr(spl, function(var, label, ref) {
  x <- d1[[var]]
  kn <- round(quantile(x, c(0.35, 0.65)), 3)
  bk <- round(quantile(x, c(0.05, 0.95)), 3)
  term <- sprintf("ns(%s, knots = c(%s), Boundary.knots = c(%s))", var, paste(kn, collapse = ", "), paste(bk, collapse = ", "))
  r <- cached(paste0("spline_", var), MIcombine(with(des, fun = function(d) svycoxph(cox_f(c(term, m2)), design = d))))
  idx <- startsWith(names(coef(r)), paste0("ns(", var))
  xs <- seq(quantile(x, 0.01), quantile(x, 0.99), length.out = 200)
  basis <- ns(x, knots = kn, Boundary.knots = bk)
  X <- predict(basis, xs) - matrix(predict(basis, ref), length(xs), ncol(basis), byrow = TRUE)
  lp <- drop(X %*% coef(r)[idx])
  se <- sqrt(rowSums((X %*% vcov(r)[idx, idx]) * X))
  tibble(variable = var, label, x = if (var == "log_acr") exp(xs) else xs, x_plot = if (var == "log_acr") xs / log(10) else xs,
         hr = exp(lp), lcl = exp(lp - 1.96 * se), ucl = exp(lp + 1.96 * se))
})
write_csv(spline_curves, file.path(out, "06_spline_curves.csv"))

spline_summary <- spline_curves %>%
  group_by(measure = sub("log10 ", "", label)) %>%
  summarise(x_lowest_risk = round(x[which.min(hr)], 1), hr_at_p1 = round(first(hr), 2), hr_at_p99 = round(last(hr), 2))
show_save(spline_summary, "06_spline_summary")

p_spline <- spline_curves %>%
  mutate(label = factor(label, spl$label)) %>%
  ggplot(aes(x_plot, hr)) +
  geom_hline(yintercept = 1, colour = "#52514e", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = pal[1], alpha = 0.15) +
  geom_line(colour = pal[1], linewidth = 0.8) +
  facet_wrap(~label, scales = "free_x", ncol = 3) +
  scale_y_log10() +
  labs(title = "Dose-response: continuous risk factors and all-cause mortality",
       subtitle = "Hazard ratio vs reference value, M2 adjustment, 1st-99th percentile of each measure",
       caption = "Reference: HbA1c 5.4%, BMI 25, SBP 120 mmHg, eGFR 95, HDL 55 mg/dL, ACR 10 mg/g (log10 = 1)",
       x = NULL, y = "Hazard ratio (log scale)")
plot_save(p_spline, "06_spline_dose_response", h = 7)

# ----------------------------------------------------------------------------
# 6. Proportional hazards check (weighted coxph on imputation 1, M3)
# With 50,000 people these tests flag even trivial departures, so the
# time-split models in section 7 show whether any departure changes the story.
# ----------------------------------------------------------------------------
header("6. PROPORTIONAL HAZARDS CHECK")

ph_fit <- coxph(cox_f(m3), data = d1, weights = wt_mec / mean(wt_mec), robust = TRUE)
ph_test <- as.data.frame(cox.zph(ph_fit)$table) %>%
  rownames_to_column("term") %>%
  as_tibble() %>%
  transmute(term, chisq = round(chisq, 1), df, p = fmt_p(p))
show_save(ph_test, "06_ph_test")

# ----------------------------------------------------------------------------
# 7. Sensitivity analyses for the M3 hazard ratios
# CVD/cancer exclusion drops anyone with either condition in any imputation, so
# the subset is identical across imputations. Only 1999-2010 cycles reach 10+
# years of follow-up, so the late-period model uses cycle as a linear term.
# ----------------------------------------------------------------------------
header("7. SENSITIVITY ANALYSES (M3)")

imp <- readRDS(file.path(inp, "05_mice.rds"))
complete_ids <- which(complete.cases(imp$data))
des_cc <- subset(des$designs[[1]], .id %in% complete_ids)
cc_fit <- cached("S_complete_case", svycoxph(cox_f(m3), design = des_cc))
cc <- tibble(term = names(coef(cc_fit)), est = unname(coef(cc_fit)), se = sqrt(diag(vcov(cc_fit)))) %>%
  mutate(lcl = est - 1.96 * se, ucl = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)))

sens <- bind_rows(
  cache[["M3"]] %>% mutate(analysis = "Main: multiple imputation"),
  cc %>% mutate(analysis = sprintf("Complete cases only (n = %s)", format(length(complete_ids), big.mark = ","))),
  cached("S_exclude_2y", fit_mi(cox_f(m3, "I(time_y - 2)"), subset(des, time_y > 2))) %>% mutate(analysis = "Excluding deaths in first 2 years"),
  cached("S_no_cvd_cancer", fit_mi(cox_f(setdiff(m3, c("cvd_history", "cancer"))), subset(des, any_cvd_cancer == 0))) %>% mutate(analysis = "Excluding baseline CVD or cancer"),
  cached("S_no_cycle", fit_mi(cox_f(setdiff(m3, "cycle")))) %>% mutate(analysis = "Without survey-cycle adjustment"),
  cached("S_first_10y", fit_mi(cox_f(m3, "t10", "d10"))) %>% mutate(analysis = "Follow-up years 0-10"),
  cached("S_after_10y", fit_mi(cox_f(c(setdiff(m3, "cycle"), "cycle_num"), "t_late"), subset(des, time_y > 10))) %>% mutate(analysis = "Follow-up beyond 10 years")) %>%
  filter(term %in% c("glycemicprediabetes", "glycemicdiabetes", "hypertension", "ckd", "bmi_catunderweight", "bmi_catoverweight", "bmi_catobese")) %>%
  mutate(term = recode(term, glycemicprediabetes = "Prediabetes", glycemicdiabetes = "Diabetes", hypertension = "Hypertension", ckd = "CKD",
                       bmi_catunderweight = "Underweight", bmi_catoverweight = "Overweight", bmi_catobese = "Obese"))

sensitivity <- sens %>% transmute(analysis, term, hr = round(exp(est), 2), lcl = round(exp(lcl), 2), ucl = round(exp(ucl), 2))
write_csv(sensitivity, file.path(out, "06_sensitivity.csv"))
sens %>%
  mutate(value = fmt_hr(est, lcl, ucl)) %>%
  select(analysis, term, value) %>%
  pivot_wider(names_from = term, values_from = value) %>%
  print(n = Inf, width = Inf)

p_sens <- sens %>%
  mutate(analysis = factor(analysis, rev(unique(sens$analysis))), term = factor(term, unique(sens$term))) %>%
  ggplot(aes(exp(est), analysis)) +
  geom_vline(xintercept = 1, colour = "#52514e", linewidth = 0.4) +
  geom_linerange(aes(xmin = exp(lcl), xmax = exp(ucl)), colour = pal[1], linewidth = 0.6) +
  geom_point(aes(colour = analysis == "Main: multiple imputation"), size = 2.2, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = pal[2], `FALSE` = pal[1])) +
  facet_wrap(~term, nrow = 2, scales = "free_x") +
  scale_x_log10() +
  labs(title = "Sensitivity of the fully adjusted hazard ratios", subtitle = "M3 model, 95% CI (log scale); orange = main analysis", x = "Hazard ratio", y = NULL) +
  theme(panel.spacing.x = unit(1.2, "lines"))
plot_save(p_sens, "06_sensitivity_forest", w = 12, h = 7)

# ----------------------------------------------------------------------------
# 8. Summary
# ----------------------------------------------------------------------------
header("8. SUMMARY")

hr_models %>%
  filter(model == "M3: + all conditions" | exposure == "Glycemic status (diagnosis)" & startsWith(model, "M2")) %>%
  transmute(exposure, level, model = sub(":.*", "", model), hr_95ci = fmt_hr(est, lcl, ucl)) %>%
  print(n = Inf, width = Inf)

writeLines(capture.output(sessionInfo()), file.path(out, "06_session_info.txt"))
cat(sprintf("\n✓ All outputs saved in: %s\n", out))
print(list.files(out))
