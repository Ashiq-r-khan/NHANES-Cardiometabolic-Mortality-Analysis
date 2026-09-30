# ============================================================================
# Script 07: Competing risks, cause-specific mortality
# CVD (heart disease + stroke), cancer, and other causes of death
# Input : outputs/R/05_survey_design_imputation/05_imputed_long.rds
# Output: outputs/R/07_competing_risks/
# ============================================================================

library(tidyverse)
library(survey)
library(mitools)
library(splines)

options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE, scipen = 999)

# Set paths
root <- "E:/NHANES Cardiometabolic Mortality Project"
inp <- file.path(root, "outputs", "R", "05_survey_design_imputation")
out <- file.path(root, "outputs", "R", "07_competing_risks")
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
causes <- c(cvd = "CVD", cancer = "Cancer", other = "Other causes")

pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
                  plot.title = element_text(face = "bold", colour = "#0b0b0b"), plot.subtitle = element_text(colour = "#52514e"),
                  axis.text = element_text(colour = "#52514e"), strip.text = element_text(face = "bold", hjust = 0, colour = "#0b0b0b"),
                  legend.position = "top", legend.justification = "left", plot.background = element_rect(fill = "#fcfcfb", colour = NA)))

# Model results are cached in 07_model_cache.rds as each one finishes
# (first run ~10-15 min). Delete the file to refit everything.
cache_file <- file.path(out, "07_model_cache.rds")
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
# 1. Load imputed data
# cause (from notebook 02): cvd = NCHS leading cause codes 1 (heart disease)
# and 5 (cerebrovascular), cancer = code 2, other = all remaining codes.
# ----------------------------------------------------------------------------
header("1. LOAD DATA")

ucod_labels <- c(`1` = "Heart disease", `2` = "Cancer", `3` = "Chronic lower respiratory disease", `4` = "Accidents", `5` = "Stroke",
                 `6` = "Alzheimer disease", `7` = "Diabetes", `8` = "Influenza and pneumonia", `9` = "Kidney disease", `10` = "All other causes")

long <- readRDS(file.path(inp, "05_imputed_long.rds")) %>%
  mutate(time_y = time_months / 12,
         cause_f = factor(as.character(cause), c("alive", names(causes))),
         cause_death = factor(if_else(death == 1, as.character(cause), NA_character_), names(causes), causes),
         ucod = factor(ucod_labels[as.character(ucod_leading)], ucod_labels),
         d_cvd = as.numeric(cause == "cvd"), d_cancer = as.numeric(cause == "cancer"), d_other = as.numeric(cause == "other"),
         glycemic4 = factor(case_when(glycemic == "diabetes" & dm_dx == 1 ~ "diagnosed diabetes", glycemic == "diabetes" ~ "undiagnosed diabetes", TRUE ~ as.character(glycemic)),
                            c("normal", "prediabetes", "undiagnosed diabetes", "diagnosed diabetes")),
         bmi_cat = factor(bmi_cat, c("normal", "underweight", "overweight", "obese")),
         hypertension_f = factor(hypertension, 0:1, c("no hypertension", "hypertension")),
         ckd_f = factor(ckd, 0:1, c("no CKD", "CKD")))

des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt_mec, nest = TRUE, data = imputationList(split(long, long$.imp)))
d1 <- filter(long, .imp == 1)
print(table(d1$cause_f))

groups <- c(glycemic4 = "Glycemic status", hypertension_f = "Hypertension", ckd_f = "Chronic kidney disease")

# ----------------------------------------------------------------------------
# 2. What people die of: weighted cause-of-death mix
# ----------------------------------------------------------------------------
header("2. CAUSE-OF-DEATH MIX (WEIGHTED % OF DEATHS)")

des_deaths <- subset(des, death == 1)

cause_mix <- map_dfr(names(groups), function(v)
  with(des_deaths, fun = function(d) svyby(~cause_death, reformulate(v), d, svymean)) %>% mi_tidy() %>%
    separate(term, c("level", "cause"), sep = ":") %>%
    mutate(group = groups[[v]], cause = sub("^cause_death", "", cause))) %>%
  transmute(group, level, cause, pct = round(100 * est, 1), lcl = round(100 * lcl, 1), ucl = round(100 * ucl, 1))
write_csv(cause_mix, file.path(out, "07_cause_mix.csv"))
cause_mix %>% select(group, level, cause, pct) %>% pivot_wider(names_from = cause, values_from = pct) %>% print(n = Inf, width = Inf)

cat("\nDetailed leading cause of death by glycemic status (weighted % of deaths):\n")
cause_detail <- with(subset(des_deaths, !is.na(ucod)), fun = function(d) svyby(~ucod, ~glycemic4, d, svymean)) %>% mi_tidy() %>%
  separate(term, c("glycemic_status", "cause"), sep = ":") %>%
  transmute(glycemic_status, cause = sub("^ucod", "", cause), pct = round(100 * est, 1))
write_csv(cause_detail, file.path(out, "07_cause_detail.csv"))
cause_detail %>% pivot_wider(names_from = glycemic_status, values_from = pct) %>% print(n = Inf, width = Inf)

cat("\nDiabetes listed anywhere on the death certificate (weighted % of deaths):\n")
dm_on_certificate <- with(subset(des_deaths, !is.na(mcod_diabetes)), fun = function(d) svyby(~mcod_diabetes, ~glycemic4, d, svymean)) %>% mi_tidy() %>%
  transmute(glycemic_status = term, pct_with_diabetes_on_certificate = round(100 * est, 1), lcl = round(100 * lcl, 1), ucl = round(100 * ucl, 1))
show_save(dm_on_certificate, "07_diabetes_on_death_certificate")

p_mix <- cause_mix %>%
  mutate(group = factor(group, groups), cause = factor(cause, rev(causes)), level = fct_rev(factor(level, unique(level)))) %>%
  ggplot(aes(pct, level, fill = cause)) +
  geom_col(width = 0.65, colour = "#fcfcfb", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_stack(vjust = 0.5), colour = "white", size = 3.4, fontface = "bold") +
  facet_wrap(~group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = setNames(pal[1:3], causes), breaks = causes) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "What people die of, by cardiometabolic status", subtitle = "Survey-weighted share of deaths by underlying cause, pooled over imputations", x = "% of deaths", y = NULL, fill = NULL)
plot_save(p_mix, "07_cause_mix", h = 7)

# ----------------------------------------------------------------------------
# 3. Cause-specific death rates per 1,000 person-years (age-standardized)
# ----------------------------------------------------------------------------
header("3. AGE-STANDARDIZED CAUSE-SPECIFIC DEATH RATES PER 1,000 PERSON-YEARS")

cause_rates <- map_dfr(names(groups), function(v) map_dfr(names(causes), function(k) {
  lv <- levels(long[[v]])
  with(des, fun = function(d) {
    r <- svyby(reformulate(paste0("d_", k)), reformulate(c(v, "age3")), denominator = ~time_y, design = d, FUN = svyratio, covmat = TRUE)
    svycontrast(r, setNames(lapply(lv, function(l) setNames(std_prop, paste(l, names(std_prop), sep = "."))), lv))
  }) %>% mi_tidy() %>% mutate(group = groups[[v]], level = term, cause = causes[[k]])
})) %>%
  transmute(group, level, cause, rate_per_1000_py = round(1000 * est, 2), lcl = round(1000 * lcl, 2), ucl = round(1000 * ucl, 2))
write_csv(cause_rates, file.path(out, "07_cause_specific_rates.csv"))
cause_rates %>% select(group, level, cause, rate_per_1000_py) %>% pivot_wider(names_from = cause, values_from = rate_per_1000_py) %>% print(n = Inf, width = Inf)

# ----------------------------------------------------------------------------
# 4. Cumulative incidence (weighted Aalen-Johansen), averaged over imputations
# 1 - Kaplan-Meier treats other deaths as censored and overstates the risk of
# a specific cause; the Aalen-Johansen estimator accounts for competing deaths.
# ----------------------------------------------------------------------------
header("4. CUMULATIVE INCIDENCE OF CAUSE-SPECIFIC DEATH")

grid <- seq(0, 18, by = 0.25)
cif <- map_dfr(names(groups), function(v) map_dfr(split(long, long$.imp), function(d) {
  s <- summary(survfit(reformulate(v, "Surv(time_y, cause_f)"), data = d, weights = wt_mec), times = grid, extend = TRUE)
  km <- map_dfr(names(causes), function(k) {
    sk <- summary(survfit(reformulate(v, sprintf("Surv(time_y, d_%s)", k)), data = d, weights = wt_mec), times = grid, extend = TRUE)
    tibble(level = sub(".*=", "", sk$strata), time = sk$time, cause = causes[[k]], naive = 1 - sk$surv)
  })
  as_tibble(s$pstate[, -1, drop = FALSE], .name_repair = ~ unname(causes)) %>%
    mutate(level = sub(".*=", "", s$strata), time = s$time) %>%
    pivot_longer(all_of(unname(causes)), names_to = "cause", values_to = "cif") %>%
    left_join(km, by = c("level", "time", "cause"))
}) %>%
  group_by(level, time, cause) %>% summarise(across(c(cif, naive), mean)) %>% ungroup() %>%
  mutate(group = groups[[v]], level = factor(level, levels(long[[v]]))))
write_csv(cif, file.path(out, "07_cumulative_incidence.csv"))

cif_15 <- cif %>%
  filter(time == 15) %>%
  transmute(group, level, cause, aalen_johansen_pct = round(100 * cif, 1), naive_1_minus_km_pct = round(100 * naive, 1),
            overestimate_pct = round(100 * (naive / cif - 1), 1)) %>%
  arrange(factor(group, groups), level, factor(cause, causes))
show_save(cif_15, "07_cumulative_incidence_15y")

walk2(names(groups), groups, function(v, g) {
  p <- cif %>%
    filter(group == g) %>%
    mutate(cause = factor(cause, causes)) %>%
    ggplot(aes(time, 100 * cif, colour = level)) +
    geom_step(linewidth = 0.8) +
    facet_wrap(~cause, ncol = 3) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(breaks = seq(0, 18, 6)) +
    labs(title = paste("Cumulative incidence of death by cause and", tolower(g)),
         subtitle = "Survey-weighted Aalen-Johansen estimates accounting for competing causes of death",
         x = "Years of follow-up", y = "Cumulative incidence (%)", colour = NULL)
  plot_save(p, paste0("07_cif_", sub("4$|_f$", "", v)), w = 11, h = 5)
})

# ----------------------------------------------------------------------------
# 5. Cause-specific Cox models (competing deaths censored), pooled over 20 imputations
# Same adjustment sets as script 06. Cause-specific HRs answer the aetiological
# question: how much faster does each cause kill people with the condition.
# ----------------------------------------------------------------------------
header("5. CAUSE-SPECIFIC COX MODELS")

age_term <- "ns(age, knots = c(35, 50, 65), Boundary.knots = c(20, 85))"
m2 <- c(age_term, "sex", "race_eth", "education", "pir", "smoking", "cycle")
m3 <- c(m2, "glycemic", "hypertension", "ckd", "bmi_cat", "cvd_history", "cancer")
cox_f <- function(terms, time = "time_y", event = "death") reformulate(unique(terms), sprintf("Surv(%s, %s)", time, event))
fit_mi <- function(f, design = des) with(design, fun = function(d) svycoxph(f, design = d)) %>% mi_tidy()
keep_terms <- c(glycemicprediabetes = "Prediabetes", glycemicdiabetes = "Diabetes", `glycemic4undiagnosed diabetes` = "Undiagnosed diabetes (M2)",
                `glycemic4diagnosed diabetes` = "Diagnosed diabetes (M2)", hypertension = "Hypertension", ckd = "CKD",
                bmi_catunderweight = "Underweight", bmi_catoverweight = "Overweight", bmi_catobese = "Obese", smokingcurrent = "Current smoking")

cs <- map_dfr(names(causes), function(k) bind_rows(
  cached(paste0("CS_M3_", k), fit_mi(cox_f(m3, event = paste0("d_", k)))),
  cached(paste0("CS_M2_glycemic4_", k), fit_mi(cox_f(c("glycemic4", m2), event = paste0("d_", k)))) %>% filter(startsWith(term, "glycemic4"))) %>%
  filter(term %in% names(keep_terms)) %>%
  mutate(cause = causes[[k]], term = keep_terms[term])) %>%
  bind_rows(readRDS(file.path(root, "outputs", "R", "06_survival_cox", "06_cox_cache.rds"))[c("M3", "M2_glycemic4")] %>%
              map_dfr(function(x) filter(x, term %in% names(keep_terms))) %>%
              filter(!(duplicated(term) & !startsWith(term, "glycemic4"))) %>%
              mutate(cause = "All causes", term = keep_terms[term]))

cs_hr <- cs %>% transmute(cause, term, hr = round(exp(est), 2), lcl = round(exp(lcl), 2), ucl = round(exp(ucl), 2), p = fmt_p(p))
write_csv(cs_hr, file.path(out, "07_cause_specific_hr.csv"))
cs %>%
  mutate(value = fmt_hr(est, lcl, ucl), cause = factor(cause, c("All causes", causes)), term = factor(term, keep_terms)) %>%
  select(term, cause, value) %>% arrange(cause) %>%
  pivot_wider(names_from = cause, values_from = value) %>% arrange(term) %>%
  print(n = Inf, width = Inf)

p_cs <- cs %>%
  mutate(cause = factor(cause, c("All causes", causes)), term = factor(term, rev(keep_terms))) %>%
  ggplot(aes(exp(est), term, colour = cause)) +
  geom_vline(xintercept = 1, colour = "#52514e", linewidth = 0.4) +
  geom_linerange(aes(xmin = exp(lcl), xmax = exp(ucl)), linewidth = 0.6, position = position_dodge(width = 0.75)) +
  geom_point(size = 2, position = position_dodge(width = 0.75)) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4)) +
  scale_colour_manual(values = c(`All causes` = "#52514e", setNames(pal[1:3], causes))) +
  labs(title = "Cause-specific hazard ratios", subtitle = "Survey-weighted Cox, M3 adjustment (diagnosis split: M2), 95% CI, 20 imputations (log scale)",
       x = "Hazard ratio", y = NULL, colour = NULL)
plot_save(p_cs, "07_cause_specific_hr_forest", h = 7.5)

# ----------------------------------------------------------------------------
# 6. Fine-Gray subdistribution hazard ratios (CVD and cancer)
# Answers the prediction question: how much more likely is death from this cause
# in the real world where other causes also kill. finegray() expands the data,
# so follow-up is rounded up to quarter-years (HRs unchanged to 3 decimals,
# rows 890k -> 350k) and the first 5 imputations are used (Rubin's rules
# remain valid with m = 5).
# ----------------------------------------------------------------------------
header("6. FINE-GRAY SUBDISTRIBUTION HAZARD RATIOS (5 IMPUTATIONS)")

fg_vars <- c("time_y", "cause_f", "age", "sex", "race_eth", "education", "pir", "smoking", "cycle", "glycemic", "hypertension", "ckd", "bmi_cat",
             "cvd_history", "cancer", "psu", "strata", "wt_mec")
fg_fit <- function(i, k) {
  fgd <- finegray(Surv(time_y, cause_f) ~ ., data = filter(long, .imp == i) %>% mutate(time_y = ceiling(time_months / 3) / 4) %>% select(all_of(fg_vars)), etype = k)
  fit <- svycoxph(cox_f(m3, "fgstart, fgstop", "fgstatus"), design = svydesign(ids = ~psu, strata = ~strata, weights = ~I(wt_mec * fgwt), nest = TRUE, data = fgd))
  res <- list(coef = coef(fit), vcov = vcov(fit))
  rm(fgd, fit); gc()
  cat(sprintf("    imputation %d done\n", i))
  res
}
fg <- map_dfr(c("cvd", "cancer"), function(k) cached(paste0("FG_", k), {
  fits <- map(1:5, fg_fit, k = k)
  MIcombine(map(fits, "coef"), map(fits, "vcov"))
}) %>%
  { tibble(term = names(coef(.)), est = unname(coef(.)), se = sqrt(unname(diag(vcov(.))))) } %>%
  mutate(lcl = est - 1.96 * se, ucl = est + 1.96 * se, cause = causes[[k]])) %>%
  filter(term %in% names(keep_terms)) %>%
  mutate(term = keep_terms[term])

cs_vs_fg <- bind_rows(cs %>% filter(cause %in% causes[c("cvd", "cancer")], !grepl("M2", term)) %>% mutate(model = "Cause-specific HR"),
                      fg %>% mutate(model = "Subdistribution HR (Fine-Gray)")) %>%
  mutate(value = fmt_hr(est, lcl, ucl)) %>%
  select(cause, term, model, value) %>%
  pivot_wider(names_from = model, values_from = value)
show_save(cs_vs_fg, "07_cause_specific_vs_finegray")

# ----------------------------------------------------------------------------
# 7. Summary
# ----------------------------------------------------------------------------
header("7. SUMMARY")

cat("\n15-year cumulative incidence of CVD death (Aalen-Johansen) by glycemic status:\n")
cif_15 %>% filter(cause == "CVD", group == "Glycemic status") %>% print(width = Inf)
cat("\nCause-specific HR for diabetes (M3) by cause:\n")
cs %>% filter(term == "Diabetes") %>% transmute(cause, hr_95ci = fmt_hr(est, lcl, ucl)) %>% print(width = Inf)

writeLines(capture.output(sessionInfo()), file.path(out, "07_session_info.txt"))
cat(sprintf("\n✓ All outputs saved in: %s\n", out))
print(list.files(out))
