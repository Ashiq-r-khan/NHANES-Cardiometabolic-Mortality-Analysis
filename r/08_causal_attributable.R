# ============================================================================
# Script 08: Causal inference and attributable burden
# IPW marginal hazard ratios, covariate balance, E-values, population
# attributable fractions (PAF) for all-cause and CVD mortality
# Input : outputs/R/05_survey_design_imputation/05_imputed_long.rds
#         outputs/R/06_survival_cox/06_cox_cache.rds (M3 hazard ratios)
#         outputs/R/07_competing_risks/07_model_cache.rds (CVD hazard ratios)
# Output: outputs/R/08_causal_attributable/
# ============================================================================

library(tidyverse)
library(survey)
library(mitools)
library(splines)

options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE, scipen = 999)

# Set paths
root <- "E:/NHANES Cardiometabolic Mortality Project"
inp <- file.path(root, "outputs", "R", "05_survey_design_imputation")
out <- file.path(root, "outputs", "R", "08_causal_attributable")
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
header <- function(x) cat("\n============================================================\n", x, "\n============================================================\n", sep = "")
fmt_hr <- function(est, lcl, ucl) sprintf("%.2f (%.2f-%.2f)", exp(est), exp(lcl), exp(ucl))
wmean <- function(x, w) sum(w * x) / sum(w)
wvar <- function(x, w) sum(w * (x - wmean(x, w))^2) / sum(w)

pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
                  plot.title = element_text(face = "bold", colour = "#0b0b0b"), plot.subtitle = element_text(colour = "#52514e"),
                  axis.text = element_text(colour = "#52514e"), strip.text = element_text(face = "bold", hjust = 0, colour = "#0b0b0b"),
                  legend.position = "top", legend.justification = "left", plot.background = element_rect(fill = "#fcfcfb", colour = NA)))

# ----------------------------------------------------------------------------
# 1. Load data and define exposures with their confounder sets
# Confounders are common causes of exposure and death. Downstream conditions
# (mediators) are left out, e.g. hypertension and CKD are not adjusted for
# when the exposure is diabetes, because diabetes causes them.
# ----------------------------------------------------------------------------
header("1. LOAD DATA AND DEFINE EXPOSURES")

long <- readRDS(file.path(inp, "05_imputed_long.rds")) %>%
  mutate(time_y = time_months / 12, bmi_cat = factor(bmi_cat, c("normal", "underweight", "overweight", "obese")))
imps <- split(long, long$.imp)

age_term <- "ns(age, knots = c(35, 50, 65), Boundary.knots = c(20, 85))"
bmi_term <- "ns(bmi, knots = c(25, 30), Boundary.knots = c(19, 42))"
base <- c(age_term, "sex", "race_eth", "education", "pir", "smoking", "cycle", bmi_term)
spec <- list(diabetes = list(label = "Diabetes", conf = base),
             hypertension = list(label = "Hypertension", conf = c(base, "diabetes")),
             ckd = list(label = "Chronic kidney disease", conf = c(base, "diabetes", "hypertension")))
walk(names(spec), function(e) cat(sprintf("%-13s confounders: %s\n", e, paste(sub(",.*", ")", spec[[e]]$conf), collapse = ", "))))

# ----------------------------------------------------------------------------
# 2. Inverse probability weights, per imputation
# Propensity score from a survey-weighted logistic model; stabilized weights
# P(A = a) / P(A = a | L) multiplied by the survey weight. Weights are not
# truncated in the main analysis: truncation left age imbalance for diabetes
# (section 4 shows the trade-off), and the largest weights stay moderate.
# ----------------------------------------------------------------------------
header("2. PROPENSITY SCORES AND STABILIZED WEIGHTS")

add_ipw <- function(d, e) {
  f <- reformulate(spec[[e]]$conf, e)
  ps <- fitted(svyglm(f, svydesign(ids = ~psu, strata = ~strata, weights = ~wt_mec, nest = TRUE, data = d), family = quasibinomial()))
  pa <- wmean(d[[e]], d$wt_mec)
  sw <- if_else(d[[e]] == 1, pa / ps, (1 - pa) / (1 - ps))
  mutate(d, ps = ps, sw = sw, w_ipw = wt_mec * sw)
}
ipw <- map(names(spec), function(e) map(imps, add_ipw, e = e)) %>% set_names(names(spec))

weight_summary <- map_dfr(names(spec), function(e) {
  d <- ipw[[e]][[1]]
  tibble(exposure = spec[[e]]$label, exposed_pct_weighted = round(100 * wmean(d[[e]], d$wt_mec), 1),
         ps_min = round(min(d$ps), 3), ps_max = round(max(d$ps), 3),
         sw_mean = round(mean(d$sw), 3), sw_p1 = round(quantile(d$sw, 0.01), 2), sw_p99 = round(quantile(d$sw, 0.99), 2), sw_max = round(max(d$sw), 1))
})
cat("Stabilized weights should average close to 1 (imputation 1):\n")
show_save(weight_summary, "08_weight_summary")

p_overlap <- map_dfr(names(spec), function(e) ipw[[e]][[1]] %>% transmute(exposure = spec[[e]]$label, group = if_else(.data[[e]] == 1, "Exposed", "Unexposed"), ps, wt_mec)) %>%
  group_by(exposure, group) %>% mutate(w = wt_mec / sum(wt_mec)) %>% ungroup() %>%
  ggplot(aes(ps, weight = w, fill = group)) +
  geom_density(alpha = 0.45, colour = NA, bw = 0.02) +
  facet_wrap(~exposure, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Exposed = pal[2], Unexposed = pal[1])) +
  labs(title = "Propensity score overlap", subtitle = "Survey-weighted distribution of the propensity score by exposure status (imputation 1)",
       x = "Propensity score", y = NULL, fill = NULL) +
  theme(axis.text.y = element_blank())
plot_save(p_overlap, "08_ps_overlap", h = 4.5)

# ----------------------------------------------------------------------------
# 3. Covariate balance: standardized mean differences before and after IPW
# |SMD| < 0.1 is the usual threshold for adequate balance.
# ----------------------------------------------------------------------------
header("3. COVARIATE BALANCE")

smd_tbl <- function(d, e, w) {
  X <- model.matrix(reformulate(c("age", "bmi", "sex", "race_eth", "education", "pir", "smoking", setdiff(spec[[e]]$conf, base))), d)[, -1]
  a <- d[[e]] == 1
  tibble(covariate = colnames(X), smd = map_dbl(seq_len(ncol(X)), function(j)
    (wmean(X[a, j], w[a]) - wmean(X[!a, j], w[!a])) / sqrt((wvar(X[a, j], w[a]) + wvar(X[!a, j], w[!a])) / 2)))
}
balance <- map_dfr(names(spec), function(e) {
  d <- ipw[[e]][[1]]
  bind_rows(before = smd_tbl(d, e, d$wt_mec), after = smd_tbl(d, e, d$w_ipw), .id = "weighting") %>% mutate(exposure = spec[[e]]$label)
})
write_csv(balance, file.path(out, "08_covariate_balance.csv"))

balance %>%
  group_by(exposure, weighting) %>%
  summarise(max_abs_smd = round(max(abs(smd)), 3), n_above_0.1 = sum(abs(smd) > 0.1)) %>%
  pivot_wider(names_from = weighting, values_from = c(max_abs_smd, n_above_0.1)) %>%
  print(width = Inf)

p_love <- balance %>%
  mutate(weighting = factor(if_else(weighting == "before", "Survey weight only", "Survey weight x IPW"), c("Survey weight only", "Survey weight x IPW")),
         covariate = fct_reorder(covariate, abs(smd), max)) %>%
  ggplot(aes(abs(smd), covariate, colour = weighting)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", colour = "#52514e", linewidth = 0.4) +
  geom_point(size = 2) +
  facet_wrap(~exposure, ncol = 3) +
  scale_colour_manual(values = pal[2:1]) +
  labs(title = "Covariate balance before and after inverse probability weighting",
       subtitle = "Absolute standardized mean difference; dashed line = 0.1 threshold (imputation 1)", x = "|SMD|", y = NULL, colour = NULL)
plot_save(p_love, "08_love_plot", w = 12, h = 7)

# ----------------------------------------------------------------------------
# 4. IPW marginal hazard ratios and adjusted survival, pooled over imputations
# The marginal HR compares the whole population as if everyone vs no one had
# the exposure. It differs from the conditional M3 HR of script 06, which also
# adjusts for mediators.
# ----------------------------------------------------------------------------
header("4. IPW MARGINAL HAZARD RATIOS AND ADJUSTED SURVIVAL")

m3 <- readRDS(file.path(root, "outputs", "R", "06_survival_cox", "06_cox_cache.rds"))[["M3"]]
truncate <- function(sw, q) if (q == 0) sw else pmin(pmax(sw, quantile(sw, q)), quantile(sw, 1 - q))

ipw_all <- map_dfr(names(spec), function(e) map_dfr(c(0, 0.001, 0.005, 0.01), function(q) {
  r <- MIcombine(map(ipw[[e]], function(d) svycoxph(reformulate(e, "Surv(time_y, death)"),
                     design = svydesign(ids = ~psu, strata = ~strata, weights = ~w, nest = TRUE, data = mutate(d, w = wt_mec * truncate(sw, q))))))
  d1 <- ipw[[e]][[1]]
  tibble(exposure = spec[[e]]$label, term = e, truncation = q, est = coef(r)[[1]], se = sqrt(vcov(r)[1, 1]),
         max_abs_smd = max(abs(smd_tbl(d1, e, d1$wt_mec * truncate(d1$sw, q))$smd)))
})) %>% mutate(lcl = est - 1.96 * se, ucl = est + 1.96 * se)

cat("Weight truncation sensitivity (truncation = share cut from each tail; max |SMD| from imputation 1):\n")
ipw_trunc <- ipw_all %>% transmute(exposure, truncation = paste0(100 * truncation, "%"), hr_95ci = fmt_hr(est, lcl, ucl), max_abs_smd = round(max_abs_smd, 3))
show_save(ipw_trunc, "08_ipw_truncation_sensitivity")

ipw_hr <- ipw_all %>%
  filter(truncation == 0) %>%
  mutate(cond = map(term, function(e) filter(m3, term == if_else(e == "diabetes", "glycemicdiabetes", e))),
         `IPW marginal HR` = fmt_hr(est, lcl, ucl), `Conditional HR (06, M3)` = map_chr(cond, function(x) fmt_hr(x$est, x$lcl, x$ucl))) %>%
  select(-cond)
cat("\nMain analysis (untruncated weights):\n")
show_save(ipw_hr %>% select(exposure, `IPW marginal HR`, `Conditional HR (06, M3)`), "08_ipw_hazard_ratios")

grid <- seq(0, 18, by = 0.1)
adj_surv <- map_dfr(names(spec), function(e) map_dfr(ipw[[e]], function(d) {
  s <- summary(survfit(reformulate(e, "Surv(time_y, death)"), data = d, weights = w_ipw), times = grid, extend = TRUE)
  tibble(status = if_else(grepl("=1", s$strata), "Exposed", "Unexposed"), time = s$time, surv = s$surv)
}) %>% group_by(status, time) %>% summarise(surv = mean(surv)) %>% ungroup() %>% mutate(exposure = spec[[e]]$label))
write_csv(adj_surv, file.path(out, "08_ipw_adjusted_survival.csv"))

risk_10y <- adj_surv %>%
  filter(time == 10) %>%
  transmute(exposure, status, risk = 100 * (1 - surv)) %>%
  pivot_wider(names_from = status, values_from = risk) %>%
  mutate(across(c(Exposed, Unexposed), function(x) round(x, 1)), risk_difference_pct_points = Exposed - Unexposed) %>%
  select(exposure, risk_unexposed_10y_pct = Unexposed, risk_exposed_10y_pct = Exposed, risk_difference_pct_points)
cat("\nIPW-adjusted 10-year risk of death (point estimates averaged over imputations):\n")
show_save(risk_10y, "08_ipw_10y_risk")

p_adj <- adj_surv %>%
  ggplot(aes(time, 100 * surv, colour = status)) +
  geom_step(linewidth = 0.8) +
  facet_wrap(~exposure, ncol = 3) +
  scale_colour_manual(values = c(Exposed = pal[2], Unexposed = pal[1])) +
  scale_x_continuous(breaks = seq(0, 18, 6)) +
  labs(title = "IPW-adjusted survival curves", subtitle = "Survey x inverse probability weighted Kaplan-Meier, averaged over imputations",
       x = "Years of follow-up", y = "Survival (%)", colour = NULL)
plot_save(p_adj, "08_ipw_adjusted_survival", w = 11, h = 4.5)

# ----------------------------------------------------------------------------
# 5. E-values: how strong would unmeasured confounding need to be?
# The E-value is the minimum strength of association (risk ratio scale) that an
# unmeasured confounder would need with both exposure and death to fully
# explain away the HR. Deaths are common (16%), so HRs are first converted to
# risk ratios with VanderWeele's approximation.
# ----------------------------------------------------------------------------
header("5. E-VALUES")

hr_to_rr <- function(hr) (1 - 0.5^sqrt(hr)) / (1 - 0.5^sqrt(1 / hr))
evalue <- function(rr) { rr <- if_else(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
ci_bound <- function(l, u) case_when(l > 1 ~ l, u < 1 ~ u, TRUE ~ 1)

evalues <- bind_rows(
  ipw_hr %>% transmute(estimate = paste(exposure, "(IPW marginal)"), hr = exp(est), lcl = exp(lcl), ucl = exp(ucl)),
  m3 %>% filter(term %in% c("glycemicdiabetes", "hypertension", "ckd", "bmi_catunderweight", "smokingcurrent")) %>%
    transmute(estimate = paste(recode(term, glycemicdiabetes = "Diabetes", hypertension = "Hypertension", ckd = "Chronic kidney disease",
                                      bmi_catunderweight = "Underweight", smokingcurrent = "Current smoking"), "(M3 conditional)"),
              hr = exp(est), lcl = exp(lcl), ucl = exp(ucl))) %>%
  mutate(e_value_point = round(evalue(hr_to_rr(hr)), 2), e_value_ci = round(evalue(hr_to_rr(ci_bound(lcl, ucl))), 2),
         across(c(hr, lcl, ucl), function(x) round(x, 2)))
show_save(evalues, "08_e_values")

# ----------------------------------------------------------------------------
# 6. Population attributable fractions (Miettinen's formula)
# PAF = p_c x (HR - 1) / HR, where p_c is the weighted share of deaths that
# occurred in exposed people and HR is the adjusted HR (M3 for all-cause,
# cause-specific M3 for CVD). PAF is the share of deaths that would not have
# occurred without the exposure, if the association were causal.
# ----------------------------------------------------------------------------
header("6. POPULATION ATTRIBUTABLE FRACTIONS")

cs_cvd <- readRDS(file.path(root, "outputs", "R", "07_competing_risks", "07_model_cache.rds"))[["CS_M3_cvd"]]
paf_exposures <- tribble(
  ~term, ~label, ~indicator,
  "glycemicdiabetes", "Diabetes", "diabetes == 1",
  "hypertension", "Hypertension", "hypertension == 1",
  "ckd", "Chronic kidney disease", "ckd == 1",
  "smokingcurrent", "Current smoking", "smoking == 'current'",
  "smokingformer", "Former smoking", "smoking == 'former'",
  "bmi_catunderweight", "Underweight", "bmi_cat == 'underweight'")

paf <- pmap_dfr(paf_exposures, function(term, label, indicator) map_dfr(list(`All-cause` = list(m3, "death"), CVD = list(cs_cvd, "d_cvd")), function(x) {
  h <- filter(x[[1]], term == !!term)
  pc <- mean(map_dbl(imps, function(d) {
    dd <- filter(mutate(d, d_cvd = as.numeric(cause == "cvd")), .data[[x[[2]]]] == 1)
    wmean(as.numeric(eval(parse(text = indicator), dd)), dd$wt_mec)
  }))
  tibble(exposure = label, hr = exp(h$est), hr_lcl = exp(h$lcl), hr_ucl = exp(h$ucl), pct_deaths_exposed = 100 * pc) %>%
    mutate(paf_pct = pct_deaths_exposed * (hr - 1) / hr, paf_lcl = pct_deaths_exposed * (hr_lcl - 1) / hr_lcl, paf_ucl = pct_deaths_exposed * (hr_ucl - 1) / hr_ucl)
}, .id = "outcome")) %>%
  mutate(across(where(is.numeric), function(x) round(x, 2)))
show_save(paf, "08_population_attributable_fractions")

p_paf <- paf %>%
  mutate(exposure = fct_reorder(exposure, paf_pct, max), outcome = factor(outcome, c("All-cause", "CVD"))) %>%
  ggplot(aes(paf_pct, exposure)) +
  geom_col(fill = pal[1], width = 0.6) +
  geom_linerange(aes(xmin = paf_lcl, xmax = paf_ucl), colour = "#0b0b0b", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", paf_pct), x = pmax(paf_ucl, 0) + 0.6), hjust = 0, size = 3.4, colour = "#52514e") +
  facet_wrap(~outcome, ncol = 2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Share of deaths attributable to each risk factor",
       subtitle = "Population attributable fraction (Miettinen) with fully adjusted HRs, US adults 20+",
       caption = "95% CI from the hazard ratio interval",
       x = "Population attributable fraction (%)", y = NULL)
plot_save(p_paf, "08_paf", h = 5)

# ----------------------------------------------------------------------------
# 7. Summary
# ----------------------------------------------------------------------------
header("7. SUMMARY")

ipw_hr %>% select(exposure, `IPW marginal HR`, `Conditional HR (06, M3)`) %>%
  left_join(risk_10y, by = "exposure") %>%
  left_join(evalues %>% filter(grepl("IPW", estimate)) %>% transmute(exposure = sub(" \\(IPW marginal\\)", "", estimate), e_value_point, e_value_ci), by = "exposure") %>%
  print(width = Inf)
paf %>% filter(outcome == "All-cause") %>% select(exposure, paf_pct, paf_lcl, paf_ucl) %>% arrange(desc(paf_pct)) %>% print(width = Inf)

writeLines(capture.output(sessionInfo()), file.path(out, "08_session_info.txt"))
cat(sprintf("\n✓ All outputs saved in: %s\n", out))
print(list.files(out))
