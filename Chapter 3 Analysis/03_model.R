############################################################
# 03_model.R
# Purpose: Compute prognostic indices, fit Cox models,
#          export model summaries, and create risk tables.
############################################################

##############################
# 1. PACKAGE SETUP
##############################
required_packages <- c("dplyr", "survival", "sjPlot", "broom", "openxlsx", "tableone")
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) {
  stop("Please install: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

##############################
# 2. CONFIGURATION
##############################
config <- list(
  processed_dir = file.path("data", "processed"),
  output_dir = "outputs",
  followup_years = 10
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

complication_path <- file.path(config$processed_dir, "analysis_dataset_complication.rds")
no_complication_path <- file.path(config$processed_dir, "analysis_dataset_no_complication.rds")
cleaned_cohort_path <- file.path(config$processed_dir, "cleaned_cohort.rds")

stopifnot(file.exists(complication_path), file.exists(no_complication_path), file.exists(cleaned_cohort_path))

##############################
# 3. HELPER FUNCTIONS
##############################
get_c_index_with_ci <- function(fit) {
  c_index <- summary(fit)$concordance[1]
  se_c <- summary(fit)$concordance[2]
  c(
    c_index = c_index,
    lower_95 = c_index - 1.96 * se_c,
    upper_95 = c_index + 1.96 * se_c
  )
}

make_risk_table <- function(data, group_var, outcome_var = "outcome") {
  group_sym <- rlang::sym(group_var)
  outcome_sym <- rlang::sym(outcome_var)

  df <- data %>% filter(!is.na(!!group_sym))
  total_N <- nrow(df)
  total_events <- sum(df[[outcome_var]] == 1, na.rm = TRUE)
  total_nonevents <- sum(df[[outcome_var]] == 0, na.rm = TRUE)

  df %>%
    group_by(!!group_sym) %>%
    summarise(
      women_n = n(),
      events_n = sum(!!outcome_sym == 1, na.rm = TRUE),
      nonevents_n = sum(!!outcome_sym == 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      women_pct = 100 * women_n / total_N,
      event_pct_within_group = 100 * events_n / women_n,
      pct_events_of_all = 100 * events_n / total_events,
      pct_nonevents_of_all = 100 * nonevents_n / total_nonevents,
      LR = (events_n / total_events) / (nonevents_n / total_nonevents),
      se_logLR = sqrt(
        ifelse(events_n > 0, 1 / events_n, NA_real_) +
          ifelse(nonevents_n > 0, 1 / nonevents_n, NA_real_) -
          ifelse(total_events > 0, 1 / total_events, 0) -
          ifelse(total_nonevents > 0, 1 / total_nonevents, 0)
      ),
      LR_lower = ifelse(events_n > 0 & nonevents_n > 0, exp(log(LR) - 1.96 * se_logLR), NA_real_),
      LR_upper = ifelse(events_n > 0 & nonevents_n > 0, exp(log(LR) + 1.96 * se_logLR), NA_real_)
    ) %>%
    rename(risk_group = !!group_sym)
}

get_H0_t0 <- function(fit, t0, newdata = NULL) {
  bh <- if (is.null(newdata)) basehaz(fit, centered = FALSE) else basehaz(fit, newdata = newdata)
  if (nrow(bh) == 0) stop("No baseline hazard values returned.")
  if (t0 < min(bh$time, na.rm = TRUE)) return(0)
  idx <- max(which(bh$time <= t0))
  bh$hazard[idx]
}

compute_pi <- function(data, beta_values, hyp_preexisting_labels, hyp_gestational_labels, hyp_preeclampsia_labels, hyp_superimposed_labels) {
  with(data,
    beta_values["maternal_age_spline2"] * age_linear +
      beta_values["maternal_age_spline3"] * age_rcs1 +
      beta_values["gest_age_spline2"] * ga_linear +
      beta_values["gest_age_spline3"] * ga_rcs1 +
      beta_values["socioecon"] * as.numeric(socioeconomic_deprivation == "yes") +
      beta_values["substance_use"] * as.numeric(substance_use == "yes") +
      beta_values["smm_now"] * as.numeric(smm_now == "yes") +
      beta_values["previous_complications"] * as.numeric(previous_complications == "Yes, multiparous") +
      beta_values["previous_complications_primip"] * as.numeric(previous_complications == "Primiparous") +
      beta_values["hyp_disorders_preexisting"] * as.numeric(hyp_disorders %in% hyp_preexisting_labels) +
      beta_values["hyp_disorders_gestational"] * as.numeric(hyp_disorders %in% hyp_gestational_labels) +
      beta_values["hyp_disorders_preeclampsia"] * as.numeric(hyp_disorders %in% hyp_preeclampsia_labels) +
      beta_values["hyp_disorders_superimposed"] * as.numeric(hyp_disorders %in% hyp_superimposed_labels) +
      beta_values["maternal_age_smm_spline2"] * age_linear * as.numeric(smm_now == "yes") +
      beta_values["maternal_age_smm_spline3"] * age_rcs1 * as.numeric(smm_now == "yes") +
      beta_values["maternal_age_hyp_preexist_spline2"] * age_linear * as.numeric(hyp_disorders %in% hyp_preexisting_labels) +
      beta_values["maternal_age_hyp_preexist_spline3"] * age_rcs1 * as.numeric(hyp_disorders %in% hyp_preexisting_labels) +
      beta_values["maternal_age_hyp_gestational_spline2"] * age_linear * as.numeric(hyp_disorders %in% hyp_gestational_labels) +
      beta_values["maternal_age_hyp_gestational_spline3"] * age_rcs1 * as.numeric(hyp_disorders %in% hyp_gestational_labels) +
      beta_values["maternal_age_hyp_preeclampsia_spline2"] * age_linear * as.numeric(hyp_disorders %in% hyp_preeclampsia_labels) +
      beta_values["maternal_age_hyp_preeclampsia_spline3"] * age_rcs1 * as.numeric(hyp_disorders %in% hyp_preeclampsia_labels) +
      beta_values["maternal_age_hyp_superimposed_spline2"] * age_linear * as.numeric(hyp_disorders %in% hyp_superimposed_labels) +
      beta_values["maternal_age_hyp_superimposed_spline3"] * age_rcs1 * as.numeric(hyp_disorders %in% hyp_superimposed_labels)
  )
}

fit_and_export_models <- function(data, prefix) {
  fit_pi <- coxph(Surv(time, outcome) ~ PI, data = data)
  fit_test <- coxph(Surv(time, outcome) ~ PI + offset(PI), data = data)
  fit_coeff <- coxph(
    Surv(time, outcome) ~
      age_linear + age_rcs1 + ga_linear + ga_rcs1 +
      socioeconomic_deprivation + substance_use + smm_now +
      previous_complications + hyp_disorders +
      age_linear * smm_now + age_rcs1 * smm_now +
      age_linear * hyp_disorders + age_rcs1 * hyp_disorders,
    data = data
  )

  model_table <- tidy(fit_coeff, conf.int = TRUE, exponentiate = FALSE)
  write.xlsx(model_table, file.path(config$output_dir, paste0(prefix, "_model_coefficients.xlsx")), rowNames = FALSE)
  write.csv(model_table, file.path(config$output_dir, paste0(prefix, "_model_coefficients.csv")), row.names = FALSE)

  list(
    fit_pi = fit_pi,
    fit_test = fit_test,
    fit_coeff = fit_coeff,
    c_index_pi = get_c_index_with_ci(fit_pi),
    c_index_coeff = get_c_index_with_ci(fit_coeff)
  )
}

##############################
# 4. LOAD DATA
##############################
final_dataset_complication <- readRDS(complication_path)
final_dataset_no_complication <- readRDS(no_complication_path)
cleaned_cohort <- readRDS(cleaned_cohort_path)

##############################
# 5. BETA COEFFICIENTS
##############################
beta_values <- c(
  maternal_age_spline2 = -0.0964,
  maternal_age_spline3 = 0.0805,
  gest_age_spline2 = 0.0521,
  gest_age_spline3 = -0.3406,
  socioecon = 0.7941,
  substance_use = 0.4916,
  smm_now = 2.7722,
  previous_complications = 0.2469,
  previous_complications_primip = 0.0447,
  hyp_disorders_preexisting = -1.6100,
  hyp_disorders_gestational = -0.0100,
  hyp_disorders_preeclampsia = 1.6897,
  hyp_disorders_superimposed = -9.7517,
  maternal_age_smm_spline2 = 0.1121,
  maternal_age_smm_spline3 = 0.1359,
  maternal_age_hyp_preexist_spline2 = 0.0820,
  maternal_age_hyp_preexist_spline3 = -0.1052,
  maternal_age_hyp_gestational_spline2 = 0.0329,
  maternal_age_hyp_gestational_spline3 = -0.0790,
  maternal_age_hyp_preeclampsia_spline2 = -0.0422,
  maternal_age_hyp_preeclampsia_spline3 = 0.0639,
  maternal_age_hyp_superimposed_spline2 = 0.3724,
  maternal_age_hyp_superimposed_spline3 = -0.2730
)

##############################
# 6. COMPUTE PI
##############################
final_dataset_complication$PI <- compute_pi(
  final_dataset_complication,
  beta_values,
  hyp_preexisting_labels = c("pre-existing or unspecified hypertension", "preexisting_unspecified"),
  hyp_gestational_labels = c("gestational hypertension", "gestational_hypertension"),
  hyp_preeclampsia_labels = c("pre-eclampsia/HELLP syndrome", "preeclampsia_hellp_eclampsia"),
  hyp_superimposed_labels = c("superimposed pre-eclampsia", "superimposed_preeclampsia")
)

final_dataset_no_complication$PI <- compute_pi(
  final_dataset_no_complication,
  beta_values,
  hyp_preexisting_labels = c("pre-existing or unspecified hypertension", "preexisting_unspecified"),
  hyp_gestational_labels = c("gestational hypertension", "gestational_hypertension"),
  hyp_preeclampsia_labels = c("pre-eclampsia/HELLP syndrome", "preeclampsia_hellp_eclampsia"),
  hyp_superimposed_labels = c("superimposed pre-eclampsia", "superimposed_preeclampsia")
)

##############################
# 7. FIT MAIN MODELS
##############################
comp_results <- fit_and_export_models(final_dataset_complication, "03_complication")
no_comp_results <- fit_and_export_models(final_dataset_no_complication, "03_no_complication")

print(comp_results$c_index_pi)
print(comp_results$c_index_coeff)
print(no_comp_results$c_index_pi)
print(no_comp_results$c_index_coeff)

##############################
# 8. RECALIBRATION EXAMPLES
##############################
final_dataset_complication <- final_dataset_complication %>% mutate(PI_recal = 0.28 * PI)
fit_comp_recal <- coxph(Surv(time, outcome) ~ PI_recal, data = final_dataset_complication)
print(get_c_index_with_ci(fit_comp_recal))

final_dataset_no_complication <- final_dataset_no_complication %>% mutate(PI_recal = 0.18 * PI)
fit_no_comp_recal <- coxph(Surv(time, outcome) ~ PI_recal, data = final_dataset_no_complication)
print(get_c_index_with_ci(fit_no_comp_recal))

##############################
# 9. OPTIONAL EXTENDED MODELS
##############################
final_dataset_complication <- final_dataset_complication %>%
  mutate(
    gen_ethnicity = relevel(factor(gen_ethnicity), ref = "White"),
    treatment_method = relevel(factor(treatment_method), ref = "lifestyle")
  )

final_dataset_no_complication <- final_dataset_no_complication %>%
  mutate(
    gen_ethnicity = relevel(factor(gen_ethnicity), ref = "White"),
    treatment_method = relevel(factor(treatment_method), ref = "lifestyle")
  )

fit_comp_add <- coxph(Surv(time, outcome) ~ PI + gen_ethnicity + treatment_method + bmi_value, data = final_dataset_complication)
fit_no_comp_add <- coxph(Surv(time, outcome) ~ PI + gen_ethnicity + treatment_method + bmi_value, data = final_dataset_no_complication)

print(get_c_index_with_ci(fit_comp_add))
print(get_c_index_with_ci(fit_no_comp_add))

##############################
# 10. RISK TABLES
##############################
t0 <- 365.25 * config$followup_years

# Complication endpoint risk groups.
mf_comp <- model.frame(comp_results$fit_pi)
rows_used_comp <- as.integer(rownames(mf_comp))
lp_used_comp <- predict(comp_results$fit_pi, newdata = mf_comp, type = "lp")
sf_comp <- survfit(comp_results$fit_pi)
S0_t0_comp <- if (t0 > max(sf_comp$time, na.rm = TRUE)) tail(sf_comp$surv, 1) else sf_comp$surv[max(which(sf_comp$time <= t0))]
final_dataset_complication$pred_risk <- NA_real_
final_dataset_complication$pred_risk[rows_used_comp] <- 1 - S0_t0_comp ^ exp(lp_used_comp)
final_dataset_complication <- final_dataset_complication %>%
  mutate(risk_group = cut(pred_risk, breaks = c(-Inf, 0.006, 0.008, 0.03, Inf), labels = c("<0.006", ">=0.006–<0.008", ">=0.008–<0.03", ">=0.03"), right = FALSE))
risk_table_comp <- make_risk_table(final_dataset_complication, "risk_group")
write.xlsx(risk_table_comp, file.path(config$output_dir, "03_risk_table_complication.xlsx"), rowNames = FALSE)
write.csv(risk_table_comp, file.path(config$output_dir, "03_risk_table_complication.csv"), row.names = FALSE)

# No-complication endpoint risk groups.
mf_no_comp <- model.frame(no_comp_results$fit_pi)
rows_used_no_comp <- as.integer(rownames(mf_no_comp))
lp_used_no_comp <- predict(no_comp_results$fit_pi, newdata = mf_no_comp, type = "lp")
sf_no_comp <- survfit(no_comp_results$fit_pi)
S0_t0_no_comp <- if (t0 > max(sf_no_comp$time, na.rm = TRUE)) tail(sf_no_comp$surv, 1) else sf_no_comp$surv[max(which(sf_no_comp$time <= t0))]
final_dataset_no_complication$pred_risk <- NA_real_
final_dataset_no_complication$pred_risk[rows_used_no_comp] <- 1 - S0_t0_no_comp ^ exp(lp_used_no_comp)
final_dataset_no_complication <- final_dataset_no_complication %>%
  mutate(risk_group = cut(pred_risk, breaks = c(-Inf, 0.006, 0.008, 0.03, Inf), labels = c("<0.006", ">=0.006–<0.008", ">=0.008–<0.03", ">=0.03"), right = FALSE))
risk_table_no_comp <- make_risk_table(final_dataset_no_complication, "risk_group")
write.xlsx(risk_table_no_comp, file.path(config$output_dir, "03_risk_table_no_complication.xlsx"), rowNames = FALSE)
write.csv(risk_table_no_comp, file.path(config$output_dir, "03_risk_table_no_complication.csv"), row.names = FALSE)

##############################
# 11. INCLUDED VS EXCLUDED CHECK
##############################
vars_required <- c(
  "age_linear", "age_rcs1", "ga_linear", "ga_rcs1",
  "socioeconomic_deprivation", "substance_use", "smm_now",
  "previous_complications", "hyp_disorders"
)

comparison_df <- final_dataset_complication %>%
  mutate(
    included = if_else(
      if_all(all_of(vars_required), ~ !is.na(.)),
      "Included (complete-case)",
      "Excluded (missing data)"
    )
  )

vars_compare <- c(
  "time", "maternal_age", "bmi_value", "gest_age", "birthweight", "numpreg",
  "parity", "delivery_method", "smm_now", "alcohol_status", "previous_complications",
  "hyp_disorders", "smoking_status", "gen_ethnicity", "treatment_method"
)

tab_missing <- CreateTableOne(vars = vars_compare, strata = "included", data = comparison_df, test = TRUE)
capture.output(
  print(tab_missing, showAllLevels = TRUE, smd = TRUE, quote = FALSE, noSpaces = TRUE, digits = 2),
  file = file.path(config$output_dir, "03_tableone_included_vs_excluded.txt")
)

##############################
# 12. SIMPLE DESCRIPTIVE TABLES
##############################
table_ses <- cleaned_cohort %>%
  mutate(
    ses_quintile = case_when(
      ses %in% c(1, 2) ~ "Q1 - Least deprived",
      ses %in% c(3, 4) ~ "Q2",
      ses %in% c(5, 6) ~ "Q3",
      ses %in% c(7, 8) ~ "Q4",
      ses %in% c(9, 10) ~ "Q5 - Most deprived",
      TRUE ~ "Not recorded"
    )
  ) %>%
  count(ses_quintile, name = "Number") %>%
  mutate(`Per cent` = round(100 * Number / sum(Number), 1))

table_region <- cleaned_cohort %>%
  count(region_name, name = "Number") %>%
  mutate(`Per cent` = round(100 * Number / sum(Number), 1))

table_ethnicity <- cleaned_cohort %>%
  count(gen_ethnicity, name = "Number") %>%
  mutate(`Per cent` = round(100 * Number / sum(Number), 1))

write.xlsx(
  list(
    ses_quintile = table_ses,
    region = table_region,
    ethnicity = table_ethnicity
  ),
  file.path(config$output_dir, "03_descriptive_tables.xlsx")
)

##############################
# 13. SHRINKAGE EXAMPLE
##############################
analysis_df <- final_dataset_complication %>%
  mutate(outcome = as.integer(outcome)) %>%
  filter(!is.na(time), !is.na(outcome), !is.na(PI)) %>%
  droplevels()

fit_offset <- coxph(Surv(time, outcome) ~ offset(PI), data = analysis_df, x = TRUE, y = TRUE)
H0_t0_offset <- get_H0_t0(fit_offset, t0, newdata = data.frame(PI = 0))
analysis_df <- analysis_df %>% mutate(pred_risk_offset = 1 - exp(-H0_t0_offset * exp(PI)))

fit_gamma <- coxph(Surv(time, outcome) ~ PI, data = analysis_df, x = TRUE, y = TRUE)
gamma_hat <- unname(coef(fit_gamma)["PI"])
gamma_ci <- confint(fit_gamma)["PI", ]
analysis_df <- analysis_df %>% mutate(PI_shrunk = gamma_hat * PI)

fit_shrunk_offset <- coxph(Surv(time, outcome) ~ offset(PI_shrunk), data = analysis_df, x = TRUE, y = TRUE)
H0_t0_shrunk <- get_H0_t0(fit_shrunk_offset, t0, newdata = data.frame(PI_shrunk = 0))
analysis_df <- analysis_df %>% mutate(pred_risk_shrunk = 1 - exp(-H0_t0_shrunk * exp(PI_shrunk)))

analysis_df <- analysis_df %>%
  mutate(
    risk_group_offset = cut(pred_risk_offset, breaks = c(-Inf, 0.006, 0.008, 0.03, Inf), labels = c("<0.006", ">=0.006–<0.008", ">=0.008–<0.03", ">=0.03"), right = FALSE),
    risk_group_shrunk = cut(pred_risk_shrunk, breaks = c(-Inf, 0.006, 0.008, 0.03, Inf), labels = c("<0.006", ">=0.006–<0.008", ">=0.008–<0.03", ">=0.03"), right = FALSE)
  )

table_offset <- make_risk_table(analysis_df, "risk_group_offset")
table_shrunk <- make_risk_table(analysis_df, "risk_group_shrunk")

write.xlsx(
  list(
    offset_only = table_offset,
    offset_plus_shrinkage = table_shrunk,
    shrinkage_summary = data.frame(gamma_hat = gamma_hat, gamma_ci_lower = gamma_ci[1], gamma_ci_upper = gamma_ci[2])
  ),
  file.path(config$output_dir, "03_recalibration_and_shrinkage.xlsx")
)

writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "03_session_info.txt"))

cat("03_model.R completed successfully.\n")
