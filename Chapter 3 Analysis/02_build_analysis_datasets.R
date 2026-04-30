############################################################
# 02_build_analysis_datasets.R
# Purpose: Build time-to-event analysis datasets, create
#          spline terms, and export descriptive tables.
############################################################

##############################
# 1. PACKAGE SETUP
##############################
required_packages <- c("dplyr", "lubridate", "rms", "Hmisc", "openxlsx", "tibble", "knitr")
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
  collection_end_date = as.Date("2023-12-20"),
  followup_years = 10
)

dir.create(config$processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(config$processed_dir, "cleaned_cohort.rds")
stopifnot(file.exists(input_path))

##############################
# 3. HELPER FUNCTIONS
##############################
build_time_to_event_dataset <- function(data, event_date_col, status_col, negative_value = "no") {
  event_date_col <- rlang::ensym(event_date_col)
  status_col <- rlang::ensym(status_col)

  data %>%
    mutate(
      delivery_date = as.Date(delivery_date),
      event_date = as.Date(!!event_date_col),
      dod = as.Date(dod),
      end_time = add_with_rollback(delivery_date, years(config$followup_years)),
      collection_end = config$collection_end_date,
      event_date = if_else(as.character(!!status_col) == negative_value, as.Date(NA), event_date),
      outcome = case_when(
        as.character(!!status_col) == negative_value ~ 0L,
        !is.na(event_date) & (is.na(dod) | event_date < dod) & (event_date <= end_time) ~ 1L,
        TRUE ~ 0L
      ),
      event_time = pmin(event_date, collection_end, dod, end_time, na.rm = TRUE),
      time = as.numeric(difftime(event_time, delivery_date, units = "days")) / 365.25
    )
}

fmt_mean_sd <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  paste0(
    format(round(mean(x), digits), nsmall = digits, big.mark = ","),
    " (",
    format(round(sd(x), digits), nsmall = digits, big.mark = ","),
    ")"
  )
}

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

smd_cont <- function(x1, x0) {
  x1 <- x1[!is.na(x1)]
  x0 <- x0[!is.na(x0)]
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  pooled_sd <- sqrt(((length(x1) - 1) * sd(x1)^2 + (length(x0) - 1) * sd(x0)^2) / (length(x1) + length(x0) - 2))
  if (pooled_sd == 0) return(NA_real_)
  abs((mean(x1) - mean(x0)) / pooled_sd)
}

smd_binary <- function(x1, x0) {
  x1 <- x1[!is.na(x1)]
  x0 <- x0[!is.na(x0)]
  if (length(x1) == 0 || length(x0) == 0) return(NA_real_)
  p1 <- mean(x1)
  p0 <- mean(x0)
  denom <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  if (denom == 0) return(NA_real_)
  abs((p1 - p0) / denom)
}

make_cont_row <- function(data, var, label, digits = 2) {
  var_name <- dplyr::ensym(var)
  x_excluded <- data %>% filter(included == "Excluded (missing data)") %>% pull(!!var_name)
  x_included <- data %>% filter(included == "Included (complete-case)") %>% pull(!!var_name)
  p_val <- tryCatch(t.test(x_included, x_excluded)$p.value, error = function(e) NA_real_)

  tibble(
    Characteristic = label,
    excluded = fmt_mean_sd(x_excluded, digits),
    included = fmt_mean_sd(x_included, digits),
    p_value = fmt_p(p_val),
    smd = ifelse(is.na(smd_cont(x_included, x_excluded)), NA_character_, sprintf("%.2f", smd_cont(x_included, x_excluded)))
  )
}

make_binary_row <- function(data, var, label) {
  var_name <- dplyr::ensym(var)
  df <- data %>%
    select(included, !!var_name) %>%
    filter(!is.na(included), !is.na(!!var_name)) %>%
    mutate(var_num = ifelse(!!var_name %in% c(1, "1", "yes", "Yes", TRUE), 1, 0))

  x_excluded <- df %>% filter(included == "Excluded (missing data)") %>% pull(var_num)
  x_included <- df %>% filter(included == "Included (complete-case)") %>% pull(var_num)
  tab <- table(df$included, df$var_num)
  p_val <- tryCatch(if (any(tab < 5)) fisher.test(tab)$p.value else chisq.test(tab)$p.value, error = function(e) NA_real_)

  tibble(
    Characteristic = label,
    excluded = sprintf("%s (%.1f%%)", format(sum(x_excluded == 1, na.rm = TRUE), big.mark = ","), 100 * mean(x_excluded, na.rm = TRUE)),
    included = sprintf("%s (%.1f%%)", format(sum(x_included == 1, na.rm = TRUE), big.mark = ","), 100 * mean(x_included, na.rm = TRUE)),
    p_value = fmt_p(p_val),
    smd = ifelse(is.na(smd_binary(x_included, x_excluded)), NA_character_, sprintf("%.2f", smd_binary(x_included, x_excluded)))
  )
}

calc_incidence <- function(df, group_var) {
  group_sym <- enquo(group_var)
  df %>%
    filter(!is.na(!!group_sym)) %>%
    group_by(!!group_sym) %>%
    summarise(
      n_women = n(),
      events = sum(outcome == 1, na.rm = TRUE),
      person_years = sum(time, na.rm = TRUE),
      incidence_per_10000 = ifelse(person_years > 0, (events / person_years) * 10000, NA_real_),
      lower_ci = ifelse(events == 0, 0, (qchisq(0.025, 2 * events) / 2) / person_years * 10000),
      upper_ci = (qchisq(0.975, 2 * (events + 1)) / 2) / person_years * 10000,
      .groups = "drop"
    ) %>%
    mutate(incidence_ci = sprintf("%.1f (%.1f–%.1f)", incidence_per_10000, lower_ci, upper_ci)) %>%
    select(group = !!group_sym, n_women, incidence_ci)
}

poisson_ci_per10000 <- function(events, py) {
  rate <- (events / py) * 10000
  lower <- if (events == 0) 0 else (qchisq(0.025, 2 * events) / 2) / py * 10000
  upper <- (qchisq(0.975, 2 * (events + 1)) / 2) / py * 10000
  list(rate = rate, lower = lower, upper = upper)
}

##############################
# 4. LOAD CLEANED COHORT
##############################
cleaned_cohort <- readRDS(input_path)

##############################
# 5. BUILD ANALYSIS DATASETS
##############################
final_dataset_complication <- build_time_to_event_dataset(
  cleaned_cohort,
  first_event_date_complication,
  complications_10yr,
  negative_value = "no"
)

final_dataset_no_complication <- build_time_to_event_dataset(
  cleaned_cohort,
  first_event_date_no_complication,
  diabetes_status_10yr,
  negative_value = "no"
)

cat("Median follow-up time for complications:", median(final_dataset_complication$time, na.rm = TRUE), "years\n")
cat("Median follow-up time for diabetes:", median(final_dataset_no_complication$time, na.rm = TRUE), "years\n")

##############################
# 6. SPLINE TERMS
##############################
knots_age <- as.numeric(quantile(cleaned_cohort$maternal_age, probs = c(0.04, 0.45, 0.96), na.rm = TRUE))
knots_ga <- as.numeric(quantile(cleaned_cohort$gest_age, probs = c(0.03, 0.10, 0.75), na.rm = TRUE))

age_spl <- rcspline.eval(cleaned_cohort$maternal_age, knots = knots_age, inclx = TRUE)
ga_spl <- rcspline.eval(cleaned_cohort$gest_age, knots = knots_ga, inclx = TRUE)

if (ncol(age_spl) != 2 || ncol(ga_spl) != 2) {
  stop("Unexpected number of spline columns returned. Review knot placement and spline output.")
}

colnames(age_spl) <- c("age_linear", "age_rcs1")
colnames(ga_spl) <- c("ga_linear", "ga_rcs1")

final_dataset_complication <- bind_cols(final_dataset_complication, as.data.frame(age_spl), as.data.frame(ga_spl))
final_dataset_no_complication <- bind_cols(final_dataset_no_complication, as.data.frame(age_spl), as.data.frame(ga_spl))

final_dataset_complication <- final_dataset_complication %>%
  mutate(
    hyp_disorders = relevel(factor(hyp_disorders), ref = "no"),
    previous_complications = relevel(factor(previous_complications), ref = "No, multiparous")
  )

final_dataset_no_complication <- final_dataset_no_complication %>%
  mutate(
    hyp_disorders = relevel(factor(hyp_disorders), ref = "no"),
    previous_complications = relevel(factor(previous_complications), ref = "No, multiparous")
  )

##############################
# 7. COMPLETE-CASE FLAG
##############################
vars_required <- c(
  "age_linear", "age_rcs1", "ga_linear", "ga_rcs1",
  "socioeconomic_deprivation", "substance_use", "smm_now",
  "previous_complications", "hyp_disorders"
)

analysis_dataset_full <- final_dataset_complication %>%
  mutate(
    included = if_else(
      if_all(all_of(vars_required), ~ !is.na(.)),
      "Included (complete-case)",
      "Excluded (missing data)"
    )
  ) %>%
  filter(!is.na(time) & time >= 0)

analysis_dataset_complete <- analysis_dataset_full %>%
  filter(included == "Included (complete-case)")

##############################
# 8. INCLUDED VS EXCLUDED TABLE
##############################
table_included_excluded <- bind_rows(
  make_cont_row(analysis_dataset_full, time, "Follow-up time, years (mean ± SD)", digits = 2),
  make_cont_row(analysis_dataset_full, maternal_age, "Maternal age, years (mean ± SD)", digits = 2),
  make_cont_row(analysis_dataset_full, bmi_value, "BMI, kg/m² (mean ± SD)", digits = 2),
  make_cont_row(analysis_dataset_full, gest_age, "Gestational age, weeks (mean ± SD)", digits = 2),
  make_cont_row(analysis_dataset_full, birthweight, "Birth weight, g (mean ± SD)", digits = 0),
  make_cont_row(analysis_dataset_full, numpreg, "Number of previous pregnancies (mean ± SD)", digits = 2),
  make_binary_row(analysis_dataset_full, outcome, "Complications within 10 years, n (%)")
)

n_excluded <- sum(analysis_dataset_full$included == "Excluded (missing data)")
n_included <- sum(analysis_dataset_full$included == "Included (complete-case)")

colnames(table_included_excluded) <- c(
  "Characteristic",
  paste0("Excluded (missing data)\n(n = ", format(n_excluded, big.mark = ","), ")"),
  paste0("Included (complete-case)\n(n = ", format(n_included, big.mark = ","), ")"),
  "p-value",
  "SMD"
)

write.xlsx(table_included_excluded, file.path(config$output_dir, "02_included_vs_excluded_table.xlsx"), rowNames = FALSE)
write.csv(table_included_excluded, file.path(config$output_dir, "02_included_vs_excluded_table.csv"), row.names = FALSE)

##############################
# 9. INCIDENCE TABLES
##############################
analysis_dataset_incidence <- final_dataset_complication %>%
  filter(!is.na(time) & time >= 0) %>%
  filter(if_all(all_of(vars_required), ~ !is.na(.))) %>%
  mutate(
    age_cat = factor(age_category, levels = c("<25", "25–29", "30–34", "≥35")),
    primiparity_label = ifelse(as.character(primiparity) %in% c("1", "yes", "Yes", "TRUE", "T"), "Yes", "No") %>% factor(levels = c("Yes", "No")),
    morbid_obese_label = case_when(
      is.na(morbid_obese) ~ NA_character_,
      morbid_obese %in% c("yes", "Yes", "1", 1, TRUE) ~ "Yes",
      TRUE ~ "No"
    ) %>% factor(levels = c("Yes", "No")),
    substance_use_label = case_when(
      is.na(substance_use) ~ NA_character_,
      substance_use %in% c("yes", "Yes", "1", 1, TRUE) ~ "Yes",
      TRUE ~ "No"
    ) %>% factor(levels = c("Yes", "No")),
    deprived_label = case_when(
      socioeconomic_deprivation %in% c("yes") ~ "Deprived",
      !is.na(socioeconomic_deprivation) ~ "Not deprived",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Deprived", "Not deprived")),
    hdp_type = case_when(
      hyp_disorders %in% c("no", "none") ~ "No HDP",
      hyp_disorders %in% c("preexisting_unspecified", "pre-existing or unspecified hypertension", "preexisting", "unspecified") ~ "Pre-existing or unspecified hypertension",
      hyp_disorders %in% c("gestational_hypertension", "gestational hypertension", "gestational") ~ "Gestational hypertension",
      hyp_disorders %in% c("preeclampsia_hellp_eclampsia", "pre-eclampsia/HELLP syndrome", "preeclampsia", "hellp", "eclampsia") ~ "Pre-eclampsia/HELLP syndrome",
      hyp_disorders %in% c("superimposed_preeclampsia", "superimposed pre-eclampsia", "superimposed") ~ "Superimposed pre-eclampsia",
      TRUE ~ as.character(hyp_disorders)
    ) %>% factor(levels = c(
      "No HDP",
      "Pre-existing or unspecified hypertension",
      "Gestational hypertension",
      "Pre-eclampsia/HELLP syndrome",
      "Superimposed pre-eclampsia"
    ))
  )

inc_age <- calc_incidence(analysis_dataset_incidence, age_cat) %>% mutate(characteristic = "Age (years)") %>% select(characteristic, group, n_women, incidence_ci)
inc_parity <- calc_incidence(analysis_dataset_incidence, primiparity_label) %>% mutate(characteristic = "Primiparity") %>% select(characteristic, group, n_women, incidence_ci)
inc_morbid <- calc_incidence(analysis_dataset_incidence, morbid_obese_label) %>% mutate(characteristic = "Morbid obesity") %>% select(characteristic, group, n_women, incidence_ci)
inc_subs <- calc_incidence(analysis_dataset_incidence, substance_use_label) %>% mutate(characteristic = "Substance use") %>% select(characteristic, group, n_women, incidence_ci)
inc_depr <- calc_incidence(analysis_dataset_incidence, deprived_label) %>% mutate(characteristic = "Socioeconomic deprivation") %>% select(characteristic, group, n_women, incidence_ci)
inc_hdp <- calc_incidence(analysis_dataset_incidence, hdp_type) %>% mutate(characteristic = "HDP type") %>% select(characteristic, group, n_women, incidence_ci)

inc_total <- analysis_dataset_incidence %>%
  summarise(group = "Total", n_women = n(), events = sum(outcome == 1, na.rm = TRUE), person_years = sum(time, na.rm = TRUE)) %>%
  mutate(
    lower_ci = ifelse(events == 0, 0, (qchisq(0.025, 2 * events) / 2) / person_years * 10000),
    upper_ci = (qchisq(0.975, 2 * (events + 1)) / 2) / person_years * 10000,
    incidence_per_10000 = (events / person_years) * 10000,
    incidence_ci = sprintf("%.1f (%.1f–%.1f)", incidence_per_10000, lower_ci, upper_ci),
    characteristic = "Total"
  ) %>%
  select(characteristic, group, n_women, incidence_ci)

final_incidence_table <- bind_rows(
  inc_age,
  tibble(characteristic = "", group = NA_character_, n_women = NA_integer_, incidence_ci = NA_character_),
  inc_parity,
  tibble(characteristic = "", group = NA_character_, n_women = NA_integer_, incidence_ci = NA_character_),
  inc_morbid,
  tibble(characteristic = "", group = NA_character_, n_women = NA_integer_, incidence_ci = NA_character_),
  inc_subs,
  tibble(characteristic = "", group = NA_character_, n_women = NA_integer_, incidence_ci = NA_character_),
  inc_depr,
  tibble(characteristic = "", group = NA_character_, n_women = NA_integer_, incidence_ci = NA_character_),
  inc_hdp,
  tibble(characteristic = "", group = NA_character_, n_women = NA_integer_, incidence_ci = NA_character_),
  inc_total
)

write.xlsx(final_incidence_table, file.path(config$output_dir, "02_incidence_characteristics_table.xlsx"), rowNames = FALSE)
write.csv(final_incidence_table, file.path(config$output_dir, "02_incidence_characteristics_table.csv"), row.names = FALSE)

# Incidence by complication type.
df <- analysis_dataset_incidence
total_py <- sum(df$time, na.rm = TRUE)
comp_codes <- c(
  "diabetes_complications_coma",
  "diabetes_complications_acidosis",
  "diabetes_complications_circulation",
  "diabetes_complications_renal",
  "diabetes_complications_ophtalmo",
  "diabetes_complications_neuro",
  "diabetes_complications_other",
  "diabetes_complications_multiple"
)
friendly_labels <- c(
  "Coma",
  "Acidosis",
  "Circulation",
  "Kidney complications",
  "Ophthalmic",
  "Neurological",
  "Other outcomes",
  "Multiple complications"
)

rows <- lapply(seq_along(comp_codes), function(i) {
  events <- df %>% filter(complications_10yr == comp_codes[i]) %>% nrow()
  ci <- poisson_ci_per10000(events, total_py)
  tibble(
    Complication = friendly_labels[i],
    `Total number of women` = events,
    `Incidence per 10,000 person-years (95% CI)` = sprintf("%.1f (%.1f–%.1f)", ci$rate, ci$lower, ci$upper)
  )
})

any_events <- df %>% filter(!is.na(complications_10yr) & complications_10yr != "no") %>% nrow()
ci_any <- poisson_ci_per10000(any_events, total_py)
any_row <- tibble(
  Complication = "Diabetes with any complication",
  `Total number of women` = any_events,
  `Incidence per 10,000 person-years (95% CI)` = sprintf("%.1f (%.1f–%.1f)", ci_any$rate, ci_any$lower, ci_any$upper)
)

total_row <- tibble(
  Complication = "Total (cohort size)",
  `Total number of women` = nrow(df),
  `Incidence per 10,000 person-years (95% CI)` = sprintf("%.1f (%.1f–%.1f)", (sum(df$complications_10yr != "no", na.rm = TRUE) / total_py) * 10000, ci_any$lower, ci_any$upper)
)

final_complication_incidence_table <- bind_rows(any_row, bind_rows(rows), total_row)
write.xlsx(final_complication_incidence_table, file.path(config$output_dir, "02_incidence_complications_table.xlsx"), rowNames = FALSE)
write.csv(final_complication_incidence_table, file.path(config$output_dir, "02_incidence_complications_table.csv"), row.names = FALSE)

##############################
# 10. SAVE ANALYSIS DATASETS
##############################
saveRDS(final_dataset_complication, file.path(config$processed_dir, "analysis_dataset_complication.rds"))
saveRDS(final_dataset_no_complication, file.path(config$processed_dir, "analysis_dataset_no_complication.rds"))
saveRDS(analysis_dataset_complete, file.path(config$processed_dir, "analysis_dataset_complete_case_complication.rds"))
writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "02_session_info.txt"))

cat("02_build_analysis_datasets.R completed successfully.\n")
