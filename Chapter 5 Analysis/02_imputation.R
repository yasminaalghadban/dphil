############################################################
# 02_imputation.R
# Purpose: Run missingness diagnostic models and perform
#          multiple imputation for the sensitivity cohort.
############################################################

##############################
# 1. PACKAGE SETUP
##############################
required_packages <- c(
  "dplyr", "forcats", "tidyr", "mice", "broom",
  "survival", "readr", "ggplot2"
)
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
  output_dir    = "outputs",
  seed          = 12345,
  test_n        = 20000,
  test_m        = 2,
  test_maxit    = 2,
  final_m       = 40,
  final_maxit   = 20
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(config$processed_dir, "sensitivity_cohort.rds")
stopifnot(file.exists(input_path))

##############################
# 3. LOAD DATA
##############################
final_dataset_no_complication <- readRDS(input_path)

##############################
# 4. HELPER FUNCTIONS
##############################
or_table <- function(fit) {
  broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
}

save_chain_plot <- function(imp_object, output_file) {
  png(output_file, width = 1400, height = 900, res = 150)
  plot(imp_object)
  dev.off()
}

##############################
# 5. MISSINGNESS DIAGNOSTIC DATASET
##############################
model1_df <- final_dataset_no_complication %>%
  mutate(
    missed = ifelse(selection == "hes", 1L, 0L),
    smoking_status = factor(smoking_status, levels = c("never", "former", "current")),
    gen_ethnicity = factor(gen_ethnicity),
    gen_ethnicity = relevel(gen_ethnicity, ref = "White"),
    region_name = factor(region_name),
    ses_quintile = factor(
      ses_quintile,
      levels = c(1, 2, 3, 4, 5),
      labels = c("1 (least deprived)", "2", "3", "4", "5 (most deprived)"),
      ordered = TRUE
    ),
    previous_complications = factor(previous_complications),
    treatment_method = factor(treatment_method),
    delivery_period = case_when(
      yodelivery >= 2010 & yodelivery <= 2012 ~ "2010-2012",
      yodelivery >= 2013 & yodelivery <= 2015 ~ "2013-2015",
      yodelivery >= 2016 & yodelivery <= 2018 ~ "2016-2018",
      yodelivery >= 2019 & yodelivery <= 2021 ~ "2019-2021",
      TRUE ~ NA_character_
    ),
    delivery_period = factor(
      delivery_period,
      levels = c("2010-2012", "2013-2015", "2016-2018", "2019-2021"),
      ordered = TRUE
    ),
    hyp_disorders = factor(hyp_disorders),
    smm_now = factor(smm_now),
    diabetes_status_10yr = factor(diabetes_status_10yr)
  ) %>%
  select(
    missed, maternal_age, bmi_value, smoking_status, gen_ethnicity,
    ses_quintile, previous_complications, treatment_method,
    region_name, delivery_period, hyp_disorders, smm_now,
    diabetes_status_10yr, gest_age, birweit
  )

##############################
# 6. MISSINGNESS SUMMARIES
##############################
missing_summary <- sapply(model1_df, function(x) mean(is.na(x)) * 100)
missing_table <- data.frame(
  variable = names(missing_summary),
  percent_missing = round(as.numeric(missing_summary), 2)
) %>%
  arrange(desc(percent_missing))

write_csv(missing_table, file.path(config$output_dir, "02_missingness_summary.csv"))

png(file.path(config$output_dir, "02_missing_data_pattern.png"), width = 1400, height = 900, res = 150)
md.pattern(model1_df, rotate.names = TRUE)
dev.off()

##############################
# 7. MISSINGNESS INDICATORS
##############################
diag_df <- model1_df %>%
  mutate(
    bmi_missing = as.integer(is.na(bmi_value)),
    smoking_missing = as.integer(is.na(smoking_status)),
    any_model_missing = as.integer(
      is.na(bmi_value) |
        is.na(smoking_status) |
        is.na(gen_ethnicity) |
        is.na(region_name) |
        is.na(ses_quintile) |
        is.na(delivery_period)
    )
  )

##############################
# 8. MISSINGNESS DIAGNOSTIC MODELS
##############################
fit_bmi_missing <- glm(
  bmi_missing ~ missed + maternal_age + smoking_status + gen_ethnicity +
    ses_quintile + previous_complications + treatment_method + region_name +
    delivery_period + hyp_disorders + smm_now + diabetes_status_10yr +
    gest_age + birweit,
  family = binomial(),
  data = diag_df
)

fit_smoking_missing <- glm(
  smoking_missing ~ missed + maternal_age + bmi_value + gen_ethnicity +
    ses_quintile + previous_complications + treatment_method + region_name +
    delivery_period + hyp_disorders + smm_now + diabetes_status_10yr +
    gest_age + birweit,
  family = binomial(),
  data = diag_df
)

fit_any_missing <- glm(
  any_model_missing ~ missed + maternal_age + gen_ethnicity + ses_quintile +
    previous_complications + treatment_method + region_name + delivery_period +
    hyp_disorders + smm_now + diabetes_status_10yr + gest_age + birweit,
  family = binomial(),
  data = diag_df
)

bmi_missing_results <- or_table(fit_bmi_missing)
smoking_missing_results <- or_table(fit_smoking_missing)
any_missing_results <- or_table(fit_any_missing)

write_csv(bmi_missing_results, file.path(config$output_dir, "02_bmi_missing_results.csv"))
write_csv(smoking_missing_results, file.path(config$output_dir, "02_smoking_missing_results.csv"))
write_csv(any_missing_results, file.path(config$output_dir, "02_any_missing_results.csv"))

##############################
# 9. BUILD IMPUTATION DATASET
##############################
imp_df <- final_dataset_no_complication %>%
  mutate(
    .id = row_number(),
    selection = factor(selection, levels = c("pcr", "hes")),
    missed = ifelse(selection == "hes", 1L, 0L),
    diabetes_status_10yr = case_when(
      diabetes_status_10yr %in% c(1, "1", "Yes", "yes", "T2DM", "type2") ~ 1L,
      diabetes_status_10yr %in% c(0, "0", "No", "no") ~ 0L,
      TRUE ~ NA_integer_
    ),
    smoking_status = factor(tolower(as.character(smoking_status)), levels = c("never", "former", "current")),
    gen_ethnicity = factor(gen_ethnicity),
    gen_ethnicity = fct_relevel(gen_ethnicity, "White"),
    region_name = factor(region_name),
    ses_quintile = factor(
      ses_quintile,
      levels = c(1, 2, 3, 4, 5),
      labels = c("1 (least deprived)", "2", "3", "4", "5 (most deprived)"),
      ordered = TRUE
    ),
    previous_complications = factor(previous_complications),
    treatment_method = factor(treatment_method),
    delivery_period = case_when(
      yodelivery >= 2010 & yodelivery <= 2012 ~ "2010-2012",
      yodelivery >= 2013 & yodelivery <= 2015 ~ "2013-2015",
      yodelivery >= 2016 & yodelivery <= 2018 ~ "2016-2018",
      yodelivery >= 2019 & yodelivery <= 2021 ~ "2019-2021",
      TRUE ~ NA_character_
    ),
    delivery_period = factor(
      delivery_period,
      levels = c("2010-2012", "2013-2015", "2016-2018", "2019-2021"),
      ordered = TRUE
    ),
    hyp_disorders = factor(hyp_disorders),
    smm_now = factor(smm_now)
  ) %>%
  select(
    .id, patid, selection, missed, time, diabetes_status_10yr,
    maternal_age, bmi_value, smoking_status, gen_ethnicity,
    ses_quintile, previous_complications, treatment_method,
    region_name, delivery_period, hyp_disorders, smm_now,
    gest_age, birweit
  )

##############################
# 10. SET UP MICE
##############################
ini  <- mice(imp_df, maxit = 0, printFlag = FALSE)
meth <- ini$method
meth[] <- ""

if (anyNA(imp_df$bmi_value))      meth["bmi_value"]      <- "pmm"
if (anyNA(imp_df$smoking_status)) meth["smoking_status"] <- "polyreg"
if (anyNA(imp_df$gen_ethnicity))  meth["gen_ethnicity"]  <- "polyreg"
if (anyNA(imp_df$ses_quintile))   meth["ses_quintile"]   <- "polr"
if (anyNA(imp_df$region_name))    meth["region_name"]    <- "polyreg"
if (anyNA(imp_df$delivery_period)) meth["delivery_period"] <- "polr"

meth[c(
  ".id", "patid", "selection", "missed", "time", "diabetes_status_10yr",
  "maternal_age", "previous_complications", "treatment_method",
  "hyp_disorders", "smm_now", "gest_age", "birweit"
)] <- ""

pred <- quickpred(
  data = imp_df,
  mincor = 0.05,
  minpuc = 0.25,
  include = c(
    "missed", "time", "diabetes_status_10yr", "maternal_age",
    "previous_complications", "treatment_method", "delivery_period",
    "hyp_disorders", "smm_now"
  )
)

pred[c(
  ".id", "patid", "selection", "missed", "time", "diabetes_status_10yr",
  "maternal_age", "previous_complications", "treatment_method",
  "delivery_period", "hyp_disorders", "smm_now", "gest_age", "birweit"
), ] <- 0

pred[, c(".id", "patid")] <- 0
pred[, "selection"] <- 0
pred[, c(
  "missed", "time", "diabetes_status_10yr", "maternal_age",
  "previous_complications", "treatment_method", "delivery_period",
  "hyp_disorders", "smm_now"
)] <- 1
pred[, c("gest_age", "birweit")] <- 0
diag(pred) <- 0

write_csv(data.frame(variable = names(meth), method = unname(meth)), file.path(config$output_dir, "02_mice_methods.csv"))
write_csv(as.data.frame(pred), file.path(config$output_dir, "02_mice_predictor_matrix.csv"))

##############################
# 11. PRACTICE RUN
##############################
set.seed(config$seed)
test_df <- imp_df %>% slice_sample(n = min(config$test_n, n()))

imp_test <- mice(
  data = test_df,
  m = config$test_m,
  maxit = config$test_maxit,
  method = meth,
  predictorMatrix = pred,
  seed = config$seed,
  printFlag = TRUE
)

logged_events_test <- imp_test$loggedEvents
if (!is.null(logged_events_test)) {
  write_csv(as.data.frame(logged_events_test), file.path(config$output_dir, "02_mice_test_logged_events.csv"))
}

save_chain_plot(imp_test, file.path(config$output_dir, "02_mice_test_chain_plot.png"))

##############################
# 12. FINAL RUN
##############################
set.seed(config$seed)
imp <- mice(
  data = imp_df,
  m = config$final_m,
  maxit = config$final_maxit,
  method = meth,
  predictorMatrix = pred,
  seed = config$seed,
  printFlag = TRUE
)

##############################
# 13. SAVE OUTPUTS
##############################
saveRDS(imp, file.path(config$processed_dir, "sensitivity_mice_imp.rds"))

logged_events_final <- imp$loggedEvents
if (!is.null(logged_events_final)) {
  write_csv(as.data.frame(logged_events_final), file.path(config$output_dir, "02_mice_final_logged_events.csv"))
}

save_chain_plot(imp, file.path(config$output_dir, "02_mice_final_chain_plot.png"))

# Save a single completed dataset for quick inspection.
complete1 <- complete(imp, 1)
write_csv(complete1, file.path(config$output_dir, "02_complete_dataset_imp1.csv"))

writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "02_session_info.txt"))
cat("02_imputation.R completed successfully.\n")
