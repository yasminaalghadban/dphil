############################################################
# 01_clean_data.R
# Purpose: Import raw files, apply exclusions, standardise
#          variables, and save a cleaned cohort for analysis.
############################################################

##############################
# 1. PACKAGE SETUP
##############################
required_packages <- c("dplyr", "readxl", "lubridate")
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) {
  stop("Please install: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

##############################
# 2. CONFIGURATION
##############################
config <- list(
  raw_dir = file.path("data", "raw"),
  processed_dir = file.path("data", "processed"),
  output_dir = "outputs",
  dataset_file = "final_dataset.xlsx",
  exclusion_file = "exclusion_final.xlsx"
)

dir.create(config$processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_path <- file.path(config$raw_dir, config$dataset_file)
exclusion_path <- file.path(config$raw_dir, config$exclusion_file)

stopifnot(file.exists(dataset_path), file.exists(exclusion_path))

##############################
# 3. IMPORT DATA
##############################
dataset <- read_excel(dataset_path)
exclusion <- read_excel(exclusion_path)

##############################
# 4. APPLY COHORT EXCLUSIONS
##############################
cleaned_cohort <- dataset %>%
  mutate(yod = year(diag_date)) %>%
  anti_join(exclusion, by = "patid")

n_after_exclusion_file <- nrow(cleaned_cohort)

# Restrict to maternal age 18-45 years.
cleaned_cohort <- cleaned_cohort %>%
  filter(maternal_age >= 18, maternal_age <= 45)

n_after_age_filter <- nrow(cleaned_cohort)

# Exclude deaths within 42 days of delivery.
cleaned_cohort <- cleaned_cohort %>%
  mutate(dod_delivery_diff = as.numeric(difftime(dod, delivery_date, units = "days"))) %>%
  filter(dod_delivery_diff > 42 | is.na(dod_delivery_diff))

n_after_mortality_filter <- nrow(cleaned_cohort)

##############################
# 5. STANDARDISE VARIABLES
##############################
cleaned_cohort <- cleaned_cohort %>%
  mutate(
    region_name = case_when(
      region == 1 ~ "North East",
      region == 2 ~ "North West",
      region == 3 ~ "Yorkshire and the Humber",
      region == 4 ~ "East Midlands",
      region == 5 ~ "West Midlands",
      region == 6 ~ "East of England",
      region == 7 ~ "London",
      region == 8 ~ "South East",
      region == 9 ~ "South West",
      TRUE ~ NA_character_
    ),
    alcohol_status = case_when(
      alcohol_status == "unknown" ~ NA_character_,
      alcohol_status %in% c("current", "former", "never") ~ alcohol_status,
      TRUE ~ NA_character_
    ),
    smoking_status = case_when(
      smoking_status == "unknown" ~ NA_character_,
      smoking_status %in% c("current", "former", "never") ~ smoking_status,
      TRUE ~ NA_character_
    ),
    smm_now = case_when(
      is.na(smm_now) ~ NA_character_,
      smm_now == 1 ~ "yes",
      smm_now == 0 ~ "no",
      TRUE ~ NA_character_
    ),
    socioeconomic_deprivation = case_when(
      ses >= 8 ~ "yes",
      ses < 8 ~ "no",
      TRUE ~ NA_character_
    ),
    ses_quintile = case_when(
      ses %in% c(1, 2) ~ "1",
      ses %in% c(3, 4) ~ "2",
      ses %in% c(5, 6) ~ "3",
      ses %in% c(7, 8) ~ "4",
      ses %in% c(9, 10) ~ "5",
      TRUE ~ NA_character_
    ),
    age_category = case_when(
      maternal_age < 25 ~ "<25",
      maternal_age >= 25 & maternal_age <= 29 ~ "25–29",
      maternal_age >= 30 & maternal_age <= 34 ~ "30–34",
      maternal_age >= 35 ~ "≥35",
      TRUE ~ NA_character_
    ),
    morbid_obese = case_when(
      bmi_category_final == "morbid obese" ~ "yes",
      !is.na(bmi_category_final) ~ "no",
      TRUE ~ NA_character_
    ),
    # Carry forward complication date if the no-complication date is missing.
    first_event_date_no_complication = if_else(
      is.na(first_event_date_no_complication) & !is.na(first_event_date_complication),
      first_event_date_complication,
      first_event_date_no_complication
    )
  ) %>%
  mutate(
    across(
      c(
        age_category, primiparity, bmi_category_final, substance_use,
        socioeconomic_deprivation, hyp_disorders, diabetes_status_10yr,
        complications_10yr, diabetes_status, complications, gen_ethnicity,
        morbid_obese
      ),
      as.factor
    )
  )

##############################
# 6. QUICK QC SUMMARY
##############################
cohort_summary <- data.frame(
  step = c(
    "After exclusion file",
    "After maternal age filter",
    "After mortality filter"
  ),
  n = c(
    n_after_exclusion_file,
    n_after_age_filter,
    n_after_mortality_filter
  )
)

practice_summary <- cleaned_cohort %>%
  summarise(n_practices = n_distinct(pracid))

##############################
# 7. SAVE OUTPUTS
##############################
saveRDS(cleaned_cohort, file.path(config$processed_dir, "cleaned_cohort.rds"))
write.csv(cohort_summary, file.path(config$output_dir, "01_cohort_summary.csv"), row.names = FALSE)
write.csv(practice_summary, file.path(config$output_dir, "01_practice_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "01_session_info.txt"))

cat("01_clean_data.R completed successfully.\n")
