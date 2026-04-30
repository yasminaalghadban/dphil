############################################################
# 03_finalize_outcomes.R
# Purpose: Derive diabetes outcomes after delivery, create
#          an exclusion list for prior diabetes before the
#          index pregnancy, and export the final datasets.
############################################################

##############################
# 1. PACKAGE SETUP
##############################
required_packages <- c("dplyr", "lubridate", "openxlsx")
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
  followup_years = 10
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

##############################
# 3. LOAD DATA
##############################
pregnancy_dataset <- readRDS(file.path(config$processed_dir, "pregnancy_dataset.rds"))
pcr_observation   <- readRDS(file.path(config$processed_dir, "pcr_observation.rds"))
hes_diagnoses     <- readRDS(file.path(config$processed_dir, "hes_diagnoses.rds"))
ons_patient       <- readRDS(file.path(config$processed_dir, "ons_patient.rds"))
all_codes         <- readRDS(file.path(config$processed_dir, "all_codes.rds"))

##############################
# 4. HELPER FUNCTIONS
##############################
pull_codes <- function(criteria, code_type, data = all_codes) {
  data %>%
    filter(Criteria == criteria, Code_Type == code_type) %>%
    pull(Code)
}

safe_pmin_date <- function(...) {
  x <- do.call(pmin, c(lapply(list(...), as.numeric), na.rm = TRUE))
  x[is.infinite(x)] <- NA_real_
  as.Date(x, origin = "1970-01-01")
}

summarise_outcome_source_hes <- function(data, delivery_df, diabetes_codes, complication_codes, source_prefix = "hes") {
  data %>%
    inner_join(delivery_df, by = "patid") %>%
    mutate(epistart = as.Date(epistart)) %>%
    filter(epistart > delivery_date) %>%
    group_by(patid) %>%
    summarise(
      first_event_date_complication = {
        idx <- ICD %in% complication_codes$all
        if (any(idx, na.rm = TRUE)) min(epistart[idx], na.rm = TRUE) else as.Date(NA)
      },
      first_event_date_no_complication = {
        idx <- ICD %in% diabetes_codes
        if (any(idx, na.rm = TRUE)) min(epistart[idx], na.rm = TRUE) else as.Date(NA)
      },
      diabetes_complications_coma        = ifelse(any(ICD %in% complication_codes$coma,        na.rm = TRUE), "yes", "no"),
      diabetes_complications_acidosis    = ifelse(any(ICD %in% complication_codes$acidosis,    na.rm = TRUE), "yes", "no"),
      diabetes_complications_renal       = ifelse(any(ICD %in% complication_codes$renal,       na.rm = TRUE), "yes", "no"),
      diabetes_complications_ophtalmo    = ifelse(any(ICD %in% complication_codes$ophtalmo,    na.rm = TRUE), "yes", "no"),
      diabetes_complications_neuro       = ifelse(any(ICD %in% complication_codes$neuro,       na.rm = TRUE), "yes", "no"),
      diabetes_complications_circulation = ifelse(any(ICD %in% complication_codes$circulation, na.rm = TRUE), "yes", "no"),
      diabetes_complications_other       = ifelse(any(ICD %in% complication_codes$other,       na.rm = TRUE), "yes", "no"),
      diabetes_complications_multiple    = ifelse(any(ICD %in% complication_codes$multiple,    na.rm = TRUE), "yes", "no"),
      diabetes_no_complications          = ifelse(any(ICD %in% diabetes_codes,                  na.rm = TRUE), "yes", "no"),
      .groups = "drop"
    )
}

summarise_outcome_source_pcr <- function(data, delivery_df, diabetes_codes, complication_codes) {
  data %>%
    inner_join(delivery_df, by = "patid") %>%
    mutate(obsdate = as.Date(obsdate)) %>%
    filter(obsdate > delivery_date) %>%
    group_by(patid) %>%
    summarise(
      first_event_date_complication = {
        idx <- medcodeid %in% complication_codes$all
        if (any(idx, na.rm = TRUE)) min(obsdate[idx], na.rm = TRUE) else as.Date(NA)
      },
      first_event_date_no_complication = {
        idx <- medcodeid %in% diabetes_codes
        if (any(idx, na.rm = TRUE)) min(obsdate[idx], na.rm = TRUE) else as.Date(NA)
      },
      diabetes_complications_coma        = ifelse(any(medcodeid %in% complication_codes$coma,        na.rm = TRUE), "yes", "no"),
      diabetes_complications_acidosis    = ifelse(any(medcodeid %in% complication_codes$acidosis,    na.rm = TRUE), "yes", "no"),
      diabetes_complications_renal       = ifelse(any(medcodeid %in% complication_codes$renal,       na.rm = TRUE), "yes", "no"),
      diabetes_complications_ophtalmo    = ifelse(any(medcodeid %in% complication_codes$ophtalmo,    na.rm = TRUE), "yes", "no"),
      diabetes_complications_neuro       = ifelse(any(medcodeid %in% complication_codes$neuro,       na.rm = TRUE), "yes", "no"),
      diabetes_complications_circulation = ifelse(any(medcodeid %in% complication_codes$circulation, na.rm = TRUE), "yes", "no"),
      diabetes_complications_other       = ifelse(any(medcodeid %in% complication_codes$other,       na.rm = TRUE), "yes", "no"),
      diabetes_complications_multiple    = ifelse(any(medcodeid %in% complication_codes$multiple,    na.rm = TRUE), "yes", "no"),
      diabetes_no_complications          = ifelse(any(medcodeid %in% diabetes_codes,                  na.rm = TRUE), "yes", "no"),
      .groups = "drop"
    )
}

summarise_outcome_source_ons <- function(data, delivery_df, diabetes_codes, complication_codes) {
  data %>%
    inner_join(delivery_df, by = "patid") %>%
    mutate(dod = as.Date(dod, origin = "1970-01-01")) %>%
    filter(dod > delivery_date) %>%
    mutate(
      any_diabetes_complication = if_any(starts_with("cause"), ~ . %in% complication_codes$all),
      any_diabetes_no_complication = if_any(starts_with("cause"), ~ . %in% diabetes_codes),
      diabetes_complications_coma_row        = if_any(starts_with("cause"), ~ . %in% complication_codes$coma),
      diabetes_complications_acidosis_row    = if_any(starts_with("cause"), ~ . %in% complication_codes$acidosis),
      diabetes_complications_renal_row       = if_any(starts_with("cause"), ~ . %in% complication_codes$renal),
      diabetes_complications_ophtalmo_row    = if_any(starts_with("cause"), ~ . %in% complication_codes$ophtalmo),
      diabetes_complications_neuro_row       = if_any(starts_with("cause"), ~ . %in% complication_codes$neuro),
      diabetes_complications_circulation_row = if_any(starts_with("cause"), ~ . %in% complication_codes$circulation),
      diabetes_complications_other_row       = if_any(starts_with("cause"), ~ . %in% complication_codes$other),
      diabetes_complications_multiple_row    = if_any(starts_with("cause"), ~ . %in% complication_codes$multiple)
    ) %>%
    group_by(patid) %>%
    summarise(
      first_event_date_complication = if (any(any_diabetes_complication, na.rm = TRUE)) min(dod[any_diabetes_complication], na.rm = TRUE) else as.Date(NA),
      first_event_date_no_complication = if (any(any_diabetes_no_complication, na.rm = TRUE)) min(dod[any_diabetes_no_complication], na.rm = TRUE) else as.Date(NA),
      diabetes_complications_coma        = ifelse(any(diabetes_complications_coma_row,        na.rm = TRUE), "yes", "no"),
      diabetes_complications_acidosis    = ifelse(any(diabetes_complications_acidosis_row,    na.rm = TRUE), "yes", "no"),
      diabetes_complications_renal       = ifelse(any(diabetes_complications_renal_row,       na.rm = TRUE), "yes", "no"),
      diabetes_complications_ophtalmo    = ifelse(any(diabetes_complications_ophtalmo_row,    na.rm = TRUE), "yes", "no"),
      diabetes_complications_neuro       = ifelse(any(diabetes_complications_neuro_row,       na.rm = TRUE), "yes", "no"),
      diabetes_complications_circulation = ifelse(any(diabetes_complications_circulation_row, na.rm = TRUE), "yes", "no"),
      diabetes_complications_other       = ifelse(any(diabetes_complications_other_row,       na.rm = TRUE), "yes", "no"),
      diabetes_complications_multiple    = ifelse(any(diabetes_complications_multiple_row,    na.rm = TRUE), "yes", "no"),
      diabetes_no_complications          = ifelse(any(any_diabetes_no_complication,           na.rm = TRUE), "yes", "no"),
      .groups = "drop"
    )
}

merge_outcome_flags <- function(base_data, hes_data, ons_data = NULL, pcr_data = NULL) {
  out <- base_data %>% left_join(hes_data, by = "patid")

  if (!is.null(ons_data)) {
    out <- out %>% left_join(ons_data, by = "patid", suffix = c("", "_ons"))
  }
  if (!is.null(pcr_data)) {
    out <- out %>% left_join(pcr_data, by = "patid", suffix = c("", "_pcr"))
  }

  out %>%
    mutate(
      first_event_date_complication = safe_pmin_date(first_event_date_complication, first_event_date_complication_ons, first_event_date_complication_pcr),
      first_event_date_no_complication = safe_pmin_date(first_event_date_no_complication, first_event_date_no_complication_ons, first_event_date_no_complication_pcr),
      diabetes_complications_coma = ifelse(coalesce(diabetes_complications_coma, "no") == "yes" | coalesce(diabetes_complications_coma_ons, "no") == "yes" | coalesce(diabetes_complications_coma_pcr, "no") == "yes", "yes", "no"),
      diabetes_complications_acidosis = ifelse(coalesce(diabetes_complications_acidosis, "no") == "yes" | coalesce(diabetes_complications_acidosis_ons, "no") == "yes" | coalesce(diabetes_complications_acidosis_pcr, "no") == "yes", "yes", "no"),
      diabetes_complications_renal = ifelse(coalesce(diabetes_complications_renal, "no") == "yes" | coalesce(diabetes_complications_renal_ons, "no") == "yes" | coalesce(diabetes_complications_renal_pcr, "no") == "yes", "yes", "no"),
      diabetes_complications_ophtalmo = ifelse(coalesce(diabetes_complications_ophtalmo, "no") == "yes" | coalesce(diabetes_complications_ophtalmo_ons, "no") == "yes" | coalesce(diabetes_complications_ophtalmo_pcr, "no") == "yes", "yes", "no"),
      diabetes_complications_neuro = ifelse(coalesce(diabetes_complications_neuro, "no") == "yes" | coalesce(diabetes_complications_neuro_ons, "no") == "yes" | coalesce(diabetes_complications_neuro_pcr, "no") == "yes", "yes", "no"),
      diabetes_complications_circulation = ifelse(coalesce(diabetes_complications_circulation, "no") == "yes" | coalesce(diabetes_complications_circulation_ons, "no") == "yes" | coalesce(diabetes_complications_circulation_pcr, "no") == "yes", "yes", "no"),
      diabetes_complications_other = ifelse(coalesce(diabetes_complications_other, "no") == "yes" | coalesce(diabetes_complications_other_ons, "no") == "yes" | coalesce(diabetes_complications_other_pcr, "no") == "yes", "yes", "no"),
      diabetes_complications_multiple = ifelse(coalesce(diabetes_complications_multiple, "no") == "yes" | coalesce(diabetes_complications_multiple_ons, "no") == "yes" | coalesce(diabetes_complications_multiple_pcr, "no") == "yes", "yes", "no"),
      diabetes_no_complications = ifelse(coalesce(diabetes_no_complications, "no") == "yes" | coalesce(diabetes_no_complications_ons, "no") == "yes" | coalesce(diabetes_no_complications_pcr, "no") == "yes", "yes", "no")
    ) %>%
    select(-matches("_ons$"), -matches("_pcr$"))
}

assign_complication_label <- function(data) {
  data %>%
    mutate(
      complications = case_when(
        diabetes_complications_coma == "yes"        ~ "diabetes_complications_coma",
        diabetes_complications_acidosis == "yes"    ~ "diabetes_complications_acidosis",
        diabetes_complications_renal == "yes"       ~ "diabetes_complications_renal",
        diabetes_complications_ophtalmo == "yes"    ~ "diabetes_complications_ophtalmo",
        diabetes_complications_neuro == "yes"       ~ "diabetes_complications_neuro",
        diabetes_complications_circulation == "yes" ~ "diabetes_complications_circulation",
        diabetes_complications_other == "yes"       ~ "diabetes_complications_other",
        diabetes_complications_multiple == "yes"    ~ "diabetes_complications_multiple",
        TRUE ~ "no"
      ),
      diabetes_status = case_when(
        diabetes_no_complications == "yes" ~ "yes",
        complications != "no"              ~ "yes",
        TRUE ~ "no"
      )
    )
}

##############################
# 5. OUTCOME CODE LISTS
##############################
diabetes_codes_icd <- pull_codes("type 2 diabetes - HES", "ICD10")
diabetes_codes_medcodeid <- pull_codes("type 2 diabetes - pcr", "MedCodeId")

complication_codes_icd <- list(
  coma        = pull_codes("diabetes_complications_coma - HES",        "ICD10"),
  acidosis    = pull_codes("diabetes_complications_acidosis - HES",    "ICD10"),
  renal       = pull_codes("diabetes_complications_renal - HES",       "ICD10"),
  ophtalmo    = pull_codes("diabetes_complications_ophtalmo - HES",    "ICD10"),
  neuro       = pull_codes("diabetes_complications_neuro - HES",       "ICD10"),
  circulation = pull_codes("diabetes_complications_circulation - HES", "ICD10"),
  other       = pull_codes("diabetes_complications_other - HES",       "ICD10"),
  multiple    = pull_codes("diabetes_complications_multiple - HES",    "ICD10")
)
complication_codes_icd$all <- unique(unlist(complication_codes_icd))

complication_codes_medcodeid <- list(
  coma        = pull_codes("diabetes_complications_coma - pcr",        "MedCodeId"),
  acidosis    = pull_codes("diabetes_complications_acidosis - pcr",    "MedCodeId"),
  renal       = pull_codes("diabetes_complications_renal - pcr",       "MedCodeId"),
  ophtalmo    = pull_codes("diabetes_complications_ophtalmo - pcr",    "MedCodeId"),
  neuro       = pull_codes("diabetes_complications_neuro - pcr",       "MedCodeId"),
  circulation = pull_codes("diabetes_complications_circulation - pcr", "MedCodeId"),
  other       = pull_codes("diabetes_complications_other - pcr",       "MedCodeId"),
  multiple    = pull_codes("diabetes_complications_multiple - pcr",    "MedCodeId")
)
complication_codes_medcodeid$all <- unique(unlist(complication_codes_medcodeid))

##############################
# 6. POST-DELIVERY OUTCOME DETERMINATION
##############################
delivery_df <- pregnancy_dataset %>% select(patid, delivery_date)

hes_complications <- summarise_outcome_source_hes(hes_diagnoses, delivery_df, diabetes_codes_icd, complication_codes_icd)
ons_complications <- summarise_outcome_source_ons(ons_patient, delivery_df, diabetes_codes_icd, complication_codes_icd)
pcr_complications <- summarise_outcome_source_pcr(pcr_observation, delivery_df, diabetes_codes_medcodeid, complication_codes_medcodeid)

pregnancy_dataset_outcome <- pregnancy_dataset %>%
  merge_outcome_flags(hes_complications, ons_complications, pcr_complications) %>%
  assign_complication_label()

##############################
# 7. PRE-DELIVERY DIABETES EXCLUSION
##############################
hes_prior <- hes_diagnoses %>%
  inner_join(delivery_df, by = "patid") %>%
  mutate(epistart = as.Date(epistart)) %>%
  filter(epistart < delivery_date) %>%
  group_by(patid) %>%
  summarise(
    first_event_date_complication = { idx <- ICD %in% complication_codes_icd$all; if (any(idx, na.rm = TRUE)) min(epistart[idx], na.rm = TRUE) else as.Date(NA) },
    first_event_date_no_complication = { idx <- ICD %in% diabetes_codes_icd; if (any(idx, na.rm = TRUE)) min(epistart[idx], na.rm = TRUE) else as.Date(NA) },
    diabetes_complications_coma        = ifelse(any(ICD %in% complication_codes_icd$coma,        na.rm = TRUE), "yes", "no"),
    diabetes_complications_acidosis    = ifelse(any(ICD %in% complication_codes_icd$acidosis,    na.rm = TRUE), "yes", "no"),
    diabetes_complications_renal       = ifelse(any(ICD %in% complication_codes_icd$renal,       na.rm = TRUE), "yes", "no"),
    diabetes_complications_ophtalmo    = ifelse(any(ICD %in% complication_codes_icd$ophtalmo,    na.rm = TRUE), "yes", "no"),
    diabetes_complications_neuro       = ifelse(any(ICD %in% complication_codes_icd$neuro,       na.rm = TRUE), "yes", "no"),
    diabetes_complications_circulation = ifelse(any(ICD %in% complication_codes_icd$circulation, na.rm = TRUE), "yes", "no"),
    diabetes_complications_other       = ifelse(any(ICD %in% complication_codes_icd$other,       na.rm = TRUE), "yes", "no"),
    diabetes_complications_multiple    = ifelse(any(ICD %in% complication_codes_icd$multiple,    na.rm = TRUE), "yes", "no"),
    diabetes_no_complications          = ifelse(any(ICD %in% diabetes_codes_icd,                  na.rm = TRUE), "yes", "no"),
    .groups = "drop"
  )

pcr_prior <- pcr_observation %>%
  inner_join(delivery_df, by = "patid") %>%
  mutate(obsdate = as.Date(obsdate)) %>%
  filter(obsdate < delivery_date) %>%
  group_by(patid) %>%
  summarise(
    first_event_date_complication = { idx <- medcodeid %in% complication_codes_medcodeid$all; if (any(idx, na.rm = TRUE)) min(obsdate[idx], na.rm = TRUE) else as.Date(NA) },
    first_event_date_no_complication = { idx <- medcodeid %in% diabetes_codes_medcodeid; if (any(idx, na.rm = TRUE)) min(obsdate[idx], na.rm = TRUE) else as.Date(NA) },
    diabetes_complications_coma        = ifelse(any(medcodeid %in% complication_codes_medcodeid$coma,        na.rm = TRUE), "yes", "no"),
    diabetes_complications_acidosis    = ifelse(any(medcodeid %in% complication_codes_medcodeid$acidosis,    na.rm = TRUE), "yes", "no"),
    diabetes_complications_renal       = ifelse(any(medcodeid %in% complication_codes_medcodeid$renal,       na.rm = TRUE), "yes", "no"),
    diabetes_complications_ophtalmo    = ifelse(any(medcodeid %in% complication_codes_medcodeid$ophtalmo,    na.rm = TRUE), "yes", "no"),
    diabetes_complications_neuro       = ifelse(any(medcodeid %in% complication_codes_medcodeid$neuro,       na.rm = TRUE), "yes", "no"),
    diabetes_complications_circulation = ifelse(any(medcodeid %in% complication_codes_medcodeid$circulation, na.rm = TRUE), "yes", "no"),
    diabetes_complications_other       = ifelse(any(medcodeid %in% complication_codes_medcodeid$other,       na.rm = TRUE), "yes", "no"),
    diabetes_complications_multiple    = ifelse(any(medcodeid %in% complication_codes_medcodeid$multiple,    na.rm = TRUE), "yes", "no"),
    diabetes_no_complications          = ifelse(any(medcodeid %in% diabetes_codes_medcodeid,                  na.rm = TRUE), "yes", "no"),
    .groups = "drop"
  )

pregnancy_dataset_exclusion <- pregnancy_dataset %>%
  left_join(hes_prior, by = "patid") %>%
  left_join(pcr_prior, by = "patid", suffix = c("", "_pcr")) %>%
  mutate(
    first_event_date_complication = safe_pmin_date(first_event_date_complication, first_event_date_complication_pcr),
    first_event_date_no_complication = safe_pmin_date(first_event_date_no_complication, first_event_date_no_complication_pcr),
    diabetes_complications_coma = ifelse(coalesce(diabetes_complications_coma, "no") == "yes" | coalesce(diabetes_complications_coma_pcr, "no") == "yes", "yes", "no"),
    diabetes_complications_acidosis = ifelse(coalesce(diabetes_complications_acidosis, "no") == "yes" | coalesce(diabetes_complications_acidosis_pcr, "no") == "yes", "yes", "no"),
    diabetes_complications_renal = ifelse(coalesce(diabetes_complications_renal, "no") == "yes" | coalesce(diabetes_complications_renal_pcr, "no") == "yes", "yes", "no"),
    diabetes_complications_ophtalmo = ifelse(coalesce(diabetes_complications_ophtalmo, "no") == "yes" | coalesce(diabetes_complications_ophtalmo_pcr, "no") == "yes", "yes", "no"),
    diabetes_complications_neuro = ifelse(coalesce(diabetes_complications_neuro, "no") == "yes" | coalesce(diabetes_complications_neuro_pcr, "no") == "yes", "yes", "no"),
    diabetes_complications_circulation = ifelse(coalesce(diabetes_complications_circulation, "no") == "yes" | coalesce(diabetes_complications_circulation_pcr, "no") == "yes", "yes", "no"),
    diabetes_complications_other = ifelse(coalesce(diabetes_complications_other, "no") == "yes" | coalesce(diabetes_complications_other_pcr, "no") == "yes", "yes", "no"),
    diabetes_complications_multiple = ifelse(coalesce(diabetes_complications_multiple, "no") == "yes" | coalesce(diabetes_complications_multiple_pcr, "no") == "yes", "yes", "no"),
    diabetes_no_complications = ifelse(coalesce(diabetes_no_complications, "no") == "yes" | coalesce(diabetes_no_complications_pcr, "no") == "yes", "yes", "no")
  ) %>%
  select(-matches("_pcr$")) %>%
  assign_complication_label()

exclusion_dataset <- pregnancy_dataset_exclusion %>%
  filter(diabetes_status == "yes") %>%
  select(patid, diabetes_status, complications) %>%
  distinct()

##############################
# 8. DEFINE 10-YEAR OUTCOMES
##############################
final_full_dataset <- pregnancy_dataset_outcome %>%
  mutate(
    delivery_date = as.Date(delivery_date),
    first_event_date_no_complication = as.Date(first_event_date_no_complication),
    first_event_date_complication = as.Date(first_event_date_complication)
  ) %>%
  mutate(
    delivery_plus_10yrs = delivery_date %m+% months(12 * config$followup_years),
    diabetes_status_10yr = case_when(
      !is.na(first_event_date_no_complication) & first_event_date_no_complication <= delivery_plus_10yrs ~ "yes",
      TRUE ~ "no"
    ),
    complications_10yr = case_when(
      !is.na(first_event_date_complication) & first_event_date_complication <= delivery_plus_10yrs ~ as.character(complications),
      TRUE ~ "no"
    )
  ) %>%
  mutate(
    diabetes_status_10yr = if_else(complications_10yr != "no", "yes", diabetes_status_10yr)
  ) %>%
  select(-delivery_plus_10yrs)

##############################
# 9. CREATE TRIMMED FINAL ANALYSIS DATASET
##############################
final_dataset <- final_full_dataset %>%
  select(
    patid, pracid, delivery_date, diag_date, regenddate, yob,
    maternal_age, ses, substance_use, gest_age, numpreg, primiparity,
    birweit, birthweight, delivery_method, smm_now, previous_complications,
    hyp_disorders, first_event_date_complication, first_event_date_no_complication,
    complications, diabetes_status, bmi_value, bmi_category_final,
    smoking_status, alcohol_status, gen_ethnicity, region,
    socioeconomic_deprivation, dod, treatment_method,
    diabetes_status_10yr, complications_10yr
  )

##############################
# 10. SAVE OUTPUTS
##############################
saveRDS(exclusion_dataset,    file.path(config$processed_dir, "exclusion_dataset.rds"))
saveRDS(final_full_dataset,   file.path(config$processed_dir, "final_full_dataset.rds"))
saveRDS(final_dataset,        file.path(config$processed_dir, "final_dataset.rds"))

write.xlsx(exclusion_dataset,  file.path(config$output_dir, "03_exclusion_prior_diabetes.xlsx"), rowNames = FALSE)
write.xlsx(final_full_dataset, file.path(config$output_dir, "03_final_full_dataset.xlsx"), rowNames = FALSE)
write.xlsx(final_dataset,      file.path(config$output_dir, "03_final_dataset.xlsx"), rowNames = FALSE)

write.csv(exclusion_dataset,  file.path(config$output_dir, "03_exclusion_prior_diabetes.csv"), row.names = FALSE)
write.csv(final_dataset,      file.path(config$output_dir, "03_final_dataset.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "03_session_info.txt"))

cat("Rows in exclusion dataset:", nrow(exclusion_dataset), "\n")
cat("Rows in final dataset:", nrow(final_dataset), "\n")
cat("03_finalize_outcomes.R completed successfully.\n")
