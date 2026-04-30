############################################################
# 02_build_pregnancy_dataset.R
# Purpose: Link GDM to an index pregnancy and derive all
#          predictor variables used in the final model.
#
# Predictors derived:
#   maternal age, IMD/SES, substance use, gestational age,
#   parity, birth weight, SMM, previous complications,
#   hypertensive disorders, BMI, ethnicity, region,
#   smoking, alcohol, treatment method, delivery method
############################################################

##############################
# 1. PACKAGE SETUP
##############################
required_packages <- c("dplyr", "lubridate", "data.table", "openxlsx", "DBI", "odbc")
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
  output_dir    = "outputs"
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

##############################
# 3. LOAD DATA
##############################
gdm_combined   <- readRDS(file.path(config$processed_dir, "gdm_combined.rds"))
pcr_observation <- readRDS(file.path(config$processed_dir, "pcr_observation.rds"))
hes_diagnoses  <- readRDS(file.path(config$processed_dir, "hes_diagnoses.rds"))
hes_maternity  <- readRDS(file.path(config$processed_dir, "hes_maternity.rds"))
hes_procedures <- readRDS(file.path(config$processed_dir, "hes_procedures.rds"))
pcr_patient    <- readRDS(file.path(config$processed_dir, "pcr_patient.rds"))
ons_patient    <- readRDS(file.path(config$processed_dir, "ons_patient.rds"))
imd_patient    <- readRDS(file.path(config$processed_dir, "imd_patient.rds"))
hes_patient    <- readRDS(file.path(config$processed_dir, "hes_patient.rds"))
pcr_practice   <- readRDS(file.path(config$processed_dir, "pcr_practice.rds"))
all_codes      <- readRDS(file.path(config$processed_dir, "all_codes.rds"))

##############################
# 4. IDENTIFY DELIVERY DATES
##############################
delivery_procedure_codes <- all_codes %>%
  filter(Criteria == "Delivery - procedure", Code_Type == "OPCS-4") %>%
  pull(Code)

delivery_diagnosis_codes <- all_codes %>%
  filter(Criteria == "Delivery - HES", Code_Type == "ICD10") %>%
  pull(Code)

delivery_from_procedures <- hes_procedures %>%
  filter(OPCS %in% delivery_procedure_codes) %>%
  transmute(patid, delivery_date = as.Date(evdate, format = "%Y-%m-%d")) %>%
  distinct()

delivery_from_diagnoses <- hes_diagnoses %>%
  filter(ICD %in% delivery_diagnosis_codes) %>%
  transmute(patid, delivery_date = as.Date(epiend, format = "%Y-%m-%d")) %>%
  distinct()

delivery_from_maternity <- hes_maternity %>%
  transmute(patid, delivery_date = as.Date(epiend, format = "%Y-%m-%d")) %>%
  distinct()

all_delivery_dates <- bind_rows(
  delivery_from_procedures,
  delivery_from_diagnoses,
  delivery_from_maternity
) %>%
  distinct() %>%
  arrange(patid, delivery_date)

rm(delivery_procedure_codes, delivery_diagnosis_codes,
   delivery_from_procedures, delivery_from_diagnoses, delivery_from_maternity)

##############################
# 5. LINK DELIVERY TO GDM INDEX PREGNANCY
##############################
pregnancy_identified <- all_delivery_dates %>%
  full_join(as.data.frame(gdm_combined), by = "patid") %>%
  mutate(
    # A delivery is considered valid if it occurred 30-240 days after GDM diagnosis.
    valid_pregnancy = as.integer(
      delivery_date >= min_diag_date + 30 &
        delivery_date <= min_diag_date + 240
    ),
    # If invalid, estimate diagnosis as delivery minus 14 weeks (98 days).
    final_diag_date = as.Date(ifelse(
      valid_pregnancy == 1, min_diag_date, delivery_date - 98
    ), origin = "1970-01-01"),
    diff_days = abs(as.numeric(final_diag_date - min_diag_date))
  ) %>%
  group_by(patid) %>%
  slice_min(order_by = diff_days, with_ties = FALSE) %>%
  ungroup()

# Estimate parity using all delivery dates.
pregnancy_identified <- pregnancy_identified %>%
  left_join(all_delivery_dates, by = "patid", suffix = c("_valid", "_all")) %>%
  group_by(patid) %>%
  mutate(
    had_previous_delivery = any(delivery_date_all < delivery_date_valid - 294, na.rm = TRUE),
    primiparity_estimated = ifelse(had_previous_delivery, "no", "yes")
  ) %>%
  ungroup() %>%
  select(-delivery_date_all) %>%
  filter(final_diag_date >= as.Date("2010-01-01")) %>%
  distinct()

##############################
# 6. INITIALISE PREGNANCY DATASET
##############################
pregnancy_dataset <- pregnancy_identified %>%
  mutate(delivery_date = as.Date(delivery_date_valid)) %>%
  select(-delivery_date_valid) %>%
  rename(diag_date = final_diag_date) %>%
  select(patid, delivery_date, diag_date, primiparity_estimated) %>%
  mutate(patid = as.character(patid))

##############################
# 7. MATERNAL AGE
##############################
pcr_patient_filtered <- pcr_patient %>%
  filter(patid %in% gdm_combined$patid)

pregnancy_dataset <- pregnancy_dataset %>%
  left_join(select(pcr_patient_filtered, patid, yob), by = "patid") %>%
  mutate(
    maternal_age = as.numeric(year(delivery_date) - yob),
    maternal_age = case_when(
      maternal_age < 15 | maternal_age > 60 ~ NA_real_,
      TRUE ~ maternal_age
    )
  )

##############################
# 8. SOCIOECONOMIC STATUS (IMD)
##############################
pregnancy_dataset <- pregnancy_dataset %>%
  left_join(select(imd_patient, patid, e2019_imd_10), by = "patid") %>%
  rename(ses = e2019_imd_10) %>%
  mutate(
    socioeconomic_deprivation = case_when(
      ses >= 8 ~ "yes",
      ses < 8  ~ "no",
      TRUE     ~ NA_character_
    )
  )

rm(imd_patient)

##############################
# 9. SUBSTANCE USE
##############################
medcodeid_substance_use <- all_codes %>%
  filter(Criteria == "Substance abuse - pcr", Code_Type == "MedCodeId", Inclusion == "yes") %>%
  pull(Code)

icd_substance_use <- all_codes %>%
  filter(Criteria == "Substance abuse - HES", Code_Type == "ICD10", Inclusion == "yes") %>%
  pull(Code)

substance_use_pcr <- pcr_observation %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  mutate(obsdate = as.Date(obsdate)) %>%
  filter(medcodeid %in% medcodeid_substance_use, obsdate < delivery_date) %>%
  distinct(patid) %>%
  mutate(substance_use_pcr = "yes")

substance_use_hes <- hes_diagnoses %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  mutate(epistart = as.Date(epistart)) %>%
  filter(ICD %in% icd_substance_use, epistart < delivery_date) %>%
  distinct(patid) %>%
  mutate(substance_use_hes = "yes")

pregnancy_dataset <- pregnancy_dataset %>%
  left_join(substance_use_pcr, by = "patid") %>%
  left_join(substance_use_hes, by = "patid") %>%
  mutate(
    substance_use = case_when(
      !is.na(substance_use_pcr) | !is.na(substance_use_hes) ~ "yes",
      TRUE ~ "no"
    )
  ) %>%
  select(-substance_use_pcr, -substance_use_hes)

rm(medcodeid_substance_use, icd_substance_use, substance_use_pcr, substance_use_hes)

##############################
# 10. GESTATIONAL AGE, PARITY & BIRTH WEIGHT
##############################
pregnancy_dataset <- pregnancy_dataset %>%
  left_join(
    select(hes_maternity, patid, gestat, numpreg, epistart, birweit, delmeth),
    by = "patid"
  ) %>%
  mutate(
    gest_age = case_when(gestat < 24 | gestat > 50 ~ NA_integer_, TRUE ~ gestat),
    birthweight = case_when(birweit < 1000 | birweit > 6000 ~ NA_real_, TRUE ~ birweit),
    primiparity = case_when(
      !is.na(numpreg) & numpreg == 0 ~ "yes",
      !is.na(numpreg) & numpreg > 0  ~ "no",
      is.na(numpreg) & primiparity_estimated == "yes" ~ "yes",
      is.na(numpreg) & primiparity_estimated == "no"  ~ "no",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(patid) %>%
  mutate(
    epistart_diff = abs(as.numeric(difftime(delivery_date, epistart, units = "days"))),
    epistart_diff = ifelse(is.na(epistart_diff), Inf, epistart_diff)
  ) %>%
  slice_min(order_by = epistart_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-epistart_diff)

##############################
# 11. SEVERE MATERNAL MORBIDITY (SMM / EMMOI)
##############################
emmoi_diag_map <- bind_rows(lapply(1:17, function(i) {
  tibble(
    code = all_codes %>%
      filter(Criteria == paste0("EMMOI Diag ", i), Code_Type == "ICD10 EMMOI") %>%
      pull(Code),
    emmoi_category = paste0("EMMOI Diag ", i)
  )
}))

emmoi_proc_map <- bind_rows(lapply(1:9, function(i) {
  tibble(
    code = all_codes %>%
      filter(Criteria == paste0("EMMOI Proc ", i), Code_Type == "OPCS-4") %>%
      pull(Code),
    emmoi_category = paste0("EMMOI Proc ", i)
  )
}))

emmoi_events <- bind_rows(
  hes_diagnoses %>%
    semi_join(distinct(select(pregnancy_dataset, patid)), by = "patid") %>%
    mutate(epistart = as.Date(epistart)) %>%
    inner_join(emmoi_diag_map, by = c("ICD" = "code")) %>%
    transmute(patid, event_date = epistart, emmoi_category),
  hes_procedures %>%
    semi_join(distinct(select(pregnancy_dataset, patid)), by = "patid") %>%
    mutate(evdate = as.Date(evdate)) %>%
    inner_join(emmoi_proc_map, by = c("OPCS" = "code")) %>%
    transmute(patid, event_date = evdate, emmoi_category)
) %>%
  distinct()

emmoi_summary <- pregnancy_dataset %>%
  select(patid, delivery_date) %>%
  distinct() %>%
  left_join(
    emmoi_events %>%
      inner_join(distinct(select(pregnancy_dataset, patid, delivery_date)), by = "patid") %>%
      group_by(patid, delivery_date) %>%
      summarise(
        EMMOI_now  = n_distinct(emmoi_category[event_date >= (delivery_date - 294) & event_date <= (delivery_date + 42)]),
        EMMOI_past = n_distinct(emmoi_category[event_date < (delivery_date - 294)]),
        .groups = "drop"
      ),
    by = c("patid", "delivery_date")
  ) %>%
  mutate(
    EMMOI_now  = coalesce(EMMOI_now,  0L),
    EMMOI_past = coalesce(EMMOI_past, 0L)
  )

pregnancy_dataset <- pregnancy_dataset %>%
  left_join(select(emmoi_summary, patid, EMMOI_now, EMMOI_past), by = "patid") %>%
  mutate(smm_now = as.integer(EMMOI_now > 0))

rm(emmoi_diag_map, emmoi_proc_map, emmoi_events, emmoi_summary)

##############################
# 12. PREVIOUS PREGNANCY COMPLICATIONS
##############################
hes_maternity_flag <- hes_maternity %>%
  mutate(epistart = as.Date(epistart)) %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  filter(
    epistart < (delivery_date - 294),
    birstat %in% c(2, 3, 4) | birweit < 2500 | neocare %in% c(1, 2, 3) | gestat < 37
  ) %>%
  distinct(patid) %>%
  mutate(previous_complication_hes = TRUE)

hes_diagnoses_flag <- hes_diagnoses %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  mutate(epistart = as.Date(epistart)) %>%
  filter(ICD == "Z35.2", epistart < (delivery_date - 294)) %>%
  distinct(patid) %>%
  mutate(previous_complication_diag = TRUE)

pregnancy_dataset <- pregnancy_dataset %>%
  left_join(hes_maternity_flag, by = "patid") %>%
  left_join(hes_diagnoses_flag, by = "patid") %>%
  mutate(
    previous_complication_hes  = coalesce(previous_complication_hes,  FALSE),
    previous_complication_diag = coalesce(previous_complication_diag, FALSE),
    previous_complications = case_when(
      is.na(primiparity) ~ NA_character_,
      primiparity == "yes" ~ "Primiparous",
      primiparity == "no" & (previous_complication_hes | previous_complication_diag) ~ "Yes, multiparous",
      primiparity == "no" & EMMOI_past > 0 ~ "Yes, multiparous",
      TRUE ~ "No, multiparous"
    )
  )

rm(hes_maternity_flag, hes_diagnoses_flag)

##############################
# 13. HYPERTENSIVE DISORDERS
##############################
# Helper to pull codes cleanly.
pull_codes <- function(criteria, code_type, data = all_codes) {
  data %>%
    filter(Criteria == criteria, Code_Type == code_type) %>%
    pull(Code)
}

# Code lists.
hyp_code_lists <- list(
  preexisting_icd         = pull_codes("Hypertension_unspecified - HES",    "ICD10"),
  gestational_icd         = pull_codes("Hypertension_gestational - HES",    "ICD10"),
  preeclampsia_icd        = pull_codes("Hypertension_preeclampsia - HES",   "ICD10"),
  superimposed_icd        = pull_codes("Hypertension_superimposed - HES",   "ICD10"),
  preexisting_medcode     = pull_codes("Hypertension_unspecified - pcr",    "MedCodeId"),
  gestational_medcode     = pull_codes("Hypertension_gestational - pcr",    "MedCodeId"),
  preeclampsia_medcode    = pull_codes("Hypertension_preeclampsia - pcr",   "MedCodeId"),
  superimposed_medcode    = pull_codes("Hypertension_superimposed - pcr",   "MedCodeId")
)

# Function to flag a hypertension type from PCR.
flag_hyp_pcr <- function(codes, flag_name) {
  pcr_observation %>%
    inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
    mutate(obsdate = as.Date(obsdate)) %>%
    filter(
      medcodeid %in% codes,
      obsdate >= (delivery_date - 294),
      obsdate <= delivery_date
    ) %>%
    distinct(patid) %>%
    mutate(!!flag_name := "yes")
}

# Function to flag a hypertension type from HES.
flag_hyp_hes <- function(codes, flag_name) {
  hes_diagnoses %>%
    inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
    mutate(epistart = as.Date(epistart)) %>%
    filter(
      ICD %in% codes,
      epistart >= (delivery_date - 294),
      epistart <= delivery_date
    ) %>%
    distinct(patid) %>%
    mutate(!!flag_name := "yes")
}

hyp_types <- c("preexisting_unspecified", "gestational_hypertension",
                "preeclampsia_hellp_eclampsia", "superimposed_preeclampsia")
medcode_keys <- c("preexisting_medcode", "gestational_medcode",
                  "preeclampsia_medcode", "superimposed_medcode")
icd_keys <- c("preexisting_icd", "gestational_icd",
               "preeclampsia_icd", "superimposed_icd")

for (i in seq_along(hyp_types)) {
  flag_pcr <- flag_hyp_pcr(hyp_code_lists[[medcode_keys[i]]], hyp_types[i])
  flag_hes <- flag_hyp_hes(hyp_code_lists[[icd_keys[i]]],    hyp_types[i])
  flag_combined <- bind_rows(flag_pcr, flag_hes) %>% distinct()
  pregnancy_dataset <- pregnancy_dataset %>%
    left_join(flag_combined, by = "patid") %>%
    mutate(!!hyp_types[i] := ifelse(is.na(.data[[hyp_types[i]]]), "no", "yes"))
}

pregnancy_dataset <- pregnancy_dataset %>%
  mutate(
    hyp_disorders = case_when(
      superimposed_preeclampsia   == "yes" ~ "superimposed_preeclampsia",
      preeclampsia_hellp_eclampsia == "yes" ~ "preeclampsia_hellp_eclampsia",
      gestational_hypertension    == "yes" ~ "gestational_hypertension",
      preexisting_unspecified     == "yes" ~ "preexisting_unspecified",
      TRUE ~ "no"
    )
  )

rm(hyp_code_lists, hyp_types, medcode_keys, icd_keys)

##############################
# 14. BMI
##############################
bmi_value_codes          <- pull_codes("bmi_value",              "MedCodeId")
bmi_centile_codes        <- pull_codes("bmi_centile",            "MedCodeId")
bmi_cat_obese_codes      <- pull_codes("bmi_category_obese",     "MedCodeId")
bmi_cat_morbid_codes     <- pull_codes("bmi_category_morbid_obese", "MedCodeId")
bmi_cat_normal_codes     <- pull_codes("bmi_category_normal",    "MedCodeId")
bmi_cat_underweight_codes <- pull_codes("bmi_category_underweight", "MedCodeId")
bmi_cat_overweight_codes <- pull_codes("bmi_category_overweight", "MedCodeId")

bmi_value_df <- pcr_observation %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  filter(medcodeid %in% bmi_value_codes, obsdate < (delivery_date - 294)) %>%
  group_by(patid) %>%
  slice_max(order_by = obsdate, n = 1, with_ties = FALSE) %>%
  select(patid, bmi_value = value) %>%
  ungroup()

bmi_centile_df <- pcr_observation %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  filter(medcodeid %in% bmi_centile_codes, obsdate < (delivery_date - 294)) %>%
  group_by(patid) %>%
  slice_max(order_by = obsdate, n = 1, with_ties = FALSE) %>%
  select(patid, bmi_centile = value) %>%
  ungroup()

bmi_category_df <- pcr_observation %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  filter(
    medcodeid %in% c(bmi_cat_obese_codes, bmi_cat_morbid_codes, bmi_cat_normal_codes,
                     bmi_cat_overweight_codes, bmi_cat_underweight_codes),
    obsdate < (delivery_date - 294)
  ) %>%
  group_by(patid) %>%
  slice_max(order_by = obsdate, n = 1, with_ties = FALSE) %>%
  mutate(
    obese       = ifelse(medcodeid %in% bmi_cat_obese_codes,       "yes", "no"),
    morbid_obese = ifelse(medcodeid %in% bmi_cat_morbid_codes,     "yes", "no"),
    normal      = ifelse(medcodeid %in% bmi_cat_normal_codes,      "yes", "no"),
    overweight  = ifelse(medcodeid %in% bmi_cat_overweight_codes,  "yes", "no"),
    underweight = ifelse(medcodeid %in% bmi_cat_underweight_codes, "yes", "no")
  ) %>%
  select(patid, obese, morbid_obese, normal, overweight, underweight) %>%
  ungroup()

bmi_final <- bmi_value_df %>%
  full_join(bmi_centile_df, by = "patid") %>%
  full_join(bmi_category_df, by = "patid") %>%
  mutate(
    bmi_value = case_when(bmi_value < 15 | bmi_value > 60 ~ NA_real_, TRUE ~ bmi_value),
    bmi_category_final = case_when(
      !is.na(bmi_value) & bmi_value < 18.5                       ~ "underweight",
      !is.na(bmi_value) & bmi_value >= 18.5 & bmi_value < 25    ~ "normal",
      !is.na(bmi_value) & bmi_value >= 25   & bmi_value < 30    ~ "overweight",
      !is.na(bmi_value) & bmi_value >= 30   & bmi_value < 40    ~ "obese",
      !is.na(bmi_value) & bmi_value >= 40                        ~ "morbid obese",
      is.na(bmi_value) & !is.na(bmi_centile) & bmi_centile < 5       ~ "underweight",
      is.na(bmi_value) & !is.na(bmi_centile) & bmi_centile < 85      ~ "normal",
      is.na(bmi_value) & !is.na(bmi_centile) & bmi_centile < 95      ~ "overweight",
      is.na(bmi_value) & !is.na(bmi_centile) & bmi_centile < 99.6    ~ "obese",
      is.na(bmi_value) & !is.na(bmi_centile) & bmi_centile >= 99.6   ~ "morbid obese",
      morbid_obese == "yes" ~ "morbid obese",
      obese       == "yes"  ~ "obese",
      overweight  == "yes"  ~ "overweight",
      normal      == "yes"  ~ "normal",
      underweight == "yes"  ~ "underweight",
      TRUE ~ NA_character_
    )
  )

pregnancy_dataset <- pregnancy_dataset %>%
  left_join(select(bmi_final, patid, bmi_category_final, bmi_value), by = "patid")

rm(bmi_value_df, bmi_centile_df, bmi_category_df, bmi_final,
   bmi_value_codes, bmi_centile_codes, bmi_cat_obese_codes, bmi_cat_morbid_codes,
   bmi_cat_normal_codes, bmi_cat_underweight_codes, bmi_cat_overweight_codes)

##############################
# 15. ETHNICITY
##############################
hes_patient_eth <- hes_patient %>%
  select(-any_of("gen_hesid")) %>%
  filter(patid %in% pregnancy_dataset$patid) %>%
  mutate(
    gen_ethnicity = case_when(
      gen_ethnicity == "White"                                                            ~ "White",
      gen_ethnicity %in% c("Bl_Afric", "Bl_Carib", "Bl_Other")                          ~ "Black",
      gen_ethnicity %in% c("Pakistani", "Indian", "Chinese", "Bangladesi", "Oth_Asian")  ~ "Asian",
      gen_ethnicity == "Mixed"                                                            ~ "Mixed",
      gen_ethnicity == "Other"                                                            ~ "Other",
      TRUE ~ NA_character_
    )
  )

eth_code_map <- list(
  White   = pull_codes("ethnicity_white",   "MedCodeId"),
  Black   = pull_codes("ethnicity_black",   "MedCodeId"),
  Asian   = pull_codes("ethnicity_asian",   "MedCodeId"),
  Mixed   = pull_codes("ethnicity_mixed",   "MedCodeId"),
  Other   = pull_codes("ethnicity_other",   "MedCodeId")
)

most_common <- function(x) {
  if (all(is.na(x))) return(NA_character_)
  names(sort(table(na.omit(x)), decreasing = TRUE))[1]
}

missing_eth_ids <- hes_patient_eth %>% filter(is.na(gen_ethnicity)) %>% select(patid)

pcr_eth_fill <- pcr_observation %>%
  semi_join(missing_eth_ids, by = "patid") %>%
  mutate(
    ethnicity = case_when(
      medcodeid %in% eth_code_map$White ~ "White",
      medcodeid %in% eth_code_map$Black ~ "Black",
      medcodeid %in% eth_code_map$Asian ~ "Asian",
      medcodeid %in% eth_code_map$Mixed ~ "Mixed",
      medcodeid %in% eth_code_map$Other ~ "Other",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(patid) %>%
  summarise(ethnicity = most_common(ethnicity), .groups = "drop")

hes_patient_eth <- hes_patient_eth %>%
  left_join(pcr_eth_fill, by = "patid") %>%
  mutate(gen_ethnicity = ifelse(is.na(gen_ethnicity), ethnicity, gen_ethnicity)) %>%
  select(-ethnicity)

pregnancy_dataset <- pregnancy_dataset %>%
  left_join(hes_patient_eth, by = "patid")

rm(hes_patient_eth, pcr_eth_fill, missing_eth_ids, eth_code_map)

##############################
# 16. REGION
##############################
pregnancy_dataset <- pregnancy_dataset %>%
  left_join(pcr_practice, by = "pracid")

##############################
# 17. DATE OF DEATH & END OF REGISTRATION
##############################
pregnancy_dataset <- pregnancy_dataset %>%
  left_join(select(ons_patient, patid, dod), by = "patid") %>%
  left_join(select(pcr_patient_filtered, patid, regenddate), by = "patid")

##############################
# 18. SMOKING STATUS
##############################
smoking_codes <- list(
  current_pcr = pull_codes("current smoker - pcr", "MedCodeId"),
  former_pcr  = pull_codes("ex-smoker - pcr",      "MedCodeId"),
  never_pcr   = pull_codes("non-smoker - pcr",     "MedCodeId"),
  current_hes = pull_codes("current smoker - HES", "ICD10"),
  former_hes  = pull_codes("ex-smoker - HES",      "ICD10"),
  never_hes   = pull_codes("non-smoker - HES",     "ICD10")
)

pcr_smoking <- pcr_observation %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  mutate(
    event_date   = as.Date(obsdate),
    smoking_type = case_when(
      medcodeid %in% smoking_codes$current_pcr ~ "current",
      medcodeid %in% smoking_codes$former_pcr  ~ "former",
      medcodeid %in% smoking_codes$never_pcr   ~ "never",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(smoking_type), event_date < delivery_date) %>%
  select(patid, event_date, smoking_type)

hes_smoking <- hes_diagnoses %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  mutate(
    event_date   = as.Date(epistart),
    smoking_type = case_when(
      ICD %in% smoking_codes$current_hes ~ "current",
      ICD %in% smoking_codes$former_hes  ~ "former",
      ICD %in% smoking_codes$never_hes   ~ "never",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(smoking_type), event_date < delivery_date) %>%
  select(patid, event_date, smoking_type)

latest_smoking <- bind_rows(pcr_smoking, hes_smoking) %>%
  arrange(patid, desc(event_date)) %>%
  group_by(patid) %>%
  slice(1) %>%
  ungroup() %>%
  rename(smoking_status = smoking_type)

pregnancy_dataset <- pregnancy_dataset %>%
  left_join(select(latest_smoking, patid, smoking_status), by = "patid") %>%
  mutate(
    smoking_status = factor(
      ifelse(is.na(smoking_status), "unknown", smoking_status),
      levels = c("current", "former", "never", "unknown")
    )
  )

rm(smoking_codes, pcr_smoking, hes_smoking, latest_smoking)

##############################
# 19. ALCOHOL STATUS
##############################
alcohol_codes <- list(
  current_pcr = pull_codes("current - pcr", "MedCodeId"),
  former_pcr  = pull_codes("ex - pcr",      "MedCodeId"),
  never_pcr   = pull_codes("non - pcr",     "MedCodeId"),
  current_hes = all_codes %>% filter(Criteria == "current - HES", Code_Type %in% c("ICD10", "ICD-10")) %>% pull(Code),
  former_hes  = all_codes %>% filter(Criteria == "ex - HES",      Code_Type %in% c("ICD10", "ICD-10")) %>% pull(Code),
  never_hes   = all_codes %>% filter(Criteria == "non - HES",     Code_Type %in% c("ICD10", "ICD-10")) %>% pull(Code)
)

pcr_alcohol <- pcr_observation %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  mutate(
    event_date   = as.Date(obsdate),
    alcohol_type = case_when(
      medcodeid %in% alcohol_codes$current_pcr ~ "current",
      medcodeid %in% alcohol_codes$former_pcr  ~ "former",
      medcodeid %in% alcohol_codes$never_pcr   ~ "never",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(alcohol_type), event_date < delivery_date) %>%
  select(patid, event_date, alcohol_type)

hes_alcohol <- hes_diagnoses %>%
  inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
  mutate(
    event_date   = as.Date(epistart),
    alcohol_type = case_when(
      ICD %in% alcohol_codes$current_hes ~ "current",
      ICD %in% alcohol_codes$former_hes  ~ "former",
      ICD %in% alcohol_codes$never_hes   ~ "never",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(alcohol_type), event_date < delivery_date) %>%
  select(patid, event_date, alcohol_type)

latest_alcohol <- bind_rows(pcr_alcohol, hes_alcohol) %>%
  arrange(patid, desc(event_date)) %>%
  group_by(patid) %>%
  slice(1) %>%
  ungroup() %>%
  rename(alcohol_status = alcohol_type)

pregnancy_dataset <- pregnancy_dataset %>%
  left_join(select(latest_alcohol, patid, alcohol_status), by = "patid") %>%
  mutate(
    alcohol_status = factor(
      ifelse(is.na(alcohol_status), "unknown", alcohol_status),
      levels = c("current", "former", "never", "unknown")
    )
  )

rm(alcohol_codes, pcr_alcohol, hes_alcohol, latest_alcohol)

##############################
# 20. TREATMENT METHOD
##############################
# Drug issue table must be downloaded from the database.
# We reuse the connection from script 01 if running interactively.
# If running standalone, reconnect using environment variables.
if (!exists("con") || !dbIsValid(con)) {
  con <- dbConnect(
    odbc::odbc(),
    Driver   = Sys.getenv("CPRD_DB_DRIVER"),
    Server   = Sys.getenv("CPRD_DB_SERVER"),
    Database = Sys.getenv("CPRD_DB_NAME"),
    UID      = Sys.getenv("CPRD_DB_UID"),
    PWD      = Sys.getenv("CPRD_DB_PWD"),
    Port     = as.integer(Sys.getenv("CPRD_DB_PORT", "1433"))
  )
  on.exit(if (exists("con") && dbIsValid(con)) dbDisconnect(con), add = TRUE)
}

patid_list   <- unique(pregnancy_dataset$patid)
patid_chunks <- split(patid_list, ceiling(seq_along(patid_list) / 10000))

drug_issue <- rbindlist(
  lapply(patid_chunks, function(chunk) {
    chunk_sql <- paste0("(", paste(chunk, collapse = ","), ")")
    dbGetQuery(con, paste0("SELECT * FROM eRAP_22_002383_Predict_GDB_Data.Drug_Issue WHERE patid IN ", chunk_sql))
  }),
  use.names = TRUE, fill = TRUE
) %>%
  mutate(
    patid      = as.character(patid),
    prodcodeid = as.character(prodcodeid),
    issuedate  = as.Date(issuedate, format = "%Y-%m-%d")
  ) %>%
  select(patid, prodcodeid, issuedate)

rm(patid_chunks, patid_list)

insulin_codes  <- all_codes %>% filter(Criteria == "Insulin",   grepl("ProdCodeId", Code_Type, ignore.case = TRUE)) %>% pull(Code) %>% unique() %>% as.character()
metformin_codes <- all_codes %>% filter(Criteria == "Metformin", grepl("ProdCodeId", Code_Type, ignore.case = TRUE)) %>% pull(Code) %>% unique() %>% as.character()

flag_drug_window <- function(drug_data, codes, flag_col) {
  drug_data %>%
    inner_join(select(pregnancy_dataset, patid, delivery_date), by = "patid") %>%
    mutate(delivery_date = as.Date(delivery_date)) %>%
    filter(
      !is.na(issuedate),
      prodcodeid %in% codes,
      issuedate >= (delivery_date - days(240)),
      issuedate <= delivery_date
    ) %>%
    distinct(patid) %>%
    mutate(!!flag_col := "yes")
}

insulin_window  <- flag_drug_window(drug_issue, insulin_codes,  "insulin")
metformin_window <- flag_drug_window(drug_issue, metformin_codes, "metformin")

pregnancy_dataset <- pregnancy_dataset %>%
  left_join(select(insulin_window,  patid, insulin),  by = "patid") %>%
  left_join(select(metformin_window, patid, metformin), by = "patid") %>%
  mutate(
    insulin   = if_else(is.na(insulin),   "no", insulin),
    metformin = if_else(is.na(metformin), "no", metformin),
    treatment_method = factor(
      case_when(
        insulin == "yes" ~ "insulin",
        metformin == "yes" ~ "metformin",
        TRUE ~ "lifestyle"
      ),
      levels = c("insulin", "metformin", "lifestyle")
    )
  )

rm(drug_issue, insulin_codes, metformin_codes, insulin_window, metformin_window)

##############################
# 21. DELIVERY METHOD
##############################
spontaneous_proc_codes <- pull_codes("Spontaneous - procedure", "OPCS-4")
instrumental_proc_codes <- pull_codes("Instrumental - procedure", "OPCS-4")
spontaneous_hes_codes   <- pull_codes("Spontaneous - HES", "ICD10")
instrumental_hes_codes  <- pull_codes("Instrumental - HES", "ICD10")

# Start from maternity delmeth field.
pregnancy_dataset <- pregnancy_dataset %>%
  mutate(
    delivery_method = case_when(
      is.na(delmeth)       ~ NA_character_,
      delmeth %in% c(0, 1) ~ "spontaneous",
      delmeth %in% 2:9     ~ "instrumental",
      TRUE ~ NA_character_
    )
  )

# Fill in missing values from HES procedures and diagnoses.
missing_ids <- pregnancy_dataset %>% filter(is.na(delivery_method)) %>% select(patid, delivery_date)

proc_method <- hes_procedures %>%
  mutate(patid = as.character(patid), evdate = as.Date(evdate)) %>%
  inner_join(missing_ids, by = "patid") %>%
  filter(!is.na(evdate), abs(as.numeric(evdate - delivery_date)) <= 7, OPCS %in% c(spontaneous_proc_codes, instrumental_proc_codes)) %>%
  mutate(
    method_from_proc = case_when(
      OPCS %in% instrumental_proc_codes ~ "instrumental",
      OPCS %in% spontaneous_proc_codes  ~ "spontaneous",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(patid, desc(method_from_proc == "instrumental")) %>%
  group_by(patid) %>% slice(1) %>% ungroup() %>%
  select(patid, method_from_proc)

diag_method <- hes_diagnoses %>%
  mutate(patid = as.character(patid), epistart = as.Date(epistart), epiend = as.Date(epiend)) %>%
  inner_join(missing_ids, by = "patid") %>%
  filter(
    ICD %in% c(spontaneous_hes_codes, instrumental_hes_codes),
    ((!is.na(epistart) & abs(as.numeric(epistart - delivery_date)) <= 7) |
       (!is.na(epiend) & abs(as.numeric(epiend - delivery_date)) <= 7))
  ) %>%
  mutate(
    method_from_diag = case_when(
      ICD %in% instrumental_hes_codes ~ "instrumental",
      ICD %in% spontaneous_hes_codes  ~ "spontaneous",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(patid, desc(method_from_diag == "instrumental")) %>%
  group_by(patid) %>% slice(1) %>% ungroup() %>%
  select(patid, method_from_diag)

pregnancy_dataset <- pregnancy_dataset %>%
  left_join(proc_method, by = "patid") %>%
  left_join(diag_method, by = "patid") %>%
  mutate(
    delivery_method = case_when(
      !is.na(delivery_method)                        ~ delivery_method,
      method_from_proc == "instrumental"             ~ "instrumental",
      method_from_diag == "instrumental"             ~ "instrumental",
      method_from_proc == "spontaneous"              ~ "spontaneous",
      method_from_diag == "spontaneous"              ~ "spontaneous",
      TRUE ~ NA_character_
    ),
    delivery_method = factor(delivery_method, levels = c("spontaneous", "instrumental"))
  ) %>%
  select(-method_from_proc, -method_from_diag)

rm(spontaneous_proc_codes, instrumental_proc_codes,
   spontaneous_hes_codes, instrumental_hes_codes,
   missing_ids, proc_method, diag_method)

##############################
# 22. SAVE
##############################
saveRDS(pregnancy_dataset, file.path(config$processed_dir, "pregnancy_dataset.rds"))
writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "02_session_info.txt"))
cat("02_build_pregnancy_dataset.R completed successfully.\n")
