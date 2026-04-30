############################################################
# 01_extract_gdm_cohort.R
# Purpose: Connect to the database, identify women with a
#          first-ever GDM diagnosis from 2010 onward, and
#          download the core source tables needed for all
#          downstream derivations.
############################################################

##############################
# 1. PACKAGE SETUP
##############################
required_packages <- c("dplyr", "readxl", "lubridate", "odbc", "DBI", "data.table", "ggplot2", "openxlsx")
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) {
  stop("Please install the following packages before running: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

##############################
# 2. CONFIGURATION
##############################
config <- list(
  raw_dir = file.path("data", "raw"),
  processed_dir = file.path("data", "processed"),
  output_dir = "outputs",
  codes_file = "ALL_CODES.xlsx",
  codes_sheet = "Sheet1",
  gdm_start_date = as.Date("2010-01-01"),
  batch_size = 10000
)

dir.create(config$processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

codes_path <- file.path(config$raw_dir, config$codes_file)
stopifnot(file.exists(codes_path))

##############################
# 3. DATABASE CONNECTION
#    Credentials are read from environment variables.
#    Set them in your .Renviron file (never commit credentials to GitHub).
##############################
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

##############################
# 4. LOAD CODE LISTS
##############################
all_codes <- read_excel(codes_path, sheet = config$codes_sheet) %>%
  mutate(
    Criteria  = trimws(Criteria),
    Code_Type = trimws(Code_Type)
  )

##############################
# 5. IDENTIFY GDM CASES
##############################

# -- Primary care (PCR) codes --
medcode_gdm_codes <- all_codes %>%
  filter(Criteria == "GDM - pcr", Code_Type == "MedCodeId") %>%
  pull(Code)

medcode_gdm_sql <- paste0("(", paste(medcode_gdm_codes, collapse = ","), ")")

gdm_pcr <- dbGetQuery(
  con,
  paste0("SELECT * FROM eRAP_22_002383_Predict_GDB_Data.Observation WHERE medcodeid IN ", medcode_gdm_sql)
) %>%
  mutate(
    patid     = as.character(patid),
    medcodeid = as.character(medcodeid),
    obsdate   = as.Date(obsdate, format = "%Y-%m-%d")
  ) %>%
  select(patid, medcodeid, obsdate)

# -- HES codes --
icd_gdm_codes <- all_codes %>%
  filter(Criteria == "GDM - HES", Code_Type == "ICD10") %>%
  pull(Code)

icd_gdm_sql <- paste0("('", paste(icd_gdm_codes, collapse = "','"), "')")

gdm_hes <- dbGetQuery(
  con,
  paste0("SELECT * FROM eRAP_22_002383_Predict_GDB_Data.hes_diagnosis_epi_DM WHERE ICD IN ", icd_gdm_sql)
) %>%
  mutate(
    patid    = as.character(patid),
    ICD      = as.character(ICD),
    epistart = as.Date(epistart, format = "%Y-%m-%d")
  ) %>%
  select(patid, ICD, epistart)

rm(medcode_gdm_sql, icd_gdm_sql, medcode_gdm_codes, icd_gdm_codes)

##############################
# 6. DERIVE FIRST-EVER GDM DATE
##############################
setDT(gdm_hes)
setDT(gdm_pcr)

gdm_hes_first <- gdm_hes[, .(first_event_date_hes = min(epistart, na.rm = TRUE)), by = patid]
gdm_pcr_first <- gdm_pcr[, .(first_event_date_pcr = min(obsdate,  na.rm = TRUE)), by = patid]

gdm_combined <- merge(gdm_hes_first, gdm_pcr_first, by = "patid", all = TRUE)

# Take the earliest date across both sources.
gdm_combined[, min_diag_date := fifelse(
  is.na(first_event_date_hes), first_event_date_pcr,
  fifelse(
    is.na(first_event_date_pcr), first_event_date_hes,
    pmin(first_event_date_hes, first_event_date_pcr)
  )
)]

# Restrict to first diagnosis from 2010 onward.
gdm_combined <- gdm_combined[min_diag_date >= config$gdm_start_date]
gdm_combined[, year_of_diagnosis := year(min_diag_date)]

cat("Unique women with GDM (first diagnosis 2010+):", uniqueN(gdm_combined$patid), "\n")

##############################
# 7. PLOT DIAGNOSIS YEAR DISTRIBUTION
##############################
p <- ggplot(gdm_combined, aes(x = year_of_diagnosis)) +
  geom_histogram(binwidth = 1, fill = "#69b3a2", color = "white", alpha = 0.8) +
  scale_x_continuous(
    breaks = seq(
      min(gdm_combined$year_of_diagnosis, na.rm = TRUE),
      max(gdm_combined$year_of_diagnosis, na.rm = TRUE),
      by = 2
    )
  ) +
  labs(
    title    = "Distribution of First GDM Diagnosis Years",
    subtitle = "Women whose first GDM diagnosis occurred from 2010 onward",
    x        = "Year of first GDM diagnosis",
    y        = "Number of women"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title   = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(config$output_dir, "01_gdm_diagnosis_year_distribution.png"), p, width = 9, height = 5)

##############################
# 8. DOWNLOAD SOURCE TABLES IN BATCHES
##############################
patid_list   <- unique(gdm_combined$patid)
patid_chunks <- split(patid_list, ceiling(seq_along(patid_list) / config$batch_size))

batch_query <- function(con, table_name, patid_chunks) {
  results <- lapply(patid_chunks, function(chunk) {
    chunk_sql <- paste0("(", paste(chunk, collapse = ","), ")")
    dbGetQuery(con, paste0("SELECT * FROM ", table_name, " WHERE patid IN ", chunk_sql))
  })
  rbindlist(results, use.names = TRUE, fill = TRUE)
}

cat("Downloading PCR observations...\n")
pcr_observation <- batch_query(con, "eRAP_22_002383_Predict_GDB_Data.Observation", patid_chunks)
pcr_observation$patid <- as.character(pcr_observation$patid)

cat("Downloading HES diagnoses...\n")
hes_diagnoses <- batch_query(con, "eRAP_22_002383_Predict_GDB_Data.hes_diagnosis_epi_DM", patid_chunks)
hes_diagnoses$patid <- as.character(hes_diagnoses$patid)

cat("Downloading HES maternity...\n")
hes_maternity <- batch_query(con, "eRAP_22_002383_Predict_GDB_Data.hes_maternity_DM", patid_chunks)
hes_maternity$patid <- as.character(hes_maternity$patid)

cat("Downloading HES procedures...\n")
hes_procedures <- batch_query(con, "eRAP_22_002383_Predict_GDB_Data.hes_procedures_epi_DM", patid_chunks)
hes_procedures$patid <- as.character(hes_procedures$patid)

rm(patid_chunks, patid_list)

##############################
# 9. DOWNLOAD SUPPLEMENTARY TABLES
##############################
cat("Downloading patient-level tables...\n")

pcr_patient <- dbGetQuery(con, "SELECT * FROM eRAP_22_002383_Predict_GDB_Data.Patient") %>%
  rename(patid = PatId) %>%
  mutate(patid = as.character(patid))

ons_patient <- dbGetQuery(con, "SELECT * FROM eRAP_22_002383_Predict_GDB_Data.Death_Patient_DM") %>%
  mutate(patid = as.character(patid))

imd_patient <- dbGetQuery(con, "SELECT * FROM eRAP_22_002383_Predict_GDB_Data.patient_2019_imd") %>%
  mutate(patid = as.character(patid))

hes_patient <- dbGetQuery(con, "SELECT * FROM eRAP_22_002383_Predict_GDB_Data.hes_patient_DM") %>%
  mutate(
    patid  = as.character(patid),
    pracid = as.character(pracid)
  )

pcr_practice <- dbGetQuery(con, "SELECT * FROM eRAP_22_002383_Predict_GDB_Data.Practice") %>%
  mutate(pracid = as.character(pracid))

##############################
# 10. CLEAN MATERNITY TABLE
##############################
hes_maternity <- hes_maternity %>%
  mutate(
    patid  = as.character(patid),
    neocare = na_if(na_if(as.numeric(neocare), 8), 9),
    birstat = na_if(as.numeric(birstat), 9),
    birweit = as.numeric(birweit),
    gestat  = na_if(as.numeric(gestat), 99),
    numpreg = case_when(
      numpreg == 99  ~ NA_real_,
      numpreg == "X" ~ NA_real_,
      TRUE           ~ suppressWarnings(as.numeric(numpreg))
    )
  )

##############################
# 11. SAVE OUTPUTS
##############################
saveRDS(gdm_combined,    file.path(config$processed_dir, "gdm_combined.rds"))
saveRDS(gdm_pcr,         file.path(config$processed_dir, "gdm_pcr.rds"))
saveRDS(gdm_hes,         file.path(config$processed_dir, "gdm_hes.rds"))
saveRDS(pcr_observation, file.path(config$processed_dir, "pcr_observation.rds"))
saveRDS(hes_diagnoses,   file.path(config$processed_dir, "hes_diagnoses.rds"))
saveRDS(hes_maternity,   file.path(config$processed_dir, "hes_maternity.rds"))
saveRDS(hes_procedures,  file.path(config$processed_dir, "hes_procedures.rds"))
saveRDS(pcr_patient,     file.path(config$processed_dir, "pcr_patient.rds"))
saveRDS(ons_patient,     file.path(config$processed_dir, "ons_patient.rds"))
saveRDS(imd_patient,     file.path(config$processed_dir, "imd_patient.rds"))
saveRDS(hes_patient,     file.path(config$processed_dir, "hes_patient.rds"))
saveRDS(pcr_practice,    file.path(config$processed_dir, "pcr_practice.rds"))
saveRDS(all_codes,       file.path(config$processed_dir, "all_codes.rds"))

writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "01_session_info.txt"))
cat("01_extract_gdm_cohort.R completed successfully.\n")
