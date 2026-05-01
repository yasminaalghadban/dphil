############################################################
# 01_prepare_sensitivity_cohort.R
# Purpose: Identify GDM ascertainment sources, merge with
#          the final cohort, flag HES-only women, build the
#          time-to-event dataset, and produce Table 1.
############################################################

##############################
# 1. PACKAGE SETUP
##############################
required_packages <- c(
  "ggplot2", "dplyr", "readxl", "lubridate", "odbc", "DBI",
  "data.table", "openxlsx", "flextable", "officer", "tibble", "tableone"
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
  raw_dir       = file.path("data", "raw"),
  processed_dir = file.path("data", "processed"),
  output_dir    = "outputs",
  codes_file    = "ALL_CODES.xlsx",
  codes_sheet   = "Sheet1",
  dataset_file  = "final_dataset.xlsx",
  exclusion_file = "exclusion_final.xlsx",
  gdm_start_date = as.Date("2010-01-01"),
  collection_end  = as.Date("2023-12-20"),
  followup_years  = 10
)

dir.create(config$processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$output_dir,    recursive = TRUE, showWarnings = FALSE)

##############################
# 3. SHARED PLOT THEME
##############################
set_journal_style <- function(mfrow = c(1, 1)) {
  par(
    mfrow    = mfrow,
    mar      = c(4.8, 6.0, 2.2, 1.2),
    las      = 1,
    bty      = "l",
    cex.axis = 1.1,
    cex.lab  = 1.2,
    cex.main = 1.4,
    family   = "serif",
    mgp      = c(3.5, 1, 0)
  )
}

##############################
# 4. DATABASE CONNECTION
#    Credentials are read from environment variables.
#    Set them in your .Renviron file — never commit to GitHub.
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
# 5. LOAD CODE LISTS
##############################
codes_path <- file.path(config$raw_dir, config$codes_file)
stopifnot(file.exists(codes_path))

all_codes <- read_excel(codes_path, sheet = config$codes_sheet) %>%
  mutate(Criteria = trimws(Criteria), Code_Type = trimws(Code_Type))

##############################
# 6. IDENTIFY GDM CASES BY SOURCE
##############################

# -- Primary care (PCR) codes --
medcode_gdm_codes <- all_codes %>%
  filter(Criteria == "GDM - pcr", Code_Type == "MedCodeId") %>%
  pull(Code)

gdm_pcr_raw <- dbGetQuery(
  con,
  paste0(
    "SELECT * FROM eRAP_22_002383_Predict_GDB_Data.Observation WHERE medcodeid IN (",
    paste(medcode_gdm_codes, collapse = ","), ")"
  )
) %>%
  mutate(
    patid     = as.character(patid),
    medcodeid = as.character(medcodeid),
    obsdate   = as.Date(obsdate, format = "%Y-%m-%d")
  ) %>%
  select(patid, medcodeid, obsdate) %>%
  mutate(selection = "pcr")

# -- Secondary care (HES) codes --
icd_gdm_codes <- all_codes %>%
  filter(Criteria == "GDM - HES", Code_Type == "ICD10") %>%
  pull(Code)

gdm_hes_raw <- dbGetQuery(
  con,
  paste0(
    "SELECT * FROM eRAP_22_002383_Predict_GDB_Data.hes_diagnosis_epi_DM WHERE ICD IN ('",
    paste(icd_gdm_codes, collapse = "','"), "')"
  )
) %>%
  mutate(
    patid    = as.character(patid),
    ICD      = as.character(ICD),
    epistart = as.Date(epistart, format = "%Y-%m-%d")
  ) %>%
  select(patid, ICD, epistart) %>%
  mutate(selection = "hes")

##############################
# 7. EARLIEST-EVER GDM DATE ACROSS SOURCES
##############################
setDT(gdm_hes_raw)
setDT(gdm_pcr_raw)

gdm_hes_first <- gdm_hes_raw[, .(first_event_date_hes = min(epistart, na.rm = TRUE)), by = patid]
gdm_pcr_first <- gdm_pcr_raw[, .(first_event_date_pcr = min(obsdate,  na.rm = TRUE)), by = patid]

gdm_combined_dt <- merge(gdm_hes_first, gdm_pcr_first, by = "patid", all = TRUE)

gdm_combined_dt[, min_diag_date := fifelse(
  is.na(first_event_date_hes), first_event_date_pcr,
  fifelse(
    is.na(first_event_date_pcr), first_event_date_hes,
    pmin(first_event_date_hes, first_event_date_pcr)
  )
)]

# Restrict to first-ever GDM from 2010 onward.
eligible_ids <- gdm_combined_dt[min_diag_date >= config$gdm_start_date]

##############################
# 8. REBUILD SOURCE COHORTS AFTER ELIGIBILITY FILTER
##############################

# PCR cohort: women with a primary care GDM code.
gdm_pcr <- eligible_ids[!is.na(first_event_date_pcr), .(
  patid,
  epistart = first_event_date_pcr
)] %>%
  as_tibble() %>%
  mutate(selection = "pcr")

# HES cohort: women with a secondary care GDM code.
gdm_hes <- eligible_ids[!is.na(first_event_date_hes), .(
  patid,
  epistart = first_event_date_hes
)] %>%
  as_tibble() %>%
  mutate(selection = "hes")

# HES-only: women identified only through secondary care.
hes_only <- gdm_hes %>%
  anti_join(gdm_pcr, by = "patid")

# Full combined source cohort (PCR takes priority for the selection label).
all_gdm <- bind_rows(gdm_pcr, hes_only)

##############################
# 9. LOAD FINAL DATASET AND EXCLUSIONS
##############################
dataset   <- read_excel(file.path(config$raw_dir, config$dataset_file))
exclusion <- read_excel(file.path(config$raw_dir, config$exclusion_file))

##############################
# 10. APPLY INCLUSION/EXCLUSION CRITERIA
##############################
final_dataset <- dataset %>%
  filter(maternal_age >= 18, maternal_age <= 45) %>%
  mutate(dod_delivery_diff = as.numeric(difftime(dod, delivery_date, units = "days"))) %>%
  filter(dod_delivery_diff > 42 | is.na(dod_delivery_diff)) %>%
  anti_join(exclusion, by = "patid") %>%
  left_join(all_gdm, by = "patid")

##############################
# 11. OVERLAP TABLE
##############################
final_ids     <- final_dataset %>% distinct(patid)
pcr_ids_final <- gdm_pcr  %>% distinct(patid) %>% semi_join(final_ids, by = "patid")
hes_ids_final <- gdm_hes  %>% distinct(patid) %>% semi_join(final_ids, by = "patid")

n_pcr  <- nrow(pcr_ids_final)
n_hes  <- nrow(hes_ids_final)
n_both <- inner_join(pcr_ids_final, hes_ids_final, by = "patid") %>% nrow()

overlap_tab <- tibble(
  label = c(
    "Women in the GDM secondary care cohort that are:",
    "  - In the GDM primary care cohort",
    "  - Not in the GDM primary care cohort",
    "Women in the GDM primary care cohort that are:",
    "  - In the GDM secondary care cohort",
    "  - Not in the GDM secondary care cohort"
  ),
  n = c(
    "",
    format(n_both, big.mark = ","),
    format(n_hes - n_both, big.mark = ","),
    "",
    format(n_both, big.mark = ","),
    format(n_pcr - n_both, big.mark = ",")
  ),
  pct = c(
    "",
    sprintf("%.1f", 100 * n_both / n_hes),
    sprintf("%.1f", 100 * (n_hes - n_both) / n_hes),
    "",
    sprintf("%.1f", 100 * n_both / n_pcr),
    sprintf("%.1f", 100 * (n_pcr - n_both) / n_pcr)
  ),
  is_header = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE)
)

ft_overlap <- flextable(overlap_tab %>% select(label, n, pct)) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(j = 1, align = "left", part = "all") %>%
  align(j = 2:3, align = "right", part = "all") %>%
  bold(i = which(overlap_tab$is_header), j = 1:3, bold = TRUE, part = "body") %>%
  width(j = 1, width = 5.6) %>%
  width(j = 2, width = 1.2) %>%
  width(j = 3, width = 0.9) %>%
  add_header_lines(values = "Overlap between GDM primary care and secondary care cohorts")

save_as_docx("Overlap table" = ft_overlap,
             path = file.path(config$output_dir, "01_overlap_table.docx"))

##############################
# 12. STANDARDISE VARIABLES
##############################
final_dataset <- final_dataset %>%
  mutate(
    region_name = case_when(
      region == 1 ~ "North East",       region == 2 ~ "North West",
      region == 3 ~ "Yorkshire and the Humber", region == 4 ~ "East Midlands",
      region == 5 ~ "West Midlands",    region == 6 ~ "East of England",
      region == 7 ~ "London",           region == 8 ~ "South East",
      region == 9 ~ "South West",       TRUE ~ NA_character_
    ),
    ses_quintile = case_when(
      ses %in% c(1, 2) ~ "1", ses %in% c(3, 4) ~ "2",
      ses %in% c(5, 6) ~ "3", ses %in% c(7, 8) ~ "4",
      ses %in% c(9, 10) ~ "5", TRUE ~ NA_character_
    ),
    alcohol_status = case_when(
      alcohol_status %in% c("current", "former", "never") ~ alcohol_status,
      TRUE ~ NA_character_
    ),
    smoking_status = case_when(
      smoking_status %in% c("current", "former", "never") ~ smoking_status,
      TRUE ~ NA_character_
    ),
    smm_now = case_when(
      smm_now == 1 ~ "yes", smm_now == 0 ~ "no", TRUE ~ NA_character_
    ),
    yodelivery = year(delivery_date)
  ) %>%
  mutate(across(
    c(primiparity, bmi_category_final, substance_use, socioeconomic_deprivation,
      hyp_disorders, diabetes_status_10yr, complications_10yr,
      diabetes_status, complications, gen_ethnicity),
    as.factor
  ))

##############################
# 13. BUILD TIME-TO-EVENT DATASET
##############################
final_dataset <- final_dataset %>%
  mutate(
    first_event_date_no_complication = if_else(
      is.na(first_event_date_no_complication) & !is.na(first_event_date_complication),
      first_event_date_complication,
      first_event_date_no_complication
    )
  )

final_dataset_no_complication <- final_dataset %>%
  mutate(
    delivery_date                    = as.Date(delivery_date),
    first_event_date_no_complication = as.Date(first_event_date_no_complication),
    dod                              = as.Date(dod),
    end_time                         = add_with_rollback(as.Date(delivery_date), years(config$followup_years)),
    collection_end                   = config$collection_end,
    outcome = case_when(
      diabetes_status_10yr == "no" ~ 0L,
      !is.na(first_event_date_no_complication) &
        (is.na(dod) | first_event_date_no_complication < dod) &
        first_event_date_no_complication <= end_time ~ 1L,
      TRUE ~ 0L
    ),
    event_time = pmin(first_event_date_no_complication, collection_end, dod, end_time, na.rm = TRUE),
    time       = as.numeric(difftime(event_time, delivery_date, units = "days")) / 365.25,
    missed     = as.integer(selection == "hes")
  )

cat("Median follow-up time:", median(final_dataset_no_complication$time, na.rm = TRUE), "years\n")

##############################
# 14. TABLE 1
##############################
final_dataset_no_complication$gen_ethnicity <-
  relevel(factor(final_dataset_no_complication$gen_ethnicity), ref = "White")

table1_df <- final_dataset_no_complication %>%
  mutate(
    ascertainment  = factor(ifelse(missed == 1, "HES only", "Primary care coded"),
                            levels = c("Primary care coded", "HES only")),
    smoking_status = factor(smoking_status, levels = c("never", "former", "current")),
    alcohol_status = factor(alcohol_status, levels = c("never", "former", "current")),
    gen_ethnicity  = factor(gen_ethnicity),
    region_name    = factor(region_name),
    primiparity    = factor(primiparity),
    previous_complications = factor(previous_complications),
    treatment_method = factor(treatment_method),
    ses_quintile = factor(ses_quintile, levels = 1:5,
                          labels = c("1 (least deprived)", "2", "3", "4", "5 (most deprived)")),
    delivery_period = factor(case_when(
      yodelivery >= 2010 & yodelivery <= 2012 ~ "2010-2012",
      yodelivery >= 2013 & yodelivery <= 2015 ~ "2013-2015",
      yodelivery >= 2016 & yodelivery <= 2018 ~ "2016-2018",
      yodelivery >= 2019 & yodelivery <= 2021 ~ "2019-2021",
      TRUE ~ NA_character_
    ), levels = c("2010-2012", "2013-2015", "2016-2018", "2019-2021"))
  ) %>%
  select(ascertainment, maternal_age, bmi_value, smoking_status, alcohol_status,
         gen_ethnicity, region_name, ses_quintile, previous_complications,
         delivery_period, treatment_method)

continuous_vars  <- c("maternal_age", "bmi_value")
categorical_vars <- c("smoking_status", "alcohol_status", "gen_ethnicity", "region_name",
                      "ses_quintile", "previous_complications", "delivery_period", "treatment_method")
all_vars <- c(continuous_vars, categorical_vars)

var_labels <- c(
  maternal_age = "Maternal age, years",
  bmi_value = "BMI, kg/m²",
  smoking_status = "Smoking status",
  alcohol_status = "Alcohol status",
  gen_ethnicity = "Ethnicity",
  region_name = "Region",
  ses_quintile = "IMD quintile",
  previous_complications = "Previous pregnancy complications",
  delivery_period = "Year of delivery",
  treatment_method = "Treatment method"
)

tab1_smd <- CreateTableOne(
  vars = all_vars, strata = "ascertainment", data = table1_df,
  factorVars = categorical_vars, includeNA = FALSE, test = FALSE
)
smd_obj <- ExtractSmd(tab1_smd)
smd_lookup <- data.frame(
  variable = rownames(smd_obj),
  smd      = as.numeric(smd_obj[, 1]),
  stringsAsFactors = FALSE
)

fmt_mean_sd <- function(x, d = 1) {
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(x), sd(x))
}
fmt_n_pct <- function(n, denom, d = 1) {
  if (is.na(denom) || denom == 0) return("")
  sprintf(paste0("%d (%.", d, "f%%)"), n, 100 * n / denom)
}
fmt_smd  <- function(x) if (is.na(x) || !length(x)) "" else sprintf("%.3f", x[1])
fmt_p    <- function(p)  if (is.na(p)) "" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
get_smd  <- function(v, lu) { o <- lu$smd[lu$variable == v]; if (!length(o)) NA_real_ else o[1] }

get_pvalue <- function(data, var, grp = "ascertainment", type = c("continuous", "categorical")) {
  type <- match.arg(type)
  dat  <- data %>% filter(!is.na(.data[[grp]]), !is.na(.data[[var]]))
  if (!nrow(dat) || dplyr::n_distinct(dat[[grp]]) < 2) return(NA_real_)
  if (type == "continuous") return(tryCatch(t.test(dat[[var]] ~ dat[[grp]])$p.value, error = function(e) NA_real_))
  tab <- table(dat[[var]], dat[[grp]])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  tryCatch(suppressWarnings(chisq.test(tab)$p.value), error = function(e) NA_real_)
}

g0 <- "Primary care coded"; g1 <- "HES only"
n_g0 <- sum(table1_df$ascertainment == g0)
n_g1 <- sum(table1_df$ascertainment == g1)

make_cont_rows <- function(data, var, label) {
  x0 <- data %>% filter(ascertainment == g0) %>% pull(.data[[var]])
  x1 <- data %>% filter(ascertainment == g1) %>% pull(.data[[var]])
  bind_rows(
    data.frame(Characteristic = label, Group0 = fmt_mean_sd(x0), Group1 = fmt_mean_sd(x1),
               SMD = fmt_smd(get_smd(var, smd_lookup)), P_value = fmt_p(get_pvalue(data, var, type = "continuous")),
               is_header = FALSE, stringsAsFactors = FALSE),
    data.frame(Characteristic = "  Missing", Group0 = fmt_n_pct(sum(is.na(x0)), length(x0)),
               Group1 = fmt_n_pct(sum(is.na(x1)), length(x1)), SMD = "", P_value = "",
               is_header = FALSE, stringsAsFactors = FALSE)
  )
}

make_cat_rows <- function(data, var, label) {
  x0   <- data %>% filter(ascertainment == g0) %>% pull(.data[[var]])
  x1   <- data %>% filter(ascertainment == g1) %>% pull(.data[[var]])
  levs <- if (is.factor(data[[var]])) levels(data[[var]]) else sort(unique(na.omit(data[[var]])))
  rows <- list(data.frame(Characteristic = label, Group0 = "", Group1 = "",
                          SMD = fmt_smd(get_smd(var, smd_lookup)),
                          P_value = fmt_p(get_pvalue(data, var, type = "categorical")),
                          is_header = TRUE, stringsAsFactors = FALSE))
  for (lv in levs) {
    rows[[length(rows) + 1]] <- data.frame(
      Characteristic = paste0("  ", lv),
      Group0 = fmt_n_pct(sum(x0 == lv, na.rm = TRUE), sum(!is.na(x0))),
      Group1 = fmt_n_pct(sum(x1 == lv, na.rm = TRUE), sum(!is.na(x1))),
      SMD = "", P_value = "", is_header = FALSE, stringsAsFactors = FALSE)
  }
  rows[[length(rows) + 1]] <- data.frame(
    Characteristic = "  Missing",
    Group0 = fmt_n_pct(sum(is.na(x0)), length(x0)),
    Group1 = fmt_n_pct(sum(is.na(x1)), length(x1)),
    SMD = "", P_value = "", is_header = FALSE, stringsAsFactors = FALSE)
  bind_rows(rows)
}

table1_rows <- c(
  lapply(continuous_vars,  function(v) make_cont_rows(table1_df, v, var_labels[[v]])),
  lapply(categorical_vars, function(v) make_cat_rows(table1_df, v, var_labels[[v]]))
)
table1_out <- bind_rows(table1_rows)
names(table1_out)[names(table1_out) == "Group0"] <- paste0("Primary care coded (N=", format(n_g0, big.mark = ","), ")")
names(table1_out)[names(table1_out) == "Group1"] <- paste0("HES only (N=", format(n_g1, big.mark = ","), ")")

ft_table1 <- flextable(table1_out %>% select(-is_header)) %>%
  theme_booktabs() %>% autofit() %>%
  bold(i = which(table1_out$is_header), bold = TRUE, part = "body") %>%
  align(j = 1, align = "left", part = "all") %>%
  align(j = 2:5, align = "center", part = "all") %>%
  width(j = 1, width = 3.4) %>% width(j = 2:5, width = 1.5) %>%
  add_header_lines("Table 1. Baseline characteristics by GDM ascertainment source") %>%
  add_footer_lines(c(
    "Values are mean (SD) for continuous variables and n (%) for categorical variables.",
    "SMD = standardised mean difference."
  ))

save_as_docx("Table 1" = ft_table1,
             path = file.path(config$output_dir, "01_table1_baseline_characteristics.docx"))

##############################
# 15. SAVE
##############################
saveRDS(final_dataset_no_complication, file.path(config$processed_dir, "sensitivity_cohort.rds"))
writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "01_session_info.txt"))
cat("01_prepare_sensitivity_cohort.R completed successfully.\n")
