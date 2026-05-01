############################################################
# 03_pooled_models_and_figures.R
# Purpose: Fit pooled models from multiply imputed data,
#          build grouped regression tables, generate
#          standardised adjusted survival curves, and
#          assess the proportional hazards assumption.
############################################################

##############################
# 1. PACKAGE SETUP
##############################
required_packages <- c(
  "dplyr", "mice", "broom", "survival", "readr", "stringr",
  "purrr", "tibble", "flextable", "officer", "ggplot2",
  "splines", "MASS"
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
  processed_dir   = file.path("data", "processed"),
  output_dir      = "outputs",
  n_imp_use       = NULL,
  interval_width  = 1,
  horizon         = 10,
  n_sim           = 500,
  ph_plots_to_show = 6
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

imp_path <- file.path(config$processed_dir, "sensitivity_mice_imp.rds")
stopifnot(file.exists(imp_path))

##############################
# 3. LOAD IMPUTATION OBJECT
##############################
imp <- readRDS(imp_path)
if (is.null(config$n_imp_use)) config$n_imp_use <- imp$m

##############################
# 4. HELPER FUNCTIONS
##############################
fmt_est <- function(x, digits = 2) ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
fmt_ci  <- function(low, high, digits = 2) {
  ifelse(is.na(low) | is.na(high), "",
         paste0(sprintf(paste0("%.", digits, "f"), low), " to ", sprintf(paste0("%.", digits, "f"), high)))
}
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "",
    p < 0.001  ~ "<0.001",
    TRUE       ~ sprintf("%.3f", p)
  )
}

clean_term_names <- function(term_vec) {
  term_vec %>%
    str_replace_all('relevel\\(factor\\(as\\.character\\(ses_quintile\\)\\),\\s*ref\\s*=\\s*"[^"]+"\\)', "ses_quintile") %>%
    str_replace_all('relevel\\(factor\\(ses_quintile\\),\\s*ref\\s*=\\s*"[^"]+"\\)', "ses_quintile") %>%
    str_replace_all('relevel\\(factor\\(as\\.character\\(delivery_period\\)\\),\\s*ref\\s*=\\s*"[^"]+"\\)', "delivery_period") %>%
    str_replace_all('relevel\\(factor\\(delivery_period\\),\\s*ref\\s*=\\s*"[^"]+"\\)', "delivery_period") %>%
    str_replace_all('relevel\\(factor\\(as\\.character\\(previous_complications\\)\\),\\s*ref\\s*=\\s*"[^"]+"\\)', "previous_complications") %>%
    str_replace_all('relevel\\(factor\\(previous_complications\\),\\s*ref\\s*=\\s*"[^"]+"\\)', "previous_complications") %>%
    str_replace_all('relevel\\(factor\\(as\\.character\\(treatment_method\\)\\),\\s*ref\\s*=\\s*"[^"]+"\\)', "treatment_method") %>%
    str_replace_all('relevel\\(factor\\(treatment_method\\),\\s*ref\\s*=\\s*"[^"]+"\\)', "treatment_method") %>%
    str_replace_all('relevel\\(factor\\(as\\.character\\(region_name\\)\\),\\s*ref\\s*=\\s*"[^"]+"\\)', "region_name") %>%
    str_replace_all('relevel\\(factor\\(region_name\\),\\s*ref\\s*=\\s*"[^"]+"\\)', "region_name")
}

clean_model_table <- function(df, REGION_REF, TREATMENT_REF, IMD_REF, DELIVERY_REF, PARITY_REF) {
  df %>%
    mutate(term_clean = clean_term_names(term)) %>%
    filter(term_clean != "(Intercept)") %>%
    mutate(
      variable_group = case_when(
        term_clean == "missed" ~ "Ascertainment source",
        term_clean == "maternal_age" ~ "Maternal age, years",
        term_clean == "bmi_value" ~ "BMI, kg/m²",
        str_detect(term_clean, "^smoking_status") ~ "Smoking status",
        str_detect(term_clean, "^gen_ethnicity") ~ "Ethnicity",
        str_detect(term_clean, "^ses_quintile") ~ "IMD quintile",
        str_detect(term_clean, "^previous_complications") ~ "Previous pregnancy complications",
        str_detect(term_clean, "^treatment_method") ~ "Treatment method",
        str_detect(term_clean, "^region_name") ~ "Region",
        str_detect(term_clean, "^delivery_period") ~ "Year of delivery",
        TRUE ~ term_clean
      ),
      comparison = case_when(
        term_clean == "missed" ~ "Secondary care only vs primary care coded",
        term_clean == "maternal_age" ~ "Per 1-year increase",
        term_clean == "bmi_value" ~ "Per 1-unit increase",
        term_clean == "smoking_statusformer" ~ "Former vs never",
        term_clean == "smoking_statuscurrent" ~ "Current vs never",
        term_clean == "gen_ethnicityAsian" ~ "Asian vs White",
        term_clean == "gen_ethnicityBlack" ~ "Black vs White",
        term_clean == "gen_ethnicityMixed" ~ "Mixed vs White",
        term_clean == "gen_ethnicityOther" ~ "Other vs White",
        term_clean == "ses_quintile2" ~ paste0("2 vs ", IMD_REF),
        term_clean == "ses_quintile3" ~ paste0("3 vs ", IMD_REF),
        term_clean == "ses_quintile4" ~ paste0("4 vs ", IMD_REF),
        term_clean == "ses_quintile5 (most deprived)" ~ paste0("5 (most deprived) vs ", IMD_REF),
        term_clean == "previous_complicationsNo, multiparous" ~ paste0("No, multiparous vs ", PARITY_REF),
        term_clean == "previous_complicationsYes, multiparous" ~ paste0("Yes, multiparous vs ", PARITY_REF),
        str_detect(term_clean, "^treatment_method") ~ paste0(str_remove(term_clean, "^treatment_method"), " vs ", TREATMENT_REF),
        str_detect(term_clean, "^region_name") ~ paste0(str_remove(term_clean, "^region_name"), " vs ", REGION_REF),
        term_clean == "delivery_period2013-2015" ~ paste0("2013-2015 vs ", DELIVERY_REF),
        term_clean == "delivery_period2016-2018" ~ paste0("2016-2018 vs ", DELIVERY_REF),
        term_clean == "delivery_period2019-2021" ~ paste0("2019-2021 vs ", DELIVERY_REF),
        TRUE ~ term_clean
      ),
      group_order = case_when(
        variable_group == "Ascertainment source" ~ 1,
        variable_group == "Maternal age, years" ~ 2,
        variable_group == "BMI, kg/m²" ~ 3,
        variable_group == "Smoking status" ~ 4,
        variable_group == "Ethnicity" ~ 5,
        variable_group == "IMD quintile" ~ 6,
        variable_group == "Previous pregnancy complications" ~ 7,
        variable_group == "Treatment method" ~ 8,
        variable_group == "Region" ~ 9,
        variable_group == "Year of delivery" ~ 10,
        TRUE ~ 99
      ),
      estimate_display = fmt_est(estimate, 2),
      `95% CI` = fmt_ci(`2.5 %`, `97.5 %`, 2),
      `P-value` = fmt_p(p.value)
    ) %>%
    arrange(group_order, variable_group, comparison) %>%
    select(variable_group, comparison, estimate_display, `95% CI`, `P-value`)
}

make_grouped_display_table <- function(clean_df, estimate_label = "Adjusted OR") {
  split_df <- split(clean_df, clean_df$variable_group)
  out <- map_dfr(split_df, function(x) {
    bind_rows(
      tibble(Characteristic = unique(x$variable_group), estimate = "", `95% CI` = "", `P-value` = "", is_header = TRUE),
      tibble(Characteristic = paste0("   ", x$comparison), estimate = x$estimate_display, `95% CI` = x$`95% CI`, `P-value` = x$`P-value`, is_header = FALSE)
    )
  })
  out %>% rename(!!estimate_label := estimate)
}

make_regression_flextable <- function(df, title, footnote_lines = NULL) {
  ft <- flextable(df %>% select(-is_header)) %>%
    theme_booktabs() %>%
    autofit() %>%
    width(j = 1, width = 5.0) %>%
    width(j = 2, width = 1.2) %>%
    width(j = 3, width = 1.8) %>%
    width(j = 4, width = 1.0) %>%
    align(j = 1, align = "left", part = "all") %>%
    align(j = 2:4, align = "center", part = "all") %>%
    add_header_lines(values = title) %>%
    fontsize(size = 10, part = "body") %>%
    fontsize(size = 10.5, part = "header")

  header_rows <- which(df$is_header)
  if (length(header_rows) > 0) {
    ft <- bold(ft, i = header_rows, j = 1:4, bold = TRUE, part = "body")
  }
  if (!is.null(footnote_lines)) {
    ft <- add_footer_lines(ft, values = footnote_lines) %>%
      fontsize(size = 9, part = "footer")
  }
  ft
}

make_person_period <- function(dat, width = 1, horizon = 10) {
  dat <- dat %>%
    mutate(time = pmin(time, horizon), event = as.integer(diabetes_status_10yr))
  out_list <- vector("list", nrow(dat))
  for (i in seq_len(nrow(dat))) {
    ti <- dat$time[i]; yi <- dat$event[i]
    starts <- seq(0, horizon - width, by = width)
    stops  <- starts + width
    keep <- starts < ti
    starts <- starts[keep]; stops <- stops[keep]
    if (!length(starts)) next
    event_interval <- as.integer(yi == 1 & ti <= stops & ti > starts)
    out_list[[i]] <- bind_cols(
      dat[rep(i, length(starts)), , drop = FALSE],
      tibble(tstart = starts, tstop = stops, event_interval = event_interval)
    )
  }
  bind_rows(out_list)
}

make_prediction_period <- function(dat, width = 1, horizon = 10) {
  tidyr::expand_grid(id = dat$id, tstart = seq(0, horizon - width, by = width)) %>%
    mutate(tstop = tstart + width) %>%
    left_join(
      dat %>% dplyr::select(id, maternal_age, bmi_value, smoking_status, gen_ethnicity,
                            ses_quintile, previous_complications, treatment_method),
      by = "id"
    )
}

predict_survival_curve <- function(fit, pred_base, exposure_value, n_sim = 500) {
  pred_dat <- pred_base %>% mutate(missed = exposure_value)
  X <- model.matrix(delete.response(terms(fit)), data = pred_dat)
  beta_hat <- coef(fit)
  V_hat    <- vcov(fit)

  p_hat <- plogis(as.vector(X %*% beta_hat))
  point_curve <- pred_dat %>%
    arrange(id, tstop) %>%
    mutate(hazard_prob = p_hat) %>%
    group_by(id) %>%
    mutate(survival = cumprod(1 - hazard_prob)) %>%
    ungroup() %>%
    group_by(tstop) %>%
    summarise(survival = mean(survival, na.rm = TRUE), .groups = "drop")

  beta_draws <- MASS::mvrnorm(n = n_sim, mu = beta_hat, Sigma = V_hat)
  sim_curves <- map_dfr(seq_len(n_sim), function(s) {
    p_sim <- plogis(as.vector(X %*% beta_draws[s, ]))
    pred_dat %>%
      arrange(id, tstop) %>%
      mutate(hazard_prob = p_sim) %>%
      group_by(id) %>%
      mutate(survival = cumprod(1 - hazard_prob)) %>%
      ungroup() %>%
      group_by(tstop) %>%
      summarise(survival = mean(survival, na.rm = TRUE), .groups = "drop") %>%
      mutate(.draw = s)
  })

  within_var <- sim_curves %>%
    group_by(tstop) %>%
    summarise(U = var(survival, na.rm = TRUE), .groups = "drop")

  point_curve %>%
    left_join(within_var, by = "tstop") %>%
    mutate(missed = exposure_value)
}

get_curve_from_imp <- function(i, imp_object, width = 1, horizon = 10, n_sim = 500) {
  dat <- complete(imp_object, i) %>%
    mutate(
      id = row_number(),
      missed = as.integer(missed),
      diabetes_status_10yr = as.integer(diabetes_status_10yr),
      smoking_status = factor(as.character(smoking_status)),
      gen_ethnicity = factor(as.character(gen_ethnicity)),
      ses_quintile = factor(as.character(ses_quintile)),
      previous_complications = factor(as.character(previous_complications)),
      treatment_method = factor(as.character(treatment_method))
    )

  pp_dat <- make_person_period(dat = dat, width = width, horizon = horizon)
  plr_fit <- glm(
    event_interval ~ missed * ns(tstop, df = 4) +
      maternal_age + bmi_value + smoking_status + gen_ethnicity +
      ses_quintile + previous_complications + treatment_method,
    family = binomial(),
    data = pp_dat
  )

  pred_base <- make_prediction_period(dat = dat, width = width, horizon = horizon)

  bind_rows(
    predict_survival_curve(plr_fit, pred_base, exposure_value = 0, n_sim = n_sim),
    predict_survival_curve(plr_fit, pred_base, exposure_value = 1, n_sim = n_sim)
  ) %>% mutate(.imp = i)
}

##############################
# 5. POOLED MISSINGNESS MODEL
##############################
fit_missed <- with(
  imp,
  {
    ses_quintile    <- factor(ses_quintile, ordered = FALSE)
    delivery_period <- factor(delivery_period, ordered = FALSE)
    glm(
      missed ~ maternal_age + bmi_value + smoking_status + gen_ethnicity +
        ses_quintile + relevel(factor(previous_complications), ref = "Primiparous") +
        relevel(factor(treatment_method), ref = "lifestyle") +
        relevel(factor(region_name), ref = "London") +
        delivery_period,
      family = binomial
    )
  }
)

pool_missed <- pool(fit_missed)
missed_results <- summary(pool_missed, conf.int = TRUE, exponentiate = TRUE)
write_csv(missed_results, file.path(config$output_dir, "03_missed_results.csv"))

##############################
# 6. POOLED COX MODEL
##############################
fit_cox <- with(
  imp,
  {
    ses_quintile <- factor(ses_quintile, ordered = FALSE)
    coxph(
      Surv(time, diabetes_status_10yr) ~
        missed + maternal_age + bmi_value + smoking_status + gen_ethnicity +
        relevel(factor(as.character(ses_quintile)), ref = "1 (least deprived)") +
        relevel(factor(previous_complications), ref = "Primiparous") +
        relevel(factor(treatment_method), ref = "lifestyle")
    )
  }
)

cox_results <- summary(pool(fit_cox), conf.int = TRUE, exponentiate = TRUE)
write_csv(cox_results, file.path(config$output_dir, "03_cox_results.csv"))

##############################
# 7. PUBLICATION-STYLE TABLES
##############################
REGION_REF <- "London"
TREATMENT_REF <- "lifestyle"
IMD_REF <- "1 (least deprived)"
DELIVERY_REF <- "2010-2012"
PARITY_REF <- "Primiparous"

missed_clean <- clean_model_table(missed_results, REGION_REF, TREATMENT_REF, IMD_REF, DELIVERY_REF, PARITY_REF)
diabetes_clean <- clean_model_table(cox_results, REGION_REF, TREATMENT_REF, IMD_REF, DELIVERY_REF, PARITY_REF)

missed_table <- make_grouped_display_table(missed_clean, estimate_label = "Adjusted OR")
diabetes_table <- make_grouped_display_table(diabetes_clean, estimate_label = "Adjusted HR")

missed_footnotes <- c(
  "OR = odds ratio; CI = confidence interval.",
  "Estimates are pooled across multiply imputed datasets using Rubin's rules.",
  "Outcome: absence of a primary care GDM code.",
  paste0("Reference categories: never smoker, White ethnicity, ", IMD_REF, ", ", PARITY_REF, ", ", TREATMENT_REF, ", ", REGION_REF, ", and ", DELIVERY_REF, ".")
)

diabetes_footnotes <- c(
  "HR = hazard ratio; CI = confidence interval.",
  "Estimates are pooled across multiply imputed datasets using Rubin's rules.",
  "Outcome: incident diabetes during follow-up.",
  "For ascertainment source, the comparison is secondary care only versus primary care coded.",
  "Participants were followed until diabetes diagnosis, death, end of data collection, or 10 years, whichever occurred first.",
  paste0("Reference categories: never smoker, White ethnicity, ", IMD_REF, ", ", PARITY_REF, ", ", TREATMENT_REF, ", ", REGION_REF, ", and ", DELIVERY_REF, ".")
)

ft_missed <- make_regression_flextable(
  df = missed_table,
  title = "Table X. Adjusted odds ratios for characteristics associated with absence of a primary care GDM code",
  footnote_lines = missed_footnotes
)

ft_diabetes <- make_regression_flextable(
  df = diabetes_table,
  title = "Table Y. Adjusted hazard ratios for the association between GDM ascertainment source and incident diabetes during follow-up",
  footnote_lines = diabetes_footnotes
)

save_as_docx(
  "Table X" = ft_missed,
  "Table Y" = ft_diabetes,
  path = file.path(config$output_dir, "03_regression_tables_grouped.docx")
)

##############################
# 8. DESCRIPTIVE SUMMARY BY ASCERTAINMENT GROUP
##############################
dat1 <- complete(imp, 1) %>%
  mutate(
    ascertainment_group = ifelse(missed == 1, "Secondary care only", "Primary care coded")
  )

summary_by_group <- dat1 %>%
  group_by(ascertainment_group) %>%
  summarise(
    N = n(),
    outcome_n = sum(diabetes_status_10yr == 1, na.rm = TRUE),
    outcome_pct = 100 * mean(diabetes_status_10yr == 1, na.rm = TRUE),
    median_followup = median(time, na.rm = TRUE),
    q1_followup = quantile(time, 0.25, na.rm = TRUE),
    q3_followup = quantile(time, 0.75, na.rm = TRUE),
    mean_followup = mean(time, na.rm = TRUE),
    sd_followup = sd(time, na.rm = TRUE),
    person_years = sum(time, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_by_group, file.path(config$output_dir, "03_summary_by_ascertainment_group.csv"))

##############################
# 9. ADJUSTED SURVIVAL CURVE EXAMPLE
##############################
fit_adj <- coxph(
  Surv(time, diabetes_status_10yr) ~ missed + maternal_age + bmi_value +
    smoking_status + gen_ethnicity + ses_quintile + previous_complications +
    treatment_method,
  data = dat1
)

new_primary <- data.frame(
  missed = 0,
  maternal_age = mean(dat1$maternal_age, na.rm = TRUE),
  bmi_value = mean(dat1$bmi_value, na.rm = TRUE),
  smoking_status = "current",
  gen_ethnicity = "Asian",
  ses_quintile = "2",
  previous_complications = "Yes, multiparous",
  treatment_method = "metformin"
)
new_secondary <- new_primary
new_secondary$missed <- 1

sf <- survfit(fit_adj, newdata = rbind(new_primary, new_secondary))
png(file.path(config$output_dir, "03_adjusted_survival_example.png"), width = 1400, height = 1000, res = 150)
plot(
  sf, col = c("black", "grey40"), lty = c(1, 2), lwd = 2,
  xlab = "Follow-up time (years)",
  ylab = "Adjusted diabetes-free survival probability",
  ylim = c(0, 1.00), mark.time = FALSE
)
legend("bottomleft", legend = c("Primary care coded", "Secondary care only"),
       col = c("black", "grey40"), lty = c(1, 2), lwd = 2, bty = "n")
dev.off()

##############################
# 10. STANDARDISED ADJUSTED SURVIVAL CURVES WITH 95% CI
##############################
all_curves <- map_dfr(
  seq_len(config$n_imp_use),
  ~ get_curve_from_imp(
    i = .x,
    imp_object = imp,
    width = config$interval_width,
    horizon = config$horizon,
    n_sim = config$n_sim
  )
)

pooled_curves <- all_curves %>%
  group_by(tstop, missed) %>%
  summarise(
    m = n(),
    Qbar = mean(survival, na.rm = TRUE),
    Ubar = mean(U, na.rm = TRUE),
    B = var(survival, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    B = ifelse(is.na(B), 0, B),
    Tvar = Ubar + (1 + 1 / m) * B,
    se = sqrt(Tvar),
    df = ifelse(B == 0, Inf, (m - 1) * (1 + Ubar / ((1 + 1 / m) * B))^2),
    lower = pmax(0, Qbar - qt(0.975, df) * se),
    upper = pmin(1, Qbar + qt(0.975, df) * se),
    ascertainment_group = factor(missed, levels = c(0, 1), labels = c("Primary care coded", "Secondary care only"))
  )

risk_10y <- pooled_curves %>%
  filter(tstop == config$horizon) %>%
  mutate(
    risk_10y = 1 - Qbar,
    risk_10y_lower = 1 - upper,
    risk_10y_upper = 1 - lower
  ) %>%
  select(
    ascertainment_group,
    survival_10y = Qbar,
    survival_10y_lower = lower,
    survival_10y_upper = upper,
    risk_10y,
    risk_10y_lower,
    risk_10y_upper
  )

write_csv(risk_10y, file.path(config$output_dir, "03_standardized_10y_risk.csv"))

surv_plot <- ggplot(pooled_curves, aes(x = tstop, y = Qbar)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = ascertainment_group), alpha = 0.18, colour = NA) +
  geom_line(aes(linetype = ascertainment_group), linewidth = 1.1, colour = "black") +
  scale_linetype_manual(values = c("solid", "dotted")) +
  scale_fill_manual(values = c("grey70", "grey40")) +
  coord_cartesian(ylim = c(0.90, 1.00)) +
  scale_y_continuous(breaks = seq(0.90, 1.00, by = 0.02)) +
  scale_x_continuous(breaks = seq(1, config$horizon, by = 1)) +
  labs(
    x = "Follow-up time (years)",
    y = "Standardized diabetes-free survival probability",
    linetype = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.23, 0.18))

ggsave(file.path(config$output_dir, "03_standardized_survival_curves.png"), surv_plot, width = 9, height = 6)

##############################
# 11. PROPORTIONAL HAZARDS CHECKS
##############################
ph_list <- lapply(fit_cox$analyses, cox.zph, transform = "km")

ph_results <- map_dfr(seq_along(ph_list), function(i) {
  as.data.frame(ph_list[[i]]$table) %>%
    rownames_to_column("term") %>%
    mutate(.imp = i)
})

ph_summary <- ph_results %>%
  group_by(term) %>%
  summarise(
    n_imp = n(),
    median_p = median(p, na.rm = TRUE),
    min_p = min(p, na.rm = TRUE),
    max_p = max(p, na.rm = TRUE),
    n_p_lt_0.05 = sum(p < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(ph_summary, file.path(config$output_dir, "03_ph_summary.csv"))

png(file.path(config$output_dir, "03_ph_plots_missed.png"), width = 1600, height = 1000, res = 150)
par(mfrow = c(2, 3))
for (i in 1:min(config$ph_plots_to_show, length(ph_list))) {
  plot(ph_list[[i]], var = "missed", main = paste("Imputation", i))
}
par(mfrow = c(1, 1))
dev.off()

##############################
# 12. SAVE SESSION INFO
##############################
writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "03_session_info.txt"))
cat("03_pooled_models_and_figures.R completed successfully.\n")
