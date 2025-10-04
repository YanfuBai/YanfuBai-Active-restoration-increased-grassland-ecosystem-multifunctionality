# ============================================================
# Restoration Recovery Analysis Batch Scenarios (v9)
# Modified model (all scenarios):
#   recovery_index ~ years_c * taxon + (1 | plot_id) + (1 | metric_name)
# ============================================================

# Input and output paths
INPUT_WIDE <- "restoration_metrics_clean_wide_v2.csv"
BASE_OUT   <- "analysis_outputs"

# Suppress package startup messages
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(ggplot2)
  library(performance)
  library(ggrepel)
})

# Ensure output directory exists
if (!dir.exists(BASE_OUT)) dir.create(BASE_OUT, recursive = TRUE)

# Helper function to write CSV files
write_csv2 <- function(x, path) { 
  readr::write_csv(x, path, na = ""); 
  message("Wrote: ", path) 
}

# Helper function to ensure a directory exists
ensure_dir <- function(path) { 
  if (!dir.exists(path)) dir.create(path, recursive = TRUE) 
}

# ----------------------------
# 1) Load master wide file and build Recovery Index (RI)
# ----------------------------
stopifnot(file.exists(INPUT_WIDE))

dat_wide <- readr::read_csv(INPUT_WIDE, show_col_types = FALSE)

id_cols   <- c("plot_id", "quadrat_rep", "treatment_group", "years")
stopifnot(all(id_cols %in% names(dat_wide)))

metric_cols <- setdiff(names(dat_wide), id_cols)
div_cols    <- metric_cols[grepl("_hill_q", metric_cols)]
func_cols   <- setdiff(metric_cols, div_cols)

# Save original wide file to output directory
write_csv2(dat_wide, file.path(BASE_OUT, "restoration_metrics_original_wide.csv"))

# Calculate CK (control) mean for each metric
ck_means <- dat_wide %>%
  filter(treatment_group == "CK") %>%
  summarise(across(all_of(c(func_cols, div_cols)), ~ mean(.x, na.rm = TRUE)))

# Check for problematic years
bad_years <- unique(dat_wide$years[is.na(readr::parse_number(dat_wide$years))])
if (length(bad_years) > 0) {
  message("Warning: Non-numeric years found:")
  print(bad_years)
}

# ----------------------------
# Helper: Build Recovery Index (RI) with Taxon Information
# ----------------------------
build_ri_with_taxon <- function(df, diversity_cols, func_cols) {
  # Diversity data with taxon extraction
  div_long <- df %>%
    select(all_of(id_cols), all_of(diversity_cols)) %>%
    pivot_longer(cols = all_of(diversity_cols), names_to = "metric_name", values_to = "value") %>%
    mutate(
      taxon = case_when(
        str_detect(metric_name, "plant_") ~ "plants",
        str_detect(metric_name, "fungi_") ~ "fungi", 
        str_detect(metric_name, "bacteria_") ~ "bacteria"
      ),
      metric_type = "diversity"
    )
  
  # Function data
  func_long <- df %>%
    select(all_of(id_cols), all_of(func_cols)) %>%
    pivot_longer(cols = all_of(func_cols), names_to = "metric_name", values_to = "value") %>%
    mutate(
      taxon = "other_functions",
      metric_type = "function"
    )
  
  # Combine and add recovery index
  combined <- bind_rows(div_long, func_long)
  
  # Add reference values
  ck_vec <- ck_means %>% 
    select(all_of(c(diversity_cols, func_cols))) %>%
    pivot_longer(everything(), names_to = "metric_name", values_to = "ck_mean")
  
  # Calculate recovery index for individual metrics
  individual_metrics <- combined %>% 
    left_join(ck_vec, by = "metric_name") %>%
    mutate(
      recovery_index = 100 * value / ck_mean,
      years = readr::parse_number(years)
    ) %>%
    filter(!is.na(recovery_index), !is.na(years))
  
  # Calculate EMF as mean recovery index of function metrics  
  emf_data <- individual_metrics %>%
    filter(metric_type == "function") %>%
    group_by(plot_id, quadrat_rep, treatment_group, years) %>%
    summarise(recovery_index = mean(recovery_index, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      metric_name = "emf",
      taxon = "function", 
      metric_type = "function"
    )
  
  # Combine everything together
  bind_rows(individual_metrics, emf_data)
}

# ----------------------------
# Define Stage Labels
# ----------------------------
stage_map <- c("DG"="DG", "RG1"="Early", "RG2"="Middle", "RG3"="Late")

# ----------------------------
# Helper: Pretty Term Labels for Model Output
# ----------------------------
pretty_term <- function(term) {
  if (term == "(Intercept)") return("Intercept (function at mean age)")
  if (term == "taxonplants") return("Plants (vs function at mean age)")
  if (term == "taxonfungi") return("Fungi (vs function at mean age)")
  if (term == "taxonbacteria") return("Bacteria (vs function at mean age)")
  if (term == "years_c") return("Years since restoration (function)")
  if (term == "years_c:taxonplants") return("Interaction: Years × Plants")
  if (term == "years_c:taxonfungi") return("Interaction: Years × Fungi")
  if (term == "years_c:taxonbacteria") return("Interaction: Years × Bacteria")
  return(term)
}

# Accumulator for combined publication table
combined_pub <- NULL
combined_plot_data <- NULL

# ----------------------------
# Helper: Run One Scenario
# ----------------------------
run_scenario <- function(div_q = 0) {
  scen <- sprintf("diversity_q%d_vs_emf", div_q)
  out_dir <- file.path(BASE_OUT, scen)
  ensure_dir(out_dir)
  
  # Define diversity metrics for this q value
  div_set <- c(sprintf("bacteria_hill_q%d", div_q),
               sprintf("fungi_hill_q%d", div_q),
               sprintf("plant_hill_q%d", div_q))
  div_set <- div_set[div_set %in% div_cols]
  
  # Build the analysis dataset with taxon information
  ri_sub <- build_ri_with_taxon(dat_wide, div_set, func_cols) %>%
    filter(treatment_group != "CK") %>%
    mutate(
      years_c = years - mean(years, na.rm = TRUE),
      plot_id = factor(plot_id),
      taxon = factor(taxon, levels = c("function", "plants", "fungi", "bacteria", "other_functions")),
      stage = recode(treatment_group, !!!stage_map),
      stage = factor(stage, levels = c("DG", "Early", "Middle", "Late"))
    ) %>%
    drop_na(years, recovery_index, plot_id, taxon)
  
  # Save the analysis dataset for the current scenario
  write_csv2(ri_sub, file.path(out_dir, "analysis_dataset.csv"))
  
  # For modeling, filter out "other_functions" and keep only EMF and diversity metrics
  model_data <- ri_sub %>% filter(taxon != "other_functions")
  
  # Fit the linear mixed model (LMM) with taxon as fixed effect
  m_lmm <- lmer(recovery_index ~ years_c * taxon + 
                  (1 | plot_id) + (1 | metric_name),
                data = model_data, REML = TRUE)
  
  # Save the fixed effects, variance components, and model summary
  coefs_lmm <- broom.mixed::tidy(m_lmm, conf.int = FALSE, effects = "fixed")
  
  # Manually calculate confidence intervals
  coefs_lmm <- coefs_lmm %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error
    )
  
  write_csv2(coefs_lmm, file.path(out_dir, "lmm_fixed_effects.csv"))
  vc <- as.data.frame(VarCorr(m_lmm))
  write_csv2(vc, file.path(out_dir, "lmm_variance_components.csv"))
  capture.output(summary(m_lmm), file = file.path(out_dir, "lmm_summary_printout.txt"))
  
  # Create publication-style table with model results
  r2 <- tryCatch(performance::r2_nakagawa(m_lmm), error = function(e) NULL)
  r2_marg <- if (!is.null(r2)) unname(r2$R2_marginal) else NA_real_
  r2_cond <- if (!is.null(r2)) unname(r2$R2_conditional) else NA_real_
  aic_v <- AIC(m_lmm); bic_v <- BIC(m_lmm); ll_v <- as.numeric(logLik(m_lmm)); n_v <- nobs(m_lmm)
  
  # Create publication table
  pub <- coefs_lmm %>%
    mutate(
      term = vapply(term, pretty_term, character(1))
    ) %>%
    transmute(
      term_raw = term,
      term,
      estimate, std_error = std.error, conf_low = conf.low, conf_high = conf.high, p_value = p.value
    ) %>%
    mutate(
      scenario = scen,
      r2_marginal = r2_marg,
      r2_conditional = r2_cond,
      AIC = aic_v, BIC = bic_v, logLik = ll_v, N = n_v
    ) %>%
    select(scenario, term, estimate, std_error, conf_low, conf_high, p_value,
           r2_marginal, r2_conditional, AIC, BIC, logLik, N)
  
  # Write the publication table to a CSV file
  write_csv2(pub, file.path(out_dir, "publication_table_fixed_effects.csv"))
  
  # Accumulate the results for later use
  if (!exists("combined_pub", envir = .GlobalEnv) || is.null(get("combined_pub", envir = .GlobalEnv))) {
    assign("combined_pub", pub, envir = .GlobalEnv)
  } else {
    old_pub <- get("combined_pub", envir = .GlobalEnv)
    assign("combined_pub", dplyr::bind_rows(old_pub, pub), envir = .GlobalEnv)
  }
  
  # Create Grouped-Stage Figure
  plot_means <- ri_sub %>%
    group_by(metric_name, taxon, plot_id, stage) %>%
    summarise(plot_mean = mean(recovery_index, na.rm = TRUE), .groups = "drop")
  
  summ <- plot_means %>%
    group_by(metric_name, taxon, stage) %>%
    summarise(mean = mean(plot_mean, na.rm = TRUE),
              std  = sd(plot_mean,  na.rm = TRUE),
              n    = dplyr::n(),
              .groups = "drop") %>%
    mutate(se = std / sqrt(pmax(n, 1)),
           conf_low = mean - 1.96 * se,
           conf_high = mean + 1.96 * se)
  
  # Add q value to the summary data
  summ <- summ %>% mutate(q_value = div_q)
  
  # Save the summary data for plotting
  write_csv2(summ, file.path(out_dir, "grouped_stage_summary.csv"))
  
  # Accumulate the plot data for later combination
  if (!exists("combined_plot_data", envir = .GlobalEnv) || is.null(get("combined_plot_data", envir = .GlobalEnv))) {
    assign("combined_plot_data", summ, envir = .GlobalEnv)
  } else {
    old_plot_data <- get("combined_plot_data", envir = .GlobalEnv)
    assign("combined_plot_data", dplyr::bind_rows(old_plot_data, summ), envir = .GlobalEnv)
  }
  
  # Prepare the plot and labels for the grouped stage figure
  order <- levels(ri_sub$stage)
  
  # Create base plot with error bars and lines for the metrics
  p <- ggplot() +
    geom_hline(yintercept = 100, color = "black", linetype = "solid", size = 0.4) +
    geom_errorbar(
      data = summ %>% filter(taxon == "other_functions"),
      aes(x = stage, ymin = mean - se, ymax = mean + se),
      width = 0.15, color = "#9E9E9E", alpha = 0.9, linewidth = 0.4
    ) +
    geom_line(
      data = summ %>% filter(taxon == "other_functions"),
      aes(x = stage, y = mean, group = metric_name),
      color = "#9E9E9E", linewidth = 0.6, linetype = "dashed"
    ) +
    geom_point(
      data = summ %>% filter(taxon == "other_functions"),
      aes(x = stage, y = mean, group = metric_name),
      color = "#9E9E9E", size = 1.5
    ) +
    geom_errorbar(
      data = summ %>% filter(taxon != "other_functions"),
      aes(x = stage, ymin = mean - se, ymax = mean + se, color = taxon),
      width = 0.15, alpha = 0.9, linewidth = 0.4, position = position_dodge(width = 0.2)
    ) +
    geom_line(
      data = summ %>% filter(taxon != "other_functions"),
      aes(x = stage, y = mean, group = interaction(taxon, metric_name), color = taxon),
      linewidth = 0.6, position = position_dodge(width = 0.2)
    ) +
    geom_point(
      data = summ %>% filter(taxon != "other_functions"),
      aes(x = stage, y = mean, group = interaction(taxon, metric_name), color = taxon),
      size = 1.5, position = position_dodge(width = 0.2)
    ) +
    scale_color_manual(
      name = NULL,
      values = c("function" = "#E41A1C", "plants" = "#377EB8", 
                 "fungi" = "#4DAF4A", "bacteria" = "#984EA3"),
      labels = c("function" = "EMF", "plants" = "Plants", 
                 "fungi" = "Fungi", "bacteria" = "Bacteria")
    )
  
  # Add labels for the last stage
  last_stage <- tail(order, 1)
  
  # Add labels for EMF and diversity metrics
  main_lab_data <- summ %>%
    filter(stage == last_stage, taxon != "other_functions") %>%
    mutate(label = case_when(
      taxon == "function" ~ "EMF",
      taxon == "plants" ~ "Plants",
      taxon == "fungi" ~ "Fungi",
      taxon == "bacteria" ~ "Bacteria"
    ))
  
  if (nrow(main_lab_data) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = main_lab_data,
        aes(x = stage, y = mean, label = label, color = taxon),
        nudge_x = 0.25, direction = "y", hjust = 0,
        segment.size = 0.2, size = 3, show.legend = FALSE, max.overlaps = Inf
      )
  }
  
  # Add labels for a few representative function metrics
  func_lab_data <- summ %>%
    filter(stage == last_stage, taxon == "other_functions") %>%
    group_by(taxon) %>%
    slice(1:3) %>%
    mutate(label = metric_name)
  
  if (nrow(func_lab_data) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = func_lab_data,
        aes(x = stage, y = mean, label = label),
        nudge_x = 0.25, direction = "y", hjust = 0,
        segment.size = 0.2, size = 2.5, color = "#9E9E9E",
        show.legend = FALSE, max.overlaps = Inf
      )
  }
  
  # Modify axis and layout
  p <- p +
    scale_x_discrete(expand = expansion(mult = c(0.02, 0.15))) +
    labs(x = "Restoration stage", y = "Recovery Index (% of CK mean)",
         title = sprintf("Recovery by Stage: Diversity (q=%d) vs Functions", div_q),
         caption = "Gray dashed lines: individual function metrics (others).") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.margin = margin(t = 5, r = 25, b = 5, l = 5)) +
    coord_cartesian(clip = "off")
  
  # Save the figure
  ggsave(filename = file.path(out_dir, "figure_grouped_stage.png"),
         plot = p, width = 8.6, height = 5.2, dpi = 300)
  
  message("Finished scenario: ", scen)
  
  return(summ)
}

# ----------------------------
# Run Multiple Scenarios
# ----------------------------

all_plot_data <- list()

# Running scenarios for different q values
q0_data <- run_scenario(div_q = 0)
q1_data <- run_scenario(div_q = 1)
q2_data <- run_scenario(div_q = 2)

# Combine all plot data
all_plot_data <- bind_rows(q0_data, q1_data, q2_data)

# Save the combined plot data to a CSV file
write_csv2(all_plot_data, file.path(BASE_OUT, "combined_plot_data.csv"))

# Save the combined publication table
if (exists("combined_pub", envir = .GlobalEnv)) {
  write_csv2(get("combined_pub", envir = .GlobalEnv),
             file.path(BASE_OUT, "combined_publication_table_fixed_effects.csv"))
}

message("\nAll batch scenarios complete. Publication tables and plot data saved.\n")
