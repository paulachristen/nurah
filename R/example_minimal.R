#' ============================================================================
#' NURAH Simulation Feasibility Study - Minimal Example
#' ============================================================================
#'
#' This script demonstrates:
#' 1. Defining a DAG with parameters
#' 2. Simulating 10 governorates x 12 months with negative binomial mortality
#' 3. Applying MAR and MNAR missingness regimes
#' 4. Fitting models with IPW missingness handling
#' 5. Computing recovery metrics vs truth
#'
#' Dependencies: brms, dplyr, magrittr (and optionally mice for imputation)
#' ============================================================================

# Load required packages
library(brms)
library(dplyr)
library(magrittr)

# Source the nurah functions (in package context, these would be loaded via library(nurah))
# source("dag.R")
# source("simulate.R")
# source("fit.R")
# source("study.R")
# source("structures.R")
# source("utils.R")

# ============================================================================
# 1. DEFINE DAG AND PARAMETERS
# ============================================================================

# Use the Checchi et al. (2017) DAG with dummy parameters
# This creates a complex causal structure modeling indirect mortality pathways
dag <- checchi_2017_dag(parameters = dummy_checchi_2017_parameters())

# View the DAG structure
cat("DAG nodes:\n")
print(unique(c(dag$parameters$from, dag$parameters$to)))

cat("\nMortality node parents:\n")
mort_parents <- dag$parameters$from[dag$parameters$to == "Population mortality"]
print(mort_parents)

# ============================================================================
# 2. SIMULATE 10 GOVERNORATES x 12 MONTHS
# ============================================================================

set.seed(12345)

# Simulate complete data (no missingness yet)
sim_data <- simulate_crisis_data(
  dag = dag,
  start_date = "2023-10-01",
  n_periods = 12,                          # 12 months
  resolution = "month",                     # Monthly aggregation
  spatial_structure = paste0("gov", 1:10), # 10 governorates
  initial_population = 100000,             # 100k population per governorate
  noise_level = 1,                         # Full stochastic noise
  mortality_node = "Population mortality", # Outcome node
  mortality_phi = 1,                       # Overdispersion
  mortality_intercept = -8                 # Baseline log-rate (~0.03 per 1000/day)
)

cat("\nSimulated data dimensions:", nrow(sim_data), "rows x", ncol(sim_data), "columns\n")
cat("Columns:", paste(head(names(sim_data), 10), collapse = ", "), "...\n")

# Summary of true mortality
cat("\nTrue mortality summary:\n")
print(summary(sim_data[["Population mortality"]]))
cat("Total true deaths:", sum(sim_data[["Population mortality"]], na.rm = TRUE), "\n")

# ============================================================================
# 3. APPLY MISSINGNESS REGIMES
# ============================================================================

# Define predictors that will have missingness
predictors_with_missingness <- c(
  "Nutritional status",
  "Burden of endemic infectious diseases",
  "Burden of NCDs"
)

# Only use predictors that exist in the simulated data
available_predictors <- intersect(predictors_with_missingness, names(sim_data))
cat("\nApplying missingness to:", paste(available_predictors, collapse = ", "), "\n")

# 3a. MAR Missingness (depends on observed conflict/displacement)
cat("\n--- MAR Regime ---\n")
mar_result <- simulate_missingness(
  data = sim_data,
  predictors = available_predictors,
  regime = "MAR",
  pi_base = 0.7,                           # 70% baseline observation probability
  mar_drivers = c("Exposure to armed attacks or mechanical force of nature"),
  mar_weights = c("Exposure to armed attacks or mechanical force of nature" = -0.5)
)
sim_data_mar <- mar_result$data
cat("MAR missingness applied. Missing rate by variable:\n")
for (p in available_predictors) {
  miss_rate <- mean(is.na(sim_data_mar[[p]]))
  cat(sprintf("  %s: %.1f%% missing\n", p, 100 * miss_rate))
}

# 3b. MNAR Missingness (depends on the unobserved value itself)
cat("\n--- MNAR Regime ---\n")
mnar_result <- simulate_missingness(
  data = sim_data,
  predictors = available_predictors,
  regime = "MNAR",
  pi_base = 0.7,
  mar_drivers = c("Exposure to armed attacks or mechanical force of nature"),
  mar_weights = c("Exposure to armed attacks or mechanical force of nature" = -0.3),
  mnar_omega = -0.6                        # Higher values -> less likely to observe
)
sim_data_mnar <- mnar_result$data
cat("MNAR missingness applied. Missing rate by variable:\n")
for (p in available_predictors) {
  miss_rate <- mean(is.na(sim_data_mnar[[p]]))
  cat(sprintf("  %s: %.1f%% missing\n", p, 100 * miss_rate))
}

# ============================================================================
# 4. OPTIONAL: APPLY GAZA-STYLE MASK
# ============================================================================

# Example: Create a mask that mimics Gaza reporting collapse
# (observations decrease over time, especially in northern governorates)
cat("\n--- Creating Gaza-style mask example ---\n")

gaza_mask <- expand.grid(
  gov = paste0("gov", 1:10),
  month = 1:12,
  var = available_predictors[1],
  stringsAsFactors = FALSE
)
# Northern governorates (gov1-3) lose reporting after month 4
# Southern governorates (gov4-10) lose reporting after month 8
gaza_mask$observed <- ifelse(
  gaza_mask$gov %in% c("gov1", "gov2", "gov3") & gaza_mask$month > 4, 0,
  ifelse(gaza_mask$month > 8, 0, 1)
)

cat("Mask preview (first 15 rows):\n")
print(head(gaza_mask, 15))

# Apply mask (force mode)
sim_data_masked <- apply_mask(sim_data, gaza_mask, mode = "force")
cat("\nMask applied. Missing rate in masked data:\n")
cat(sprintf("  %s: %.1f%% missing\n",
            available_predictors[1],
            100 * mean(is.na(sim_data_masked[[available_predictors[1]]]))))

# ============================================================================
# 5. FIT MODELS
# ============================================================================

cat("\n========================================\n")
cat("FITTING MODELS\n")
cat("========================================\n")

# True values for comparison
true_total <- sum(sim_data[["Population mortality"]], na.rm = TRUE)
cat("\nTrue total deaths:", true_total, "\n")

# 5a. Complete case analysis (for comparison)
cat("\n--- Model 1: Complete Case Analysis ---\n")
fit_cc <- fit_dag(
  data = sim_data_mar,
  dag = dag,
  outcome = "Population mortality",
  spatial_levels = "region",
  complete_case = TRUE,                    # Drop rows with missing data
  family = negbinomial(),
  offset_col = "population",
  chains = 2,
  iter = 1000,
  cores = 2
)

# 5b. IPW analysis
cat("\n--- Model 2: IPW Analysis ---\n")
fit_ipw <- fit_dag(
  data = sim_data_mar,
  dag = dag,
  outcome = "Population mortality",
  spatial_levels = "region",
  complete_case = FALSE,                   # Handle missingness
  missing_method = "ipw",                  # Use inverse probability weighting
  family = negbinomial(),
  offset_col = "population",
  chains = 2,
  iter = 1000,
  cores = 2
)

# 5c. Oracle model (full data, no missingness)
cat("\n--- Model 3: Oracle (Full Data) ---\n")
fit_oracle <- fit_dag(
  data = sim_data,
  dag = dag,
  outcome = "Population mortality",
  spatial_levels = "region",
  complete_case = TRUE,
  family = negbinomial(),
  offset_col = "population",
  chains = 2,
  iter = 1000,
  cores = 2
)

# ============================================================================
# 6. COMPUTE RECOVERY METRICS
# ============================================================================

cat("\n========================================\n")
cat("RECOVERY METRICS\n")
cat("========================================\n")

# Function to get predictions and metrics
get_model_metrics <- function(fit, fit_name) {
  preds <- posterior_predict(fit)
  est_total <- mean(rowSums(preds))
  est_sd <- sd(rowSums(preds))
  ci <- quantile(rowSums(preds), c(0.025, 0.975))

  metrics <- compute_recovery_metrics(
    true_value = true_total,
    est_value = est_total,
    est_sd = est_sd,
    ci_lower = ci[1],
    ci_upper = ci[2]
  )

  cat(sprintf("\n%s:\n", fit_name))
  cat(sprintf("  Estimated total: %.0f (SD: %.0f)\n", est_total, est_sd))
  cat(sprintf("  95%% CI: [%.0f, %.0f]\n", ci[1], ci[2]))
  cat(sprintf("  Bias: %.0f (%.1f%%)\n", metrics$bias, 100 * metrics$rel_bias))
  cat(sprintf("  Coverage: %s\n", ifelse(metrics$coverage == 1, "Yes", "No")))

  return(data.frame(
    model = fit_name,
    true = true_total,
    estimated = est_total,
    bias = metrics$bias,
    rel_bias = metrics$rel_bias,
    coverage = metrics$coverage
  ))
}

results_cc <- get_model_metrics(fit_cc, "Complete Case")
results_ipw <- get_model_metrics(fit_ipw, "IPW")
results_oracle <- get_model_metrics(fit_oracle, "Oracle")

# Combine results
results_table <- rbind(results_cc, results_ipw, results_oracle)
cat("\n--- Summary Table ---\n")
print(results_table)

# ============================================================================
# 7. RUN FULL SIMULATION STUDY (OPTIONAL)
# ============================================================================

cat("\n========================================\n")
cat("SIMULATION STUDY EXAMPLE\n")
cat("========================================\n")
cat("(This would run multiple replicates across scenarios)\n")

# Create a small scenario grid for demonstration
demo_scenarios <- create_scenario_grid(
  n_govs = 5,                              # Smaller for demo
  n_months = 6,                            # Shorter for demo
  noise_levels = 1,
  regimes = c("MCAR", "MAR"),              # Skip MNAR for speed
  pi_bases = c(0.7),
  mnar_omegas = 0
)

cat("\nDemo scenarios:\n")
print(demo_scenarios)

# To run the full study (commented out for speed):
# study_results <- run_simulation_study(
#   dag = dag,
#   scenarios = demo_scenarios,
#   n_replicates = 5,                       # Use 5 for demo, 100+ for real study
#   mortality_node = "Population mortality",
#   predictors = available_predictors,
#   chains = 2,
#   iter = 500,
#   verbose = TRUE
# )
#
# # Summarize results
# summary_results <- summarize_study_results(study_results)
# print(summary_results)

cat("\n========================================\n")
cat("EXAMPLE COMPLETE\n")
cat("========================================\n")
