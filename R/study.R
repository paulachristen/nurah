#' Run Simulation Study for Mortality Estimation
#'
#' Main function for running a simulation study to evaluate mortality estimators
#' under various scenarios (effect strengths, noise levels, missingness regimes).
#'
#' @param dag A nurah_dag object defining the causal structure
#' @param scenarios Data frame defining simulation scenarios. Each row is a scenario with columns:
#'   - scenario_id: Unique identifier
#'   - n_govs: Number of governorates (spatial units)
#'   - n_months: Number of months (time periods)
#'   - initial_pop: Initial population per governorate
#'   - noise_level: Noise scaling factor (0-1)
#'   - missingness_regime: "MCAR", "MAR", or "MNAR"
#'   - pi_base: Baseline observation probability
#'   - mnar_omega: MNAR sensitivity parameter (if regime = "MNAR")
#'   - Any additional columns for scenario-specific parameters
#' @param models Named list of model specifications to fit. Each element should be a list with:
#'   - formula: Model formula (or NULL to use DAG-derived formula)
#'   - family: brms family (default: negbinomial())
#'   - complete_case: Whether to use complete case analysis (default: FALSE)
#'   - missing_method: "ipw" or "imputation" (if complete_case = FALSE)
#'   Default models include "complete_case", "ipw", and "full_obs" (oracle with no missingness)
#' @param n_replicates Number of simulation replicates per scenario (default: 100)
#' @param mortality_node Name of the mortality outcome node in the DAG
#' @param predictors Character vector of predictor variables to include in missingness
#' @param mar_drivers Character vector of variables driving MAR missingness
#' @param seed Random seed for reproducibility
#' @param chains Number of MCMC chains per model fit
#' @param iter Number of MCMC iterations per chain
#' @param cores Number of CPU cores to use
#' @param verbose Whether to print progress messages
#'
#' @return A data frame with one row per scenario x model x replicate containing:
#'   - scenario_id, model_id, replicate
#'   - true_total_deaths: True total deaths from simulation
#'   - est_total_deaths: Estimated total deaths from model
#'   - bias: Absolute bias (estimated - true)
#'   - rel_bias: Relative bias (bias / true)
#'   - rmse: Root mean squared error
#'   - coverage: Whether 95% CI contains true value
#'   - Additional governorate-level metrics
#'
#' @examples
#' \dontrun{
#' # Define scenarios
#' scenarios <- expand.grid(
#'   n_govs = 10,
#'   n_months = 12,
#'   initial_pop = 100000,
#'   noise_level = c(0.5, 1.0),
#'   missingness_regime = c("MCAR", "MAR", "MNAR"),
#'   pi_base = c(0.5, 0.8),
#'   mnar_omega = c(0, -0.5)
#' )
#' scenarios$scenario_id <- seq_len(nrow(scenarios))
#'
#' # Run study
#' results <- run_simulation_study(
#'   dag = checchi_2017_dag(dummy_checchi_2017_parameters()),
#'   scenarios = scenarios,
#'   n_replicates = 10,
#'   mortality_node = "Population mortality"
#' )
#' }
#'
#' @export
run_simulation_study <- function(dag,
                                  scenarios,
                                  models = NULL,
                                  n_replicates = 100,
                                  mortality_node = "Population mortality",
                                  predictors = NULL,
                                  mar_drivers = NULL,
                                  seed = NULL,
                                  chains = 2,
                                  iter = 1000,
                                  cores = 1,
                                  verbose = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  # Validate inputs
  if (!inherits(dag, "nurah_dag")) {
    stop("dag must be a nurah_dag object")
  }

  required_cols <- c("scenario_id", "n_govs", "n_months", "missingness_regime", "pi_base")
  missing_cols <- setdiff(required_cols, names(scenarios))
  if (length(missing_cols) > 0) {
    stop("scenarios data frame missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Set default models if not provided
  if (is.null(models)) {
    models <- list(
      complete_case = list(
        complete_case = TRUE,
        family = brms::negbinomial()
      ),
      ipw = list(
        complete_case = FALSE,
        missing_method = "ipw",
        family = brms::negbinomial()
      )
    )
  }

  # Set default predictors if not provided
  if (is.null(predictors)) {
    # Use parents of mortality node from DAG
    params <- dag$parameters
    if (!is.null(params)) {
      predictors <- unique(as.character(params$from[params$to == mortality_node]))
      predictors <- predictors[predictors %in% unique(c(params$from, params$to))]
    }
  }

  # Set default MAR drivers
  if (is.null(mar_drivers)) {
    mar_drivers <- c(
      "Exposure to armed attacks or mechanical force of nature",
      "Forced displacement"
    )
  }

  # Initialize results storage
  results_list <- list()
  result_idx <- 1

  n_scenarios <- nrow(scenarios)
  n_models <- length(models)
  total_iterations <- n_scenarios * n_replicates

  if (verbose) {
    message(sprintf("Starting simulation study: %d scenarios x %d replicates x %d models",
                    n_scenarios, n_replicates, n_models))
  }

  # Loop over scenarios
  for (s in seq_len(n_scenarios)) {
    scenario <- scenarios[s, ]
    scenario_id <- scenario$scenario_id

    if (verbose) {
      message(sprintf("\n=== Scenario %d/%d (ID: %s) ===", s, n_scenarios, scenario_id))
      message(sprintf("  n_govs=%d, n_months=%d, regime=%s, pi_base=%.2f",
                      scenario$n_govs, scenario$n_months,
                      scenario$missingness_regime, scenario$pi_base))
    }

    # Loop over replicates
    for (r in seq_len(n_replicates)) {
      if (verbose && r %% 10 == 0) {
        message(sprintf("  Replicate %d/%d", r, n_replicates))
      }

      # Set replicate-specific seed for reproducibility
      if (!is.null(seed)) {
        set.seed(seed + s * 1000 + r)
      }

      # ========== SIMULATE DATA ==========
      tryCatch({
        sim_data <- simulate_crisis_data(
          dag = dag,
          start_date = "2023-10-01",
          n_periods = scenario$n_months,
          resolution = "month",
          spatial_structure = paste0("gov", seq_len(scenario$n_govs)),
          initial_population = scenario$initial_pop %||% 100000,
          noise_level = scenario$noise_level %||% 1,
          mortality_node = mortality_node,
          mortality_phi = scenario$mortality_phi %||% 1,
          mortality_intercept = scenario$mortality_intercept %||% -8
        )

        # Store true values
        true_total_deaths <- sum(sim_data[[mortality_node]], na.rm = TRUE)
        true_deaths_by_gov <- stats::aggregate(
          sim_data[[mortality_node]],
          by = list(gov = sim_data[[names(sim_data)[1]]]),
          FUN = sum, na.rm = TRUE
        )
        names(true_deaths_by_gov)[2] <- "true_deaths"

        # ========== APPLY MISSINGNESS ==========

        # Identify which predictors exist in simulated data
        available_predictors <- intersect(predictors, names(sim_data))
        available_drivers <- intersect(mar_drivers, names(sim_data))

        if (length(available_predictors) > 0) {
          miss_result <- simulate_missingness(
            data = sim_data,
            predictors = available_predictors,
            regime = scenario$missingness_regime,
            pi_base = scenario$pi_base,
            mar_drivers = available_drivers,
            mnar_omega = scenario$mnar_omega %||% 0
          )
          sim_data_missing <- miss_result$data

          # Calculate actual missingness rate
          actual_miss_rate <- mean(is.na(sim_data_missing[, available_predictors]))
        } else {
          sim_data_missing <- sim_data
          actual_miss_rate <- 0
        }

        # ========== FIT MODELS ==========

        for (m_name in names(models)) {
          m_spec <- models[[m_name]]

          # Skip oracle model on missing data
          if (m_name == "full_obs") {
            fit_data <- sim_data  # Use complete data for oracle
          } else {
            fit_data <- sim_data_missing
          }

          fit_result <- tryCatch({
            fit <- fit_dag(
              data = fit_data,
              dag = dag,
              outcome = mortality_node,
              spatial_levels = names(sim_data)[1],  # governorate column
              complete_case = m_spec$complete_case %||% FALSE,
              missing_method = m_spec$missing_method %||% "ipw",
              family = m_spec$family %||% brms::negbinomial(),
              offset_col = "population",
              chains = chains,
              iter = iter,
              cores = cores,
              seed = seed,
              control = list(adapt_delta = 0.9),
              refresh = 0  # Suppress Stan output
            )

            # Get predictions
            preds <- brms::posterior_predict(fit)
            est_total_deaths <- mean(rowSums(preds))
            est_total_deaths_sd <- stats::sd(rowSums(preds))
            ci_lower <- stats::quantile(rowSums(preds), 0.025)
            ci_upper <- stats::quantile(rowSums(preds), 0.975)

            list(
              success = TRUE,
              est_total_deaths = est_total_deaths,
              est_total_deaths_sd = est_total_deaths_sd,
              ci_lower = ci_lower,
              ci_upper = ci_upper
            )
          }, error = function(e) {
            list(
              success = FALSE,
              error_msg = e$message,
              est_total_deaths = NA,
              est_total_deaths_sd = NA,
              ci_lower = NA,
              ci_upper = NA
            )
          })

          # ========== COMPUTE METRICS ==========

          metrics <- compute_recovery_metrics(
            true_value = true_total_deaths,
            est_value = fit_result$est_total_deaths,
            est_sd = fit_result$est_total_deaths_sd,
            ci_lower = fit_result$ci_lower,
            ci_upper = fit_result$ci_upper
          )

          # Store results
          results_list[[result_idx]] <- data.frame(
            scenario_id = scenario_id,
            replicate = r,
            model_id = m_name,
            n_govs = scenario$n_govs,
            n_months = scenario$n_months,
            noise_level = scenario$noise_level %||% 1,
            missingness_regime = scenario$missingness_regime,
            pi_base = scenario$pi_base,
            mnar_omega = scenario$mnar_omega %||% 0,
            actual_miss_rate = actual_miss_rate,
            true_total_deaths = true_total_deaths,
            est_total_deaths = fit_result$est_total_deaths,
            bias = metrics$bias,
            rel_bias = metrics$rel_bias,
            rmse = metrics$rmse,
            coverage = metrics$coverage,
            fit_success = fit_result$success,
            stringsAsFactors = FALSE
          )
          result_idx <- result_idx + 1
        }

      }, error = function(e) {
        if (verbose) {
          message(sprintf("  ERROR in replicate %d: %s", r, e$message))
        }
        # Store error result for all models
        for (m_name in names(models)) {
          results_list[[result_idx]] <- data.frame(
            scenario_id = scenario_id,
            replicate = r,
            model_id = m_name,
            n_govs = scenario$n_govs,
            n_months = scenario$n_months,
            noise_level = scenario$noise_level %||% 1,
            missingness_regime = scenario$missingness_regime,
            pi_base = scenario$pi_base,
            mnar_omega = scenario$mnar_omega %||% 0,
            actual_miss_rate = NA,
            true_total_deaths = NA,
            est_total_deaths = NA,
            bias = NA,
            rel_bias = NA,
            rmse = NA,
            coverage = NA,
            fit_success = FALSE,
            stringsAsFactors = FALSE
          )
          result_idx <- result_idx + 1
        }
      })
    }
  }

  # Combine results
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

  if (verbose) {
    message("\n=== Study Complete ===")
    message(sprintf("Total results: %d rows", nrow(results)))
    message(sprintf("Success rate: %.1f%%", 100 * mean(results$fit_success, na.rm = TRUE)))
  }

  return(results)
}


#' Compute Recovery Metrics for Simulation Study
#'
#' Computes bias, RMSE, and coverage metrics comparing estimated to true values.
#'
#' @param true_value True value from simulation
#' @param est_value Estimated value from model
#' @param est_sd Standard deviation of estimate
#' @param ci_lower Lower bound of confidence interval
#' @param ci_upper Upper bound of confidence interval
#'
#' @return List with bias, rel_bias, rmse, coverage
#' @export
compute_recovery_metrics <- function(true_value, est_value, est_sd = NULL,
                                      ci_lower = NULL, ci_upper = NULL) {

  # Handle NA/NULL inputs
  if (is.na(true_value) || is.null(true_value) ||
      is.na(est_value) || is.null(est_value)) {
    return(list(
      bias = NA_real_,
      rel_bias = NA_real_,
      rmse = NA_real_,
      coverage = NA
    ))
  }

  # Bias
  bias <- est_value - true_value

  # Relative bias (as proportion)
  rel_bias <- if (abs(true_value) > 1e-6) bias / true_value else NA_real_


  # RMSE (for single estimate, this equals |bias|; more meaningful across replicates)
  rmse <- sqrt(bias^2)

  # Coverage: does CI contain true value?
  coverage <- if (!is.null(ci_lower) && !is.null(ci_upper) &&
                   !is.na(ci_lower) && !is.na(ci_upper)) {
    as.integer(ci_lower <= true_value && true_value <= ci_upper)
  } else {
    NA_integer_
  }

  return(list(
    bias = bias,
    rel_bias = rel_bias,
    rmse = rmse,
    coverage = coverage
  ))
}


#' Summarize Simulation Study Results
#'
#' Aggregates results across replicates to compute mean bias, RMSE, and coverage.
#'
#' @param results Data frame from run_simulation_study
#' @param by Character vector of grouping variables (default: scenario_id, model_id)
#'
#' @return Data frame with summary statistics by group
#' @export
summarize_study_results <- function(results, by = c("scenario_id", "model_id")) {

  # Filter to successful fits
  results_success <- results[results$fit_success == TRUE, ]

  if (nrow(results_success) == 0) {
    warning("No successful fits to summarize")
    return(NULL)
  }

  # Aggregate
  summary_df <- stats::aggregate(
    cbind(bias, rel_bias, coverage, true_total_deaths, est_total_deaths) ~
      scenario_id + model_id + missingness_regime + pi_base + mnar_omega,
    data = results_success,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      n = sum(!is.na(x))
    )
  )

  # Flatten the matrix columns
  summary_flat <- data.frame(
    scenario_id = summary_df$scenario_id,
    model_id = summary_df$model_id,
    missingness_regime = summary_df$missingness_regime,
    pi_base = summary_df$pi_base,
    mnar_omega = summary_df$mnar_omega,
    mean_bias = summary_df$bias[, "mean"],
    sd_bias = summary_df$bias[, "sd"],
    mean_rel_bias = summary_df$rel_bias[, "mean"],
    coverage_rate = summary_df$coverage[, "mean"],
    n_successful = summary_df$coverage[, "n"],
    mean_true_deaths = summary_df$true_total_deaths[, "mean"],
    mean_est_deaths = summary_df$est_total_deaths[, "mean"],
    stringsAsFactors = FALSE
  )

  # Compute RMSE across replicates
  for (i in seq_len(nrow(summary_flat))) {
    idx <- results_success$scenario_id == summary_flat$scenario_id[i] &
           results_success$model_id == summary_flat$model_id[i]
    biases <- results_success$bias[idx]
    summary_flat$rmse[i] <- sqrt(mean(biases^2, na.rm = TRUE))
  }

  return(summary_flat)
}


#' Create Default Scenario Grid
#'
#' Helper function to create a standard scenario grid for simulation studies.
#'
#' @param n_govs Vector of governorate counts to test
#' @param n_months Vector of month counts to test
#' @param noise_levels Vector of noise levels to test
#' @param regimes Missingness regimes to test
#' @param pi_bases Baseline observation probabilities to test
#' @param mnar_omegas MNAR sensitivity parameters to test
#'
#' @return Data frame of scenarios with scenario_id
#' @export
create_scenario_grid <- function(n_govs = 10,
                                  n_months = 12,
                                  noise_levels = c(0.5, 1),
                                  regimes = c("MCAR", "MAR", "MNAR"),
                                  pi_bases = c(0.5, 0.7, 0.9),
                                  mnar_omegas = c(0, -0.3, -0.6)) {

  # Create base grid
  grid <- expand.grid(
    n_govs = n_govs,
    n_months = n_months,
    initial_pop = 100000,
    noise_level = noise_levels,
    missingness_regime = regimes,
    pi_base = pi_bases,
    mnar_omega = mnar_omegas,
    stringsAsFactors = FALSE
  )

  # Remove invalid combinations: mnar_omega only matters for MNAR regime
  grid <- grid[!(grid$missingness_regime != "MNAR" & grid$mnar_omega != 0), ]

  # Add scenario ID
  grid$scenario_id <- seq_len(nrow(grid))

  # Reorder columns
  grid <- grid[, c("scenario_id", "n_govs", "n_months", "initial_pop",
                    "noise_level", "missingness_regime", "pi_base", "mnar_omega")]

  return(grid)
}
