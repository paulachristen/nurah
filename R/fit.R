#' Fit DAG-based Bayesian Hierarchical Model
#'
#' This function constructs and fits a Bayesian hierarchical linear model (via **brms** and Stan)
#' to estimate "indirect mortality" given a user-specified causal DAG and observed data. It uses
#' the DAG to determine which predictors (direct causes of the outcome) to include, and adds
#' specified spatial hierarchy levels as random intercept effects.
#'
#' @param data A data frame of observed data. Must contain the `outcome` variable and all predictor and grouping variables.
#' @param dag A `nurah_dag` object defining the causal structure. This DAG is used to extract the direct causes of `outcome` as fixed-effect predictors.
#' @param outcome Name (string) of the dependent variable in `data` to model.
#' @param spatial_levels Character vector of spatial grouping variables for hierarchical random effects.
#' @param complete_case Logical. If TRUE, drop rows with any missing values (default behavior pre-update).
#'   If FALSE (default), handle missingness using IPW or imputation.
#' @param missing_method Method for handling missing data when complete_case = FALSE.
#'   Options: "ipw" (inverse probability weighting), "imputation" (multiple imputation via mice).
#'   Default is "ipw".
#' @param ipw_formula Formula for the observation/missingness model (logistic regression predicting
#'   whether each row is fully observed). If NULL and missing_method = "ipw", a default formula
#'   using available non-outcome columns will be constructed.
#' @param ipw_stabilize Logical. If TRUE, use stabilized IPW weights (recommended). Default TRUE.
#' @param n_imputations Number of imputations if missing_method = "imputation". Default 5.
#' @param priors Optional brms prior specification. If `NULL`, defaults to weakly informative priors.
#' @param chains Number of MCMC chains to run (defaults to 4).
#' @param cores Number of CPU cores to use for parallel sampling.
#' @param iter Number of iterations per chain (including warmup; default 2000).
#' @param warmup Number of warmup (burn-in) iterations per chain.
#' @param seed Random seed for reproducibility.
#' @param control A list of control parameters passed to Stan.
#' @param formula_override Optional formula object to use as the model formula.
#' @param family The likelihood family for the outcome. Defaults to `gaussian()` for continuous
#'   outcomes. For mortality counts, use `negbinomial()` with offset.
#' @param offset_col Name of the population/offset column for count models. If provided, log(offset)
#'   will be added to the model.
#' @param ... Additional arguments passed to `brms::brm()`.
#'
#' @details
#' **Missingness Handling (NEW):**
#' When `complete_case = FALSE`, the function implements one of two approaches:
#'
#' 1. **IPW (Inverse Probability Weighting)**: Estimates P(fully observed | covariates) using
#'    logistic regression, then weights each observation by 1/P in the outcome model. This
#'    corrects for selection bias under MAR (missing at random) assumptions. Stabilized weights
#'    (recommended) multiply by marginal P(observed) to reduce variance.
#'
#' 2. **Multiple Imputation**: Uses the `mice` package to generate multiple imputed datasets,
#'    fits the model to each, and pools results using Rubin's rules.
#'
#' @return A `brmsfit` object (or list of fits for imputation) containing the fitted Bayesian model.
#'
#' @examples
#' \dontrun{
#' # Fit with IPW to handle missingness
#' fit <- fit_dag(
#'   data = sim_data_missing,
#'   dag = dag,
#'   outcome = "Population mortality",
#'   spatial_levels = "region",
#'   complete_case = FALSE,
#'   missing_method = "ipw",
#'   family = negbinomial(),
#'   offset_col = "population"
#' )
#' }
#'
#' @export
fit_dag <- function(data, dag, outcome,
                    spatial_levels = NULL,
                    complete_case = FALSE,
                    missing_method = c("ipw", "imputation"),
                    ipw_formula = NULL,
                    ipw_stabilize = TRUE,
                    n_imputations = 5,
                    priors = NULL,
                    chains = 4,
                    cores = getOption("mc.cores", 1),
                    iter = 2000,
                    warmup = floor(iter / 2),
                    seed = NULL,
                    control = list(),
                    formula_override = NULL,
                    family = gaussian(),
                    offset_col = NULL,
                    ...) {

  missing_method <- match.arg(missing_method)

  # Basic input validations
  if (!is.character(outcome) || length(outcome) != 1) {
    stop("`outcome` must be a single string (name of the outcome variable).")
  }
  if (!(outcome %in% names(data))) {
    stop("Outcome variable '", outcome, "' not found in `data`.")
  }


  # Helper function to make syntactically valid names
  make_safe_name <- function(x) {
    gsub("[^a-zA-Z0-9]", "_", x)
  }

  # Create column name mapping (original -> safe)
  original_names <- names(data)
  safe_names <- make_safe_name(original_names)
  name_map <- setNames(original_names, safe_names)

  # Rename columns in data
  names(data) <- safe_names

  # Update outcome and spatial_levels to safe names
  outcome_safe <- make_safe_name(outcome)
  if (!is.null(spatial_levels)) {
    spatial_levels <- make_safe_name(spatial_levels)
  }
  if (!is.null(offset_col)) {
    offset_col <- make_safe_name(offset_col)
  }

  # Construct fixed-effects formula from DAG (unless overridden)
  if (is.null(formula_override)) {
    if (!inherits(dag, "nurah_dag")) {
      stop("`dag` must be a 'nurah_dag' object created with define_dag().")
    }
    # Extract direct causes (parents) of the outcome from the DAG
    direct_causes <- NULL
    if ("dagitty" %in% class(dag) || !is.null(dag$dag)) {
      # Underlying dagitty object may be the object itself or stored in dag$dag
      dag_obj <- if ("dagitty" %in% class(dag)) dag else dag$dag
      direct_causes <- tryCatch(
        dagitty::parents(dag_obj, outcome),
        error = function(e) NULL
      )
    }

    # Also try to get parents from parameters dataframe
    if (is.null(direct_causes) && !is.null(dag$parameters)) {
      params_to_outcome <- dag$parameters[dag$parameters$to == outcome, ]
      if (nrow(params_to_outcome) > 0) {
        direct_causes <- unique(as.character(params_to_outcome$from))
      }
    }

    if (is.null(direct_causes) || length(direct_causes) == 0) {
      warning("No parents of '", outcome, "' defined in DAG; using an intercept-only model.")
      direct_causes <- character(0)
    }

    # Convert direct causes to safe names and filter to those in data
    direct_causes_safe <- make_safe_name(direct_causes)
    direct_causes_safe <- intersect(direct_causes_safe, names(data))

    # Build the fixed effects formula string
    if (length(direct_causes_safe) == 0) {
      fixed_formula_str <- paste(outcome_safe, "~ 1")
    } else {
      fixed_formula_str <- paste(outcome_safe, "~", paste(direct_causes_safe, collapse = " + "))
    }
    formula_fixed <- stats::as.formula(fixed_formula_str)
  } else {
    if (!inherits(formula_override, "formula")) {
      stop("`formula_override` must be a formula object (e.g., outcome ~ predictors) if provided.")
    }
    formula_fixed <- formula_override
    direct_causes_safe <- all.vars(formula_fixed[[3]])
  }

  # Construct random-effects formula for spatial hierarchy if specified
  random_formula_str <- NULL
  if (!is.null(spatial_levels) && length(spatial_levels) > 0) {
    # Ensure grouping variables exist in data
    missing_groups <- setdiff(spatial_levels, names(data))
    if (length(missing_groups) > 0) {
      stop("Grouping variable(s) not found in data: ", paste(missing_groups, collapse = ", "))
    }
    if (length(spatial_levels) > 1) {
      # Assume nested hierarchy in the given order
      random_formula_str <- paste0("(1|", paste(spatial_levels, collapse = "/"), ")")
    } else {
      random_formula_str <- paste0("(1|", spatial_levels, ")")
    }
  }

  # Add offset for count models
  offset_str <- NULL
  if (!is.null(offset_col) && offset_col %in% names(data)) {
    # Create safe column name for log offset
    log_offset_name <- paste0("log_", offset_col)
    # Ensure positive values for log offset
    data[[log_offset_name]] <- log(pmax(data[[offset_col]], 1))
    offset_str <- paste0("offset(", log_offset_name, ")")
  }

  # Combine fixed, random, and offset parts into the final formula
  formula_str <- deparse(formula_fixed, width.cutoff = 500)
  if (length(formula_str) > 1) {
    formula_str <- paste(formula_str, collapse = " ")
  }
  if (!is.null(offset_str)) {
    formula_str <- paste(formula_str, "+", offset_str)
  }
  if (!is.null(random_formula_str)) {
    formula_str <- paste(formula_str, "+", random_formula_str)
  }
  final_formula <- stats::as.formula(formula_str)

  # Determine variables needed for model
  vars_needed <- c(outcome_safe, direct_causes_safe)
  if (!is.null(spatial_levels)) {
    vars_needed <- c(vars_needed, spatial_levels)
  }
  if (!is.null(offset_col)) {
    vars_needed <- c(vars_needed, offset_col)
  }
  vars_needed <- unique(vars_needed)

  # Handle missingness
  if (complete_case) {
    # Original behavior: drop rows with any missing values
    data_complete <- data[stats::complete.cases(data[, vars_needed, drop = FALSE]), ]
    if (nrow(data_complete) < nrow(data)) {
      message("Dropped ", nrow(data) - nrow(data_complete),
              " rows due to missing values in outcome or predictors (complete_case = TRUE).")
    }
    weights_vec <- NULL

  } else if (missing_method == "ipw") {
    # IPW approach
    ipw_result <- compute_ipw_weights(
      data = data,
      vars_needed = vars_needed,
      outcome = outcome,
      ipw_formula = ipw_formula,
      stabilize = ipw_stabilize
    )
    data_complete <- ipw_result$data
    weights_vec <- ipw_result$weights
    message("Using IPW with ", sum(!is.na(weights_vec)), " complete observations.")

  } else if (missing_method == "imputation") {
    # Multiple imputation approach
    return(fit_with_imputation(
      data = data,
      final_formula = final_formula,
      vars_needed = vars_needed,
      n_imputations = n_imputations,
      priors = priors,
      chains = chains,
      cores = cores,
      iter = iter,
      warmup = warmup,
      seed = seed,
      control = control,
      family = family,
      ...
    ))
  }

  # Set default priors if none provided
  if (is.null(priors)) {
    priors <- c(
      brms::prior("normal(0, 1)", class = "b"),        # Coefficients ~ Normal(0, 1)
      brms::prior("normal(0, 1)", class = "sd"),      # Random-effect SDs ~ Half-Normal(0, 1)
      brms::prior("normal(0, 1)", class = "Intercept")
    )
    # Add sigma prior only for Gaussian family
    if (inherits(family, "family") && family$family == "gaussian") {
      priors <- c(priors, brms::prior("normal(0, 1)", class = "sigma"))
    }
  }

  # Fit the Bayesian model using brms (HMC via Stan)
  if (!is.null(weights_vec)) {
    # Add weights as a column in data_complete
    data_complete$ipw_weights <- weights_vec

    # Reconstruct formula with weights on LHS: outcome | weights(w) ~ predictors
    formula_char <- deparse(final_formula, width.cutoff = 500)
    if (length(formula_char) > 1) formula_char <- paste(formula_char, collapse = " ")

    # Split at ~ and add weights
    formula_parts <- strsplit(formula_char, "~")[[1]]
    lhs <- trimws(formula_parts[1])
    rhs <- trimws(paste(formula_parts[-1], collapse = "~"))
    weighted_formula_str <- paste0(lhs, " | weights(ipw_weights) ~ ", rhs)
    weighted_formula <- stats::as.formula(weighted_formula_str)

    fit <- brms::brm(
      formula = weighted_formula,
      data    = data_complete,
      family  = family,
      prior   = priors,
      chains  = chains,
      cores   = cores,
      iter    = iter,
      warmup  = warmup,
      seed    = seed,
      control = control,
      ...
    )
  } else {
    fit <- brms::brm(
      formula = final_formula,
      data    = data_complete,
      family  = family,
      prior   = priors,
      chains  = chains,
      cores   = cores,
      iter    = iter,
      warmup  = warmup,
      seed    = seed,
      control = control,
      ...
    )
  }

  return(fit)
}


#' Compute Inverse Probability Weights for Missingness
#'
#' Fits a logistic regression model for P(fully observed | covariates) and
#' returns inverse probability weights.
#'
#' @param data Data frame
#' @param vars_needed Variables needed for the outcome model
#' @param outcome Name of outcome variable
#' @param ipw_formula Optional formula for observation model
#' @param stabilize Whether to use stabilized weights
#'
#' @return List with data (complete cases only) and weights vector
#' @keywords internal
compute_ipw_weights <- function(data, vars_needed, outcome, ipw_formula = NULL, stabilize = TRUE) {

  # Create indicator for complete observation
  data$is_complete <- stats::complete.cases(data[, vars_needed, drop = FALSE])

  # Identify predictors for missingness model
  # Use variables that are fully observed to predict missingness
  potential_predictors <- setdiff(names(data), c(outcome, "is_complete"))

  # Find columns with no missing values
  complete_cols <- names(data)[sapply(data, function(x) !any(is.na(x)))]
  ipw_predictors <- intersect(potential_predictors, complete_cols)

  # Keep only numeric or factor predictors
  ipw_predictors <- ipw_predictors[sapply(data[, ipw_predictors, drop = FALSE],
                                          function(x) is.numeric(x) || is.factor(x))]

  if (length(ipw_predictors) == 0) {
    warning("No complete predictors available for IPW model. Using marginal weights.")
    # Fall back to marginal probability
    p_obs <- mean(data$is_complete)
    weights <- ifelse(data$is_complete, 1/p_obs, NA)
    return(list(
      data = data[data$is_complete, ],
      weights = weights[data$is_complete]
    ))
  }

  # Build IPW formula
  if (is.null(ipw_formula)) {
    # Limit number of predictors to avoid overfitting
    if (length(ipw_predictors) > 10) {
      ipw_predictors <- ipw_predictors[1:10]
    }
    ipw_formula <- stats::as.formula(
      paste("is_complete ~", paste(ipw_predictors, collapse = " + "))
    )
  }

  # Fit logistic regression for P(observed)
  ipw_model <- tryCatch(
    stats::glm(ipw_formula, data = data, family = stats::binomial()),
    error = function(e) {
      warning("IPW model fitting failed: ", e$message, ". Using marginal weights.")
      NULL
    }
  )

  if (is.null(ipw_model)) {
    p_obs <- mean(data$is_complete)
    weights <- ifelse(data$is_complete, 1/p_obs, NA)
    return(list(
      data = data[data$is_complete, ],
      weights = weights[data$is_complete]
    ))
  }

  # Get predicted probabilities
  p_hat <- stats::predict(ipw_model, type = "response")

  # Compute weights (only for complete cases)
  if (stabilize) {
    # Stabilized weights: multiply by marginal P(observed)
    p_marginal <- mean(data$is_complete)
    weights <- ifelse(data$is_complete, p_marginal / p_hat, NA)
  } else {
    weights <- ifelse(data$is_complete, 1 / p_hat, NA)
  }

  # Truncate extreme weights (at 1st and 99th percentile)
  weight_bounds <- stats::quantile(weights[!is.na(weights)], c(0.01, 0.99), na.rm = TRUE)
  weights <- pmin(pmax(weights, weight_bounds[1]), weight_bounds[2])

  # Return complete cases with weights
  complete_idx <- which(data$is_complete)

  return(list(
    data = data[complete_idx, ],
    weights = weights[complete_idx]
  ))
}


#' Fit Model with Multiple Imputation
#'
#' Uses mice package for multiple imputation, fits model to each imputed dataset,
#' and pools results.
#'
#' @param data Data frame with missing values
#' @param final_formula Model formula
#' @param vars_needed Variables needed for model
#' @param n_imputations Number of imputations
#' @param priors brms priors
#' @param chains Number of MCMC chains
#' @param cores Number of CPU cores
#' @param iter MCMC iterations
#' @param warmup Warmup iterations
#' @param seed Random seed
#' @param control Stan control parameters
#' @param family Model family
#' @param ... Additional arguments for brm
#'
#' @return List of brmsfit objects (one per imputation)
#' @keywords internal
fit_with_imputation <- function(data, final_formula, vars_needed, n_imputations = 5,
                                priors = NULL, chains = 4, cores = 1, iter = 2000,
                                warmup = floor(iter/2), seed = NULL, control = list(),
                                family = gaussian(), ...) {

  # Check if mice is available
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required for multiple imputation. ",
         "Install it with: install.packages('mice')")
  }

  # Subset to relevant variables for imputation
  data_subset <- data[, intersect(vars_needed, names(data)), drop = FALSE]

  # Run multiple imputation
  message("Running multiple imputation with ", n_imputations, " imputations...")
  imp <- mice::mice(data_subset, m = n_imputations, printFlag = FALSE, seed = seed)

  # Set default priors if none provided
  if (is.null(priors)) {
    priors <- c(
      brms::prior("normal(0, 1)", class = "b"),
      brms::prior("normal(0, 1)", class = "sd"),
      brms::prior("normal(0, 1)", class = "Intercept")
    )
    if (inherits(family, "family") && family$family == "gaussian") {
      priors <- c(priors, brms::prior("normal(0, 1)", class = "sigma"))
    }
  }

  # Fit model to each imputed dataset
  fits <- list()
  for (i in seq_len(n_imputations)) {
    message("Fitting model to imputation ", i, " of ", n_imputations)
    data_imp <- mice::complete(imp, i)

    fits[[i]] <- brms::brm(
      formula = final_formula,
      data    = data_imp,
      family  = family,
      prior   = priors,
      chains  = chains,
      cores   = cores,
      iter    = iter,
      warmup  = warmup,
      seed    = if (!is.null(seed)) seed + i else NULL,
      control = control,
      ...
    )
  }

  # Add pooling summary
  attr(fits, "n_imputations") <- n_imputations
  class(fits) <- c("brmsfit_mi", class(fits))

  return(fits)
}


#' Pool Results from Multiple Imputation Fits
#'
#' Combines posterior samples from multiple brmsfit objects using Rubin's rules.
#'
#' @param fits List of brmsfit objects from fit_with_imputation
#'
#' @return Data frame with pooled coefficient estimates, SE, and confidence intervals
#' @export
pool_imputation_fits <- function(fits) {
  if (!inherits(fits, "brmsfit_mi")) {
    stop("Input must be result from fit_with_imputation (class 'brmsfit_mi')")
  }

  m <- length(fits)

  # Extract fixed effects from each fit
  fixed_list <- lapply(fits, function(f) {
    fe <- brms::fixef(f)
    data.frame(
      parameter = rownames(fe),
      estimate = fe[, "Estimate"],
      se = fe[, "Est.Error"],
      stringsAsFactors = FALSE
    )
  })

  # Get parameter names
  params <- fixed_list[[1]]$parameter

  # Pool estimates using Rubin's rules
  pooled <- data.frame(
    parameter = params,
    estimate = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    stringsAsFactors = FALSE
  )

  for (p in params) {
    # Extract estimates and variances for this parameter across imputations
    q_bar <- mean(sapply(fixed_list, function(x) x$estimate[x$parameter == p]))
    u_bar <- mean(sapply(fixed_list, function(x) x$se[x$parameter == p]^2))  # within-imputation variance
    b <- stats::var(sapply(fixed_list, function(x) x$estimate[x$parameter == p]))  # between-imputation variance

    # Total variance (Rubin's rules)
    t_var <- u_bar + (1 + 1/m) * b

    pooled$estimate[pooled$parameter == p] <- q_bar
    pooled$se[pooled$parameter == p] <- sqrt(t_var)
    pooled$ci_lower[pooled$parameter == p] <- q_bar - 1.96 * sqrt(t_var)
    pooled$ci_upper[pooled$parameter == p] <- q_bar + 1.96 * sqrt(t_var)
  }

  return(pooled)
}


#' Extract Predicted Deaths from Fitted Model
#'
#' Convenience function to get predicted mortality counts from a fitted model.
#'
#' @param fit A brmsfit object
#' @param newdata Data frame for prediction (if NULL, uses original data)
#' @param summary Whether to return summary statistics (default TRUE)
#'
#' @return Data frame with predictions
#' @export
predict_deaths <- function(fit, newdata = NULL, summary = TRUE) {
  if (!inherits(fit, "brmsfit")) {
    stop("fit must be a brmsfit object")
  }

  preds <- brms::posterior_predict(fit, newdata = newdata)

  if (summary) {
    pred_summary <- data.frame(
      mean = colMeans(preds),
      median = apply(preds, 2, stats::median),
      sd = apply(preds, 2, stats::sd),
      q025 = apply(preds, 2, stats::quantile, 0.025),
      q975 = apply(preds, 2, stats::quantile, 0.975)
    )
    return(pred_summary)
  } else {
    return(preds)
  }
}

