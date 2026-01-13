#' Simulate Crisis Data Based on a DAG Model
#'
#' Generates a synthetic crisis dataset over time and across regions, using a
#' directed acyclic graph (DAG) of relationships between indicators. This high-level
#' function orchestrates the simulation by setting up time and space dimensions,
#' simulating each indicator (node) in the DAG, handling dynamic population changes
#' due to IDP flows, and aggregating results to the desired time resolution.
#'
#' @param dag A DAG object of class \code{nurah_dag} defining the nodes (indicators)
#' and their relationships. Typically created by \code{define_dag()} or a specific
#' DAG constructor (e.g., \code{checchi_2017_dag()}). If \code{NULL}, a default DAG
#' will be used (if available).
#' @param start_date Start date of the simulation (as a \code{Date} or string in
#' "YYYY-MM-DD" format).
#' @param end_date End date of the simulation (as a \code{Date} or string). Either
#' \code{end_date} or \code{n_periods} must be provided (but not both).
#' @param n_periods Number of time periods to simulate at the specified output
#' resolution (ignored if \code{end_date} is provided).
#' @param resolution Time resolution for the output data. One of \code{"day"},
#' \code{"week"}, \code{"month"}, \code{"quarter"}, or \code{"year"}. Defaults to
#' \code{"month"}.
#' @param spatial_structure Specification of the spatial layout (regions and optional
#' hierarchy). Can be:
#'   \itemize{
#'     \item A data frame with one row per lowest-level region (and optional higher
#'           level columns for hierarchy).
#'     \item A named numeric vector defining counts at each hierarchy level (top to
#'           bottom). For example, \code{c(country = 1, region = 2, district = 4)}
#'           means 1 country, 2 regions total, 4 districts total (distributed as
#'           evenly as possible under the regions).
#'     \item A single numeric value (number of lowest-level regions; one top-level
#'           grouping will be added if >1 region is implied).
#'     \item A character vector of region names (defining the lowest level; one
#'           top-level grouping will be added since no higher levels are specified).
#'   }
#'   Defaults to \code{1} (a single region with one top-level grouping if needed).
#' @param initial_population Initial population at the start date for each lowest-level
#' region. Can be a single number (applied to all regions) or a vector. If a vector,
#' it can be:
#'   \itemize{
#'     \item A named vector, where names correspond to region identifiers (at the
#'           lowest level).
#'     \item An unnamed vector with length equal to the number of lowest-level
#'           regions (in the same order as the regions in the spatial structure).
#'   }
#'   Defaults to 10000 for all regions.
#' @param noise_level Numeric factor to scale stochastic noise in node simulations.
#' \code{noise_level = 0} produces fully deterministic output, while \code{noise_level = 1}
#' uses the full specified randomness. Defaults to 1.
#' @param mortality_node Character name of the mortality outcome node. If specified,
#' this node will be simulated using negative binomial distribution with population offset.
#' @param mortality_phi Overdispersion parameter for negative binomial mortality. Default 1.
#' @param mortality_intercept Baseline log-rate for mortality (on log scale). Default -8
#' (approximately 0.03 deaths per 1000 per day).
#'
#' @return A data frame of simulated indicators aggregated at the requested time
#' resolution (default: monthly). Contains one row per governorate per month with columns
#' for each simulated indicator including mortality counts (Y_it) and population (N_it).
#'
#' @examples
#' # Simulate 10 governorates x 12 months with mortality outcome
#' dag <- checchi_2017_dag(parameters = dummy_checchi_2017_parameters())
#' sim_data <- simulate_crisis_data(
#'   dag,
#'   start_date = "2023-10-01",
#'   n_periods = 12,
#'   resolution = "month",
#'   spatial_structure = paste0("gov", 1:10),
#'   initial_population = 100000,
#'   noise_level = 1,
#'   mortality_node = "Population mortality"
#' )
#' head(sim_data)
#' @export
simulate_crisis_data <- function(dag = NULL,
                                 start_date,
                                 end_date = NULL,
                                 n_periods = NULL,
                                 resolution = "month",
                                 spatial_structure = 1,
                                 initial_population = 10000,
                                 noise_level = 1,
                                 mortality_node = NULL,
                                 mortality_phi = 1,
                                 mortality_intercept = -8) {
  # If no DAG provided, use a default if available (e.g., checchi_2017_dag)
  if (is.null(dag)) {
    if (exists("checchi_2017_dag", mode = "function")) {
      dag <- checchi_2017_dag()
    } else {
      stop("No DAG provided and no default DAG available.")
    }
  }
  if (!inherits(dag, "nurah_dag")) {
    stop("The 'dag' must be a 'nurah_dag' object (use define_dag() to create one).")
  }

  # Set up temporal structure (daily dates and output period start dates)
  temporal <- setup_temporal_structure(start_date = start_date, end_date = end_date,
                                       n_periods = n_periods, resolution = resolution)
  time_seq_daily <- temporal$time_seq_daily
  time_seq_out   <- temporal$time_seq_out

  # Set up spatial structure (region hierarchy data frame)
  region_df <- setup_spatial_structure(spatial_structure)
  n_regions <- nrow(region_df)

  # Create a base daily data frame: one row per region per day
  n_days <- length(time_seq_daily)
  # Repeat each region entry for each day (region_df rows repeated for each day)
  daily_data <- region_df[rep(seq_len(n_regions), each = n_days), , drop = FALSE]
  # Add a date column with the sequence of daily dates for each region
  daily_data$date <- rep(time_seq_daily, times = n_regions)
  # Ensure daily_data is ordered by region then date (by construction it is in that order)

  # Add governorate index for convenience
  region_id_col <- names(region_df)[ncol(region_df)]
  daily_data$gov_id <- as.numeric(factor(daily_data[[region_id_col]]))

  # Add time index (month within simulation)
  daily_data$time_id <- as.numeric(factor(format(daily_data$date, "%Y-%m")))

  # Determine simulation order for nodes (topologically sorted)
  node_order <- topological_sort_dag(dag)
  # Identify if a dynamic population node exists in DAG (case-insensitive match for "population")
  pop_node <- node_order[tolower(node_order) %in% "population"]
  if (length(pop_node) > 0) {
    # Exclude the population node from static simulation (it will be handled separately)
    node_order <- setdiff(node_order, pop_node)
  }

  # Identify mortality node for special handling
  mortality_in_order <- NULL
  if (!is.null(mortality_node) && mortality_node %in% node_order) {
    mortality_in_order <- mortality_node
    node_order <- setdiff(node_order, mortality_node)
  }

  # Simulate each node in topological order (except dynamic population and mortality)
  for (node in node_order) {
    # Determine this node's parent nodes from the DAG structure
    parents_df <- dag$parameters[dag$parameters$to == node, , drop = FALSE]
    parent_names <- as.character(parents_df$from)
    # Prepare parent data for simulation
    if (length(parent_names) == 0) {
      # No parents: use an empty data frame with correct number of rows
      parent_data <- daily_data[, FALSE]
    } else {
      parent_data <- daily_data[, parent_names, drop = FALSE]
      # Apply time lags: shift each parent's values by its lag (in days) within each region
      if ("lag" %in% names(parents_df)) {
        # Ensure daily_data is grouped by region for lagging
        # We know each region's data is in blocks of n_days rows
        for (i in seq_len(nrow(parents_df))) {
          parent <- parents_df$from[i]
          lag_days <- ifelse(is.na(parents_df$lag[i]), 0, parents_df$lag[i])
          if (lag_days > 0) {
            # Extract the parent vector and reshape into matrix [n_days x n_regions]
            parent_vec <- daily_data[[parent]]
            parent_mat <- matrix(parent_vec, nrow = n_days, ncol = n_regions, byrow = FALSE)
            if (lag_days >= n_days) {
              # If lag is larger or equal to the time series length, all values become 0 (no previous data within period)
              parent_mat[] <- 0
            } else {
              # Shift each column (region) down by lag_days, pad the top with zeros
              parent_mat[(lag_days + 1):n_days, ] <- parent_mat[1:(n_days - lag_days), ]
              parent_mat[1:lag_days, ] <- 0
            }
            # Flatten matrix back to vector and assign to parent_data
            parent_data[[as.character(parent)]] <- as.vector(parent_mat)
          }
        }
      }
    }
    # Assemble node parameter list for simulation
    node_params <- list(parents = parent_names)
    # Include effect sizes corresponding to each parent (for use in simulate_node)
    if ("effect_size" %in% names(parents_df)) {
      node_params$effect_size <- parents_df$effect_size
    }
    # Include distribution or other parameters if they exist for this node (not provided in nurah_dag by default)
    # (No additional distribution parameters are extracted here, so defaults in simulate_node will apply)

    # Simulate this node's values across all regions and days
    values <- simulate_node(node = node, parents_data = parent_data,
                            dag_params = node_params, noise_level = noise_level)
    # Append the simulated values as a new column in daily_data
    daily_data[[node]] <- values
  }

  # Handle dynamic population simulation (if applicable)
  if (length(pop_node) > 0) {
    # If a population node exists, simulate population over time using IDP flows
    pop_name <- pop_node[1]  # (Assume at most one population node)
    daily_data <- simulate_population_dynamics(daily_data, initial_population = initial_population)
    # (If the DAG had dependents on population, those could be simulated after population, but typically none do.)
  } else {
    # If no population node in DAG, set a constant population (if initial_population given)
    daily_data$population <- NA  # initialize column
    # Assign initial_population to each lowest-level region (repeating daily)
    daily_data$population <- unlist(lapply(split(daily_data, daily_data[, names(region_df)]), function(df) {
      # Each group (region) gets its initial population value
      region_id <- df[1, ncol(region_df)]  # lowest-level region identifier value
      init_val <- if (length(initial_population) == 1) {
        initial_population
      } else if (!is.null(names(initial_population))) {
        initial_population[as.character(region_id)]
      } else {
        # If unnamed vector of correct length, assume ordering matches regions in region_df
        initial_population
      }
      rep(init_val, nrow(df))
    }))
  }

  # Simulate mortality node using negative binomial if specified
  if (!is.null(mortality_in_order)) {
    daily_data <- simulate_mortality_negbin(
      daily_data = daily_data,
      dag = dag,
      mortality_node = mortality_in_order,
      phi = mortality_phi,
      intercept = mortality_intercept,
      noise_level = noise_level,
      n_days = n_days,
      n_regions = n_regions
    )
  }

  # Aggregate daily data to the desired output resolution
  if (tolower(resolution) == "day") {
    # If output resolution is daily, return the daily_data (ordered by date and region)
    output_data <- daily_data[order(daily_data$date, daily_data[, ncol(region_df)]), ]
    rownames(output_data) <- NULL
  } else {
    output_data <- aggregate_simulated_data(daily_data, time_seq_out = time_seq_out)
  }

  # Store true parameters as attributes for later recovery evaluation
  attr(output_data, "true_params") <- list(
    mortality_phi = mortality_phi,
    mortality_intercept = mortality_intercept,
    dag_parameters = dag$parameters
  )

  return(output_data)
}


#' Simulate Mortality Using Negative Binomial Distribution
#'
#' Simulates mortality counts Y_it ~ NegBin(mu_it, phi) with:
#' log(mu_it) = log(N_it) + eta_it
#' eta_it = intercept + sum_j beta_j * parent_j(it) + random effects (optional)
#'
#' @param daily_data Daily simulation data with population column
#' @param dag DAG object containing parameters
#' @param mortality_node Name of mortality outcome node
#' @param phi Overdispersion parameter (higher = less overdispersion; phi -> Inf gives Poisson)
#' @param intercept Baseline log-rate intercept
#' @param noise_level Noise scaling factor
#' @param n_days Number of days in simulation
#' @param n_regions Number of regions
#'
#' @return daily_data with mortality column added
#' @keywords internal
simulate_mortality_negbin <- function(daily_data, dag, mortality_node,
                                      phi = 1, intercept = -8,
                                      noise_level = 1, n_days, n_regions) {

  # Get parent nodes and effect sizes for mortality
  parents_df <- dag$parameters[dag$parameters$to == mortality_node, , drop = FALSE]
  parent_names <- as.character(parents_df$from)

  # Compute linear predictor eta_it
  n <- nrow(daily_data)
  eta <- rep(intercept, n)

  if (length(parent_names) > 0 && nrow(parents_df) > 0) {
    # Collect parent data with lags
    parent_data <- daily_data[, parent_names, drop = FALSE]

    # Apply lags to parent data
    if ("lag" %in% names(parents_df)) {
      for (i in seq_len(nrow(parents_df))) {
        parent <- parents_df$from[i]
        lag_days <- ifelse(is.na(parents_df$lag[i]), 0, parents_df$lag[i])
        if (lag_days > 0 && parent %in% names(parent_data)) {
          parent_vec <- daily_data[[parent]]
          parent_mat <- matrix(parent_vec, nrow = n_days, ncol = n_regions, byrow = FALSE)
          if (lag_days >= n_days) {
            parent_mat[] <- 0
          } else {
            parent_mat[(lag_days + 1):n_days, ] <- parent_mat[1:(n_days - lag_days), ]
            parent_mat[1:lag_days, ] <- 0
          }
          parent_data[[as.character(parent)]] <- as.vector(parent_mat)
        }
      }
    }

    # Compute linear combination
    effect_sizes <- parents_df$effect_size
    parent_mat <- as.matrix(parent_data)
    # Scale effect sizes for mortality (effects are typically much smaller on log-rate scale)
    scaled_effects <- effect_sizes * 0.01  # Scale down for reasonable mortality rates
    eta <- eta + as.vector(parent_mat %*% scaled_effects)
  }

  # Add optional random effects for governorate (spatial heterogeneity)
  if (noise_level > 0) {
    gov_re <- stats::rnorm(n_regions, 0, 0.2 * noise_level)
    eta <- eta + gov_re[daily_data$gov_id]
  }

  # Compute mu_it = N_it * exp(eta_it)
  # Use population as offset
  N_it <- pmax(daily_data$population, 1)  # Ensure positive
  mu_it <- N_it * exp(eta)

  # Simulate Y_it ~ NegBin(mu_it, phi)
  # Using size = phi, mu = mu_it parameterization
  # Var(Y) = mu + mu^2/phi
  if (noise_level == 0) {
    # Deterministic: return expected value
    daily_data[[mortality_node]] <- round(mu_it)
  } else {
    # Ensure phi is positive
    phi_adj <- max(phi, 0.01)
    daily_data[[mortality_node]] <- stats::rnbinom(n, size = phi_adj, mu = mu_it)
  }

  return(daily_data)
}


#' Simulate a Single DAG Node
#'
#' Generates simulated values for one node in the DAG, given its parent values and
#' predefined parameters for its distribution and relationship.
#'
#' @param node Name of the node (used for reference in error messages or debugging).
#' @param parents_data A data frame (or list) containing the values of the parent
#' nodes required to simulate this node. Each column (or list element) corresponds
#' to one parent, with each row representing one observation (e.g., one time step for
#' one region). If the node has no parents, this can be an empty data frame or list;
#' in that case, the function will simulate based on baseline parameters.
#' @param dag_params A list of parameters defining how to simulate the node. Typically
#' this would be the element of a DAG definition corresponding to the node. It should
#' include at least:
#'   \itemize{
#'     \item \code{parents}: Character vector of parent node names (can be empty or
#'           \code{NULL} if none).
#'     \item \code{dist}: Distribution type for the node's stochastic simulation. For
#'           example, \code{"normal"} for a continuous variable, \code{"poisson"} for
#'           counts, \code{"negbin"} for negative binomial, \code{"binary"} (or
#'           \code{"bernoulli"}) for 0/1 outcomes, or \code{"none"} for deterministic.
#'     \item Additional parameters specific to the distribution and relationship.
#'   }
#' @param noise_level Numeric factor to scale the stochastic noise. \code{noise_level = 0}
#' produces deterministic output (no randomness), while \code{noise_level = 1} uses the
#' full noise as specified in \code{dag_params}. Default is 1.
#'
#' @return A numeric vector of simulated values for the node.
#'
#' @export
simulate_node <- function(node, parents_data, dag_params, noise_level = 1) {
  # Determine number of observations to simulate
  n <- if (is.data.frame(parents_data)) {
    nrow(parents_data)
  } else if (is.list(parents_data)) {
    # If list of parent vectors, use length of first element (or 1 if empty list)
    if (length(parents_data) == 0) 1 else length(parents_data[[1]])
  } else {
    1
  }

  # Ensure all expected parent columns are present in parents_data (if any parent names given)
  parent_names <- dag_params$parents
  if (!is.null(parent_names)) {
    missing_parents <- setdiff(parent_names, names(parents_data))
    if (length(missing_parents) > 0) {
      stop(sprintf("simulate_node: Parent values for node '%s' not provided: %s",
                   node, paste(missing_parents, collapse = ", ")))
    }
  }

  # Compute the deterministic component of the node's value
  result_det <- NULL
  if (!is.null(dag_params$formula) && is.function(dag_params$formula)) {
    # Use the provided formula function with parent_data values
    # If parents_data is a data frame, convert it to a list of vectors for do.call
    parent_args <- if (is.data.frame(parents_data)) as.list(parents_data) else parents_data
    result_det <- do.call(dag_params$formula, parent_args)
  } else if (!is.null(dag_params$mean)) {
    # Use a provided constant mean/base value if available
    result_det <- dag_params$mean
  } else if (!is.null(dag_params$base)) {
    result_det <- dag_params$base
  }

  # If no deterministic component has been calculated yet, compute from effect sizes if available
  if (is.null(result_det)) {
    if (!is.null(dag_params$effect_size)) {
      # Compute linear combination of parent values using effect_size
      # Ensure parents_data is in matrix form
      if (is.data.frame(parents_data)) {
        parent_mat <- as.matrix(parents_data)
      } else if (is.list(parents_data)) {
        parent_mat <- as.matrix(data.frame(parents_data))
      } else {
        parent_mat <- matrix(parents_data, nrow = n)
      }
      # Make sure effect_size is a numeric vector
      effect_vec <- as.numeric(dag_params$effect_size)
      # Compute the deterministic result as row-wise dot product of parent values and effect sizes
      result_det <- as.vector(parent_mat %*% effect_vec)
    } else {
      # Default deterministic value if no formula, no mean/base, and no effect sizes
      result_det <- 0
    }
  }

  # If result_det is a single value and multiple observations needed, replicate it
  if (length(result_det) == 1 && n > 1) {
    result_det <- rep(result_det, n)
  }
  # If lengths don't match (not 1 and not n), throw an error
  if (!(length(result_det) == 1 || length(result_det) == n)) {
    stop(sprintf("simulate_node: Formula or computation for node '%s' returned length %d, but expected 1 or %d",
                 node, length(result_det), n))
  }

  # Determine distribution and add stochastic noise if applicable
  dist <- if (!is.null(dag_params$dist)) tolower(dag_params$dist) else "none"
  result <- result_det  # start with deterministic component
  if (dist %in% c("normal", "gaussian")) {
    sd_val <- dag_params$sd %||% 1  # use provided sd or default to 1
    result <- result_det + stats::rnorm(n, mean = 0, sd = sd_val * noise_level)
  } else if (dist == "poisson") {
    # Treat deterministic result as rate (lambda) for Poisson
    lambda <- pmax(result_det, 0)  # ensure lambda is non-negative
    result <- stats::rpois(n, lambda = lambda)
  } else if (dist == "negbin") {
    # Negative binomial: treat result_det as mu (mean), use provided phi (size)
    mu <- pmax(result_det, 0.001)
    phi <- dag_params$phi %||% 1
    result <- stats::rnbinom(n, size = phi, mu = mu)
  } else if (dist %in% c("binary", "bernoulli")) {
    # Treat deterministic result as probability for binary outcome
    prob <- pmin(pmax(result_det, 0), 1)  # clamp between 0 and 1
    result <- stats::rbinom(n, size = 1, prob = prob)
  } else if (dist %in% c("none", "deterministic")) {
    # No noise: use deterministic result directly
    result <- result_det
  } else {
    warning(sprintf("simulate_node: Unknown distribution '%s' for node '%s'. Returning deterministic component.",
                    dist, node))
    result <- result_det
  }

  return(result)
}


#' Simulate Missingness in Predictor Variables
#'
#' Generates observation indicators S_itk ~ Bernoulli(pi_itk) for each predictor column k,
#' where logit(pi_itk) depends on the specified missingness regime.
#'
#' @param data Data frame containing simulated data at governorate-month level
#' @param predictors Character vector of predictor column names to apply missingness to
#' @param regime Missingness regime: "MCAR", "MAR", or "MNAR"
#' @param pi_base Baseline observation probability (for MCAR or intercept)
#' @param mar_drivers Character vector of column names driving MAR missingness (e.g., conflict indicators)
#' @param mar_weights Named numeric vector of weights for MAR drivers on logit scale
#' @param mnar_omega MNAR sensitivity parameter: weight on the (unobserved) value itself
#' @param outcome_missing Logical: also apply missingness to outcome variable? Default FALSE.
#' @param outcome_col Name of outcome column if outcome_missing = TRUE
#'
#' @return A list containing:
#'   - data: Original data with values set to NA where not observed
#'   - observed_indicators: Data frame of observation indicators S_itk (1=observed, 0=missing)
#'   - missingness_params: List of parameters used
#'
#' @details
#' Three missingness regimes:
#'   (a) MCAR: S_itk ~ Bernoulli(pi_base) independently
#'   (b) MAR-access: logit(pi_itk) = alpha + w_A*A_it + w_D*D_it + ...
#'   (c) MNAR-severity: logit(pi_itk) = alpha + w_A*A_it + omega*X_itk (depends on unobserved value)
#'
#' @examples
#' # Apply MAR missingness driven by conflict intensity
#' result <- simulate_missingness(
#'   data = sim_data,
#'   predictors = c("Nutritional status", "Burden of endemic infectious diseases"),
#'   regime = "MAR",
#'   pi_base = 0.7,
#'   mar_drivers = c("Exposure to armed attacks or mechanical force of nature"),
#'   mar_weights = c("Exposure to armed attacks or mechanical force of nature" = -0.5)
#' )
#' @export
simulate_missingness <- function(data,
                                 predictors,
                                 regime = c("MCAR", "MAR", "MNAR"),
                                 pi_base = 0.8,
                                 mar_drivers = NULL,
                                 mar_weights = NULL,
                                 mnar_omega = 0,
                                 outcome_missing = FALSE,
                                 outcome_col = NULL) {

  regime <- match.arg(regime)
  n <- nrow(data)

  # Validate predictors exist
  missing_preds <- setdiff(predictors, names(data))
  if (length(missing_preds) > 0) {
    stop("Predictors not found in data: ", paste(missing_preds, collapse = ", "))
  }

  # Initialize observation indicators (1 = observed)
  obs_indicators <- data.frame(matrix(1L, nrow = n, ncol = length(predictors)))
  names(obs_indicators) <- predictors

  # Compute logit(pi) baseline
  alpha <- stats::qlogis(pi_base)  # logit(pi_base)

  for (k in predictors) {
    # Compute logit(pi_itk) based on regime
    logit_pi <- rep(alpha, n)

    if (regime == "MAR" && !is.null(mar_drivers)) {
      # MAR: add contributions from observed drivers
      for (driver in mar_drivers) {
        if (driver %in% names(data)) {
          w <- if (!is.null(mar_weights) && driver %in% names(mar_weights)) {
            mar_weights[[driver]]
          } else {
            -0.3  # default negative weight (more conflict -> less observation)
          }
          # Standardize driver for stable coefficients
          driver_vals <- data[[driver]]
          driver_std <- (driver_vals - mean(driver_vals, na.rm = TRUE)) /
            (sd(driver_vals, na.rm = TRUE) + 1e-8)
          logit_pi <- logit_pi + w * driver_std
        }
      }
    }

    if (regime == "MNAR" && mnar_omega != 0) {
      # MNAR: probability depends on the value itself
      x_vals <- data[[k]]
      x_std <- (x_vals - mean(x_vals, na.rm = TRUE)) / (sd(x_vals, na.rm = TRUE) + 1e-8)
      logit_pi <- logit_pi + mnar_omega * x_std

      # Also add MAR drivers if specified
      if (!is.null(mar_drivers)) {
        for (driver in mar_drivers) {
          if (driver %in% names(data) && driver != k) {
            w <- if (!is.null(mar_weights) && driver %in% names(mar_weights)) {
              mar_weights[[driver]]
            } else {
              -0.3
            }
            driver_vals <- data[[driver]]
            driver_std <- (driver_vals - mean(driver_vals, na.rm = TRUE)) /
              (sd(driver_vals, na.rm = TRUE) + 1e-8)
            logit_pi <- logit_pi + w * driver_std
          }
        }
      }
    }

    # Convert to probability and simulate observation indicator
    pi_itk <- stats::plogis(logit_pi)
    obs_indicators[[k]] <- stats::rbinom(n, size = 1, prob = pi_itk)
  }

  # Apply missingness to data
  data_missing <- data
  for (k in predictors) {
    data_missing[[k]][obs_indicators[[k]] == 0] <- NA
  }

  # Optionally apply to outcome
  if (outcome_missing && !is.null(outcome_col) && outcome_col %in% names(data)) {
    logit_pi_y <- rep(alpha, n)
    if (regime != "MCAR" && !is.null(mar_drivers)) {
      for (driver in mar_drivers) {
        if (driver %in% names(data)) {
          w <- mar_weights[[driver]] %||% -0.3
          driver_vals <- data[[driver]]
          driver_std <- (driver_vals - mean(driver_vals, na.rm = TRUE)) /
            (sd(driver_vals, na.rm = TRUE) + 1e-8)
          logit_pi_y <- logit_pi_y + w * driver_std
        }
      }
    }
    pi_y <- stats::plogis(logit_pi_y)
    obs_y <- stats::rbinom(n, size = 1, prob = pi_y)
    obs_indicators[[outcome_col]] <- obs_y
    data_missing[[outcome_col]][obs_y == 0] <- NA
  }

  return(list(
    data = data_missing,
    observed_indicators = obs_indicators,
    missingness_params = list(
      regime = regime,
      pi_base = pi_base,
      mar_drivers = mar_drivers,
      mar_weights = mar_weights,
      mnar_omega = mnar_omega
    )
  ))
}


#' Apply Missingness Mask to Simulated Data
#'
#' Forces the simulated dataset to follow an external missingness pattern (e.g., from Gaza data).
#'
#' @param data Data frame of simulated data at governorate-month level
#' @param mask Data frame specifying missingness pattern with columns:
#'   - gov: governorate identifier (must match data)
#'   - month: month identifier or date
#'   - var: variable name
#'   - observed: 1 if observed, 0 if missing
#' @param gov_col Name of governorate column in data (default: auto-detected)
#' @param date_col Name of date/month column in data (default: "date")
#' @param mode How to apply mask:
#'   - "force": directly apply mask (set to NA where observed=0)
#'   - "calibrate": estimate missingness parameters to match mask proportions, then simulate
#'
#' @return Data frame with missingness applied according to mask
#'
#' @examples
#' # Create a mask matching Gaza reporting patterns
#' gaza_mask <- data.frame(
#'   gov = rep(c("North Gaza", "Gaza City"), each = 6),
#'   month = rep(1:6, 2),
#'   var = "Nutritional status",
#'   observed = c(1,1,0,0,0,0, 1,1,1,0,0,0)  # reporting stops mid-crisis
#' )
#' data_masked <- apply_mask(sim_data, gaza_mask, mode = "force")
#'
#' @export
apply_mask <- function(data,
                       mask,
                       gov_col = NULL,
                       date_col = "date",
                       mode = c("force", "calibrate")) {

  mode <- match.arg(mode)

  # Auto-detect governorate column if not specified
  if (is.null(gov_col)) {
    # Look for common names
    candidates <- c("gov", "governorate", "region", "district", "gov_id")
    gov_col <- intersect(candidates, names(data))[1]
    if (is.na(gov_col)) {
      # Use last non-date, non-numeric column
      non_num_cols <- names(data)[!sapply(data, is.numeric)]
      gov_col <- setdiff(non_num_cols, c(date_col, "period_start"))[1]
    }
  }

  # Validate mask structure
  required_mask_cols <- c("gov", "var", "observed")
  if (!all(required_mask_cols %in% names(mask))) {
    stop("Mask must contain columns: gov, var, observed (and optionally month)")
  }

  # Add month identifier to data if not present
  if (!"month" %in% names(data) && date_col %in% names(data)) {
    data$month <- as.numeric(format(as.Date(data[[date_col]]), "%m")) +
      (as.numeric(format(as.Date(data[[date_col]]), "%Y")) -
         min(as.numeric(format(as.Date(data[[date_col]]), "%Y")))) * 12
  }

  # Standardize governorate names in mask to match data
  data_govs <- unique(data[[gov_col]])
  mask_govs <- unique(mask$gov)

  if (mode == "force") {
    # Directly apply mask
    data_masked <- data

    for (i in seq_len(nrow(mask))) {
      gov_match <- mask$gov[i]
      var_name <- mask$var[i]
      obs_val <- mask$observed[i]

      if (var_name %in% names(data_masked)) {
        # Find matching rows
        if ("month" %in% names(mask)) {
          month_match <- mask$month[i]
          row_idx <- which(data_masked[[gov_col]] == gov_match &
                             data_masked$month == month_match)
        } else {
          row_idx <- which(data_masked[[gov_col]] == gov_match)
        }

        # Apply missingness
        if (length(row_idx) > 0 && obs_val == 0) {
          data_masked[[var_name]][row_idx] <- NA
        }
      }
    }

    return(data_masked)

  } else if (mode == "calibrate") {
    # Calibrate missingness parameters to match mask proportions
    # Then simulate missingness with those parameters

    # Compute target missingness proportion by variable
    target_props <- stats::aggregate(observed ~ var, data = mask, FUN = mean)

    # For each variable, find pi_base that matches observed proportion
    data_masked <- data
    for (i in seq_len(nrow(target_props))) {
      var_name <- target_props$var[i]
      target_pi <- target_props$observed[i]

      if (var_name %in% names(data_masked)) {
        # Simulate with calibrated pi_base
        obs_ind <- stats::rbinom(nrow(data_masked), size = 1, prob = target_pi)
        data_masked[[var_name]][obs_ind == 0] <- NA
      }
    }

    return(data_masked)
  }
}


#' Simulate Population Dynamics with IDP Flows
#'
#' Simulates changes in population over time for each region, driven by daily inward
#' and outward movements of internally displaced persons (IDPs). This function updates
#' the population for each day based on the previous day's population and net IDP flow.
#'
#' @param daily_data A data frame of daily simulation data by region. It must contain
#' at least a \code{date} column and one or more region identifier columns (e.g.,
#' region names or codes). It should also contain columns for IDP inflows and outflows
#' (or a single net flow column), unless no displacement is being modeled.
#' @param initial_population Numeric value or vector giving the initial population at
#' the start of the simulation for each lowest-level region. If a single number is
#' provided, all regions start with that population. If a vector is provided, it
#' should be either:
#'   \itemize{
#'     \item A named vector, where names correspond to the region identifiers in
#'           \code{daily_data} (matching the lowest-level region names or IDs).
#'     \item An unnamed vector with length equal to the number of lowest-level regions
#'           (in the same order as the regions appear when grouping \code{daily_data}
#'           by region).
#'   }
#' @param idp_in_col Name of the column in \code{daily_data} representing IDP inflows
#' to a region. Defaults to \code{"idp_in"}. If this column is not found, the function
#' will interpret the data differently (see Details).
#' @param idp_out_col Name of the column in \code{daily_data} representing IDP outflows
#' from a region. Defaults to \code{"idp_out"}. If this column is not found or set to
#' \code{NULL}, and an inflow column is present, that inflow is assumed to represent
#' net population change (positive for influx, negative for outflux).
#'
#' @return The \code{daily_data} data frame with an additional column \code{population},
#' giving the simulated population for each region on each day.
#'
#' @export
simulate_population_dynamics <- function(daily_data, initial_population,
                                         idp_in_col = "idp_in", idp_out_col = "idp_out") {
  # Ensure required columns exist
  if (!("date" %in% names(daily_data))) {
    stop("daily_data must contain a 'date' column.")
  }
  # Identify region identifier columns (all non-date, non-flow, non-population columns)
  region_cols <- names(daily_data)[!names(daily_data) %in% c("date", idp_in_col, idp_out_col,
                                                             "population", "gov_id", "time_id")]
  if (length(region_cols) < 1) {
    stop("daily_data must contain at least one region identifier column.")
  }
  # Use the lowest-level region column (assumed to be the last in hierarchical order)
  region_id_col <- region_cols[length(region_cols)]

  # Prepare initial population vector by region
  region_ids <- unique(daily_data[[region_id_col]])
  if (length(initial_population) == 1) {
    init_pop_vec <- setNames(rep(initial_population, length(region_ids)), region_ids)
  } else if (!is.null(names(initial_population))) {
    init_pop_vec <- initial_population
    if (!all(region_ids %in% names(init_pop_vec))) {
      stop("Names of initial_population vector do not match region identifiers in data.")
    }
  } else if (length(initial_population) == length(region_ids)) {
    init_pop_vec <- setNames(initial_population, region_ids)
  } else {
    stop("initial_population must be a single number or a vector (named or with length equal to number of regions).")
  }

  # Sort data by region and date to ensure sequential updates per region
  daily_data <- daily_data[order(daily_data[[region_id_col]], daily_data$date), ]
  rownames(daily_data) <- NULL  # reset row indices

  # Determine how to interpret IDP flows
  have_in  <- idp_in_col %in% names(daily_data)
  have_out <- !is.null(idp_out_col) && idp_out_col %in% names(daily_data)
  net_flow_mode <- FALSE
  if (!have_in && have_out) {
    # If only outflows are present, treat outflows as negative inflows (net out-migration)
    idp_in_col <- idp_out_col
    have_in <- TRUE
    net_flow_mode <- TRUE
  } else if (have_in && !have_out) {
    # Only inflow present (and no separate outflow) -> interpret inflow as net change
    net_flow_mode <- TRUE
  } else if (!have_in && !have_out) {
    # No IDP flow data: assign constant population for all days
    daily_data$population <- init_pop_vec[ daily_data[[region_id_col]] ]
    return(daily_data)
  }

  # Initialize population column and set starting populations
  daily_data$population <- NA_real_
  # Track current population per region in a named vector
  current_pop <- init_pop_vec

  # Iterate through each record in chronological order
  current_region <- daily_data[[region_id_col]][1]
  for (i in seq_len(nrow(daily_data))) {
    region <- daily_data[[region_id_col]][i]
    # If we move to a new region block in the sorted data, reset current_region
    if (region != current_region) {
      current_region <- region
      # (current_pop for this region is already stored or initialized)
    }
    if (i == 1 || daily_data[[region_id_col]][i - 1] != region) {
      # First day for this region: initialize population
      current_pop[region] <- init_pop_vec[region]
      daily_data$population[i] <- current_pop[region]
    } else {
      # Subsequent days: update population based on flows
      if (net_flow_mode) {
        # Single net flow column: treat idp_in_col as net change (positive = influx, negative = outflux)
        net_change <- daily_data[i, idp_in_col]
        current_pop[region] <- current_pop[region] + as.numeric(net_change)
      } else {
        # Separate inflow and outflow columns
        inflow  <- daily_data[i, idp_in_col]
        outflow <- daily_data[i, idp_out_col]
        current_pop[region] <- current_pop[region] + as.numeric(inflow) - as.numeric(outflow)
      }
      daily_data$population[i] <- current_pop[region]
    }
  }

  return(daily_data)
}


#' Aggregate Daily Simulation Data to Periods
#'
#' Aggregates detailed daily simulation outputs into the specified temporal resolution
#' (e.g., weekly, monthly) by region. This function sums flow variables over each
#' period and carries forward stock variables like population as of the end of each period.
#'
#' @param daily_data A data frame of daily data by region, containing at least a
#' \code{date} column and various numeric indicator columns to aggregate.
#' @param time_seq_out (Optional) A Date vector of period start dates to aggregate by.
#' If provided, these define the boundaries for aggregation (with an additional boundary
#' at the end of the simulation). If not provided, \code{resolution} must be given.
#' @param resolution (Optional) The target time resolution for aggregation, as a string
#' (e.g., \code{"week"}, \code{"month"}). Ignored if \code{time_seq_out} is supplied.
#' If \code{time_seq_out} is not given, the function will derive period breaks from
#' the minimum date in \code{daily_data} at the specified resolution.
#' @param start_date (Optional) If \code{time_seq_out} is not given, an anchor start date
#' for the periods (defaults to the first date in \code{daily_data}). This date will
#' be used as the start of the first period.
#'
#' @return A data frame aggregated by region and time period. It contains one row per
#' region per period. The output includes all region identifier columns, a \code{date}
#' column indicating the start date of each period, and all numeric indicator columns
#' aggregated. Flow metrics (e.g., daily counts) are summed over the period, while the
#' \code{population} column (if present) is taken as the last day's value in the period.
#'
#' @export
aggregate_simulated_data <- function(daily_data, time_seq_out = NULL, resolution = NULL, start_date = NULL) {
  if (!("date" %in% names(daily_data))) {
    stop("daily_data must contain a 'date' column.")
  }
  # Determine period breaks based on provided output sequence or resolution
  if (is.null(time_seq_out)) {
    if (is.null(resolution)) {
      stop("Either time_seq_out or resolution must be provided.")
    }
    # Determine anchor start date for periods
    anchor <- if (!is.null(start_date)) as.Date(start_date) else min(as.Date(daily_data$date))
    # Generate sequence of period start dates from anchor to cover all data
    max_date <- max(as.Date(daily_data$date))
    time_seq_out <- seq.Date(from = anchor, to = max_date, by = resolution)
    # Ensure the sequence covers the final partial period (add one more period start if needed)
    if (tail(time_seq_out, 1) < max_date) {
      time_seq_out <- c(time_seq_out, seq.Date(from = tail(time_seq_out, 1), by = resolution, length.out = 2)[2])
    }
  } else {
    # Ensure provided time_seq_out is sorted and in Date class
    time_seq_out <- sort(as.Date(time_seq_out))
  }
  # Define period breakpoints (end of each period is just before the next period start)
  period_breaks <- c(time_seq_out, max(as.Date(daily_data$date)) + 1)
  # Assign each daily date to a period start label
  daily_data$period_start <- cut(as.Date(daily_data$date), breaks = period_breaks, right = FALSE, labels = time_seq_out)
  daily_data$period_start <- as.Date(as.character(daily_data$period_start))

  # Determine grouping columns (all non-numeric columns, replacing 'date' with 'period_start')
  grouping_cols <- setdiff(names(daily_data)[!sapply(daily_data, is.numeric)], "date")
  grouping_cols <- c(grouping_cols, "period_start")
  grouping_cols <- intersect(grouping_cols, names(daily_data))

  # Aggregate by region identifiers and period_start
  agg_data <- daily_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
    dplyr::summarise(
      # Sum all numeric columns, except treat population specially
      dplyr::across(dplyr::where(is.numeric) & !dplyr::matches("^population$"), ~ sum(.x, na.rm = TRUE)),
      population = if ("population" %in% names(daily_data)) dplyr::last(.data$population) else NULL,
      .groups = "drop"
    )

  # Rename period_start to date for output (indicating period start date)
  names(agg_data)[names(agg_data) == "period_start"] <- "date"
  return(as.data.frame(agg_data))
}
