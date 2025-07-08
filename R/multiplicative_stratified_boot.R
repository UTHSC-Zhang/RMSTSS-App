#' @title Internal Factory for RMST Simulation Functions (Bootstrap)
#' @description This is a non-exported helper function that sets up and returns
#'   another function responsible for running the core bootstrap simulation logic. This
#'   factory pattern helps reduce code duplication between the public-facing
#'   power and sample size functions. It prepares the model formula and pilot data
#'   and returns the `run_power_sim` function which is then called repeatedly.
#' @return A function that takes `n_per_stratum` and runs the simulation.
#' @keywords internal
.get_internal_simulation_runner <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                                            linear_terms, tau, alpha, n_sim, parallel.cores) {

   # --- Data Preparation & Model Formula ---
   core_vars <- c(time_var, status_var, arm_var, strata_var)
   all_vars <- c(core_vars, linear_terms)
   cleaned_pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
   pilot_groups <- split(cleaned_pilot_data, cleaned_pilot_data[[strata_var]])

   # --- Model Formula Construction ---
   # Model with stratum-specific intercepts and stratum-specific treatment effects
   all_terms <- c(strata_var, paste0(strata_var, ":", arm_var), linear_terms)
   model_rhs <- paste(all_terms[!sapply(all_terms, is.null)], collapse = " + ")
   model_formula <- as.formula(paste("log_pseudo_obs ~", model_rhs))
   message("Model: log(pseudo_obs) ~ ", model_rhs)
   # Robust pattern to find interaction terms like 'stratumA:arm1', 'stratumB:arm1'
   test_term_pattern <- paste0(":", arm_var, "1$")

   # --- Helper to calculate jackknife pseudo-observations ---
   get_pseudo_obs <- function(time, status, tau_val) {
      n <- length(time)
      if (n == 0) return(numeric(0))
      km_fit_full <- survival::survfit(survival::Surv(time, status) ~ 1)
      km_step_full <- stats::stepfun(km_fit_full$time, c(1, km_fit_full$surv))
      rmst_full <- tryCatch(stats::integrate(km_step_full, 0, tau_val, subdivisions = 2000, stop.on.error = FALSE)$value, error = function(e) 0)

      rmst_loo <- vapply(seq_len(n), function(i) {
         if (n > 1) {
            km_fit_loo <- survival::survfit(survival::Surv(time[-i], status[-i]) ~ 1)
            km_step_loo <- stats::stepfun(km_fit_loo$time, c(1, km_fit_loo$surv))
            tryCatch(stats::integrate(km_step_loo, 0, tau_val, subdivisions = 2000, stop.on.error = FALSE)$value, error = function(e) 0)
         } else { 0 }
      }, FUN.VALUE = numeric(1))

      return(n * rmst_full - (n - 1) * rmst_loo)
   }

   # --- The actual simulation function that will be returned ---
   run_power_sim <- function(n_per_stratum) {
      # This function is designed to be called with future_lapply
      single_iteration <- function(iter) {
         p_val <- NA_real_
         estimate <- NA_real_

         boot_list <- lapply(pilot_groups, function(group_df) {
            if(nrow(group_df) > 0) {
               group_df[sample(seq_len(nrow(group_df)), size = n_per_stratum, replace = TRUE), ]
            } else { NULL }
         })
         boot_data <- do.call(rbind, boot_list)
         if(is.null(boot_data) || nrow(boot_data) == 0) return(list(p_value=p_val, estimate=estimate))

         # Explicitly set arm to factor for consistent coefficient naming
         boot_data[[arm_var]] <- factor(boot_data[[arm_var]], levels = c(0, 1))

         pseudo_obs_list <- by(boot_data, boot_data[[strata_var]], function(sub_data) {
            sub_data$pseudo_obs <- get_pseudo_obs(sub_data[[time_var]], sub_data[[status_var]], tau)
            sub_data
         })
         boot_data <- do.call(rbind, pseudo_obs_list)

         boot_data <- boot_data[which(boot_data$pseudo_obs > 0), ]

         if (nrow(boot_data) > (length(all_vars) + length(pilot_groups))) {
            boot_data$log_pseudo_obs <- log(boot_data$pseudo_obs)

            # Check for sufficient data across all arms/strata combinations
            if (all(table(boot_data[[arm_var]], boot_data[[strata_var]]) > 0)) {
               fit <- tryCatch(stats::lm(model_formula, data = boot_data), error = function(e) NULL)
               if (!is.null(fit)) {
                  sfit <- summary(fit)
                  coeffs <- sfit$coefficients
                  matching_rows <- grep(test_term_pattern, rownames(coeffs), value = TRUE)

                  if (length(matching_rows) > 0) {
                     # For multiple strata, take the minimum p-value
                     p_val <- min(coeffs[matching_rows, "Pr(>|t|)"], na.rm = TRUE)
                     estimate <- mean(coeffs[matching_rows, "Estimate"], na.rm = TRUE)
                  }
               }
            }
         }
         return(list(p_value = p_val, estimate = estimate))
      }

      if (parallel.cores > 1) {
         old_plan <- future::plan(future::multisession, workers = parallel.cores)
         on.exit(future::plan(old_plan), add = TRUE)
         all_sim_outputs <- future.apply::future_lapply(1:n_sim, single_iteration, future.seed = TRUE)
      } else {
         all_sim_outputs <- lapply(1:n_sim, single_iteration)
      }

      p_values <- vapply(all_sim_outputs, `[[`, "p_value", FUN.VALUE = numeric(1))
      estimates <- vapply(all_sim_outputs, `[[`, "estimate", FUN.VALUE = numeric(1))

      power <- mean(p_values < alpha, na.rm = TRUE)
      valid_estimates <- estimates[is.finite(estimates)]
      # Return estimates as RMST Ratio (exp(beta))
      return(list(power = power, estimates = exp(valid_estimates)))
   }

   # Return the fully configured simulation function
   return(run_power_sim)
}


# Power Calculation -------------------------------------------------------

#' @title Analyze Power for a Multiplicative Stratified RMST Model via Simulation
#' @description Performs power analysis based on a multiplicative model for RMST
#'   for stratified trials, using a bootstrap simulation approach.
#'
#' @details This function estimates power through bootstrap simulation, resampling
#'   within each stratum defined by `strata_var`. In each of the `n_sim` iterations:
#'   1. Jackknife pseudo-observations for the RMST are calculated.
#'   2. A log-linear model (`stats::lm`) is fitted to the `log(pseudo_obs)`.
#'      This models the multiplicative relationship on the original RMST scale,
#'      i.e., \eqn{\mu_{ij} = \mu_{0j} \exp\{\beta'Z_i\}}. The model formula
#'      includes stratum-specific intercepts and interactions with the treatment arm.
#'   3. The p-value for the treatment effect is extracted from the model summary.
#'
#'   Power is determined as the proportion of simulations where the p-value is
#'   less than `alpha`.
#'
#' @note `status_var` should be `1`/`0`. `arm_var` should be `1`/`0`. `strata_var`
#'   is a mandatory argument.
#'
#' @param pilot_data A `data.frame` with pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable.
#' @param arm_var A character string for the treatment arm variable.
#' @param strata_var A character string for the stratification variable.
#' @param sample_sizes A numeric vector of sample sizes *per stratum* to calculate power for.
#' @param linear_terms Optional character vector of covariates for the model.
#' @param tau The numeric truncation time for RMST.
#' @param n_sim Number of bootstrap simulations.
#' @param alpha The significance level.
#' @param parallel.cores Number of cores for parallel processing.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` of sample sizes and corresponding powers.}
#' \item{results_plot}{A `ggplot` object visualizing the power curve.}
#' \item{results_summary}{A `data.frame` with the estimated treatment effect (RMST Ratio).}
#'
#' @importFrom survival survfit Surv
#' @importFrom stats lm as.formula complete.cases integrate na.omit quantile stepfun sd
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs theme_minimal ylim
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom knitr kable
#' @export
#' @examples
#' \dontrun{
#' pilot_df_strat <- data.frame(
#'  time = rexp(120, 0.15),
#'  status = rbinom(120, 1, 0.6),
#'  arm = rep(0:1, each = 60),
#'  region = factor(rep(c("A", "B", "C"), each = 40))
#' )
#' pilot_df_strat$time[pilot_df_strat$arm == 1] <- pilot_df_strat$time[pilot_df_strat$arm == 1] * 1.4
#'
#' power_results <- MS.power.boot(
#'  pilot_data = pilot_df_strat,
#'  time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
#'  sample_sizes = c(50, 75),
#'  tau = 10,
#'  n_sim = 100 # Low n_sim for example
#' )
#' print(power_results$results_data)
#' }
MS.power.boot <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                          sample_sizes, linear_terms = NULL, tau, n_sim = 1000,
                          alpha = 0.05, parallel.cores = 1) {
   # --- Input Validation ---
   if (is.null(sample_sizes)) stop("You must provide the 'sample_sizes' argument.")
   if (parallel.cores > 1) {
      if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
         stop("Packages 'future' and 'future.apply' are required for parallel processing.")
      }
   }

   # --- Setup & Internal Functions ---
   sim_function <- .get_internal_simulation_runner(
      pilot_data = pilot_data,
      time_var = time_var, status_var = status_var, arm_var = arm_var,
      strata_var = strata_var, linear_terms = linear_terms, tau = tau,
      alpha = alpha, n_sim = n_sim, parallel.cores = parallel.cores
   )

   # --- Main Execution ---
   cat("--- Calculating Power (Method: Multiplicative Stratified RMST Model) ---\n")
   start_time <- proc.time()

   sim_outputs <- lapply(sample_sizes, function(n) {
      cat(paste0("Simulating for n = ", n, "/stratum...\n"))
      sim_function(n)
   })

   power_values <- sapply(sim_outputs, `[[`, "power")
   results_df <- data.frame(N_per_Stratum = sample_sizes, Power = power_values)

   # --- Summarize Estimated Effect ---
   results_summary <- NULL
   # Use estimates from the largest sample size simulation for the summary
   if (length(sim_outputs) > 0) {
      est <- na.omit(sim_outputs[[length(sim_outputs)]]$estimates)
      if (length(est) > 1) {
         results_summary <- data.frame(
            Statistic = c("Mean RMST Ratio", "95% CI Lower", "95% CI Upper"),
            Value = c(mean(est), stats::quantile(est, 0.025), stats::quantile(est, 0.975))
         )
      }
   }

   # --- Create Plot ---
   p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Stratum, y = Power)) +
      ggplot2::geom_line(color = "#E69F00", linewidth = 1) +
      ggplot2::geom_point(color = "#E69F00", size = 3) +
      ggplot2::labs(
         title = "Power Curve: Multiplicative Stratified RMST Model",
         x = "Sample Size Per Stratum", y = "Estimated Power"
      ) +
      ggplot2::ylim(0, 1) + ggplot2::theme_minimal()

   # --- Final Output ---
   end_time <- proc.time()
   elapsed_time <- round((end_time - start_time)["elapsed"], 2)
   message(paste("Total simulation time:", elapsed_time, "seconds"))
   cat("\n--- Simulation Summary ---\n")
   if (!is.null(results_summary)) {
      print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Ratio)"))
   } else {
      cat("No valid estimates were generated to create a summary.\n")
   }

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}


# Sample Size Search ------------------------------------------------------

#' @title Estimate Sample Size for a Multiplicative Stratified RMST Model via Simulation
#' @description Performs sample size estimation based on a multiplicative model for RMST
#'   for stratified trials, using iterative bootstrap simulations.
#'
#' @details
#' This function iteratively searches for the sample size required to
#' achieve a `target_power`. At each step of the search, it runs a full bootstrap simulation
#' (as described in `MS.power.boot`) to estimate the power for the
#' current sample size. The search proceeds until the target power is met or
#' other stopping criteria are satisfied. This process can be very computationally
#' intensive.
#'
#' @note `status_var` should be `1`/`0`. `arm_var` should be `1`/`0`. `strata_var`
#'   is a mandatory argument.
#'
#' @param pilot_data A `data.frame` with pilot study data.
#' @param time_var A character string for the time-to-event variable.
#' @param status_var A character string for the event status variable.
#' @param arm_var A character string for the treatment arm variable.
#' @param strata_var A character string for the stratification variable.
#' @param target_power A single numeric value for the target power (e.g., 0.80).
#' @param linear_terms Optional vector of covariates for the model.
#' @param tau The numeric truncation time for RMST.
#' @param n_sim Number of bootstrap simulations per search step.
#' @param alpha The significance level.
#' @param parallel.cores Number of cores for parallel processing.
#' @param patience Number of consecutive non-improving steps in the search before terminating.
#' @param n_start Starting sample size per stratum for the search.
#' @param n_step Increment for the sample size search.
#' @param max_n_per_arm Maximum sample size per stratum to try.
#'
#' @return A `list` containing:
#' \item{results_data}{A `data.frame` with the target power and required N.}
#' \item{results_plot}{A `ggplot` object showing the search path.}
#' \item{results_summary}{A `data.frame` with the estimated treatment effect.}
#'
#' @importFrom survival survfit Surv
#' @importFrom stats lm as.formula complete.cases integrate na.omit quantile stepfun sd
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline geom_hline labs theme_minimal
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom knitr kable
#' @export
#' @examples
#' \dontrun{
#' pilot_df_strat_effect <- data.frame(
#'  time = c(rexp(60, 0.15), rexp(60, 0.08)), # Effect
#'  status = rbinom(120, 1, 0.7),
#'  arm = rep(0:1, each = 60),
#'  region = factor(rep(c("A", "B", "C"), each = 40))
#' )
#' ss_results <- MS.ss.boot(
#'  pilot_data = pilot_df_strat_effect,
#'  time_var = "time", status_var = "status", arm_var = "arm", strata_var = "region",
#'  target_power = 0.80, tau = 10,
#'  n_sim = 100, # Low n_sim for example
#'  n_start = 100,
#'  n_step = 50, patience = 2
#' )
#' print(ss_results$results_data)
#' }
MS.ss.boot <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                       target_power, linear_terms = NULL, tau, n_sim = 1000,
                       alpha = 0.05, parallel.cores = 1, patience = 5,
                       n_start = 50, n_step = 25, max_n_per_arm = 2000) {
   # --- Input Validation ---
   if (is.null(target_power) || length(target_power) != 1 || !is.numeric(target_power)) {
      stop("You must provide a single numeric value for 'target_power'.")
   }
   if (parallel.cores > 1) {
      if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
         stop("Packages 'future' and 'future.apply' are required for parallel processing.")
      }
   }

   # --- Setup & Internal Functions ---
   sim_function <- .get_internal_simulation_runner(
      pilot_data = pilot_data,
      time_var = time_var, status_var = status_var, arm_var = arm_var,
      strata_var = strata_var, linear_terms = linear_terms, tau = tau,
      alpha = alpha, n_sim = n_sim, parallel.cores = parallel.cores
   )

   # --- Main Execution ---
   cat("--- Searching for Sample Size (Method: Multiplicative Stratified RMST Model) ---\n")
   start_time <- proc.time()

   search_path <- list()
   final_n <- NA_integer_
   current_n <- n_start
   max_power_so_far <- -1
   n_at_max_power <- n_start # Correctly track N at best power
   stagnation_counter <- 0
   best_sim_output <- NULL

   while (current_n <= max_n_per_arm) {
      cat(paste0("  N = ", current_n, "/stratum, Calculating Power..."))
      sim_output <- sim_function(current_n)
      calculated_power <- if(is.finite(sim_output$power)) sim_output$power else 0

      search_path[[as.character(current_n)]] <- calculated_power
      cat(paste0(" Power = ", round(calculated_power, 3), "\n"))

      if (calculated_power >= target_power) {
         message("Success: Target power reached at N = ", current_n, "/stratum.")
         best_sim_output <- sim_output
         final_n <- current_n
         break
      }

      if (calculated_power > max_power_so_far) {
         max_power_so_far <- calculated_power
         n_at_max_power <- current_n # Correctly update N
         best_sim_output <- sim_output
         stagnation_counter <- 0
      } else {
         stagnation_counter <- stagnation_counter + 1
      }

      if (stagnation_counter >= patience) {
         warning(paste("Search terminated due to stagnation. Best N =", n_at_max_power,
                       "with power =", round(max_power_so_far, 3)), call. = FALSE)
         final_n <- n_at_max_power # Correctly assign the best N found
         break
      }
      current_n <- current_n + n_step
   }

   if (is.na(final_n)) {
      warning(paste("Target power", target_power, "not achieved by max N of", max_n_per_arm), call. = FALSE)
      final_n <- n_at_max_power
   }

   # --- Finalize Summary and Results ---
   results_summary <- NULL
   if (!is.null(best_sim_output)) {
      est <- na.omit(best_sim_output$estimates)
      if (length(est) > 1) {
         results_summary <- data.frame(
            Statistic = c("Mean RMST Ratio", "95% CI Lower", "95% CI Upper"),
            Value = c(mean(est), stats::quantile(est, 0.025), stats::quantile(est, 0.975))
         )
      }
   }

   results_df <- data.frame(Target_Power = target_power, Required_N_per_Stratum = final_n)
   search_path_df <- data.frame(N_per_Stratum = as.integer(names(search_path)), Power = unlist(search_path))

   # --- Create Plot ---
   p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Stratum, y = Power)) +
      ggplot2::geom_line(color = "#009E73", linewidth = 1) +
      ggplot2::geom_point(color = "#009E73", size = 3) +
      ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
      ggplot2::labs(
         title = "Sample Size Search Path: Multiplicative Stratified RMST Model",
         x = "Sample Size Per Stratum", y = "Calculated Power"
      ) +
      ggplot2::theme_minimal()

   # --- Final Output ---
   end_time <- proc.time()
   elapsed_time <- round((end_time - start_time)["elapsed"], 2)
   message(paste("Total simulation time:", elapsed_time, "seconds"))
   cat("\n--- Simulation Summary ---\n")
   if (!is.null(results_summary)) {
      print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Ratio)"))
   } else {
      cat("No valid estimates were generated to create a summary.\n")
   }

   return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}
