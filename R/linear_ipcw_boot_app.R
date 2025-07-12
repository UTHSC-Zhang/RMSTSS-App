
#' @title Internal Factory for Linear IPCW RMST Simulation
#' @description This internal function prepares and returns a configured simulation function.
#' @return A function that takes `n_per_arm` and runs the simulation.
#' @keywords internal
#' @export
.get_linear_ipcw_simulation_runner <- function(pilot_data, time_var, status_var, arm_var,
                                               linear_terms, L, alpha, n_sim, parallel.cores) {
  
  # --- 1. One-Time Setup ---
  core_vars <- c(time_var, status_var, arm_var)
  all_vars <- c(core_vars, linear_terms)
  cleaned_pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
  pilot_groups <- split(cleaned_pilot_data, cleaned_pilot_data[[arm_var]])
  
  model_rhs <- paste(c(arm_var, linear_terms), collapse = " + ")
  model_formula <- as.formula(paste("Y_rmst ~", model_rhs))
  message("Model: Y_rmst ~ ", model_rhs)
  
  # --- 2. The Returned Simulation Function ---
  run_simulation <- function(n_per_arm) {
    p_values <- rep(NA_real_, n_sim)
    estimates <- rep(NA_real_, n_sim)
    std_errors <- rep(NA_real_, n_sim)
    
    for (j in seq_len(n_sim)) {
      boot_list <- lapply(pilot_groups, function(df) df[sample(seq_len(nrow(df)), size = n_per_arm, replace = TRUE), ])
      boot_data <- do.call(rbind, boot_list)
      boot_data[[arm_var]] <- factor(boot_data[[arm_var]], levels = c(0, 1))
      
      is_censored <- boot_data[[status_var]] == 0
      cens_fit <- tryCatch(survival::survfit(Surv(boot_data[[time_var]], is_censored) ~ 1), error = function(e) NULL)
      if (is.null(cens_fit)) next
      surv_summary <- tryCatch(summary(cens_fit, times = pmin(boot_data[[time_var]], L), extend = TRUE), error = function(e) NULL)
      if (is.null(surv_summary)) next
      
      weights <- 1 / surv_summary$surv
      finite_weights <- weights[is.finite(weights)]
      if (length(finite_weights) > 0) {
        weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
        weights[weights > weight_cap] <- weight_cap
      }
      weights[!is.finite(weights)] <- NA
      
      boot_data$Y_rmst <- pmin(boot_data[[time_var]], L)
      fit_data <- boot_data[boot_data[[status_var]] == 1 & is.finite(weights), ]
      fit_weights <- weights[boot_data[[status_var]] == 1 & is.finite(weights)]
      
      if (nrow(fit_data) > (length(all_vars) + 1)) {
        fit <- tryCatch(lm(model_formula, data = fit_data, weights = fit_weights), error = function(e) NULL)
        if (!is.null(fit)) {
          sfit <- summary(fit)
          test_term_pattern <- paste0("^", arm_var, "1$")
          test_term <- grep(test_term_pattern, rownames(sfit$coefficients), value = TRUE)[1]
          if (!is.na(test_term)) {
            p_values[j] <- tryCatch(sfit$coefficients[test_term, "Pr(>|t|)"], error = function(e) NA)
            estimates[j] <- tryCatch(sfit$coefficients[test_term, "Estimate"], error = function(e) NA)
            std_errors[j] <- tryCatch(sfit$coefficients[test_term, "Std. Error"], error = function(e) NA)
          }
        }
      }
    } # End of simulation loop
    
    return(list(
      power = mean(p_values < alpha, na.rm = TRUE),
      estimates = estimates,
      std_errors = std_errors
    ))
  }
  
  # Return the fully configured simulation function
  return(run_simulation)
}


# Power Calculation -------------------------------------------------------
#' @title Analyze Power for a Linear RMST Model via Simulation
#' @description Performs a power analysis for given sample sizes based on the direct
#'   linear regression model for RMST, using a bootstrap simulation approach.
#' @inheritParams .get_linear_ipcw_simulation_runner
#' @param sample_sizes A numeric vector of sample sizes *per arm* to calculate power for.
#' @return A `list` containing results.
#' @keywords internal
#' @export
linear.power.boot.app <- function(pilot_data, time_var, status_var, arm_var,
                              sample_sizes, linear_terms = NULL, L, n_sim = 1000, alpha = 0.05, parallel.cores) {
  
  start_time <- proc.time()
  if (is.null(sample_sizes)) stop("You must provide a numeric vector for 'sample_sizes'.")
  
  # Get the configured simulation runner
  simulation_runner <- .get_linear_ipcw_simulation_runner(pilot_data, time_var, status_var, arm_var, linear_terms, L, alpha, n_sim)
  
  cat("--- Calculating Power (Method: Linear RMST with IPCW) ---\n")
  all_sim_outputs <- vector("list", length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    n_per_arm <- sample_sizes[i]
    cat("Simulating for n =", n_per_arm, "per arm...\n")
    all_sim_outputs[[i]] <- simulation_runner(n_per_arm)
  }
  
  power_values <- sapply(all_sim_outputs, `[[`, "power")
  results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)
  
  # Summarize results from the largest sample size simulation
  best_idx <- which.max(sample_sizes)
  est <- na.omit(all_sim_outputs[[best_idx]]$estimates)
  se  <- na.omit(all_sim_outputs[[best_idx]]$std_errors)
  results_summary <- NULL
  if (length(est) > 1) {
    results_summary <- data.frame(
      Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
      Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * sd(est), mean(est) + 1.96 * sd(est))
    )
  }
  
  p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
    ggplot2::geom_line(color = "#D55E00", linewidth = 1) +
    ggplot2::geom_point(color = "#D55E00", size = 3) +
    ggplot2::labs(title = "Power Curve: Linear IPCW RMST Model",
                  x = "Sample Size Per Arm", y = "Estimated Power") +
    ggplot2::ylim(0, 1) + ggplot2::theme_minimal()
  
  end_time <- proc.time()
  elapsed_time <- round((end_time - start_time)["elapsed"], 2)
  message(paste("Total simulation time:", elapsed_time, "seconds"))
  
  cat("\n--- Simulation Summary ---\n")
  if (!is.null(results_summary)) {
    print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Difference)"))
  } else {
    cat("No valid estimates were generated to create a summary.\n")
  }
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

# Sample_Size_Search ------------------------------------------------------
#' @title Find Sample Size for a Linear RMST Model via Simulation
#' @description Performs an iterative sample size search to achieve a target power.
#' @inheritParams .get_linear_ipcw_simulation_runner
#' @param target_power A single numeric value for the target power (e.g., 0.80).
#' @param patience The number of consecutive non-improving steps before terminating.
#' @param n_start The starting sample size *per arm* for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size *per arm* to search up to.
#' @return A `list` containing results.
#' @keywords internal
#' @export
linear.ss.boot.app <- function(pilot_data, time_var, status_var, arm_var,
                           target_power,
                           linear_terms = NULL, L, n_sim = 1000, alpha = 0.05,
                           patience = 5,
                           n_start = 50, n_step = 25, max_n_per_arm = 2000, parallel.cores) {
  
  start_time <- proc.time()
  if (is.null(target_power)) stop("You must provide a single numeric value for 'target_power'.")
  
  # Get the configured simulation runner
  simulation_runner <- .get_linear_ipcw_simulation_runner(pilot_data, time_var, status_var, arm_var, linear_terms, L, alpha, n_sim)
  
  # --- Main Search Loop ---
  cat("--- Searching for Sample Size (Method: Linear RMST with IPCW) ---\n")
  current_n <- n_start
  max_power_so_far <- -1
  stagnation_counter <- 0
  n_at_max_power <- n_start
  best_sim_output <- NULL
  search_results <- list()
  final_n <- NA
  
  while (current_n <= max_n_per_arm) {
    cat(paste0("  N = ", current_n, "/arm, Calculating Power..."))
    
    sim_output <- simulation_runner(current_n)
    
    calculated_power <- if(is.finite(sim_output$power)) sim_output$power else 0
    search_results[[as.character(current_n)]] <- calculated_power
    cat(paste0(" Power = ", round(calculated_power, 3), "\n"))
    
    if (calculated_power >= target_power) {
      message("Success: Target power reached at N = ", current_n, "/arm.")
      best_sim_output <- sim_output
      final_n <- current_n
      break
    }
    if (calculated_power > max_power_so_far) {
      max_power_so_far <- calculated_power
      n_at_max_power <- current_n
      best_sim_output <- sim_output
      stagnation_counter <- 0
    } else {
      stagnation_counter <- stagnation_counter + 1
    }
    if (stagnation_counter >= patience) {
      warning(paste0("Search terminated due to stagnation. Best N = ", n_at_max_power,
                     " with power = ", round(max_power_so_far, 3)), call. = FALSE)
      final_n <- n_at_max_power
      break
    }
    current_n <- current_n + n_step
  }
  
  if (is.na(final_n) && current_n > max_n_per_arm) {
    warning(paste("Target power", target_power, "not achieved by max N of", max_n_per_arm))
    final_n <- n_at_max_power
  }
  
  # --- Finalize Summary, Plot, and Results ---
  results_summary <- NULL
  if (!is.null(best_sim_output)) {
    est <- na.omit(best_sim_output$estimates)
    se  <- na.omit(best_sim_output$std_errors)
    if(length(est) > 1) {
      results_summary <- data.frame(
        Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
        Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * sd(est), mean(est) + 1.96 * sd(est))
      )
    }
  }
  
  results_df <- data.frame(Target_Power = target_power, Required_N_per_Arm = final_n)
  search_path_df <- data.frame(N_per_Arm = as.integer(names(search_results)),
                               Power = unlist(search_results))
  
  p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Arm, y = Power)) +
    ggplot2::geom_line(color = "#009E73", linewidth = 1) +
    ggplot2::geom_point(color = "#009E73", size = 3) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
    ggplot2::labs(title = "Sample Size Search Path: Linear IPCW RMST Model",
                  x = "Sample Size Per Arm", y = "Calculated Power") +
    ggplot2::theme_minimal()
  
  end_time <- proc.time()
  elapsed_time <- round((end_time - start_time)["elapsed"], 2)
  message(paste("Total simulation time:", elapsed_time, "seconds"))
  
  cat("\n--- Simulation Summary ---\n")
  if (!is.null(results_summary)) {
    print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Difference)"))
  } else {
    cat("No valid estimates were generated to create a summary.\n")
  }
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}