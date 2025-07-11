
#' @description An internal function for the Shiny app.
#' @keywords internal
#' @export
.get_gam_simulation_runner <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                                       linear_terms, smooth_terms, L, alpha, n_sim, parallel.cores) {
  
  # --- 1. One-Time Setup ---
  is_stratified <- !is.null(strata_var)
  
  # Prepare data
  core_vars <- c(time_var, status_var, arm_var, strata_var)
  all_vars <- c(core_vars, linear_terms, smooth_terms)
  cleaned_pilot_data <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
  
  # Model Formula Construction
  smooth_part <- if (!is.null(smooth_terms)) paste0("s(", smooth_terms, ")") else NULL
  if (is_stratified) {
    pilot_groups <- split(cleaned_pilot_data, cleaned_pilot_data[[strata_var]])
    all_terms <- c(strata_var, paste0(arm_var, ":", strata_var), smooth_part, linear_terms)
    test_term_pattern <- paste0("^", arm_var, "[^:]*:")
  } else {
    pilot_groups <- list(split(cleaned_pilot_data, cleaned_pilot_data[[arm_var]]))
    all_terms <- c(arm_var, smooth_part, linear_terms)
    test_term_pattern <- paste0("^", arm_var)
  }
  model_rhs <- paste(all_terms[!sapply(all_terms, is.null)], collapse = " + ")
  model_formula <- stats::as.formula(paste("pseudo_obs ~", model_rhs))
  message("Model: pseudo_obs ~ ", model_rhs)
  
  # Helper to calculate jackknife pseudo-observations
  get_pseudo_obs <- function(time, status, L_val) {
    n <- length(time)
    if (n == 0) return(numeric(0))
    km_fit_full <- survival::survfit(survival::Surv(time, status) ~ 1)
    km_step_full <- stats::stepfun(km_fit_full$time, c(1, km_fit_full$surv))
    rmst_full <- tryCatch(stats::integrate(km_step_full, 0, L_val, subdivisions=2000, stop.on.error = FALSE)$value, error = function(e) 0)
    rmst_loo <- vapply(seq_len(n), function(i) {
      if(n > 1) {
        km_fit_loo <- survival::survfit(survival::Surv(time[-i], status[-i]) ~ 1)
        km_step_loo <- stats::stepfun(km_fit_loo$time, c(1, km_fit_loo$surv))
        tryCatch(stats::integrate(km_step_loo, 0, L_val, subdivisions=2000, stop.on.error = FALSE)$value, error = function(e) 0)
      } else { 0 }
    }, FUN.VALUE = numeric(1))
    return(n * rmst_full - (n - 1) * rmst_loo)
  }
  
  # --- 2. The Returned Simulation Function ---
  run_simulation <- function(n_per_group) {
    
    sim_results_list <- future.apply::future_lapply(1:n_sim, function(j) {
      p_val <- NA_real_; est_val <- NA_real_; se_val <- NA_real_
      
      if (is_stratified) {
        boot_list <- lapply(pilot_groups, function(group_df) group_df[sample(1:nrow(group_df), size = n_per_group, replace = TRUE), ])
      } else {
        boot_list <- lapply(pilot_groups[[1]], function(arm_df) arm_df[sample(1:nrow(arm_df), size = n_per_group, replace = TRUE), ])
      }
      boot_data <- do.call(rbind, boot_list)
      
      grouping_var <- if(is_stratified) strata_var else arm_var
      pseudo_obs_list <- by(boot_data, boot_data[[grouping_var]], function(sub_data) {
        sub_data$pseudo_obs <- get_pseudo_obs(sub_data[[time_var]], sub_data[[status_var]], L)
        sub_data
      })
      boot_data <- do.call(rbind, pseudo_obs_list)
      
      fit <- tryCatch(mgcv::gam(model_formula, data = boot_data), error = function(e) NULL)
      if (!is.null(fit)) {
        sfit <- mgcv::summary.gam(fit)
        p_table <- sfit$p.table
        matching_rows <- grep(test_term_pattern, rownames(p_table), fixed = FALSE)
        if (length(matching_rows) > 0) {
          p_val <- min(p_table[matching_rows, "Pr(>|t|)"], na.rm = TRUE)
          est_val <- mean(p_table[matching_rows, "Estimate"], na.rm = TRUE)
          se_val <- mean(p_table[matching_rows, "Std. Error"], na.rm = TRUE)
        }
      }
      return(list(p_value = p_val, estimate = est_val, std_error = se_val))
    }, future.seed = TRUE)
    
    p_values <- vapply(sim_results_list, `[[`, "p_value", FUN.VALUE = numeric(1))
    estimates <- vapply(sim_results_list, `[[`, "estimate", FUN.VALUE = numeric(1))
    std_errors <- vapply(sim_results_list, `[[`, "std_error", FUN.VALUE = numeric(1))
    
    return(list(
      power = mean(p_values < alpha, na.rm = TRUE),
      estimates = estimates,
      std_errors = std_errors
    ))
  }
  
  # Return the fully configured simulation function
  return(run_simulation)
}


# Power Calculations -------------------------------------------------------
#' @description An internal function for the Shiny app.
#' @keywords internal
#' @export
additive.power.boot.app <- function(pilot_data, time_var, status_var, arm_var, strata_var = NULL,
                           sample_sizes, linear_terms = NULL, smooth_terms = NULL,
                           L, n_sim = 1000, alpha = 0.05,
                           parallel.cores = 1) {
  start_time <- proc.time()
  if (is.null(sample_sizes)) stop("You must provide a numeric vector for 'sample_sizes'.")
  
  # Configure parallel processing
  if (parallel.cores > 1) {
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Packages 'future' and 'future.apply' are required for parallel processing.")
    }
    future::plan(future::multisession, workers = parallel.cores)
    on.exit(future::plan(future::sequential), add = TRUE)
  }
  
  # Get the configured simulation runner
  simulation_runner <- .get_gam_simulation_runner(pilot_data, time_var, status_var, arm_var, strata_var,
                                                  linear_terms, smooth_terms, L, alpha, n_sim, parallel.cores)
  
  all_sim_outputs <- vector("list", length(sample_sizes))
  group_label <- if(!is.null(strata_var)) "/stratum" else "/arm"
  cat("--- Calculating Power (Method: Additive GAM for RMST) ---\n")
  
  for (i in seq_along(sample_sizes)) {
    n_per_group <- sample_sizes[i]
    cat("Simulating for n =", n_per_group, group_label, "...\n")
    all_sim_outputs[[i]] <- simulation_runner(n_per_group)
  }
  
  power_values <- sapply(all_sim_outputs, `[[`, "power")
  results_df <- data.frame(N_per_Group = sample_sizes, Power = power_values)
  
  # Summarize results from the largest sample size simulation
  best_idx <- which.max(sample_sizes)
  est <- stats::na.omit(all_sim_outputs[[best_idx]]$estimates)
  se  <- stats::na.omit(all_sim_outputs[[best_idx]]$std_errors)
  results_summary <- NULL
  if (length(est) > 1) {
    results_summary <- data.frame(
      Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
      Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * stats::sd(est), mean(est) + 1.96 * stats::sd(est))
    )
  }
  
  p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Group, y = Power)) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
    ggplot2::geom_point(color = "#0072B2", size = 3) +
    ggplot2::labs(title = "Power Curve: Additive GAM RMST Model",
                  x = paste0("Sample Size Per ", if(!is.null(strata_var)) "Stratum" else "Arm"),
                  y = "Estimated Power") +
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


# Sample Size Search ------------------------------------------------------


#' @description An internal function for the Shiny app.
#' @keywords internal
#' @export
additive.ss.boot.app <- function(pilot_data, time_var, status_var, arm_var, strata_var = NULL,
                        target_power,
                        linear_terms = NULL, smooth_terms = NULL,
                        L, n_sim = 1000, alpha = 0.05,
                        parallel.cores = 1, patience = 5,
                        n_start = 50, n_step = 25, max_n_per_arm = 2000) {
  
  start_time <- proc.time()
  if (is.null(target_power)) stop("You must provide a single numeric value for 'target_power'.")
  
  # Configure parallel processing
  if (parallel.cores > 1) {
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Packages 'future' and 'future.apply' are required for parallel processing.")
    }
    future::plan(future::multisession, workers = parallel.cores)
    on.exit(future::plan(future::sequential), add = TRUE)
  }
  
  # Get the configured simulation runner
  simulation_runner <- .get_gam_simulation_runner(pilot_data, time_var, status_var, arm_var, strata_var,
                                                  linear_terms, smooth_terms, L, alpha, n_sim, parallel.cores)
  
  # --- Main Search Loop ---
  group_label <- if(!is.null(strata_var)) "/stratum" else "/arm"
  cat("--- Searching for Sample Size (Method: Additive GAM for RMST) ---\n")
  
  current_n <- n_start
  max_power_so_far <- -1
  stagnation_counter <- 0
  n_at_max_power <- n_start
  best_sim_output <- NULL
  search_results <- list()
  final_n <- NA
  
  while (current_n <= max_n_per_arm) {
    cat(paste0("  N = ", current_n, group_label, ", Calculating Power..."))
    
    sim_output <- simulation_runner(current_n)
    
    calculated_power <- if(is.finite(sim_output$power)) sim_output$power else 0
    search_results[[as.character(current_n)]] <- calculated_power
    cat(paste0(" Power = ", round(calculated_power, 3), "\n"))
    
    if (calculated_power >= target_power) {
      message("Success: Target power reached at N = ", current_n, group_label, ".")
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
    est <- stats::na.omit(best_sim_output$estimates)
    se  <- stats::na.omit(best_sim_output$std_errors)
    if(length(est) > 1) {
      results_summary <- data.frame(
        Statistic = c("Mean RMST Difference", "Mean Standard Error", "95% CI Lower", "95% CI Upper"),
        Value = c(mean(est), mean(se, na.rm=TRUE), mean(est) - 1.96 * stats::sd(est), mean(est) + 1.96 * stats::sd(est))
      )
    }
  }
  
  results_df <- data.frame(Target_Power = target_power, Required_N_per_Group = final_n)
  search_path_df <- data.frame(N_per_Group = as.integer(names(search_results)),
                               Power = unlist(search_results))
  
  p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Group, y = Power)) +
    ggplot2::geom_line(color = "#009E73", linewidth = 1) +
    ggplot2::geom_point(color = "#009E73", size = 3) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
    ggplot2::labs(title = "Sample Size Search Path: Additive GAM RMST Model",
                  x = paste0("Sample Size Per ", if(!is.null(strata_var)) "Stratum" else "Arm"),
                  y = "Calculated Power") +
    ggplot2::theme_minimal()
  
  end_time <- proc.time()
  elapsed_time <- round((end_time - start_time)["elapsed"], 2)
  message(paste("Total simulation time:", elapsed_time, "seconds"))
  
  cat("\n--- Simulation Summary ----\n")
  if (!is.null(results_summary)) {
    print(knitr::kable(results_summary, caption = "Estimated Treatment Effect (RMST Difference)"))
  } else {
    cat("No valid estimates were generated to create a summary.\n")
  }
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}