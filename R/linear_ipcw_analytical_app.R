
.estimate_linear_ipcw_params <- function(pilot_data, time_var, status_var, arm_var, linear_terms, L) {
  
  cat("--- Estimating parameters from pilot data for analytic calculation... ---\n")
  
  core_vars <- c(time_var, status_var, arm_var)
  all_vars <- c(core_vars, linear_terms)
  df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
  n_pilot <- nrow(df)
  
  # Define model formulas
  factor_arm_str <- paste0("factor(", arm_var, ")")
  model_rhs <- paste(c(factor_arm_str, linear_terms), collapse = " + ")
  model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))
  message("Model: Y_rmst ~ ", model_rhs)
  
  # Prepare data for IPCW
  df$Y_rmst <- pmin(df[[time_var]], L)
  df$is_censored <- df[[status_var]] == 0
  df$is_event <- df[[status_var]] == 1
  
  # Fit censoring model (Kaplan-Meier for G(t))
  cens_fit <- survival::survfit(survival::Surv(Y_rmst, is_censored) ~ 1, data = df)
  cens_surv_prob <- stats::stepfun(cens_fit$time, c(1, cens_fit$surv))(df$Y_rmst)
  df$weights <- df$is_event / cens_surv_prob
  
  # Stabilize weights
  finite_weights <- df$weights[is.finite(df$weights) & df$weights > 0]
  if (length(finite_weights) > 0) {
    weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
    df$weights[df$weights > weight_cap] <- weight_cap
  }
  df$weights[!is.finite(df$weights)] <- 0
  
  # Filter for model fitting
  fit_data <- df[df$weights > 0, ]
  fit_weights <- fit_data$weights
  
  if (length(unique(fit_data[[arm_var]])) < 2) {
    stop("Pilot data contains events in only one arm after filtering.", call. = FALSE)
  }
  
  # Fit the primary linear model
  fit_lm <- stats::lm(model_formula, data = fit_data, weights = fit_weights)
  beta_hat <- stats::coef(fit_lm)
  
  arm_pattern <- paste0("^factor\\(", arm_var, "\\)1$")
  arm_coeff_name <- names(beta_hat)[grep(arm_pattern, names(beta_hat))]
  if (length(arm_coeff_name) == 0) stop("Could not find treatment effect coefficient.")
  beta_effect <- beta_hat[arm_coeff_name]
  
  # --- Calculate Asymptotic Sandwich Variance Components ---
  cat("--- Calculating asymptotic variance... ---\n")
  X <- stats::model.matrix(model_formula, data = df)
  
  A_hat <- crossprod(X * sqrt(df$weights), X * sqrt(df$weights)) / n_pilot
  
  A_hat_inv <- tryCatch({
    solve(A_hat)
  }, error = function(e) {
    stop("The design matrix (A_hat) is singular. This can happen with small pilot datasets or perfect separation.", call. = FALSE)
  })
  
  df$predicted_rmst <- predict(fit_lm, newdata = df)
  residuals <- df$Y_rmst - df$predicted_rmst
  
  epsilon <- X * residuals * df$weights
  epsilon[is.na(epsilon)] <- 0
  B_hat <- crossprod(epsilon) / n_pilot
  
  V_hat_n <- A_hat_inv %*% B_hat %*% t(A_hat_inv)
  
  # Scale to get the variance for n=1
  var_beta_n1 <- V_hat_n[arm_coeff_name, arm_coeff_name]
  se_beta_n1 <- sqrt(var_beta_n1)
  
  # Return a list containing the key estimated parameters
  return(list(
    beta_effect = beta_effect,
    se_beta_n1 = se_beta_n1,
    fit_lm = fit_lm
  ))
}



# Power Calculation -------------------------------------------------------------------

linear.power.analytical.app <- function(pilot_data, time_var, status_var, arm_var,
                                    sample_sizes, linear_terms = NULL, L, alpha = 0.05) {
  
  # 1. Estimate nuisance parameters from pilot data using the helper function
  params <- .estimate_linear_ipcw_params(pilot_data, time_var, status_var, arm_var, linear_terms, L)
  
  # 2. Calculate Power for Each Sample Size
  cat("--- Calculating power for specified sample sizes... ---\n")
  z_alpha <- stats::qnorm(1 - alpha / 2)
  power_values <- sapply(sample_sizes, function(n_per_arm) {
    total_n <- n_per_arm * 2
    se_final <- params$se_beta_n1 / sqrt(total_n)
    power <- stats::pnorm( (abs(params$beta_effect) / se_final) - z_alpha )
    return(power)
  })
  
  results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)
  
  # 3. Create Summary and Plot
  results_summary <- data.frame(
    Statistic = "Assumed RMST Difference (from pilot)",
    Value = params$beta_effect
  )
  
  p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
    ggplot2::geom_line(color = "#D55E00", linewidth = 1) +
    ggplot2::geom_point(color = "#D55E00", size = 3) +
    ggplot2::labs(
      title = "Analytic Power Curve: Linear IPCW RMST Model",
      subtitle = "Based on the asymptotic variance from Tian et al. (2014).",
      x = "Sample Size Per Arm", y = "Estimated Power"
    ) +
    ggplot2::ylim(0, 1) + ggplot2::theme_minimal()
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

# Sample Size Search ------------------------------------------------------

linear.ss.analytical.app <- function(pilot_data, time_var, status_var, arm_var,
                                 target_power, linear_terms = NULL, L, alpha = 0.05,
                                 n_start = 50, n_step = 25, max_n_per_arm = 2000) {
  
  # 1. Estimate parameters from pilot data using the helper function
  params <- .estimate_linear_ipcw_params(pilot_data, time_var, status_var, arm_var, linear_terms, L)
  
  # 2. Iterative Search for Sample Size using Analytic Formula
  cat("--- Searching for Sample Size (Method: Analytic) ---\n")
  current_n <- n_start
  search_path <- list()
  final_n <- NA_integer_
  z_alpha <- stats::qnorm(1 - alpha / 2)
  
  while (current_n <= max_n_per_arm) {
    total_n <- current_n * 2
    se_final <- params$se_beta_n1 / sqrt(total_n)
    calculated_power <- stats::pnorm((abs(params$beta_effect) / se_final) - z_alpha)
    if (!is.finite(calculated_power)) calculated_power <- 0
    
    search_path[[as.character(current_n)]] <- calculated_power
    cat(paste0("  N = ", current_n, "/arm, Calculated Power = ", round(calculated_power, 3), "\n"))
    
    if (calculated_power >= target_power) {
      final_n <- current_n
      break
    }
    current_n <- current_n + n_step
  }
  
  if (is.na(final_n)) {
    warning(paste("Target power", target_power, "not achieved by max N of", max_n_per_arm), call. = FALSE)
    final_n <- max_n_per_arm
  }
  
  # 3. Finalize and Return Results
  results_summary <- data.frame(
    Statistic = "Assumed RMST Difference (from pilot)",
    Value = params$beta_effect
  )
  results_df <- data.frame(Target_Power = target_power, Required_N_per_Arm = final_n)
  search_path_df <- data.frame(N_per_Arm = as.integer(names(search_path)), Power = unlist(search_path))
  
  p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Arm, y = Power)) +
    ggplot2::geom_line(color = "#009E73", linewidth = 1) +
    ggplot2::geom_point(color = "#009E73", size = 3) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
    ggplot2::labs(
      title = "Analytic Sample Size Search: Linear IPCW RMST Model",
      subtitle = "Power calculated from formula at each step.",
      x = "Sample Size Per Arm", y = "Calculated Power"
    ) + ggplot2::theme_minimal()
  
  cat("\n--- Calculation Summary ---\n")
  print(knitr::kable(results_df, caption = "Required Sample Size"))
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}