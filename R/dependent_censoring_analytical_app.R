#' @keywords internal
#' @export
.estimate_dependent_censoring_params <- function(pilot_data, time_var, status_var, arm_var, dep_cens_status_var, linear_terms, L) {
  
  # --- 1. Estimate Nuisance Parameters from Pilot Data ---
  cat("--- Estimating parameters from pilot data... ---\n")
  
  core_vars <- c(time_var, status_var, arm_var, dep_cens_status_var)
  covariates <- c(arm_var, linear_terms)
  all_vars <- c(core_vars, covariates)
  df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
  n_pilot <- nrow(df)
  
  # Define model formulas
  model_rhs <- paste(covariates, collapse = " + ")
  cens_formula <- stats::as.formula(paste("survival::Surv(", time_var, ", cens_ind) ~", model_rhs))
  dep_cens_formula <- stats::as.formula(paste("survival::Surv(", time_var, ", dep_cens_ind) ~", model_rhs))
  model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))
  
  message("Model: Y_rmst ~ ", model_rhs)
  
  # Define censoring indicators and truncated time
  df$cens_ind <- df[[status_var]] == 0 & df[[dep_cens_status_var]] == 0
  df$dep_cens_ind <- df[[dep_cens_status_var]] == 1
  df$Y_rmst <- pmin(df[[time_var]], L)
  
  # Fit censoring models
  fit_cens <- tryCatch(survival::coxph(cens_formula, data = df, ties = "breslow"), error = function(e) NULL)
  fit_dep_cens <- tryCatch(survival::coxph(dep_cens_formula, data = df, ties = "breslow"), error = function(e) NULL)
  
  if(is.null(fit_cens) || is.null(fit_dep_cens)) {
    stop("One or both of the censoring models failed to fit on the pilot data.", call. = FALSE)
  }
  
  # Calculate IPCW weights
  bh_cens <- survival::basehaz(fit_cens, centered = FALSE)
  H_cens <- stats::stepfun(bh_cens$time, c(0, bh_cens$hazard))(df$Y_rmst) * exp(predict(fit_cens, newdata = df, type = "lp"))
  
  bh_dep_cens <- survival::basehaz(fit_dep_cens, centered = FALSE)
  H_dep_cens <- stats::stepfun(bh_dep_cens$time, c(0, bh_dep_cens$hazard))(df$Y_rmst) * exp(predict(fit_dep_cens, newdata = df, type = "lp"))
  
  df$weights <- exp(H_cens + H_dep_cens)
  df$weights[df[[status_var]] != 1] <- 0 # Weights are only applied to uncensored subjects
  
  # Stabilize weights
  finite_weights <- df$weights[is.finite(df$weights) & df$weights > 0]
  if (length(finite_weights) > 0) {
    weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
    df$weights[df$weights > weight_cap] <- weight_cap
  }
  df$weights[!is.finite(df$weights)] <- 0
  
  # Filter data for final model
  fit_data <- df[df$weights > 0, ]
  fit_weights <- fit_data$weights
  
  if (length(unique(fit_data[[arm_var]])) < 2) {
    stop("Pilot data contains primary events in only one arm after filtering. Cannot estimate treatment effect.", call. = FALSE)
  }
  
  # Fit the primary linear model
  fit_lm <- stats::lm(model_formula, data = fit_data, weights = fit_weights)
  beta_hat <- stats::coef(fit_lm)
  arm_coeff_name <- arm_var
  beta_effect <- beta_hat[arm_coeff_name]
  
  # --- 2. Calculate Asymptotic Sandwich Variance Components ---
  cat("--- Calculating asymptotic variance... ---\n")
  X <- stats::model.matrix(model_formula, data = df)
  
  A_hat <- crossprod(X) / n_pilot
  
  df$predicted_rmst <- predict(fit_lm, newdata = df)
  residuals <- df$Y_rmst - df$predicted_rmst
  epsilon <- X * residuals * df$weights
  epsilon[is.na(epsilon)] <- 0
  B_hat <- crossprod(epsilon) / n_pilot
  
  A_hat_inv <- solve(A_hat)
  V_hat_n <- A_hat_inv %*% B_hat %*% t(A_hat_inv)
  
  var_beta_pilot <- V_hat_n[arm_coeff_name, arm_coeff_name] / n_pilot
  var_beta_n1 <- var_beta_pilot * n_pilot
  se_beta_n1 <- sqrt(var_beta_n1)
  
  return(list(
    beta_effect = beta_effect,
    se_beta_n1 = se_beta_n1
  ))
}

# Power Calculation -------------------------------------------------------

#' @keywords internal
#' @export
DC.power.analytical.app <- function(pilot_data,
                                time_var,
                                status_var,
                                arm_var,
                                dep_cens_status_var,
                                sample_sizes,
                                linear_terms = NULL,
                                L,
                                alpha = 0.05) {
  
  # 1. Estimate parameters from pilot data
  params <- .estimate_dependent_censoring_params(pilot_data, time_var, status_var, arm_var, dep_cens_status_var, linear_terms, L)
  
  # 2. Calculate Power for Each Sample Size
  cat("--- Calculating power for specified sample sizes... ---\n")
  z_alpha <- stats::qnorm(1 - alpha / 2)
  power_values <- sapply(sample_sizes, function(n_per_arm) {
    total_n <- n_per_arm * 2 # Assuming two arms
    se_final <- params$se_beta_n1 / sqrt(total_n)
    power <- stats::pnorm( (abs(params$beta_effect) / se_final) - z_alpha )
    return(power)
  })
  
  results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)
  
  # 3. Create Plot and Return
  p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
    ggplot2::geom_point(color = "#0072B2", size = 3) +
    ggplot2::labs(
      title = "Analytic Power Curve: Dependent Censoring RMST Model",
      subtitle = "Based on the asymptotic variance of the estimator.",
      x = "Sample Size Per Arm", y = "Estimated Power"
    ) +
    ggplot2::ylim(0, 1) + ggplot2::theme_minimal()
  
  results_summary <- data.frame(
    Statistic = "Assumed RMST Difference (from pilot)",
    Value = params$beta_effect
  )
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

# Sample Size Search ------------------------------------------------------
#' @keywords internal
#' @export
DC.ss.analytical.app <- function(pilot_data,
                             time_var,
                             status_var,
                             arm_var,
                             dep_cens_status_var,
                             target_power,
                             linear_terms = NULL,
                             L,
                             alpha = 0.05,
                             n_start = 50,
                             n_step = 25,
                             max_n_per_arm = 2000) {
  
  # 1. Estimate parameters from pilot data
  params <- .estimate_dependent_censoring_params(pilot_data, time_var, status_var, arm_var, dep_cens_status_var, linear_terms, L)
  
  # 2. Iterative Search for Sample Size
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
  
  search_path_df <- data.frame(
    N_per_Arm = as.integer(names(search_path)),
    Power = as.numeric(unlist(search_path))
  )
  
  p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Arm, y = Power)) +
    ggplot2::geom_line(color = "#009E73", linewidth = 1) +
    ggplot2::geom_point(color = "#009E73", size = 3) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
    ggplot2::labs(
      title = "Analytic Sample Size Search Path",
      subtitle = "Power is calculated directly from formula at each step.",
      x = "Sample Size Per Arm", y = "Calculated Power"
    ) +
    ggplot2::theme_minimal()
  
  cat("\n--- Calculation Summary ---\n")
  print(knitr::kable(results_df, caption = "Required Sample Size"))
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}