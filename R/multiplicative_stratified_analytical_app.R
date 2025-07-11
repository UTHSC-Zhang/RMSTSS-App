#' @importFrom stats aggregate
#' @keywords internal
#' @export
.estimate_multiplicative_stratified_params <- function(pilot_data, time_var, status_var, arm_var, strata_var, linear_terms, L) {
  
  # --- 1. Estimate Parameters from Pilot Data (log-linear approximation) ---
  cat("--- Estimating parameters from pilot data (log-linear approximation)... ---\n")
  
  core_vars <- c(time_var, status_var, arm_var, strata_var)
  covariates <- c(arm_var, linear_terms)
  all_vars <- c(core_vars, linear_terms)
  df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
  n_pilot <- nrow(df)
  min_events_per_stratum <- 2 # Define a minimum threshold
  event_counts <- aggregate(as.formula(paste(status_var, "~", strata_var)), 
                            data = df, 
                            FUN = sum)
  
  problem_strata <- event_counts[event_counts[[status_var]] < min_events_per_stratum, ]
  
  if (nrow(problem_strata) > 0) {
    strata_names <- paste(problem_strata[[strata_var]], collapse = ", ")
    error_msg <- paste(
      "Insufficient events in one or more strata. Each stratum must have at least",
      min_events_per_stratum, "events.",
      "Problematic strata:", strata_names
    )
    stop(error_msg, call. = FALSE)
  }
  # Prepare data for IPCW
  df$Y_rmst <- pmin(df[[time_var]], L)
  df$is_event <- df[[status_var]] == 1
  
  # Model for censoring (stratified by strata_var)
  cens_formula <- stats::as.formula(paste0("survival::Surv(Y_rmst, is_event == 0) ~ ",
                                           paste(covariates, collapse = " + "),
                                           " + strata(", strata_var, ")"))
  fit_cens <- survival::coxph(cens_formula, data = df, ties = "breslow")
  
  # Calculate IPCW weights
  bh_cens <- survival::basehaz(fit_cens, centered = FALSE)
  df$H_cens <- 0
  unique_strata_from_bh <- unique(bh_cens$strata)
  for(st in unique(df[[strata_var]])){
    st_label <- paste0(strata_var, "=", st)
    is_stratum <- df[[strata_var]] == st
    if (st_label %in% unique_strata_from_bh) {
      is_bh_stratum <- bh_cens$strata == st_label
      if(sum(is_bh_stratum) > 0){
        H_st <- stats::stepfun(bh_cens$time[is_bh_stratum], c(0, bh_cens$hazard[is_bh_stratum]))(df$Y_rmst[is_stratum])
        df$H_cens[is_stratum] <- H_st
      }
    }
  }
  df$weights <- exp(df$H_cens * exp(predict(fit_cens, newdata = df, type="lp", reference="zero")))
  
  # Stabilize weights
  finite_weights <- df$weights[is.finite(df$weights)]
  if (length(finite_weights) > 0) {
    weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
    df$weights[df$weights > weight_cap] <- weight_cap
  }
  df$weights[!is.finite(df$weights)] <- 0
  df$w_delta <- df$weights * df$is_event
  df$w_delta[is.na(df$w_delta)] <- 0
  
  # Use a log-linear model as a practical approximation to the multiplicative model
  fit_data <- df[df$w_delta > 0 & df$Y_rmst > 0, ]
  fit_weights <- fit_data$w_delta
  
  log_model_formula <- stats::as.formula(paste("log(Y_rmst) ~", paste(c(covariates, strata_var), collapse=" + ")))
  message("Approximation Model: ", deparse(log_model_formula))
  
  if(nrow(fit_data) < (length(covariates) + length(unique(fit_data[[strata_var]])))) {
    stop("Not enough data points to fit the approximation model after filtering.", call. = FALSE)
  }
  
  fit_log_lm <- lm(log_model_formula, data = fit_data, weights = fit_weights)
  
  beta_summary <- coef(summary(fit_log_lm))
  beta_effect <- beta_summary[arm_var, "Estimate"]
  
  # --- 2. Calculate Asymptotic Variance from the approximate model ---
  V_hat <- vcov(fit_log_lm) * nrow(fit_data) # Scale to Var(sqrt(n)*beta)
  var_beta_n1 <- V_hat[arm_var, arm_var]
  se_beta_n1 <- sqrt(var_beta_n1)
  
  return(list(
    beta_effect = beta_effect,
    se_beta_n1 = se_beta_n1,
    n_strata = length(unique(df[[strata_var]])),
    fit_cens = fit_cens
  ))
}


# Power Calculation -------------------------------------------------------
#' @keywords internal
#' @export
MS.power.analytical.app <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                                sample_sizes, linear_terms = NULL, L, alpha = 0.05) {
  
  # 1. Estimate parameters from pilot data
  params <- .estimate_multiplicative_stratified_params(pilot_data, time_var, status_var, arm_var, strata_var, linear_terms, L)
  
  # 2. Calculate Power
  cat("--- Calculating power for specified sample sizes... ---\n")
  z_alpha <- stats::qnorm(1 - alpha / 2)
  power_values <- sapply(sample_sizes, function(n_per_stratum) {
    total_n <- n_per_stratum * params$n_strata
    se_final <- params$se_beta_n1 / sqrt(total_n)
    power <- stats::pnorm( (abs(params$beta_effect) / se_final) - z_alpha )
    return(power)
  })
  
  results_df <- data.frame(N_per_Stratum = sample_sizes, Power = power_values)
  
  # 3. Create Plot and Return
  p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Stratum, y = Power)) +
    ggplot2::geom_line(color = "#E69F00", linewidth = 1) +
    ggplot2::geom_point(color = "#E69F00", size = 3) +
    ggplot2::labs(
      title = "Analytic Power Curve: Multiplicative Stratified RMST Model",
      subtitle = "Approximate method based on Wang et al. (2019).",
      x = "Sample Size Per Stratum", y = "Estimated Power"
    ) +
    ggplot2::ylim(0, 1) + ggplot2::theme_minimal()
  
  results_summary <- data.frame(
    Statistic = "Assumed log(RMST Ratio) (from pilot)",
    Value = params$beta_effect
  )
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}

# Sample Size Search ------------------------------------------------------
#' @keywords internal
#' @export
MS.ss.analytical.app <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                             target_power, linear_terms = NULL, L, alpha = 0.05,
                             n_start = 50, n_step = 25, max_n_per_arm = 2000) {
  
  # 1. Estimate parameters from pilot data
  params <- .estimate_multiplicative_stratified_params(pilot_data, time_var, status_var, arm_var, strata_var, linear_terms, L)
  
  # 2. Iterative Search for Sample Size
  cat("--- Searching for Sample Size (Method: Analytic/Approximation) ---\n")
  current_n <- n_start
  search_path <- list()
  final_n <- NA_integer_
  z_alpha <- stats::qnorm(1 - alpha / 2)
  
  while (current_n <= max_n_per_arm) {
    total_n <- current_n * params$n_strata
    se_final <- params$se_beta_n1 / sqrt(total_n)
    calculated_power <- stats::pnorm((abs(params$beta_effect) / se_final) - z_alpha)
    if (!is.finite(calculated_power)) calculated_power <- 0
    
    search_path[[as.character(current_n)]] <- calculated_power
    cat(paste0("  N = ", current_n, "/stratum, Calculated Power = ", round(calculated_power, 3), "\n"))
    
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
  results_summary <- data.frame(Statistic = "Assumed log(RMST Ratio) (from pilot)", Value = params$beta_effect)
  results_df <- data.frame(Target_Power = target_power, Required_N_per_Stratum = final_n)
  search_path_df <- data.frame(N_per_Stratum = as.integer(names(search_path)), Power = unlist(search_path))
  
  p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Stratum, y = Power)) +
    ggplot2::geom_line(color = "#009E73", linewidth = 1) +
    ggplot2::geom_point(color = "#009E73", size = 3) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
    ggplot2::labs(
      title = "Analytic Sample Size Search: Multiplicative Stratified RMST Model",
      subtitle = "Power calculated from formula at each step (approximate method).",
      x = "Sample Size Per Stratum", y = "Calculated Power"
    ) + ggplot2::theme_minimal()
  
  cat("\n--- Calculation Summary ---\n")
  print(knitr::kable(results_df, caption = "Required Sample Size"))
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}