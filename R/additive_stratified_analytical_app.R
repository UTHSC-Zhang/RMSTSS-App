#' @title Internal Helper to Estimate Additive Stratified Model Parameters
#' @description This internal function contains the common logic for estimating the
#'   treatment effect and its variance from pilot data for the additive stratified model.
#'   It is called by both the power and sample size calculation functions to avoid
#'   code duplication.
#' @return A list containing `beta_effect` (the estimated treatment effect) and
#'   `se_beta_n1` (the standard error for a sample size of 1).
#' @keywords internal
#' @export
.estimate_additive_stratified_params <- function(pilot_data, time_var, status_var, arm_var, strata_var, linear_terms, L) {
  
  # --- 1. Prepare Data and Calculate IPCW Weights ---
  cat("--- Estimating parameters from pilot data... ---\n")
  covariates <- c(arm_var, linear_terms)
  all_vars <- c(time_var, status_var, strata_var, covariates)
  df <- pilot_data[stats::complete.cases(pilot_data[, all_vars]), ]
  n_pilot <- nrow(df)
  min_events_per_stratum <- 2 # Define a minimum threshold
  event_counts <- stats::aggregate(as.formula(paste(status_var, "~", strata_var)), 
                            data = df, 
                            FUN = sum)
  
  problem_strata <- event_counts[event_counts[[status_var]] < min_events_per_stratum, ]
  
  if (nrow(problem_strata) > 0) {
    strata_names <- paste(problem_strata[[strata_var]], collapse = ", ")
    error_msg <- paste(
      "Insufficient events in one or more strata. Each stratum must have at least ",
      min_events_per_stratum, "events.",
      "Problematic strata:", strata_names
    )
    stop(error_msg, call. = FALSE)
  }
  # Validation for covariates
  if (!is.null(linear_terms)) {
    for (cov in linear_terms) {
      if (!is.numeric(df[[cov]])) {
        stop(paste0("Covariate '", cov, "' must be numeric, but its type is '", class(df[[cov]]), "'. Please check the column mapping in the UI."), call. = FALSE)
      }
    }
  }
  
  df$Y_rmst <- pmin(df[[time_var]], L)
  df$is_event <- df[[status_var]] == 1
  
  cens_formula <- stats::as.formula(paste0("survival::Surv(Y_rmst, is_event == 0) ~ ",
                                           paste(covariates, collapse = " + "),
                                           " + strata(", strata_var, ")"))
  fit_cens <- survival::coxph(cens_formula, data = df, ties = "breslow")
  
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
  df$weights <- exp(df$H_cens * exp(stats::predict(fit_cens, newdata=df, type="lp", reference="zero")))
  df$weights[!df$is_event] <- 0
  finite_weights <- df$weights[is.finite(df$weights) & df$weights > 0]
  if (length(finite_weights) > 0) {
    weight_cap <- stats::quantile(finite_weights, probs = 0.99, na.rm = TRUE)
    df$weights[df$weights > weight_cap] <- weight_cap
  }
  df$weights[!is.finite(df$weights)] <- 0
  
  # --- 2. Estimate Beta via Stratum-Centering ---
  cat("--- Estimating additive effect via stratum-centering... ---\n")
  vars_to_center <- c("Y_rmst", covariates)
  
  stratum_means <- df %>%
    dplyr::group_by(.data[[strata_var]]) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(vars_to_center),
        ~ weighted.mean(.x, w = weights, na.rm = TRUE)
      ),
      .groups = 'drop'
    )
  names(stratum_means) <- c(strata_var, paste0(vars_to_center, "_mean"))
  
  df_centered <- df %>%
    dplyr::left_join(stratum_means, by = strata_var)
  
  for (cov in vars_to_center) {
    df_centered[[paste0(cov, "_tilde")]] <- df_centered[[cov]] - df_centered[[paste0(cov, "_mean")]]
  }
  
  Z_tilde <- as.matrix(df_centered[, paste0(covariates, "_tilde")])
  Y_tilde <- df_centered[["Y_rmst_tilde"]]
  W <- df_centered$weights
  
  A_hat_num <- crossprod(Z_tilde * sqrt(W))
  dimnames(A_hat_num) <- list(covariates, covariates)
  A_hat <- A_hat_num / n_pilot
  
  A_hat_inv <- tryCatch({
    solve(A_hat)
  }, error = function(e) {
    stop("The covariate matrix (A_hat) is singular and cannot be inverted.\nThis may be caused by a lack of variation in covariates among subjects with an event within one or more strata.", call. = FALSE)
  })
  
  beta_hat <- (A_hat_inv / n_pilot) %*% (t(Z_tilde * W) %*% Y_tilde)
  rownames(beta_hat) <- covariates
  beta_effect <- beta_hat[arm_var, 1]
  
  # --- 3. Calculate Asymptotic Sandwich Variance ---
  cat("--- Calculating asymptotic variance... ---\n")
  mu0_hats <- stratum_means %>%
    dplyr::mutate(
      Z_matrix = as.matrix(dplyr::select(., dplyr::all_of(paste0(covariates, "_mean")))),
      mu0_hat = .data[["Y_rmst_mean"]] - Z_matrix %*% beta_hat
    ) %>%
    dplyr::select(dplyr::all_of(strata_var), mu0_hat)
  
  df_final <- df_centered %>% dplyr::left_join(mu0_hats, by = strata_var)
  Z_matrix <- as.matrix(df_final[, covariates])
  
  df_final$residuals <- df_final$Y_rmst - (df_final$mu0_hat + as.vector(Z_matrix %*% beta_hat))
  
  epsilon <- apply(Z_tilde, 2, function(z_col) z_col * W * df_final$residuals)
  B_hat <- crossprod(epsilon) / n_pilot
  dimnames(B_hat) <- list(covariates, covariates)
  
  V_hat_n <- A_hat_inv %*% B_hat %*% t(A_hat_inv)
  
  var_beta_pilot <- V_hat_n[arm_var, arm_var] / n_pilot
  var_beta_n1 <- var_beta_pilot * n_pilot
  se_beta_n1 <- sqrt(var_beta_n1)
  
  # Return the key parameters
  return(list(
    beta_effect = beta_effect,
    se_beta_n1 = se_beta_n1,
    n_strata = length(unique(df[[strata_var]]))
  ))
}

#' @title Analyze Power for a Stratified Additive RMST Model (Analytic)
#' @description Performs power analysis for a stratified, additive RMST model.
#' @inheritParams .estimate_additive_stratified_params
#' @param sample_sizes A numeric vector of sample sizes *per stratum* to calculate power for.
#' @param alpha The significance level (Type I error rate).
#' @return A list containing the results.
#' @keywords internal
#' @export
additive.power.analytical.app <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                                      sample_sizes, linear_terms = NULL, L, alpha = 0.05) {
  
  # 1. Estimate parameters from pilot data using the helper function
  params <- .estimate_additive_stratified_params(pilot_data, time_var, status_var, arm_var, strata_var, linear_terms, L)
  
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
  
  # 3. Create plot and summary
  results_summary <- data.frame(
    Statistic = "Assumed RMST Difference (from pilot)",
    Value = params$beta_effect
  )
  
  p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Stratum, y = Power)) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
    ggplot2::geom_point(color = "#0072B2", size = 3) +
    ggplot2::labs(
      title = "Analytic Power Curve: Additive Stratified RMST Model",
      subtitle = "Based on stratum-centered estimating equations.",
      x = "Sample Size Per Stratum", y = "Estimated Power"
    ) +
    ggplot2::ylim(0, 1) + ggplot2::theme_minimal()
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}


#' @title Find Sample Size for a Stratified Additive RMST Model (Analytic)
#' @description Calculates the required sample size for a target power.
#' @inheritParams .estimate_additive_stratified_params
#' @param target_power A single numeric value for the desired power.
#' @param n_start The starting sample size *per stratum* for the search.
#' @param n_step The increment in sample size at each step of the search.
#' @param max_n_per_arm The maximum sample size *per stratum* to search up to.
#' @param alpha The significance level (Type I error rate).
#' @return A list containing the results.
#' @keywords internal
#' @export
additive.ss.analytical.app <- function(pilot_data, time_var, status_var, arm_var, strata_var,
                                   target_power, linear_terms = NULL, L, alpha = 0.05,
                                   n_start = 50, n_step = 25, max_n_per_arm = 2000) {
  
  # 1. Estimate parameters from pilot data using the helper function
  params <- .estimate_additive_stratified_params(pilot_data, time_var, status_var, arm_var, strata_var, linear_terms, L)
  
  # 2. Iterative Search for Sample Size
  cat("--- Searching for Sample Size (Method: Additive Analytic) ---\n")
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
  results_summary <- data.frame(Statistic = "Assumed RMST Difference (from pilot)", Value = params$beta_effect)
  results_df <- data.frame(Target_Power = target_power, Required_N_per_Stratum = final_n)
  search_path_df <- data.frame(N_per_Stratum = as.integer(names(search_path)), Power = unlist(search_path))
  
  p <- ggplot2::ggplot(na.omit(search_path_df), ggplot2::aes(x = N_per_Stratum, y = Power)) +
    ggplot2::geom_line(color = "#009E73", linewidth = 1) +
    ggplot2::geom_point(color = "#009E73", size = 3) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = final_n, linetype = "dotted", color = "blue") +
    ggplot2::labs(
      title = "Analytic Sample Size Search: Additive Stratified RMST Model",
      subtitle = "Power calculated from formula at each step.",
      x = "Sample Size Per Stratum", y = "Calculated Power"
    ) + ggplot2::theme_minimal()
  
  cat("\n--- Calculation Summary ---\n")
  print(knitr::kable(results_df, caption = "Required Sample Size"))
  
  return(list(results_data = results_df, results_plot = p, results_summary = results_summary))
}