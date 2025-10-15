#' @title Internal Helper to Estimate Covariate-Dependent Censoring Model Parameters
#' @description Internal helper for estimating the treatment effect and its variance
#'   from pilot data under a single censoring mechanism that depends on covariates.
#'
#' This version uses one Cox model for censoring:
#' \code{Surv(time, status==0) ~ linear_terms} (treatment excluded by default),
#' builds IPCW \eqn{w_i = 1/\hat G(Y_i\mid X_i)} at \eqn{Y_i=\min(T_i,L)}, fits a weighted
#' RMST regression, and computes a sandwich variance that ignores uncertainty in \eqn{\hat G}.
#'
#' @note \code{dep_cens_status_var} is accepted for API compatibility but ignored here.
#' @return A list with \code{beta_effect} (arm effect) and \code{se_beta_n1} (SE scaled to N=1).
#' @keywords internal
#' @export
.estimate_dependent_censoring_params <- function(pilot_data, time_var, status_var, arm_var, dep_cens_status_var, linear_terms, L) {

   # --- 1) Prep & outcome ---
   covariates <- c(arm_var, linear_terms)
   # Ignore dep_cens_status_var in filtering (kept only for API compatibility)
   all_vars <- unique(c(time_var, status_var, covariates))
   df <- pilot_data[stats::complete.cases(pilot_data[, all_vars, drop = FALSE]), ]
   n_pilot <- nrow(df)
   if (n_pilot < 10) stop("Too few complete cases after filtering.", call. = FALSE)

   df$Y_rmst <- pmin(df[[time_var]], L)

   # --- 2) Censoring model: covariate-dependent only (no arm, no competing risks) ---
   cens_rhs <- if (is.null(linear_terms) || length(linear_terms) == 0) "1" else paste(linear_terms, collapse = " + ")
   cens_formula <- stats::as.formula(
      paste0("survival::Surv(", time_var, ", ", status_var, "==0) ~ ", cens_rhs)
   )
   fit_cens <- survival::coxph(cens_formula, data = df, ties = "breslow")

   bh <- survival::basehaz(fit_cens, centered = FALSE)
   H0_step <- stats::stepfun(bh$time, c(0, bh$hazard))
   lp <- if (cens_rhs == "1") 0 else stats::predict(fit_cens, newdata = df, type = "lp")
   Hc <- H0_step(df$Y_rmst) * exp(lp)
   Ghat <- exp(-Hc)

   # IPCW for all subjects (stabilize and cap heavy tails)
   eps <- 1e-6
   w <- 1 / pmax(Ghat, eps)
   w[!is.finite(w)] <- 0
   if (any(is.finite(w) & w > 0)) {
      cap <- stats::quantile(w[is.finite(w) & w > 0], 0.99, na.rm = TRUE)
      w[w > cap] <- cap
   }
   df$weights <- w

   # --- 3) Weighted RMST regression: Y_rmst ~ arm + covariates ---
   model_rhs <- paste(c(arm_var, linear_terms), collapse = " + ")
   model_formula <- stats::as.formula(paste("Y_rmst ~", model_rhs))
   fit_lm <- stats::lm(model_formula, data = df, weights = df$weights)
   beta_hat <- stats::coef(fit_lm)
   arm_coeff_name <- arm_var
   if (!(arm_coeff_name %in% names(beta_hat))) stop("Arm coefficient not found in model.", call. = FALSE)
   beta_effect <- beta_hat[arm_coeff_name]

   # --- 4) Sandwich variance: A = X' W X / n, B = X' diag((w r)^2) X / n ---
   X <- stats::model.matrix(model_formula, data = df)
   r <- df$Y_rmst - as.numeric(X %*% beta_hat)
   A_hat <- crossprod(X, X * df$weights) / n_pilot
   meat <- X * (df$weights * r)
   B_hat <- crossprod(meat) / n_pilot

   A_hat_inv <- solve(A_hat)
   V_hat_n <- A_hat_inv %*% B_hat %*% A_hat_inv

   var_beta_pilot <- V_hat_n[colnames(X) == arm_coeff_name, colnames(X) == arm_coeff_name] / n_pilot
   var_beta_n1 <- as.numeric(var_beta_pilot * n_pilot)
   se_beta_n1 <- sqrt(var_beta_n1)

   list(beta_effect = beta_effect, se_beta_n1 = se_beta_n1)
}
# Power Calculation -------------------------------------------------------

#' @title Analyze Power for RMST with Covariate-Dependent Censoring (Analytic)
#' @description Analytic power for RMST regression when censoring depends on covariates
#'   (single mechanism; no competing risks). Uses IPCW and a sandwich variance.
#' @inheritParams .estimate_dependent_censoring_params
#' @param sample_sizes Numeric vector of per-arm sample sizes for power.
#' @param alpha Two-sided Type I error (default 0.05).
#' @return A list: \code{results_data} (data.frame), \code{results_plot} (ggplot),
#'   \code{results_summary} (data.frame with pilot effect).
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

   # 1) Estimate parameters from pilot data
   params <- .estimate_dependent_censoring_params(
      pilot_data, time_var, status_var, arm_var, dep_cens_status_var, linear_terms, L
   )

   # 2) Power across N
   cat("--- Calculating power for specified sample sizes... ---\n")
   z_alpha <- stats::qnorm(1 - alpha / 2)
   power_values <- sapply(sample_sizes, function(n_per_arm) {
      total_n <- n_per_arm * 2
      se_final <- params$se_beta_n1 / sqrt(total_n)
      stats::pnorm((abs(params$beta_effect) / se_final) - z_alpha)
   })

   results_df <- data.frame(N_per_Arm = sample_sizes, Power = power_values)

   # 3) Plot + return
   p <- ggplot2::ggplot(results_df, ggplot2::aes(x = N_per_Arm, y = Power)) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
      ggplot2::geom_point(color = "#0072B2", size = 3) +
      ggplot2::labs(
         title = "Analytic Power Curve: RMST with Covariate-Dependent Censoring (IPCW)",
         subtitle = "Single censoring mechanism; variance ignores uncertainty in Ĝ.",
         x = "Sample Size Per Arm", y = "Estimated Power"
      ) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_minimal()

   results_summary <- data.frame(
      Statistic = "Assumed RMST Difference (from pilot)",
      Value = params$beta_effect
   )

   list(results_data = results_df, results_plot = p, results_summary = results_summary)
}
# Sample Size Search ------------------------------------------------------

#' @title Find Sample Size for RMST with Covariate-Dependent Censoring (Analytic)
#' @description Iterative per-arm N to reach target power using IPCW-based analytic variance.
#' @inheritParams .estimate_dependent_censoring_params
#' @param target_power Desired power.
#' @param alpha Two-sided Type I error (default 0.05).
#' @param n_start Starting per-arm N.
#' @param n_step Step size for search.
#' @param max_n_per_arm Maximum per-arm N to search.
#' @return A list: \code{results_data}, \code{results_plot}, \code{results_summary}.
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

   # 1) One-time estimation from pilot data
   params <- .estimate_dependent_censoring_params(
      pilot_data, time_var, status_var, arm_var, dep_cens_status_var, linear_terms, L
   )

   # 2) Iterative search
   cat("--- Searching for Sample Size (Method: Analytic) ---\n")
   z_alpha <- stats::qnorm(1 - alpha / 2)
   current_n <- n_start
   search_path <- list()
   final_n <- NA_integer_

   while (current_n <= max_n_per_arm) {
      total_n <- 2 * current_n
      se_final <- params$se_beta_n1 / sqrt(total_n)
      calculated_power <- stats::pnorm((abs(params$beta_effect) / se_final) - z_alpha)
      if (!is.finite(calculated_power)) calculated_power <- 0

      search_path[[as.character(current_n)]] <- calculated_power
      cat(sprintf("  N = %d/arm, Calculated Power = %.3f\n", current_n, calculated_power))

      if (calculated_power >= target_power) { final_n <- current_n; break }
      current_n <- current_n + n_step
   }

   if (is.na(final_n)) {
      warning(paste("Target power", target_power, "not achieved by max N of", max_n_per_arm), call. = FALSE)
      final_n <- max_n_per_arm
   }

   # 3) Finalize + return
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
         title = "Analytic Sample Size Search (RMST; Covariate-Dependent Censoring)",
         subtitle = "IPCW-based; single censoring mechanism.",
         x = "Sample Size Per Arm", y = "Calculated Power"
      ) +
      ggplot2::theme_minimal()

   cat("\n--- Calculation Summary ---\n")
   print(knitr::kable(results_df, caption = "Required Sample Size"))

   results_summary <- data.frame(
      Statistic = "Assumed RMST Difference (from pilot)",
      Value = params$beta_effect
   )

   list(results_data = results_df, results_plot = p, results_summary = results_summary)
}
