#' @importFrom survival survfit survdiff Surv
#' @importFrom stats pchisq
#' @importFrom scales pvalue
#' @importFrom ggplot2 ggplot aes fortify geom_step geom_ribbon labs coord_cartesian annotate theme_minimal facet_wrap
.run_survival_diagnostics <- function(pilot_data, time_var, status_var, arm_var, strata_var = NULL) {
  
  df <- pilot_data
  df[[arm_var]] <- as.factor(df[[arm_var]])
  if (!is.null(strata_var)) {
    df[[strata_var]] <- as.factor(df[[strata_var]])
  }
  
  # --- 1. Perform Log-Rank Test ---
  logrank_summary <- NULL
  p_value <- NULL
  
  if (length(unique(df[[arm_var]])) >= 2) {
    surv_formula_logrank <- if (is.null(strata_var)) {
      as.formula(paste("Surv(", time_var, ",", status_var, ") ~", arm_var))
    } else {
      as.formula(paste("Surv(", time_var, ",", status_var, ") ~", arm_var, "+ strata(", strata_var, ")"))
    }
    logrank_test <- survdiff(surv_formula_logrank, data = df)
    p_value <- pchisq(logrank_test$chisq, length(logrank_test$n) - 1, lower.tail = FALSE)
    
    logrank_summary <- data.frame(
      Statistic = "Log-Rank Test Chi-Square", Value = round(logrank_test$chisq, 3)
    )
    logrank_summary <- rbind(logrank_summary, data.frame(Statistic = "Degrees of Freedom", Value = length(logrank_test$n) - 1))
    logrank_summary <- rbind(logrank_summary, data.frame(Statistic = "P-Value", Value = scales::pvalue(p_value)))
  } else {
    logrank_summary <- data.frame(
      Statistic = "Log-Rank Test Status", Value = "Not Applicable (only one treatment arm present in data)"
    )
  }
  
  # --- 2. Generate Kaplan-Meier Plot ---
  surv_formula <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~", arm_var))
  fit <- survfit(surv_formula, data = df)
  fit_fortified <- ggplot2::fortify(fit)
  
  km_plot <- ggplot(fit_fortified, aes(x = time, y = surv, color = strata, fill = strata)) +
    geom_step(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
    labs(
      title = "Kaplan-Meier Curve by Treatment Arm",
      x = "Time", y = "Survival Probability", color = "Arm", fill = "Arm"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal()
  
  if (!is.null(p_value) && is.null(strata_var)) {
    km_plot <- km_plot + annotate("text", x = 0, y = 0, hjust = -0.1, vjust = -0.5, label = paste("Log-Rank p =", scales::pvalue(p_value)))
  }
  
  # If a strata variable is present, create a faceted plot
  if (!is.null(strata_var)) {
    km_plot <- km_plot + facet_wrap(as.formula(paste("~", strata_var))) +
      labs(title = "Kaplan-Meier Curves by Stratum and Treatment Arm")
  }
  
  return(list(
    km_plot = km_plot,
    logrank_summary = logrank_summary
  ))
}