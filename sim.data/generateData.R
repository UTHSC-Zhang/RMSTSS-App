# Load necessary libraries
library(dplyr)

# --- 1. Dataset for Linear IPCW Model ---
# A standard two-arm trial with one covariate.
set.seed(123)
n_linear <- 200
linear_pilot_data <- data.frame(
  id = 1:n_linear,
  age = rnorm(n_linear, mean = 60, sd = 10),
  arm = rep(0:1, each = n_linear / 2)
)

# Generate survival time with a treatment effect
base_rate <- 0.1
treatment_effect <- -0.04 # Additive effect on log-hazard
true_time <- rexp(n_linear, rate = base_rate * exp(linear_pilot_data$arm * treatment_effect))

# Generate independent administrative censoring
censoring_time <- runif(n_linear, min = 10, max = 30)

# Final observed data
linear_pilot_data$time <- pmin(true_time, censoring_time)
linear_pilot_data$status <- as.numeric(true_time <= censoring_time)

# Save the dataset
write.csv(linear_pilot_data, "sim.data/linear_pilot_data.csv", row.names = FALSE)
cat("1. 'linear_pilot_data.csv' created successfully.\n")


# --- 2. Dataset for Stratified Models (Additive/Multiplicative) ---
# A multi-center trial with three regions (strata).
set.seed(456)
n_strat <- 300
stratified_pilot_data <- data.frame(
  id = 1:n_strat,
  region = factor(rep(c("North", "South", "West"), each = n_strat / 3)),
  arm = rep(0:1, each = n_strat / 2)
)

# Generate survival time with a treatment effect
# The base rate varies by region
region_effect <- c("North" = 0.1, "South" = 0.15, "West" = 0.08)
base_rate_strat <- region_effect[stratified_pilot_data$region]
treatment_effect_strat <- -0.05
true_time_strat <- rexp(n_strat, rate = base_rate_strat * exp(stratified_pilot_data$arm * treatment_effect_strat))

# Generate censoring
censoring_time_strat <- runif(n_strat, min = 15, max = 40)

# Final observed data
stratified_pilot_data$time <- pmin(true_time_strat, censoring_time_strat)
stratified_pilot_data$status <- as.numeric(true_time_strat <= censoring_time_strat)

# Save the dataset
write.csv(stratified_pilot_data, "sim.data/stratified_pilot_data.csv", row.names = FALSE)
cat("2. 'stratified_pilot_data.csv' created successfully.\n")


# --- 3. Dataset for Semiparametric (GAM) Model ---
# A trial where a biomarker has a non-linear effect on survival.
set.seed(789)
n_gam <- 250
gam_pilot_data <- data.frame(
  id = 1:n_gam,
  biomarker = rnorm(n_gam, mean = 50, sd = 25),
  arm = rep(0:1, each = n_gam / 2)
)
# Ensure biomarker is non-negative
gam_pilot_data$biomarker[gam_pilot_data$biomarker < 0] <- 0

# Generate survival time with a non-linear effect from the biomarker
# Effect is U-shaped: worse survival at low and high biomarker levels
non_linear_effect <- -0.0001 * (gam_pilot_data$biomarker - 50)^2
treatment_effect_gam <- -0.3
true_time_gam <- rexp(n_gam, rate = 0.05 * exp(gam_pilot_data$arm * treatment_effect_gam + non_linear_effect))

# Generate censoring
censoring_time_gam <- runif(n_gam, min = 20, max = 50)

# Final observed data
gam_pilot_data$time <- pmin(true_time_gam, censoring_time_gam)
gam_pilot_data$status <- as.numeric(true_time_gam <= censoring_time_gam)

# Save the dataset
write.csv(gam_pilot_data, "sim.data/gam_pilot_data.csv", row.names = FALSE)
cat("3. 'gam_pilot_data.csv' created successfully.\n")


# --- 4. Dataset for Dependent Censoring Model ---
# A trial where there is a primary event (e.g., disease progression)
# and a competing, dependent event (e.g., non-disease death).
set.seed(101)
n_dep <- 400
dep_cens_pilot_data <- data.frame(
  id = 1:n_dep,
  age = rnorm(n_dep, mean = 65, sd = 8),
  arm = rep(0:1, each = n_dep / 2)
)

# Generate three potential event times:
# T1: Time to primary event (e.g., progression)
# T2: Time to dependent event (e.g., death)
# T3: Time to administrative censoring
treatment_effect_dc <- -0.4 # Reduces hazard of primary event
T1 <- rexp(n_dep, rate = 0.08 * exp(dep_cens_pilot_data$arm * treatment_effect_dc))
T2 <- rexp(n_dep, rate = 0.05 * exp(dep_cens_pilot_data$age * 0.01)) # Death hazard increases with age
T3 <- runif(n_dep, min = 24, max = 60) # Administrative censoring at 2-5 years

# Determine the observed time and event type
final_time <- pmin(T1, T2, T3)
dep_cens_pilot_data$time <- final_time
# status: 1 for primary event, 0 otherwise
dep_cens_pilot_data$status <- as.numeric(final_time == T1)
# dep_cens_status: 1 for dependent event, 0 otherwise
dep_cens_pilot_data$dep_cens_status <- as.numeric(final_time == T2)

# Save the dataset
write.csv(dep_cens_pilot_data, "sim.data/dep_cens_pilot_data.csv", row.names = FALSE)
cat("4. 'dep_cens_pilot_data.csv' created successfully.\n")