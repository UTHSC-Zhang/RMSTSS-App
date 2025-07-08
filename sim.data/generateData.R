# Load library for data manipulation
library(dplyr)

# --- Dataset 1: Standard Survival Data for Linear IPCW Model ---
set.seed(42)
n1 <- 200
linear_ipcw_data <- data.frame(
  id = 1:n1,
  time = rexp(n1, rate = 0.08),
  status = rbinom(n1, 1, 0.75),
  arm = rep(0:1, each = n1 / 2),
  age = round(rnorm(n1, 65, 8)),
  biomarker = rnorm(n1, 1.5, 0.4),
  sex = sample(c("Male", "Female"), n1, replace = TRUE, prob = c(0.55, 0.45)),
  performance_status = sample(c("Good", "Poor"), n1, replace = TRUE, prob = c(0.7, 0.3))
)

# Introduce a treatment effect (additive)
linear_ipcw_data <- linear_ipcw_data %>%
  mutate(time = if_else(arm == 1, time + 5, time))

write.csv(linear_ipcw_data, "sim.data/linear_ipcw_data.csv", row.names = FALSE)


# --- Dataset 2: Stratified Survival Data ---
set.seed(123)
n2 <- 240
stratified_data <- data.frame(
  id = 1:n2,
  time = rexp(n2, rate = 0.1),
  status = rbinom(n2, 1, 0.8),
  arm = rep(0:1, each = n2 / 2),
  age = round(rnorm(n2, 68, 10)),
  comorbidities = rpois(n2, 1.2),
  prior_therapy = sample(c("Yes", "No"), n2, replace = TRUE, prob = c(0.4, 0.6)),
  region = rep(c("North America", "Europe", "Asia"), each = n2 / 3)
)

# Introduce a region-dependent treatment effect (multiplicative)
stratified_data <- stratified_data %>%
  mutate(
    effect_multiplier = case_when(
      arm == 1 & region == "North America" ~ 1.5,
      arm == 1 & region == "Europe" ~ 1.3,
      arm == 1 & region == "Asia" ~ 1.2,
      TRUE ~ 1.0
    ),
    time = time * effect_multiplier
  ) %>%
  select(-effect_multiplier)

write.csv(stratified_data, "sim.data/stratified_data.csv", row.names = FALSE)


# --- Dataset 3: Data with Dependent Censoring (Competing Risks) ---
set.seed(789)
n3 <- 300
# Simulate times for three potential events
time_to_event <- rexp(n3, rate = 0.05)
time_to_dep_cens <- rexp(n3, rate = 0.02)
time_to_admin_cens <- runif(n3, 50, 80) # Administrative censoring time

dependent_censoring_data <- data.frame(
  id = 1:n3,
  arm = rep(0:1, each = n3 / 2),
  age = round(rnorm(n3, 62, 5)),
  frailty_score = rnorm(n3, 5, 1.5),
  treatment_line = sample(c("First", "Second"), n3, replace = TRUE),
  ecog_score = sample(0:2, n3, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
  # Observed time is the first event to occur
  time = pmin(time_to_event, time_to_dep_cens, time_to_admin_cens),
  # Status for the primary event
  status = as.integer(time_to_event <= pmin(time_to_dep_cens, time_to_admin_cens)),
  # Status for the dependent censoring event
  dep_cens_status = as.integer(time_to_dep_cens < pmin(time_to_event, time_to_admin_cens))
)

write.csv(dependent_censoring_data, "sim.data/dependent_censoring_data.csv", row.names = FALSE)

cat("Successfully generated linear_ipcw_data.csv, stratified_data.csv, and dependent_censoring_data.csv\n")