# app.R — RMSTpowerBoost (single file, no external Rmd required)

# ------------------ Packages ------------------
packages <- c(
  "shiny","shinyjs","bslib","DT","ggplot2","plotly","survival","survminer",
  "kableExtra","magrittr","rmarkdown","dplyr","tidyr","purrr","stringr","tibble"
)
invisible(lapply(packages, function(pkg){
  # if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE)
}))

# ------------------ Source your R/ scripts (simulation engine etc.) ------------------
if (dir.exists("R")) {
  r_files <- list.files("R", pattern="\\.R$", full.names = TRUE)
  sapply(r_files, source)
  cat("All R scripts in the 'R/' directory have been sourced.\n")
} else {
  cat("No R/ directory found; make sure simulation helpers are available.\n")
}

`%||%` <- function(x,y) if (is.null(x)) y else x

# ------------------ Helpers ------------------
reset_cov_builder_ui <- function(session) {
  updateTextInput(session, "cov_name", value = "")
  updateSelectInput(session, "cov_type", selected = "continuous")
  updateSelectInput(session, "cont_dist", selected = "normal")
  updateNumericInput(session, "cont_beta", value = 0)
  updateNumericInput(session, "tf_center", value = 0)
  updateNumericInput(session, "tf_scale", value = 1)
}

reset_cat_entry_ui <- function(session) {
  updateTextInput(session, "cat_add_name", value = "")
  updateNumericInput(session, "cat_add_prob", value = NA)
  updateNumericInput(session, "cat_add_coef", value = 0)
}




# --- Transparency helpers ---
transparent_theme <- function(base_size = 11) {
  theme_light(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = NA, color = NA),
      plot.background  = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA),
      legend.box.background = element_rect(fill = NA, color = NA)
    )
}

to_plotly_clear <- function(p) {
  plotly::ggplotly(p) %>%
    plotly::layout(
      paper_bgcolor = 'rgba(0,0,0,0)',
      plot_bgcolor  = 'rgba(0,0,0,0)'
    )
}
# Put a title above every subplot domain (for plotly::subplot)
add_subplot_titles <- function(p, titles, vgap = 0.04, fontsize = 14){
  xaxs <- grep("^xaxis", names(p$x$layout), value = TRUE)
  yaxs <- grep("^yaxis", names(p$x$layout), value = TRUE)
  anns <- vector("list", length(titles))
  for (i in seq_along(titles)) {
    xd <- p$x$layout[[xaxs[i]]]$domain
    yd <- p$x$layout[[yaxs[i]]]$domain
    anns[[i]] <- list(
      text = titles[i],
      x = mean(xd), y = yd[2] + vgap,
      xref = "paper", yref = "paper",
      xanchor = "center", yanchor = "bottom",
      showarrow = FALSE,
      font = list(size = fontsize)
    )
  }
  plotly::layout(p, annotations = c(p$x$layout$annotations, anns))
}

DT_25 <- function(df) DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)

as_factor_safe <- function(x) { if (is.factor(x)) x else factor(x) }

# Summaries that never error
covariate_summary <- function(df, arm_var = NULL) {
  if (is.null(df) || !nrow(df)) return(list())
  ignore <- c("time","status",arm_var)
  vars <- setdiff(names(df), ignore)
  if (!length(vars)) return(list())
  is_num <- sapply(df[vars], is.numeric)
  cont_vars <- vars[is_num]
  cat_vars  <- vars[!is_num]
  out <- list()
  if (length(cont_vars)) {
    cont_tab <- purrr::map_df(cont_vars, function(v) {
      x <- suppressWarnings(as.numeric(df[[v]]))
      tibble::tibble(
        Variable = v,
        Mean = round(mean(x, na.rm = TRUE), 3),
        SD   = round(stats::sd(x, na.rm = TRUE), 3),
        Min  = round(min(x, na.rm = TRUE), 3),
        Q1   = round(stats::quantile(x, 0.25, na.rm = TRUE), 3),
        Median  = round(stats::median(x, na.rm = TRUE), 3),
        Q3   = round(stats::quantile(x, 0.75, na.rm = TRUE), 3),
        Max  = round(max(x, na.rm = TRUE), 3),
        N_Missing = sum(!is.finite(x))
      )
    })
    out$continuous <- cont_tab
  }
  if (length(cat_vars)) {
    cat_tab <- purrr::map_df(cat_vars, function(v) {
      x <- as_factor_safe(df[[v]])
      if (all(is.na(x))) {
        tibble::tibble(Variable = v, Level = NA_character_, Count = 0, Percent = 0)
      } else {
        tt <- as.data.frame(table(x, useNA = "ifany"))
        names(tt) <- c("Level","Count")
        tt$Percent <- round(100 * tt$Count / sum(tt$Count), 1)
        tt$Variable <- v
        tt[, c("Variable","Level","Count","Percent")]
      }
    })
    out$categorical <- cat_tab
  }
  out
}

# ---------- THEMED PLOTTING HELPERS ----------

# Covariate plots that accept a palette from the app theme
covariate_plots <- function(df, arm_var = NULL, palette = c("#0d6efd","#dc3545","#198754","#0dcaf0")) {
  if (is.null(df) || !nrow(df)) return(list())
  vars <- setdiff(names(df), c("time","status",arm_var))
  if (!length(vars)) return(list())
  plots <- list()
  for (v in vars) {
    x <- df[[v]]
    if (is.numeric(x)) {
      p <- ggplot(df, aes(x = .data[[v]])) +
        geom_histogram(bins = 30, alpha = 0.9, fill = palette[1], color = NA) +
        labs(title = paste("Histogram of", v), x = v, y = "Count") +
        transparent_theme()
      plots[[v]] <- to_plotly_clear(p)
    } else {
      levs <- unique(as.factor(df[[v]]))
      p <- ggplot(df, aes(x = as.factor(.data[[v]]), fill = as.factor(.data[[v]]))) +
        geom_bar(alpha = 0.9, color = NA) +
        scale_fill_manual(values = rep(palette, length.out = length(levs))) +
        guides(fill = "none") +
        labs(title = paste("Bar chart of", v), x = v, y = "Count") +
        transparent_theme()
      plots[[v]] <- to_plotly_clear(p)
    }
  }
  plots
}



# Build a model.matrix column order from covariate definitions
build_mm_columns <- function(cov_defs, include_intercept = TRUE) {
  if (!length(cov_defs)) return(character(0))
  
  # Determine the longest categorical level count
  rows <- 1
  cols <- list()
  
  for (d in cov_defs) {
    if (d$type == "continuous") {
      cols[[d$name]] <- 0
    } else if (d$type == "categorical") {
      lv <- d$params$labels
      if (is.null(lv) || !length(lv)) {
        k <- length(d$params$prob %||% c(0,1))
        lv <- paste0(d$name, seq_len(k))
      }
      cols[[d$name]] <- factor(lv, levels = lv)
      rows <- max(rows, length(lv))
    }
  }
  
  # Rep every column to the SAME length
  df <- as.data.frame(
    lapply(cols, function(x) rep(x, length.out = max(1, rows))),
    stringsAsFactors = FALSE
  )
  
  form <- as.formula(paste0(
    if (include_intercept) "~ 1 +" else "~ -1 +",
    paste(vapply(cov_defs, function(d) d$name, character(1)), collapse = " + ")
  ))
  
  mm <- model.matrix(form, data = df)
  colnames(mm)
}


# Coeff vector in model.matrix column order
assemble_beta <- function(cov_defs, user_betas, include_intercept = TRUE, intercept_value = 0) {
  cols <- build_mm_columns(cov_defs, include_intercept = include_intercept)
  out <- numeric(0)
  if (include_intercept) out <- c(out, user_betas[["(Intercept)"]] %||% intercept_value)
  for (d in cov_defs) {
    if (d$type == "continuous") {
      b <- as.numeric(user_betas[[d$name]])
      if (length(b) != 1 || !is.finite(b)) stop("Coefficient for continuous '", d$name, "' must be a single finite number.")
      out <- c(out, b)
    } else if (d$type == "categorical") {
      lv <- d$params$labels
      if (is.null(lv) || !length(lv)) {
        k <- length(d$params$prob)
        lv <- paste0(d$name, seq_len(k))
      }
      if (include_intercept) {
        need <- length(lv) - 1
        b <- as.numeric(user_betas[[d$name]])
        if (length(b) < need) stop("Categorical '", d$name, "' requires ", need, " coefficients (intercept included).")
        out <- c(out, b[seq_len(need)])
      } else {
        need <- length(lv)
        b <- as.numeric(user_betas[[d$name]])
        if (length(b) < need) stop("Categorical '", d$name, "' requires ", need, " coefficients (no intercept).")
        out <- c(out, b[seq_len(need)])
      }
    }
  }
  out
}

# ------------------ REPORT TEMPLATE (NEW SECTIONING) ------------------
make_inline_template <- function() {
  tf <- tempfile(fileext = ".Rmd")
  txt <- c(
    "---",
    "title: \"RMSTpowerBoost Report\"",
    "output: pdf_document",
    "params:",
    "  inputs: NA",
    "  results: NA",
    "  log: NA",
    "  data_provenance: NA",
    "  data: NA",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "suppressPackageStartupMessages({",
    "  library(ggplot2); library(survival); library(survminer);",
    "  library(kableExtra); library(dplyr); library(tidyr); library(tibble)",
    "})",
    "to_kv <- function(x){",
    "  if (is.null(x)) return(\"\")",
    "  keys <- names(x); if (is.null(keys)) return(paste(x, collapse=\", \"))",
    "  paste(paste(keys, sapply(x, function(z){",
    "    if (is.null(z)) return(\"NULL\")",
    "    if (is.atomic(z)) return(paste(z, collapse=\", \"))",
    "    if (is.function(z)) return(\"<function>\")",
    "    if (is.list(z)) return(paste(unlist(z), collapse=\", \"))",
    "    as.character(z)",
    "  }), sep = \" = \"), collapse = \"; \")",
    "}",
    "or_else <- function(x, y) if (is.null(x)) y else x",
    "```",
    "",
    "# 1) Data Generating Mechanism",
    "",
    "## 1.1) Covariate generation",
    "```{r covariate-generation, echo=FALSE}",
    "cov_list <- or_else(params$data_provenance$Covariates_defined, list())",
    "if (length(cov_list)) {",
    "  cov_tbl <- do.call(rbind, lapply(cov_list, function(d) {",
    "    data.frame(",
    "      Variable     = or_else(d$name, \"\"),",
    "      Type         = or_else(d$type, \"\"),",
    "      Distribution = or_else(d$dist, \"\"),",
    "      Parameters   = to_kv(d$params),",
    "      Transform    = paste(or_else(d$transform, character(0)), collapse = \", \"),",
    "      stringsAsFactors = FALSE",
    "    )",
    "  }))",
    "  kbl(cov_tbl, booktabs = TRUE, caption = \"Covariate definitions\", row.names = FALSE) %>%",
    "    kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "} else {",
    "  cat(\"No covariates were defined or the data were uploaded.\\n\")",
    "}",
    "```",
    "",
    "## 1.2) Event time and censoring",
    "```{r event-time-censoring, echo=FALSE}",
    "prov <- params$data_provenance",
    "et <- or_else(prov$Event_time, list())",
    "cz <- or_else(prov$Censoring, list())",
    "tbl <- data.frame(",
    "  Field = c(\"Event-time model\", \"Baseline parameters\", \"Censoring mode\", \"Target overall censoring\", \"Administrative time\"),",
    "  Value = c(",
    "    as.character(or_else(et$model, or_else(et$Model, \"\"))),",
    "    to_kv(or_else(et$baseline, or_else(et$Baseline, list()))),",
    "    as.character(or_else(cz$mode, or_else(cz$Mode, \"\"))),",
    "    as.character(or_else(cz$target, or_else(cz$Target, \"\"))),",
    "    as.character(or_else(cz$admin_time, or_else(cz$Admin_time, \"\")))",
    "  ), stringsAsFactors = FALSE",
    ")",
    "kbl(tbl, booktabs = TRUE, caption = \"Event-time and censoring\", row.names = FALSE) %>%",
    "  kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "```",
    "",
    "# 2) Analysis Configuration",
    "```{r input-parameters, echo=FALSE}",
    "inputs_df <- data.frame(",
    "  Parameter = names(params$inputs),",
    "  Value     = unlist(lapply(params$inputs, function(x) paste(x, collapse = \", \"))),",
    "  stringsAsFactors = FALSE",
    ")",
    "kbl(inputs_df, booktabs = TRUE, caption = \"Input Parameters\",",
    "    row.names = FALSE, col.names = c(\"Parameter\",\"Value\")) %>%",
    "  kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "```",
    "",
    "# 3) Analysis Results",
    "",
    "## 3.1) Survival Analysis",
    "```{r km-plot, eval = !is.null(params$results$km_note)}",
    "cat(params$results$km_note)",
    "```",
    "```{r logrank-summary, echo=FALSE, eval = !is.null(params$results$logrank_summary)}",
    "cap <- if (!is.null(params$inputs$strata_var) && nzchar(params$inputs$strata_var)) {",
    "  paste0('Stratified log–rank test results (strata: ', params$inputs$strata_var, ')')",
    "} else {",
    "  'Log–rank test results'",
    "}",
    "kbl(params$results$logrank_summary, booktabs = TRUE,",
    "    caption = cap, row.names = FALSE) %>%",
    "  kable_styling(latex_options = 'hold_position')",
    "```",
    "",
    "## 3.2) Power and Sample Size",
    "```{r power-curve, eval = !is.null(params$results$results_plot)}",
    "params$results$results_plot",
    "```",
    "```{r results-table, eval = !is.null(params$results$results_data)}",
    "kbl(params$results$results_data, booktabs = TRUE,",
    "    caption = \"Power and sample size results\", row.names = FALSE) %>%",
    "  kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "```",
    "```{r effect-size, eval = !is.null(params$results$results_summary)}",
    "kbl(params$results$results_summary, booktabs = TRUE,",
    "    caption = \"Summary measures derived from the data\", row.names = FALSE) %>%",
    "  kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "```",
    "",
    "# 4) Console Output",
    "```{r console-log, results='asis', echo=FALSE, eval = !is.null(params$log)}",
    "cat(params$log)",
    "```"
  )
  writeLines(txt, tf)
  tf
}

report_inputs_builder <- function(input) {
  list(
    model_selection = input$model_selection,
    analysis_type   = input$analysis_type,
    time_var        = input$time_var,
    status_var      = input$status_var,
    arm_var         = input$arm_var,
    strata_var      = input$strata_var %||% "",
    dep_cens_var    = input$dep_cens_var %||% "",
    calc_method     = input$calc_method %||% "Analytical",
    L               = input$L,
    alpha           = input$alpha,
    sample_sizes    = input$sample_sizes,
    target_power    = input$target_power
  )
}

# ------------------ UI ------------------
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  useShinyjs(),
  titlePanel("RMSTpowerBoost: Power and Sample Size Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      # Step 1: Data (Upload or Generate)
      wellPanel(
        h4("Step 1. Upload/Generate Data"),
        radioButtons("data_mode", "Choose data source:", choices = c("Upload", "Generate"), inline = TRUE),
        shinyjs::hidden(div(id = "upload_panel", fileInput("pilot_data_upload", "Upload Pilot Data (.csv)", accept = ".csv"))),
        shinyjs::hidden(div(
          id = "simulate_panel",
          h5("1a. Covariate Builder"),
          # Row 1: name/type
          fluidRow(
            column(6, textInput("cov_name", "Variable name", value = "", placeholder = "x1, x2 … auto if empty")),
            column(6, selectInput("cov_type", "Type", choices = c("continuous","categorical")))
          ),
          # Continuous vs Categorical UI
          uiOutput("cov_details_ui"),
          # Transform (continuous only)
          shinyjs::hidden(div(id = "transform_block",
                              tags$hr(),
                              h5("Transform (continuous only)"),
                              fluidRow(
                                column(6, numericInput("tf_center", "Location (center a)", value = 0)),
                                column(6, numericInput("tf_scale",  "Scale (divide by b)", value = 1, min = 0.0001, step = 0.1))
                              ),
                              helpText("Applied after generation: (x - a) / b")
          )),
          fluidRow(
            column(6, actionButton("add_cov", "Add covariate", icon = icon("plus"), class = "btn btn-success")),
            column(6, actionButton("reset_cov_builder", "Reset builder", icon = icon("trash")))
          ),
          br(), DTOutput("cov_table"),
          div(style="margin-top:8px;",
              actionButton("remove_cov", "Remove selected covariate", icon = icon("minus"))
          ),
          tags$hr(),
          h5("1b. Simulation Settings"),
          fluidRow(
            column(4, numericInput("sim_n", "Sample size", value = 300, min = 10)),
            column(4, textInput("sim_allocation", "Allocation (a:b)", value = "1:1")),
            column(4, numericInput("sim_treat_eff", "Treatment β (arm)", value = -0.2, step = 0.05))
          ),
          fluidRow(
            column(6, checkboxInput("intercept_in_mm", "Include intercept in model.matrix (β0 inside β)", value = TRUE)),
            column(6, numericInput("user_intercept", "β0 (used if no intercept in model.matrix)", value = 0))
          ),
          fluidRow(
            column(6, selectInput("sim_model", "Event-time model",
                                  choices = c("aft_lognormal","aft_weibull","ph_exponential","ph_weibull","ph_pwexp"))),
            column(6, sliderInput("sim_cens", "Target censoring", min = 0, max = 0.9, value = 0.25, step = 0.01))
          ),
          fluidRow(column(12, uiOutput("sim_baseline_ui"))),
          fluidRow(
            column(4, numericInput("sim_seed", "Seed (optional)", value = NA)),
            column(4, actionButton("generate_sim", "Generate Pilot Dataset", icon = icon("gears"), class = "btn btn-primary")),
            column(4, actionButton("reset_generate", "Reset data", icon = icon("trash")))
          )
        ))
      ),
      
      # Step 2+3: Model + Analysis (hidden until data ready)
      shinyjs::hidden(
        wellPanel(
          id = "model_analysis_panel",
          h4("Step 2. Model & Step 3. Analysis"),
          uiOutput("col_mapping_ui"),
          fluidRow(
            column(4, radioButtons("analysis_type", "Target Quantity", choices = c("Power", "Sample Size"), selected = "Power")),
            column(4, selectInput("model_selection", "Select RMST Model",
                                  choices = c("Linear IPCW Model","Additive Stratified Model",
                                              "Multiplicative Stratified Model","Semiparametric (GAM) Model",
                                              "Dependent Censoring Model"),
                                  selected = "Linear IPCW Model")),
            column(4, numericInput("L", "RMST L (τ)", value = 365, min = 1))
          ),
          # NEW: Analytical vs Repeated
          radioButtons("calc_method", "Calculation Method", choices = c("Analytical", "Repeated"), selected = "Analytical", inline = TRUE),
          uiOutput("analysis_inputs_ui"),
          conditionalPanel(
            condition = "input.calc_method == 'Repeated'",
            fluidRow(
              column(6, numericInput("R_reps", "Replications", value = 500, min = 100, step = 100)),
              column(6, numericInput("seed_reps", "Seed (optional)", value = NA))
            )
          ),
          sliderInput("alpha", "Significance Level (α)", min = 0.01, max = 0.1, value = 0.05, step = 0.01),
          tags$hr(),
          fluidRow(
            column(4, actionButton("run_analysis", "Run Analysis", icon = icon("play"), class = "btn-primary btn-lg")),
            column(8, shinyjs::hidden(
              div(id="download_reset_row",
                  downloadButton("download_report_pdf", "Download PDF"),
                  downloadButton("download_report_html", "Download HTML"),
                  actionButton("reset_all", "Reset All", icon = icon("trash"))
              )
            ))
          )
        )
      )
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Instructions",
                 h3("Welcome to RMSTpowerBoost"),
                 p("This tool provides power and sample size analysis using Restricted Mean Survival Time."),
                 tags$ol(
                   tags$li("Step 1: Choose the data source."),
                   tags$ul(
                     tags$li("1a. Upload: Select a CSV file with pilot data."),
                     tags$li("1b. Generate: Build covariates, choose model and censoring, set coefficients and intercept, then simulate.")
                   ),
                   tags$li("Step 2 and Step 3: After data is available, map columns and select analysis settings."),
                   tags$ul(
                     tags$li("2a. Map the time, status, and treatment columns."),
                     tags$li("2b. Choose the target (Power or Sample Size) and analysis parameters."),
                     tags$li("2c. Set the analysis horizon (τ) and the significance level (α).")
                   ),
                   tags$li("Run the analysis; download a report (PDF/HTML).")
                 ),
                 hr(),
                 h4("License Information"),
                 verbatimTextOutput("license_display")
        ),
        tabPanel("Data Preview", DT::dataTableOutput("data_preview_table")),
        tabPanel("Plot Output",
                 h4("Covariate Distributions"),
                 uiOutput("cov_plots_ui"),
                 hr(),
                 h4("Kaplan–Meier Survival Plots"),
                 htmlOutput("km_note"),
                 plotlyOutput("survival_plotly_output", height = "650px"),
                 hr(),
                 h4("Power vs. Sample Size"),
                 plotlyOutput("results_plot", height = "520px")
        ),
        tabPanel("Summary",
                 h4("Analysis Results (Tables Only)"),
                 uiOutput("results_table_ui"),
                 hr(),
                 h4("Analysis Summary"),
                 uiOutput("summary_table_ui"),
                 hr(),
                 h4("Data Summary"),
                 uiOutput("data_summary_ui")
        ),
        tabPanel("Console Log", verbatimTextOutput("console_log_output"))
      )
    )
  )
)

# ------------------ Repeated Power (no 'bootstrap' wording) ------------------
repeated_power_from_pilot <- function(pilot_df, time_var, status_var, arm_var,
                                      n_per_arm_vec, alpha = 0.05, R = 500,
                                      strata_var = NULL, seed = NULL) {
  stopifnot(is.data.frame(pilot_df))
  needed <- c(time_var, status_var, arm_var, strata_var)
  needed <- needed[!is.null(needed)]
  if (!all(needed %in% names(pilot_df))) stop("Required columns not found in pilot data.")
  df <- pilot_df[, unique(c(time_var, status_var, arm_var, strata_var)), drop = FALSE]
  names(df)[match(c(time_var, status_var, arm_var), names(df))] <- c("time","status","arm")
  df$arm <- as.factor(df$arm)
  if (nlevels(df$arm) != 2L) stop("Repeated method currently expects exactly 2 arms.")
  if (!is.null(seed) && is.finite(seed)) set.seed(as.integer(seed))
  
  if (!is.null(strata_var)) {
    names(df)[names(df) == strata_var] <- "stratum"
    df$stratum <- as.factor(df$stratum)
    split_stratum <- split(df, df$stratum, drop = TRUE)
    split_by <- lapply(split_stratum, function(d) split(d, d$arm, drop = TRUE))
  } else {
    split_by <- split(df, df$arm, drop = TRUE)
  }
  
  out <- lapply(n_per_arm_vec, function(n_arm){
    rej <- logical(R)
    for (r in seq_len(R)) {
      if (is.null(strata_var)) {
        s0 <- split_by[[1]][sample.int(nrow(split_by[[1]]), n_arm, replace = TRUE), , drop = FALSE]
        s1 <- split_by[[2]][sample.int(nrow(split_by[[2]]), n_arm, replace = TRUE), , drop = FALSE]
        dat <- rbind(s0, s1)
        fml <- Surv(time, status) ~ arm
      } else {
        blocks <- lapply(split_by, function(by_arm) {
          a0 <- by_arm[[1]]; a1 <- by_arm[[2]]
          s0 <- a0[sample.int(nrow(a0), n_arm, replace = TRUE), , drop = FALSE]
          s1 <- a1[sample.int(nrow(a1), n_arm, replace = TRUE), , drop = FALSE]
          rbind(s0, s1)
        })
        dat <- do.call(rbind, blocks)
        dat$stratum <- factor(dat$stratum)
        fml <- Surv(time, status) ~ arm + strata(stratum)
      }
      dat$arm <- factor(dat$arm)
      lr <- tryCatch(survdiff(fml, data = dat), error = function(e) NULL)
      if (is.null(lr) || is.null(lr$chisq)) {
        rej[r] <- NA
      } else {
        df_chi <- length(lr$n) - 1
        p <- 1 - pchisq(lr$chisq, df = df_chi)
        rej[r] <- is.finite(p) && (p < alpha)
      }
    }
    m <- mean(rej, na.rm = TRUE); k <- sum(is.finite(rej))
    se <- if (k > 0) sqrt(m*(1-m)/k) else NA_real_
    data.frame(N_per_arm = n_arm, Power = m, SE = se, Reps = k)
  })
  do.call(rbind, out)
}

# ------------------ Server ------------------
server <- function(input, output, session) {
  bslib::bs_themer()
  license_content <- tryCatch(paste(readLines("LICENSE"), collapse = "\n"), error = function(e) "LICENSE not found.")
  output$license_display <- renderPrint({ cat(license_content) })
  
  theme_palette <- reactive({
    th <- bslib::bs_current_theme()
    # grab a few bootstrap variables; fall back if missing
    vars <- bslib::bs_get_variables(th, c("primary", "danger", "success", "info"))
    cols <- unname(unlist(vars))
    # ensure we always have at least two colors
    if (length(cols) < 2) cols <- c("#0d6efd", "#dc3545", "#198754", "#0dcaf0")
    cols
  })
  
  rv <- reactiveValues(
    covariates = list(),
    cat_rows = tibble::tibble(cat = character(), prob = numeric(), coef = numeric()),
    data_mode = "Upload",
    data_df = NULL,
    data_source = NULL,
    console_buf = character(0)
  )
  
  # Toggle Upload vs Generate
  observeEvent(input$data_mode, {
    rv$data_mode <- input$data_mode
    shinyjs::toggle(id = "upload_panel", condition = input$data_mode == "Upload")
    shinyjs::toggle(id = "simulate_panel", condition = input$data_mode == "Generate")
  }, ignoreInit = FALSE)
  
  # ---------- Covariate details UI ----------
  output$cov_details_ui <- renderUI({
    if ((input$cov_type %||% "continuous") == "continuous") {
      shinyjs::show("transform_block")
      tagList(
        fluidRow(
          column(6, selectInput("cont_dist", "Distribution", choices = c("normal","lognormal","gamma","weibull","uniform","t","beta"))),
          column(6, numericInput("cont_beta", "Coefficient β", value = 0))
        ),
        uiOutput("cont_param_ui")
      )
    } else {
      shinyjs::hide("transform_block")
      tagList(
        fluidRow(
          column(6, textInput("cat_add_name", "Add category name", placeholder = "auto if blank")),
          column(3, numericInput("cat_add_prob", "Probability", value = NA, min = 0, max = 1, step = 0.01)),
          column(3, numericInput("cat_add_coef", "Coefficient β", value = 0))
        ),
        fluidRow(
          column(4, actionButton("add_cat_row", "Add category", icon=icon("plus"))),
          column(4, actionButton("reset_cat_rows", "Reset categories", icon=icon("trash"))),
          column(4, actionButton("remove_cat_row", "Remove selected", icon=icon("minus")))
        ),
        br(),
        DTOutput("cat_table"),
        helpText("Tip: If you include intercept in model.matrix, only K−1 coefficients are used (last level’s β is ignored).")
      )
    }
  })
  
  # Continuous parameter UI
  output$cont_param_ui <- renderUI({
    switch(input$cont_dist %||% "normal",
           normal = tagList(
             numericInput("p_mean", "mean", value = 0),
             numericInput("p_sd",   "sd", value = 1, min = 0)
           ),
           lognormal = tagList(
             numericInput("p_meanlog", "meanlog", value = 0),
             numericInput("p_sdlog",   "sdlog", value = 1, min = 0)
           ),
           gamma = tagList(
             numericInput("p_shape", "shape", value = 2, min = 0.001),
             numericInput("p_scale", "scale", value = 1, min = 0.0001)
           ),
           weibull = tagList(
             numericInput("p_wshape", "shape", value = 1.5, min = 0.0001),
             numericInput("p_wscale", "scale", value = 1, min = 0.0001)
           ),
           uniform = tagList(
             numericInput("p_min", "min", value = 0),
             numericInput("p_max", "max", value = 1)
           ),
           t = tagList(
             numericInput("p_df", "df", value = 5, min = 1)
           ),
           beta = tagList(
             numericInput("p_shape1", "shape1", value = 2, min = 0.0001),
             numericInput("p_shape2", "shape2", value = 2, min = 0.0001)
           )
    )
  })
  
  # Category rows table & actions
  observeEvent(input$add_cat_row, {
    nm <- input$cat_add_name %||% ""
    if (!nzchar(nm)) nm <- paste0("L", nrow(rv$cat_rows) + 1)
    pr <- input$cat_add_prob
    if (!is.na(pr) && (pr < 0 || pr > 1)) { showNotification("Probability must be between 0 and 1.", type = "warning"); return() }
    cf <- input$cat_add_coef
    if (!is.finite(cf)) { showNotification("Coefficient must be numeric.", type = "warning"); return() }
    
    rv$cat_rows <- dplyr::bind_rows(rv$cat_rows, tibble::tibble(cat = nm, prob = pr, coef = cf))
    
    # NEW: clear the per-row inputs immediately after adding
    reset_cat_entry_ui(session)
  })
  
  observeEvent(input$reset_cat_rows, { rv$cat_rows <- tibble::tibble(cat = character(), prob = numeric(), coef = numeric()) })
  output$cat_table <- renderDT({
    if (!nrow(rv$cat_rows)) {
      DT::datatable(
        data.frame(Message="No categories yet — add rows above."),
        options = list(dom='t'),
        rownames = FALSE,
        selection = "none"
      )
    } else {
      DT::datatable(rv$cat_rows, options = list(pageLength = 25, scrollX = TRUE),
                    rownames = FALSE, selection = "single")
    }
  })
  observeEvent(input$remove_cat_row, {
    idx <- input$cat_table_rows_selected
    if (length(idx) && idx >= 1 && idx <= nrow(rv$cat_rows)) {
      rv$cat_rows <- rv$cat_rows[-idx, , drop = FALSE]
    } else {
      showNotification("Select a category row to remove.", type = "warning")
    }
  })
  
  
  # Transform visibility (continuous only)
  observe({ shinyjs::toggle(id = "transform_block", condition = (input$cov_type %||% "") == "continuous") })
  
  # Reset builder
  observeEvent(input$reset_cov_builder, {
    updateTextInput(session, "cov_name", value = "")
    updateSelectInput(session, "cov_type", selected = "continuous")
    updateSelectInput(session, "cont_dist", selected = "normal")
    updateNumericInput(session, "cont_beta", value = 0)
    updateNumericInput(session, "tf_center", value = 0)
    updateNumericInput(session, "tf_scale", value = 1)
    rv$cat_rows <- tibble::tibble(cat = character(), prob = numeric(), coef = numeric())
    reset_cat_entry_ui(session)   # ← add this
  })
  
  observeEvent(input$reset_cat_rows, {
    rv$cat_rows <- tibble::tibble(cat = character(), prob = numeric(), coef = numeric())
    reset_cat_entry_ui(session)   # ← nice touch
  })
  # Helper: auto covariate name x1, x2…
  next_cov_name <- reactive({
    nm <- input$cov_name
    if (nzchar(nm)) return(nm)
    existing <- vapply(rv$covariates, function(d) d$name, character(1))
    i <- 1
    repeat {
      cand <- paste0("x", i)
      if (!(cand %in% existing)) return(cand)
      i <- i + 1
    }
  })
  
  # Add covariate to list
  observeEvent(input$add_cov, {
    req(input$cov_type)
    vname <- next_cov_name()
    
    if (input$cov_type == "continuous") {
      # ------- CONTINUOUS -------
      pars <- switch(input$cont_dist,
                     normal    = list(mean = input$p_mean, sd = input$p_sd),
                     lognormal = list(meanlog = input$p_meanlog, sdlog = input$p_sdlog),
                     gamma     = list(shape = input$p_shape,   scale = input$p_scale),
                     weibull   = list(shape = input$p_wshape,  scale = input$p_wscale),
                     uniform   = list(min   = input$p_min,     max   = input$p_max),
                     t         = list(df    = input$p_df),
                     beta      = list(shape1 = input$p_shape1, shape2 = input$p_shape2),
                     list())
      tf <- c(sprintf("center(%s)", input$tf_center), sprintf("scale(%s)", input$tf_scale))
      if (!is.finite(as.numeric(input$cont_beta))) {
        showNotification("Continuous coefficient must be a single number.", type = "error"); return()
      }
      
      rv$covariates <- c(rv$covariates, list(list(
        name = vname, type = "continuous", dist = input$cont_dist, params = pars,
        transform = tf, beta = as.numeric(input$cont_beta)
      )))
      
      # >>> RESET UI after adding a continuous covariate
      reset_cov_builder_ui(session)
      rv$cat_rows <- tibble::tibble(cat = character(), prob = numeric(), coef = numeric())
      reset_cat_entry_ui(session)
      
    } else {
      # ------- CATEGORICAL -------
      if (!nrow(rv$cat_rows)) { showNotification("Add at least one category.", type = "error"); return() }
      
      cats <- rv$cat_rows$cat
      prob <- rv$cat_rows$prob
      coef <- rv$cat_rows$coef
      
      # Validate coefs
      if (any(!is.finite(coef))) { showNotification("All category coefficients must be numeric.", type="error"); return() }
      
      # Fill NA probs with equal shares of the remainder
      p_known <- prob
      p_known[is.na(p_known)] <- 0
      rem <- 1 - sum(p_known)
      if (rem < -1e-8) { showNotification("Sum of specified probabilities exceeds 1.", type="error"); return() }
      if (any(is.na(prob))) {
        k_na <- sum(is.na(prob))
        add_each <- if (k_na > 0) rem / k_na else 0
        prob[is.na(prob)] <- add_each
        rem <- 1 - sum(prob)
      }
      
      # If still < 1, append a remainder category
      if (rem > 1e-8) {
        new_name <- paste0("category-", length(cats) + 1)
        cats <- c(cats, new_name)
        prob <- c(prob, rem)
        coef <- c(coef, 0)  # neutral effect for auto-added level
      }
      
      # Special case: only one level entered with prob < 1
      if (length(cats) == 1 && abs(1 - prob[1]) > 1e-8) {
        cats <- c(cats, paste0("category-", 2))
        prob <- c(prob, 1 - prob[1])
        coef <- c(coef, 0)
      }
      
      # Final check
      if (any(prob < 0) || abs(sum(prob) - 1) > 1e-6) {
        showNotification("Category probabilities must be ≥ 0 and sum to 1 (after auto-completion).", type="error"); return()
      }
      
      # Bernoulli vs multiclass
      if (length(cats) == 2) {
        pars <- list(p = prob[2])
        rv$covariates <- c(rv$covariates, list(list(
          name = vname, type = "categorical", dist = "bernoulli",
          params = pars, transform = NULL, beta = coef
        )))
      } else {
        pars <- list(prob = prob, labels = cats)
        rv$covariates <- c(rv$covariates, list(list(
          name = vname, type = "categorical", dist = "categorical",
          params = pars, transform = NULL, beta = coef
        )))
      }
      
      # >>> RESET UI after adding a categorical covariate
      reset_cov_builder_ui(session)
      rv$cat_rows <- tibble::tibble(cat = character(), prob = numeric(), coef = numeric())
      reset_cat_entry_ui(session)
    }
  })
  
  
  # Covariate list table
  output$cov_table <- renderDT({
    if (!length(rv$covariates)) return(DT::datatable(data.frame(), selection = "none"))
    show <- purrr::map_df(rv$covariates, function(d) {
      beta_txt <- if (length(d$beta)>1) paste(d$beta, collapse=", ") else as.character(d$beta)
      tibble::tibble(
        name = d$name, type = d$type, dist = d$dist,
        params = paste(names(d$params), unlist(d$params), sep="=", collapse="; "),
        transform = paste(d$transform %||% character(0), collapse = ", "),
        beta = beta_txt
      )
    })
    DT::datatable(show, options = list(pageLength = 25, scrollX = TRUE),
                  rownames = FALSE, selection = "single")
  })
  # Remove selected covariate
  observeEvent(input$remove_cov, {
    idx <- input$cov_table_rows_selected
    if (length(idx) && idx >= 1) {
      # remove the idx-th covariate from the list
      rv$covariates <- rv$covariates[-idx]
    } else {
      showNotification("Select a covariate to remove.", type = "warning")
    }
  })
  
  
  # Baseline UI (grouped)
  output$sim_baseline_ui <- renderUI({
    pad <- function(x) div(style = "display:inline-block; margin-right:12px;", x)
    switch(input$sim_model %||% "aft_lognormal",
           aft_lognormal = div(pad(numericInput("b_mu", "mu", value = 2.3)),
                               pad(numericInput("b_sigma", "sigma", value = 0.5, min = 0.0001))),
           aft_weibull   = div(pad(numericInput("b_shape", "shape", value = 1.5, min = 0.0001)),
                               pad(numericInput("b_scale", "scale", value = 10,  min = 0.0001))),
           ph_exponential= div(pad(numericInput("b_rate", "rate", value = 0.05, min = 0.000001))),
           ph_weibull    = div(pad(numericInput("b_wshape2", "shape", value = 1.3, min = 0.0001)),
                               pad(numericInput("b_wscale2", "scale", value = 8,   min = 0.0001))),
           ph_pwexp      = div(pad(textInput("b_rates", "rates (comma)", value = "0.05,0.02")),
                               pad(textInput("b_cuts",  "cuts (comma)",  value = "5")))
    )
  })
  
  # Upload
  observeEvent(input$pilot_data_upload, {
    req(input$pilot_data_upload)
    df <- tryCatch(read.csv(input$pilot_data_upload$datapath, check.names = FALSE), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) { showNotification("Error reading CSV or empty data.", type = "error"); return() }
    rv$data_df <- df
    rv$data_source <- "uploaded"
    shinyjs::show(id = "model_analysis_panel")
    updateTabsetPanel(session, "main_tabs", selected = "Data Preview")
  })
  
  # Reset data (generation card)
  observeEvent(input$reset_generate, {
    rv$data_df <- NULL; rv$data_source <- NULL
    showNotification("Data reset.", type="message")
  })
  
  # Generate
  observeEvent(input$generate_sim, {
    if (!length(rv$covariates)) { showNotification("Please add at least one covariate before simulating.", type = "warning"); return() }
    include_intercept <- isTRUE(input$intercept_in_mm)
    cov_defs <- lapply(rv$covariates, function(d){
      list(name = d$name, type = d$type, dist = d$dist, params = d$params, transform = d$transform)
    })
    user_betas <- list()
    if (include_intercept) user_betas[["(Intercept)"]] <- 0
    for (d in rv$covariates) {
      if (d$type == "continuous") user_betas[[d$name]] <- as.numeric(d$beta) else {
        lv <- if (d$dist == "bernoulli") c("0","1") else (d$params$labels %||% paste0(d$name, seq_along(d$params$prob)))
        need <- if (include_intercept) length(lv)-1 else length(lv)
        if (length(d$beta) < need) { showNotification(sprintf("'%s' needs %d coefficients but %d provided.", d$name, need, length(d$beta)), type="error"); return() }
        user_betas[[d$name]] <- as.numeric(d$beta)
      }
    }
    beta_vec <- assemble_beta(cov_defs, user_betas, include_intercept = include_intercept, intercept_value = input$user_intercept)
    mm_cols  <- build_mm_columns(cov_defs, include_intercept = include_intercept)
    baseline <- switch(input$sim_model,
                       "aft_lognormal" = list(mu = input$b_mu, sigma = input$b_sigma),
                       "aft_weibull"   = list(shape = input$b_shape, scale = input$b_scale),
                       "ph_exponential"= list(rate = input$b_rate),
                       "ph_weibull"    = list(shape = input$b_wshape2, scale = input$b_wscale2),
                       "ph_pwexp"      = {
                         rates <- suppressWarnings(as.numeric(trimws(strsplit(input$b_rates, ",")[[1]])))
                         cuts  <- trimws(strsplit(input$b_cuts, ",")[[1]])
                         cuts  <- if (length(cuts) == 1 && cuts == "") numeric(0) else suppressWarnings(as.numeric(cuts))
                         list(rates = rates, cuts = cuts)
                       })
    form <- as.formula(paste0(if (include_intercept) "~ 1 +" else "~ -1 +",
                              paste(vapply(cov_defs, function(d) d$name, character(1)), collapse = " + ")))
    effects_list <- list(
      intercept = if (include_intercept) 0 else input$user_intercept,
      treatment = input$sim_treat_eff,
      formula   = deparse(form),
      beta      = beta_vec
    )
    rec <- list(
      n = as.integer(input$sim_n),
      covariates = list(defs = cov_defs),
      treatment = list(assignment = "randomization", allocation = input$sim_allocation),
      event_time = list(model = input$sim_model, baseline = baseline, effects = effects_list),
      censoring = list(mode = "target_overall", target = input$sim_cens, admin_time = Inf),
      seed = if (is.na(input$sim_seed)) NULL else as.integer(input$sim_seed)
    )
    buf <- capture.output({
      cat("---- Data Simulation ----\n"); print(str(rec))
      cat("Model matrix columns order:\n"); print(mm_cols)
      cat("Beta vector:\n"); print(beta_vec)
    })
    rv$console_buf <- c(rv$console_buf, buf)
    dat <- tryCatch({ simulate_from_recipe(rec, seed = rec$seed) }, error = function(e) { showNotification(paste("Simulation failed:", e$message), type = "error"); NULL })
    if (is.null(dat)) return()
    rv$data_df <- dat
    rv$data_source <- "simulated"
    showNotification("Simulation complete.", type = "message")
    shinyjs::hide(id = "simulate_panel")
    shinyjs::show(id = "model_analysis_panel")
    updateTabsetPanel(session, "main_tabs", selected = "Data Preview")
    rv$provenance <- list(
      Source = "simulated",
      Number_of_rows = nrow(dat),
      Variable_names = names(dat),
      Covariates_defined = lapply(rv$covariates, function(d) d),
      Event_time = list(model = input$sim_model, baseline = baseline),
      Treatment  = list(assignment = "randomization", allocation = input$sim_allocation),
      Effects    = list(treatment = input$sim_treat_eff,
                        intercept_report = if (include_intercept) "(in model.matrix β)" else input$user_intercept,
                        intercept_in_mm = include_intercept,
                        formula = deparse(form),
                        mm_cols = mm_cols),
      Censoring  = list(mode = "target_overall", target = input$sim_cens, admin_time = Inf)
    )
  })
  
  # Column mapping (now includes strata & dep. censoring status when needed)
  output$col_mapping_ui <- renderUI({
    df <- rv$data_df; req(df)
    cn <- names(df)
    needs_strata <- (input$model_selection %in% c("Additive Stratified Model","Multiplicative Stratified Model"))
    needs_dep    <- (input$model_selection %in% c("Dependent Censoring Model"))
    tagList(
      fluidRow(
        column(4, selectInput("time_var", "Time-to-Event", choices = cn,
                              selected = if ("time" %in% cn) "time" else cn[1])),
        column(4, selectInput("status_var", "Status (1=event)", choices = cn,
                              selected = if ("status" %in% cn) "status" else cn[min(2,length(cn))])),
        column(4, selectInput("arm_var", "Treatment Arm (1=treat)", choices = cn,
                              selected = if ("arm" %in% cn) "arm" else cn[min(3,length(cn))]))
      ),
      if (needs_strata) {
        fluidRow(
          column(6, selectInput("strata_var", "Stratification Variable", choices = cn,
                                selected = setdiff(cn, c(input$time_var, input$status_var, input$arm_var))[1]))
        )
      },
      if (needs_dep) {
        fluidRow(
          column(6, selectInput("dep_cens_var", "Dependent Censoring Status", choices = cn,
                                selected = setdiff(cn, c(input$time_var, input$status_var, input$arm_var, input$strata_var))[1]))
        )
      }
    )
  })
  
  # Analysis inputs
  output$analysis_inputs_ui <- renderUI({
    if (input$analysis_type == "Power") {
      textInput("sample_sizes", "Sample Sizes (per arm/stratum, comma-separated)", value = "100, 150, 200")
    } else {
      sliderInput("target_power", "Target Power", min = 0.1, max = 1, value = 0.8, step = 0.01)
    }
  })
  
  # Data Preview
  output$data_preview_table <- DT::renderDataTable({ req(rv$data_df); DT_25(rv$data_df) })
  
  # Covariate plots
  # UI container
  output$cov_plots_ui <- renderUI({
    req(rv$data_df)
    plots <- covariate_plots(rv$data_df, arm_var = input$arm_var, palette = theme_palette())
    if (!length(plots)) return(p("No covariate plots are available."))
    tagList(lapply(names(plots), function(nm) plotlyOutput(paste0("cov_plot_", nm), height = "300px")))
  })
  
  # Render each plot
  observe({
    req(rv$data_df)
    plots <- covariate_plots(rv$data_df, arm_var = input$arm_var, palette = theme_palette())
    lapply(names(plots), function(nm) {
      local({
        id <- paste0("cov_plot_", nm)
        p  <- plots[[nm]]
        output[[id]] <- renderPlotly({ p })
      })
    })
  })
  
  
  # Run analysis
  run_output <- reactiveVal(list(results = NULL, log = "Analysis has not been run yet."))
  console_log <- reactiveVal("")
  
  run_analysis_results <- eventReactive(input$run_analysis, {
    validate(need(rv$data_df, "Please upload or simulate data first."))
    validate(need(input$time_var, "Please map Time-to-Event column."))
    validate(need(input$status_var, "Please map Status column."))
    validate(need(input$arm_var, "Please map Treatment Arm column."))
    
    analysis_results <- NULL
    log_text <- capture.output({
      withProgress(message = 'Running Analysis', value = 0, {
        setProgress(0.2, detail = "Preparing analysis data...")
        analysis_data <- data.frame(
          time = rv$data_df[[input$time_var]],
          status = as.numeric(rv$data_df[[input$status_var]]),
          arm = as.factor(rv$data_df[[input$arm_var]])
        )
        if (!is.null(input$strata_var) && nzchar(input$strata_var) &&
            (input$model_selection %in% c("Additive Stratified Model","Multiplicative Stratified Model"))) {
          analysis_data$stratum <- as.factor(rv$data_df[[input$strata_var]])
        }
        # Dependent censoring variable available but not directly used in log-rank; retained for models that need it
        if (!is.null(input$dep_cens_var) && nzchar(input$dep_cens_var) &&
            (input$model_selection %in% c("Dependent Censoring Model"))) {
          analysis_data$dep_cens_status <- rv$data_df[[input$dep_cens_var]]
        }
        
        # ----- Log-rank (stratified if applicable) -----
        setProgress(0.5, detail = "Log-rank test...")
        logrank_summary_df <- NULL
        analysis_data_for_plot <- NULL
        km_note_text <- NULL
        try({
          cat("\n--- Survival Analysis ---\n")
          if ("stratum" %in% names(analysis_data)) {
            fixed_formula <- as.formula("Surv(time, status) ~ arm + strata(stratum)")
            km_note_text  <- sprintf("<em>Showing KM curves by arm within each stratum of <b>%s</b>.</em>", input$strata_var)
          } else {
            fixed_formula <- as.formula("Surv(time, status) ~ arm")
            km_note_text  <- "<em>KM curves by arm.</em>"
          }
          logrank_test <- survdiff(fixed_formula, data = analysis_data)
          print(logrank_test)
          p_value <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
          logrank_summary_df <- data.frame(
            Statistic = "Chi-Square",
            Value = round(logrank_test$chisq, 3),
            DF = length(logrank_test$n) - 1,
            `P-Value` = format.pval(p_value, eps = .001, digits = 3)
          )
          analysis_data_for_plot <- analysis_data
        })
        
        # ----- Power / Sample size -----
        setProgress(0.8, detail = if (input$calc_method == "Analytical") "Computing (analytical) …" else "Computing (repeated) …")
        make_power_plot <- function(df_power, title_suffix = "") {
          pal <- theme_palette()
          ggplot(df_power, aes(x = N_per_arm, y = Power, group = 1)) +
            geom_line(size = 1.2, alpha = 0.9, color = pal[1]) +
            geom_point(size = 3.5, color = pal[2]) +
            scale_y_continuous(limits = c(0,1)) +
            labs(x = "Sample size per arm", y = "Power", title = paste("Method:", title_suffix)) +
            transparent_theme(base_size = 13)
        }
        
        
        # parse Ns
        if (input$analysis_type == "Power") {
          n_vec <- as.numeric(trimws(strsplit(input$sample_sizes, ",")[[1]]))
          n_vec <- n_vec[is.finite(n_vec) & n_vec > 0]
          if (!length(n_vec)) n_vec <- c(100,150,200)
        } else {
          # grid for hunting minimum N meeting target_power
          grid <- seq(30, 1000, by = 10)
        }
        
        if (input$calc_method == "Repeated") {
          R <- input$R_reps %||% 500
          if (input$analysis_type == "Power") {
            power_df <- repeated_power_from_pilot(
              rv$data_df, input$time_var, input$status_var, input$arm_var,
              n_per_arm_vec = n_vec, alpha = input$alpha, R = R,
              strata_var = if ("stratum" %in% names(analysis_data)) input$strata_var else NULL,
              seed = if (is.na(input$seed_reps)) NULL else as.integer(input$seed_reps)
            )
            results_plot <- make_power_plot(power_df, "repeated")
            results_data <- power_df
          } else {
            # search minimal N achieving target power
            power_df <- repeated_power_from_pilot(
              rv$data_df, input$time_var, input$status_var, input$arm_var,
              n_per_arm_vec = grid, alpha = input$alpha, R = R,
              strata_var = if ("stratum" %in% names(analysis_data)) input$strata_var else NULL,
              seed = if (is.na(input$seed_reps)) NULL else as.integer(input$seed_reps)
            )
            meet <- subset(power_df, Power >= input$target_power)
            n_star <- if (nrow(meet)) min(meet$N_per_arm) else NA
            results_plot <- make_power_plot(power_df, "repeated") +
              geom_hline(yintercept = input$target_power, linetype = 2) +
              ggplot2::annotate("text", x = max(power_df$N_per_arm), y = input$target_power,
                                label = sprintf("Target %.2f", input$target_power), hjust = 1, vjust = -0.5, size = 3.5)
            results_data <- power_df
            if (!is.na(n_star)) {
              results_plot <- results_plot + ggplot2::annotate("point", x = n_star, y = input$target_power, size = 4)
            }
          }
        } else {
          # Analytical placeholder (keep your original behavior but prettier points/lines)
          if (input$analysis_type == "Power") {
            dfp <- data.frame(N_per_arm = c(100,150,200), Power = c(0.72,0.81,0.88))
            results_plot <- make_power_plot(dfp, "analytical")
            results_data <- dfp
          } else {
            grid <- data.frame(N_per_arm = c(120, 160, 200), Power = c(0.70, 0.80, 0.87))
            results_plot <- make_power_plot(grid, "analytical") +
              geom_hline(yintercept = input$target_power, linetype = 2)
            results_data <- grid
          }
        }
        
        results_summary <- data.frame(
          Arm_1 = sum(analysis_data$arm == levels(analysis_data$arm)[1]),
          Arm_2 = sum(analysis_data$arm == levels(analysis_data$arm)[2]),
          Events = sum(analysis_data$status == 1), Censored = sum(analysis_data$status == 0)
        )
        analysis_results <- list(
          results_plot = results_plot,
          results_data = results_data,
          results_summary = results_summary,
          logrank_summary = logrank_summary_df,
          analysis_data_for_plot = analysis_data_for_plot,
          km_note = km_note_text
        )
      })
    }, type = c("output","message"))
    
    console_log(paste(log_text, collapse = "\n"))
    list(results = analysis_results, log = paste(log_text, collapse = "\n"))
  })
  
  observeEvent(run_analysis_results(), {
    run_output(run_analysis_results())
    shinyjs::show(id = "download_reset_row")
    updateTabsetPanel(session, "main_tabs", selected = "Summary")
  })
  
  # KM plots (faceted by stratum; at most 2 per row)
  output$km_note <- renderUI({
    req(run_output()$results$km_note)
    HTML(run_output()$results$km_note)
  })
  pretty_arm_labels <- function(f) {
    lv <- levels(f)
    if (identical(lv, c("0","1")) || identical(lv, c("0", "1")) || identical(lv, c(0,1))) {
      c("Control","Treatment")
    } else lv
  }
  
  output$survival_plotly_output <- renderPlotly({
    req(run_output()$results$analysis_data_for_plot, input$alpha)
    plot_data <- run_output()$results$analysis_data_for_plot
    pal <- theme_palette()
    
    if ("stratum" %in% names(plot_data)) {
      str_levels <- levels(as.factor(plot_data$stratum))
      plots <- lapply(seq_along(str_levels), function(i) {
        st <- str_levels[i]
        d  <- subset(plot_data, stratum == st)
        d$arm <- factor(d$arm)
        fit <- survfit(Surv(time, status) ~ arm, data = d)
        
        gp <- ggsurvplot(
          fit, data = d,
          conf.int = TRUE, conf.int.alpha = 0.3, conf.int.style = "ribbon",
          conf.level = 1 - input$alpha,
          palette = pal[1:2],
          legend.title = input$arm_var,
          legend.labs  = pretty_arm_labels(d$arm),
          xlab = paste("Time in the units of", input$time_var),
          ylab = "Survival probability",
          ggtheme = transparent_theme()
        )$plot
        
        pl <- to_plotly_clear(gp)
        if (i > 1) pl <- plotly::layout(pl, showlegend = FALSE)
        pl
      })
      
      ncols <- 2
      nrows <- ceiling(length(plots) / ncols)
      sp <- do.call(plotly::subplot, c(plots, nrows = nrows,
                                       shareX = TRUE, shareY = TRUE,
                                       titleX = TRUE, titleY = TRUE,
                                       margin = 0.04))
      
      # add per-panel titles
      sp <- add_subplot_titles(sp, paste(input$strata_var, "=", str_levels))
      
      # add an overall title
      sp <- plotly::layout(
        sp,
        title = list(
          text = sprintf("Kaplan–Meier — %s by %s", input$arm_var, input$strata_var),
          x = 0.02, xanchor = "left"
        ),
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor  = "rgba(0,0,0,0)"
      )
      sp
      
    } else {
      plot_data$arm <- factor(plot_data$arm)
      fit <- survfit(Surv(time, status) ~ arm, data = plot_data)
      gp <- ggsurvplot(
        fit, data = plot_data,
        conf.int = TRUE, conf.int.alpha = 0.3, conf.int.style = "ribbon",
        conf.level = 1 - input$alpha,
        palette = pal[1:2],
        legend.title = input$arm_var,
        legend.labs  = pretty_arm_labels(plot_data$arm),
        xlab = paste("Time in the units of", input$time_var),
        ylab = "Survival probability",
        ggtheme = transparent_theme()
      )$plot
      to_plotly_clear(gp) %>%
        plotly::layout(
          title = list(
            text = sprintf("Kaplan–Meier — %s", input$arm_var),
            x = 0.02, xanchor = "left"
          )
        )
    }
  })
  
  
  
  
  
  # Power plot (already colorful and connected via make_power_plot)
  output$results_plot <- renderPlotly({
    req(run_output()$results$results_plot)
    to_plotly_clear(run_output()$results$results_plot)
    
  })
  
  # Summary (tables only)
  output$results_table_ui <- renderUI({
    req(run_output()$results$results_data)
    run_output()$results$results_data %>%
      kbl("html", caption = "Power and sample size results") %>%
      kable_styling(bootstrap_options = c("striped","hover","condensed"), full_width = FALSE) %>% HTML()
  })
  output$summary_table_ui <- renderUI({
    req(run_output()$results$results_summary)
    run_output()$results$results_summary %>%
      kbl("html", caption = "Summary measures derived from the data") %>%
      kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE) %>% HTML()
  })
  output$data_summary_ui <- renderUI({
    req(rv$data_df)
    arm_var <- isolate(input$arm_var %||% NULL)
    sm <- covariate_summary(rv$data_df, arm_var = arm_var)
    if (!length(sm)) return(HTML("<em>No covariate tables are available.</em>"))
    ui <- tagList()
    if (!is.null(sm$continuous) && nrow(sm$continuous)) {
      ui <- tagAppendChildren(ui,
                              h5("Continuous covariates"),
                              sm$continuous %>% kbl("html") %>% kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE) %>% HTML()
      )
    }
    if (!is.null(sm$categorical) && nrow(sm$categorical)) {
      ui <- tagAppendChildren(ui,
                              h5("Categorical covariates"),
                              sm$categorical %>% kbl("html") %>% kable_styling(bootstrap_options = c("striped","hover"), full_width = FALSE) %>% HTML()
      )
    }
    ui
  })
  
  # Console log (text only)
  output$console_log_output <- renderText({ paste(console_log(), collapse = "\n") })
  
  # ------------------ Downloads (PDF & HTML) ------------------
  data_provenance <- reactive({
    prov <- rv$provenance
    if (is.null(prov) && !is.null(rv$data_df)) {
      prov <- list(
        Source = rv$data_source %||% "uploaded",
        Number_of_rows = nrow(rv$data_df),
        Variable_names = names(rv$data_df),
        Covariates_defined = if (rv$data_source == "simulated") lapply(rv$covariates, function(d) d) else list()
      )
    }
    prov
  })
  get_pilot_data <- reactive({ rv$data_df })
  
  output$download_report_pdf <- downloadHandler(
    filename = function() paste0("RMSTpowerBoost_report_", Sys.Date(), ".pdf"),
    contentType = "application/pdf",
    content = function(file) {
      req(run_output()$results)
      id <- showNotification("Generating PDF report…", type="message", duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(id), add = TRUE)
      tpl <- make_inline_template()
      rmarkdown::render(
        input         = tpl,
        output_format = rmarkdown::pdf_document(),
        output_file   = basename(file),
        output_dir    = dirname(file),
        params        = list(
          inputs          = report_inputs_builder(input),
          results         = run_output()$results,
          log             = paste(console_log(), run_output()$log, sep = "\n"),
          data_provenance = data_provenance(),
          data            = get_pilot_data()
        ),
        envir         = new.env(parent = globalenv()),
        clean         = TRUE
      )
    }
  )
  output$download_report_html <- downloadHandler(
    filename = function() paste0("RMSTpowerBoost_report_", Sys.Date(), ".html"),
    contentType = "text/html",
    content = function(file) {
      req(run_output()$results)
      id <- showNotification("Generating HTML report…", type="message", duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(id), add = TRUE)
      tpl <- make_inline_template()
      rmarkdown::render(
        input         = tpl,
        output_format = rmarkdown::html_document(theme = "flatly", toc = TRUE, toc_depth = 3),
        output_file   = basename(file),
        output_dir    = dirname(file),
        params        = list(
          inputs          = report_inputs_builder(input),
          results         = run_output()$results,
          log             = paste(console_log(), run_output()$log, sep = "\n"),
          data_provenance = data_provenance(),
          data            = get_pilot_data()
        ),
        envir         = new.env(parent = globalenv()),
        clean         = TRUE
      )
    }
  )
  
  # Reveal Step2+3 when data is present; hide simulate after success
  observe({
    shinyjs::toggle(id = "model_analysis_panel", condition = !is.null(rv$data_df))
    shinyjs::toggle(id = "simulate_panel", condition = is.null(rv$data_df) && rv$data_mode == "Generate")
    shinyjs::toggle(id = "upload_panel",   condition = is.null(rv$data_df) && rv$data_mode == "Upload")
  })
  
  # Reset all (appears after analysis)
  observeEvent(input$reset_all, {
    rv$covariates <- list()
    rv$cat_rows <- tibble::tibble(cat = character(), prob = numeric(), coef = numeric())
    rv$data_df <- NULL; rv$data_source <- NULL; rv$console_buf <- character(0); rv$provenance <- NULL
    shinyjs::hide("download_reset_row")
    shinyjs::show("simulate_panel")
    updateTabsetPanel(session, "main_tabs", selected = "Instructions")
    showNotification("All inputs reset.", type="message")
  })
}

options(shiny.launch.browser = TRUE)
shinyApp(ui, server)
