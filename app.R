# app.R â€” RMSTpowerBoost (single file, no external Rmd required)

# ------------------ Packages ------------------
packages <- c(
  "shiny","shinyjs","bslib","DT","ggplot2","plotly","survival","survminer",
  "kableExtra","magrittr","rmarkdown","dplyr","tidyr","purrr","stringr","tibble"
)
invisible(lapply(packages, function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
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

DT_25 <- function(df) DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)

is_categorical_like <- function(x) is.factor(x) || is.character(x) || is.logical(x)
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

covariate_plots <- function(df, arm_var = NULL) {
  if (is.null(df) || !nrow(df)) return(list())
  vars <- setdiff(names(df), c("time","status",arm_var))
  if (!length(vars)) return(list())
  plots <- list()
  for (v in vars) {
    x <- df[[v]]
    if (is.numeric(x)) {
      p <- ggplot(df, aes(x = .data[[v]])) +
        geom_histogram(bins = 30, alpha = 0.9, fill = "#17a2b8") +
        labs(title = paste("Histogram of", v), x = v, y = "Count") +
        theme_light()
      plots[[v]] <- plotly::ggplotly(p)
    } else {
      p <- ggplot(df, aes(x = as.factor(.data[[v]]), fill = as.factor(.data[[v]]))) +
        geom_bar(alpha = 0.9) +
        guides(fill = "none") +
        labs(title = paste("Bar chart of", v), x = v, y = "Count") +
        theme_light()
      plots[[v]] <- plotly::ggplotly(p)
    }
  }
  plots
}

# Build a model.matrix column order from covariate definitions
build_mm_columns <- function(cov_defs, include_intercept = TRUE) {
  if (!length(cov_defs)) return(character(0))
  # Construct a tiny data.frame that contains all factor levels (so model.matrix creates all columns)
  rows <- 0
  cols <- list()
  for (d in cov_defs) {
    if (d$type == "continuous") {
      cols[[d$name]] <- 0
    } else if (d$type == "categorical") {
      lv <- d$params$labels
      if (is.null(lv) || !length(lv)) {
        # auto labels if not provided
        k <- length(d$params$prob %||% c(0,1))
        lv <- paste0(d$name, seq_len(k))
      }
      cols[[d$name]] <- factor(lv, levels = lv)
      rows <- max(rows, length(lv))
    }
  }
  # Recycle short columns to 'rows'
  df <- as.data.frame(lapply(cols, function(x){
    if (length(x) == 1) rep(x, max(1, rows)) else x
  }), stringsAsFactors = FALSE)
  form <- as.formula(paste0(if (include_intercept) "~ 1 +" else "~ -1 +",
                            paste(vapply(cov_defs, function(d) d$name, character(1)), collapse = " + ")))
  mm <- model.matrix(form, data = df)
  colnames(mm)
}

# Coeff vector in model.matrix column order
assemble_beta <- function(cov_defs, user_betas, include_intercept = TRUE, intercept_value = 0) {
  # user_betas: named list: for continuous a single number; for categorical a numeric vector aligned to levels (K or K-1 depending on include_intercept)
  cols <- build_mm_columns(cov_defs, include_intercept = include_intercept)
  out <- numeric(0)
  if (include_intercept) {
    # the first column is "(Intercept)"
    out <- c(out, user_betas[["(Intercept)"]] %||% intercept_value)
  }
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
        # K-1 betas, baseline is level 1
        need <- length(lv) - 1
        b <- as.numeric(user_betas[[d$name]])
        if (length(b) < need) stop("Categorical '", d$name, "' requires ", need, " coefficients (intercept included).")
        b <- b[seq_len(need)]
        out <- c(out, b)
      } else {
        # K betas
        need <- length(lv)
        b <- as.numeric(user_betas[[d$name]])
        if (length(b) < need) stop("Categorical '", d$name, "' requires ", need, " coefficients (no intercept).")
        b <- b[seq_len(need)]
        out <- c(out, b)
      }
    }
  }
  out
}

# Create an inline Rmd for both PDF and HTML â€” includes coefficients & provenance
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
    # Add this below the setup chunk
    "",
    "This document summarizes the simulation design, analysis configuration, and results for a scenario generated using the RMSTpowerBoost tool.",
    "",
    "# Covariate Generation Mechanism",
    "",
    "```{r covariate-generation, echo=FALSE}",
    "cov_list <- or_else(params$data_provenance$Covariates_defined, list())",
    "if (length(cov_list)) {",
    "  cov_tbl <- do.call(rbind, lapply(cov_list, function(d) {",
    "    data.frame(",
    "      Variable     = or_else(d$name, \"\"),",
    "      Type         = or_else(d$type, \"\"),",
    "      Distribution = or_else(d$dist, \"\"),",
    "      Parameters   = to_kv(d$params),",
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
    "# Event Time Generation",
    "",
    "```{r event-time, echo=FALSE}",
    "et <- or_else(params$data_provenance$Event_time, list())",
    "et_tbl <- data.frame(",
    "  Field = c(\"Model\", \"Baseline parameters\"),",
    "  Value = c(",
    "    as.character(or_else(et$model, or_else(et$Model, \"\"))),",
    "    to_kv(or_else(et$baseline, or_else(et$Baseline, list())))",
    "  ),",
    "  stringsAsFactors = FALSE",
    ")",
    "kbl(et_tbl, booktabs = TRUE, caption = \"Event-time model and baseline\",",
    "    row.names = FALSE, col.names = c(\"Field\",\"Value\")) %>%",
    "  kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "```",
    
    "",
    "# Data Generating Mechanism",
    "",
    "```{r provenance, echo=FALSE}",
    "prov <- params$data_provenance",
    "if (!is.null(prov) && identical(prov$Source, \"simulated\")) {",
    "  prov_tbl <- data.frame(",
    "    Field = c(\"Source\", \"Number of rows\", \"Variable names\"),",
    "    Value = c(as.character(prov$Source), as.character(prov$Number_of_rows),",
    "              paste(prov$Variable_names, collapse = \", \")),",
    "    stringsAsFactors = FALSE",
    "  )",
    "  kbl(prov_tbl, booktabs = TRUE, caption = \"Data source and provenance\",",
    "      row.names = FALSE, col.names = c(\"Field\",\"Value\")) %>%",
    "    kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "",
    "  cov_list <- prov$Covariates_defined",
    "  if (!is.null(cov_list) && length(cov_list)) {",
    "    cov_tbl <- do.call(rbind, lapply(cov_list, function(d){",
    "      data.frame(",
    "        Variable     = or_else(d$name, \"\"),",
    "        Type         = or_else(d$type, \"\"),",
    "        Distribution = or_else(d$dist, \"\"),",
    "        Parameters   = to_kv(d$params),",
    "        Transform    = paste(or_else(d$transform, character(0)), collapse = \", \"),",
    "        stringsAsFactors = FALSE",
    "      )",
    "    }))",
    "    kbl(cov_tbl, booktabs = TRUE, caption = \"Covariate definitions used for simulation\",",
    "        row.names = FALSE) %>%",
    "      kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "  }",
    "",
    "  if (!is.null(prov$Event_time)) {",
    "    et <- prov$Event_time",
    "    et_tbl <- data.frame(",
    "      Field = c(\"Model\", \"Baseline parameters\"),",
    "      Value = c(as.character(or_else(et$model, or_else(et$Model, \"\"))), to_kv(or_else(et$baseline, et$Baseline))),",
    "      stringsAsFactors = FALSE",
    "    )",
    "    kbl(et_tbl, booktabs = TRUE, caption = \"Event-time model and baseline\",",
    "        row.names = FALSE, col.names = c(\"Field\",\"Value\")) %>%",
    "      kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "  }",
    "",
    "  if (!is.null(prov$Treatment) || !is.null(prov$Effects) || !is.null(prov$Strata)) {",
    "    tr  <- prov$Treatment",
    "    eff <- prov$Effects",
    "    st  <- prov$Strata",
    "    rows <- list(",
    "      c(\"Assignment\", as.character(or_else(tr$assignment, or_else(tr$Assignment, \"\")))),",
    "      c(\"Allocation\", as.character(or_else(tr$allocation, or_else(tr$Allocation, \"\")))),",
    "      c(\"Treatment effect\", as.character(or_else(eff$treatment, or_else(eff$Treatment, \"\")))),",
    "      c(\"Intercept\", as.character(or_else(eff$intercept, or_else(eff$Intercept, \"\")))),",
    "      c(\"Formula\", as.character(or_else(eff$formula, or_else(eff$Formula, \"\")))),",
    "      c(\"Beta\", to_kv(eff$beta))",
    "    )",
    "    if (!is.null(st)) rows <- c(rows, list(c(\"Stratify by\", paste(st, collapse = \", \"))))",
    "    tr_tbl <- data.frame(Field = vapply(rows, `[[`, \"\" ,1),",
    "                         Value = vapply(rows, `[[`, \"\" ,2),",
    "                         stringsAsFactors = FALSE)",
    "    tr_tbl <- tr_tbl[nchar(trimws(tr_tbl$Value))>0, , drop=FALSE]",
    "    if (nrow(tr_tbl)) {",
    "      kbl(tr_tbl, booktabs = TRUE, caption = \"Treatment assignment and effects\",",
    "          row.names = FALSE, col.names = c(\"Field\",\"Value\")) %>%",
    "        kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "    }",
    "  }",
    "",
    "  if (!is.null(prov$Censoring) || !is.null(prov$Frailty)) {",
    "    cz <- prov$Censoring; fr <- prov$Frailty",
    "    rows <- list(",
    "      c(\"Censoring mode\", as.character(or_else(cz$mode, or_else(cz$Mode, \"\")))),",
    "      c(\"Target overall censoring\", as.character(or_else(cz$target, or_else(cz$Target, \"\")))),",
    "      c(\"Administrative time\", as.character(or_else(cz$admin_time, or_else(cz$Admin_time, \"\"))))",
    "    )",
    "    if (!is.null(fr)) {",
    "      rows <- c(rows, list(",
    "        c(\"Frailty type\", as.character(or_else(fr$type, or_else(fr$Type, \"\")))),",
    "        c(\"Frailty variance\", as.character(or_else(fr$var, or_else(fr$Var, \"\")))),",
    "        c(\"Frailty group\", as.character(or_else(fr$group, or_else(fr$Group, \"\"))))",
    "      ))",
    "    }",
    "    cz_tbl <- data.frame(Field = vapply(rows, `[[`, \"\", 1),",
    "                         Value = vapply(rows, `[[`, \"\", 2),",
    "                         stringsAsFactors = FALSE)",
    "    cz_tbl <- cz_tbl[nchar(trimws(cz_tbl$Value))>0, , drop=FALSE]",
    "    if (nrow(cz_tbl)) {",
    "      kbl(cz_tbl, booktabs = TRUE, caption = \"Censoring and frailty settings\",",
    "          row.names = FALSE, col.names = c(\"Field\",\"Value\")) %>%",
    "        kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "    }",
    "  }",
    "} else {",
    "  cat(\"Data were uploaded; a simulation mechanism was not used.\\n\")",
    "}",
    "```",
    "",
    "# Analysis Configuration",
    "",
    "```{r input-parameters, echo=FALSE}",
    "inputs_df <- data.frame(",
    "  Parameter = names(params$inputs),",
    "  Value     = unlist(lapply(params$inputs, function(x) paste(x, collapse = \", \"))),",
    "  stringsAsFactors = FALSE",
    ")",
    "kbl(inputs_df, booktabs = TRUE, caption = \"Analysis Input Parameters\",",
    "    row.names = FALSE, col.names = c(\"Parameter\",\"Value\")) %>%",
    "  kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "```",
    "",
    "# Survival Analysis of the Data",
    "",
    "```{r km-plot, eval = !is.null(params$results$analysis_data_for_plot)}",
    "plot_data <- params$results$analysis_data_for_plot",
    "fit <- survfit(Surv(time, status) ~ arm, data = plot_data)",
    "p <- ggsurvplot(",
    "  fit, data = plot_data,",
    "  palette = c(\"#007BFF\", \"#D9534F\"),",
    "  legend.title = params$inputs$arm_var,",
    "  xlab = paste(\"Time in the units of\", params$inputs$time_var),",
    "  ylab = \"Survival probability\",",
    "  ggtheme = theme_light()",
    ")",
    "p$plot",
    "```",
    "",
    "```{r logrank-summary, eval = !is.null(params$results$logrank_summary)}",
    "kbl(params$results$logrank_summary, booktabs = TRUE,",
    "    caption = \"Logâ€“rank test results\", row.names = FALSE) %>%",
    "  kable_styling(latex_options = \"hold_position\")",
    "```",
    "",
    "# Power and Sample Size Analysis",
    "",
    "```{r power-curve, eval = !is.null(params$results$results_plot)}",
    "params$results$results_plot",
    "```",
    "",
    "```{r results-table, eval = !is.null(params$results$results_data)}",
    "kbl(params$results$results_data, booktabs = TRUE,",
    "    caption = \"Power and sample size results\", row.names = FALSE) %>%",
    "  kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "```",
    "",
    "```{r effect-size, eval = !is.null(params$results$results_summary)}",
    "kbl(params$results$results_summary, booktabs = TRUE,",
    "    caption = \"Summary measures derived from the data\", row.names = FALSE) %>%",
    "  kable_styling(latex_options = c(\"striped\", \"hold_position\"))",
    "```",
    "",
    "# Console Output",
    "",
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
            column(6, textInput("cov_name", "Variable name", value = "", placeholder = "x1, x2 â€¦ auto if empty")),
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
          tags$hr(),
          h5("1b. Simulation Settings"),
          fluidRow(
            column(4, numericInput("sim_n", "Sample size", value = 300, min = 10)),
            column(4, textInput("sim_allocation", "Allocation (a:b)", value = "1:1")),
            column(4, numericInput("sim_treat_eff", "Treatment Î² (arm)", value = -0.2, step = 0.05))
          ),
          fluidRow(
            column(6, checkboxInput("intercept_in_mm", "Include intercept in model.matrix (Î²0 inside Î²)", value = TRUE)),
            column(6, numericInput("user_intercept", "Î²0 (used if no intercept in model.matrix)", value = 0))
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
            column(4, numericInput("L", "RMST L (Ï„)", value = 365, min = 1))
          ),
          uiOutput("analysis_inputs_ui"),
          sliderInput("alpha", "Significance Level (Î±)", min = 0.01, max = 0.1, value = 0.05, step = 0.01),
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
                     tags$li("2c. Set the analysis horizon (Ï„) and the significance level (Î±).")
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
                 h4("Kaplanâ€“Meier Survival Plot"),
                 plotlyOutput("survival_plotly_output", height = "500px"),
                 hr(),
                 h4("Power vs. Sample Size"),
                 plotlyOutput("results_plot", height = "500px")
        ),
        tabPanel("Summary",
                 h4("Analysis Results (Tables Only)"),
                 uiOutput("results_table_ui"),
                 hr(),
                 h4("Effect Size Summary"),
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

# ------------------ Server ------------------
server <- function(input, output, session) {
  bslib::bs_themer()
  license_content <- tryCatch(paste(readLines("LICENSE"), collapse = "\n"), error = function(e) "LICENSE not found.")
  output$license_display <- renderPrint({ cat(license_content) })
  
  rv <- reactiveValues(
    covariates = list(),          # defs with params+transform+beta (per var)
    cat_rows = tibble::tibble(cat = character(), prob = numeric(), coef = numeric()), # current builder rows
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
          column(6, numericInput("cont_beta", "Coefficient Î²", value = 0))
        ),
        uiOutput("cont_param_ui")
      )
    } else {
      shinyjs::hide("transform_block")
      tagList(
        fluidRow(
          column(6, textInput("cat_add_name", "Add category name", placeholder = "auto if blank")),
          column(3, numericInput("cat_add_prob", "Probability", value = NA, min = 0, max = 1, step = 0.01)),
          column(3, numericInput("cat_add_coef", "Coefficient Î²", value = 0))
        ),
        fluidRow(
          column(6, actionButton("add_cat_row", "Add category", icon=icon("plus"))),
          column(6, actionButton("reset_cat_rows", "Reset categories", icon=icon("trash")))
        ),
        br(),
        DTOutput("cat_table"),
        helpText("Tip: If you include intercept in model.matrix, only Kâˆ’1 coefficients are used (last levelâ€™s Î² is ignored).")
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
    if (!nzchar(nm)) {
      nm <- paste0("L", nrow(rv$cat_rows) + 1)
    }
    pr <- input$cat_add_prob
    if (!is.na(pr) && (pr < 0 || pr > 1)) {
      showNotification("Probability must be between 0 and 1.", type = "warning"); return()
    }
    cf <- input$cat_add_coef
    if (!is.finite(cf)) {
      showNotification("Coefficient must be numeric.", type = "warning"); return()
    }
    rv$cat_rows <- bind_rows(rv$cat_rows, tibble::tibble(cat = nm, prob = pr, coef = cf))
  })
  observeEvent(input$reset_cat_rows, { rv$cat_rows <- tibble::tibble(cat = character(), prob = numeric(), coef = numeric()) })
  output$cat_table <- renderDT({
    if (!nrow(rv$cat_rows)) {
      DT::datatable(data.frame(Message="No categories yet â€” add rows above."), options = list(dom='t'), rownames = FALSE)
    } else {
      DT_25(rv$cat_rows)
    }
  })
  
  # Transform visibility (continuous only)
  observe({
    shinyjs::toggle(id = "transform_block", condition = (input$cov_type %||% "") == "continuous")
  })
  
  # Reset builder
  observeEvent(input$reset_cov_builder, {
    updateTextInput(session, "cov_name", value = "")
    updateSelectInput(session, "cov_type", selected = "continuous")
    updateSelectInput(session, "cont_dist", selected = "normal")
    updateNumericInput(session, "cont_beta", value = 0)
    updateNumericInput(session, "tf_center", value = 0)
    updateNumericInput(session, "tf_scale", value = 1)
    rv$cat_rows <- tibble::tibble(cat = character(), prob = numeric(), coef = numeric())
  })
  
  # Helper: auto covariate name x1, x2â€¦
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
      # params
      pars <- switch(input$cont_dist,
                     normal    = list(mean = input$p_mean, sd = input$p_sd),
                     lognormal = list(meanlog = input$p_meanlog, sdlog = input$p_sdlog),
                     gamma     = list(shape = input$p_shape,   scale = input$p_scale),
                     weibull   = list(shape = input$p_wshape,  scale = input$p_wscale),
                     uniform   = list(min   = input$p_min,     max   = input$p_max),
                     t         = list(df    = input$p_df),
                     beta      = list(shape1 = input$p_shape1, shape2 = input$p_shape2),
                     list()
      )
      # transform
      tf <- c(sprintf("center(%s)", input$tf_center), sprintf("scale(%s)", input$tf_scale))
      # beta
      if (!is.finite(as.numeric(input$cont_beta))) {
        showNotification("Continuous coefficient must be a single number.", type = "error"); return()
      }
      rv$covariates <- c(rv$covariates, list(list(
        name = vname, type = "continuous", dist = input$cont_dist, params = pars,
        transform = tf, beta = as.numeric(input$cont_beta)
      )))
      
    } else { # categorical
      if (!nrow(rv$cat_rows)) { showNotification("Add at least one category.", type = "error"); return() }
      cats <- rv$cat_rows$cat
      # equal probs if any NA; otherwise use provided
      prob <- rv$cat_rows$prob
      if (any(is.na(prob))) prob[is.na(prob)] <- 1/length(prob)
      if (any(prob < 0) || abs(sum(prob) - 1) > 1e-6) {
        showNotification("Probabilities must be non-negative and sum to 1.", type = "error"); return()
      }
      coef <- rv$cat_rows$coef
      if (any(!is.finite(coef))) { showNotification("All category coefficients must be numeric.", type="error"); return() }
      
      # dist decision: 2 cats -> bernoulli
      if (length(cats) == 2) {
        # define bernoulli for level2 probability
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
    }
  })
  
  # Covariate list table
  output$cov_table <- renderDT({
    if (!length(rv$covariates)) return(DT_25(data.frame()))
    show <- purrr::map_df(rv$covariates, function(d) {
      beta_txt <- if (length(d$beta)>1) paste(d$beta, collapse=", ") else as.character(d$beta)
      tibble::tibble(
        name = d$name, type = d$type, dist = d$dist,
        params = paste(names(d$params), unlist(d$params), sep="=", collapse="; "),
        transform = paste(d$transform %||% character(0), collapse = ", "),
        beta = beta_txt
      )
    })
    DT_25(show)
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
    
    # Validate coefficients lengths vs encoding
    include_intercept <- isTRUE(input$intercept_in_mm)
    # Build a clean copy of cov_defs for formula
    cov_defs <- lapply(rv$covariates, function(d){
      list(name = d$name, type = d$type, dist = d$dist, params = d$params, transform = d$transform)
    })
    
    # Build user_betas
    user_betas <- list()
    if (include_intercept) user_betas[["(Intercept)"]] <- 0 # actual intercept comes from Î² vector of user per var, not b0
    for (d in rv$covariates) {
      if (d$type == "continuous") {
        user_betas[[d$name]] <- as.numeric(d$beta)
      } else {
        # categorical
        lv <- if (d$dist == "bernoulli") c("0","1") else (d$params$labels %||% paste0(d$name, seq_along(d$params$prob)))
        need <- if (include_intercept) length(lv)-1 else length(lv)
        if (length(d$beta) < need) {
          showNotification(sprintf("'%s' needs %d coefficients but %d provided.", d$name, need, length(d$beta)), type="error")
          return()
        }
        user_betas[[d$name]] <- as.numeric(d$beta)
      }
    }
    # Assemble beta in mm column order
    beta_vec <- assemble_beta(cov_defs, user_betas, include_intercept = include_intercept, intercept_value = input$user_intercept)
    mm_cols  <- build_mm_columns(cov_defs, include_intercept = include_intercept)
    
    # Baseline
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
    
    # Effects: formula includes or excludes intercept
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
      cat("---- Data Simulation ----\n")
      print(str(rec))
      cat("Model matrix columns order:\n")
      print(mm_cols)
      cat("Beta vector:\n")
      print(beta_vec)
    })
    rv$console_buf <- c(rv$console_buf, buf)
    
    dat <- tryCatch({
      simulate_from_recipe(rec, seed = rec$seed)
    }, error = function(e) {
      showNotification(paste("Simulation failed:", e$message), type = "error")
      NULL
    })
    if (is.null(dat)) return()
    
    rv$data_df <- dat
    rv$data_source <- "simulated"
    showNotification("Simulation complete.", type = "message")
    
    shinyjs::hide(id = "simulate_panel")
    shinyjs::show(id = "model_analysis_panel")
    updateTabsetPanel(session, "main_tabs", selected = "Data Preview")
    
    # keep provenance for report
    rv$provenance <- list(
      Source = "simulated",
      Number_of_rows = nrow(dat),
      Variable_names = names(dat),
      Covariates_defined = lapply(rv$covariates, function(d) d),
      Event_time = list(model = input$sim_model, baseline = baseline),
      Treatment  = list(assignment = "randomization", allocation = input$sim_allocation),
      Effects    = list(treatment = input$sim_treat_eff,
                        intercept_report = if (include_intercept) "(in model.matrix Î²)" else input$user_intercept,
                        intercept_in_mm = include_intercept,
                        formula = deparse(form),
                        mm_cols = mm_cols),
      Censoring  = list(mode = "target_overall", target = input$sim_cens, admin_time = Inf)
    )
  })
  
  # Column mapping
  output$col_mapping_ui <- renderUI({
    df <- rv$data_df; req(df)
    cn <- names(df)
    tagList(
      fluidRow(
        column(4, selectInput("time_var", "Time-to-Event", choices = cn, selected = "time" %||% cn[1])),
        column(4, selectInput("status_var", "Status (1=event)", choices = cn, selected = "status" %||% cn[min(2,length(cn))])),
        column(4, selectInput("arm_var", "Treatment Arm (1=treat)", choices = cn,
                              selected = if ("arm" %in% cn) "arm" else cn[min(3,length(cn))]))
      )
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
  output$cov_plots_ui <- renderUI({
    req(rv$data_df)
    plots <- covariate_plots(rv$data_df, arm_var = input$arm_var)
    if (!length(plots)) return(p("No covariate plots are available."))
    tagList(lapply(names(plots), function(nm) {
      plotlyOutput(paste0("cov_plot_", nm), height = "300px")
    }))
  })
  observe({
    req(rv$data_df)
    plots <- covariate_plots(rv$data_df, arm_var = input$arm_var)
    lapply(names(plots), function(nm) {
      local({ id <- paste0("cov_plot_", nm); p <- plots[[nm]]; output[[id]] <- renderPlotly({ p }) })
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
        setProgress(0.5, detail = "Log-rank test...")
        logrank_summary_df <- NULL
        analysis_data_for_plot <- NULL
        try({
          cat("\n--- Survival Analysis ---\n")
          fixed_formula <- as.formula("Surv(time, status) ~ arm")
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
        setProgress(0.8, detail = "Computing power curve (placeholder)â€¦")
        results_plot <- ggplot(data.frame(n = c(100,150,200), power = c(0.72,0.81,0.88)),
                               aes(n, power)) + geom_line(color="#17a2b8") + geom_point(color="#d9534f") +
          labs(x = "Sample size per arm", y = "Power") + theme_light()
        results_data <- data.frame(N_per_arm = c(100,150,200), Power = c(0.72,0.81,0.88))
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
          analysis_data_for_plot = analysis_data_for_plot
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
  
  # Plots
  output$survival_plotly_output <- renderPlotly({
    req(run_output()$results$analysis_data_for_plot, input$alpha)
    plot_data <- run_output()$results$analysis_data_for_plot
    fit <- survfit(Surv(time, status) ~ arm, data = plot_data)
    p <- ggsurvplot(
      fit, data = plot_data,
      conf.int = TRUE, conf.int.alpha = 0.3, conf.int.style = "ribbon",
      conf.level = 1 - input$alpha,
      palette = c("#007BFF", "#D9534F"),
      legend.title = input$arm_var,
      xlab = paste("Time in the units of", input$time_var),
      ylab = "Survival probability",
      ggtheme = theme_light()
    )
    ggplotly(p$plot)
  })
  output$results_plot <- renderPlotly({ req(run_output()$results$results_plot); plotly::ggplotly(run_output()$results$results_plot, tooltip = c("x","y")) })
  
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
      id <- showNotification("Generating PDF reportâ€¦", type="message", duration = NULL, closeButton = FALSE)
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
      id <- showNotification("Generating HTML reportâ€¦", type="message", duration = NULL, closeButton = FALSE)
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
    # wipe everything except sourced functions
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
