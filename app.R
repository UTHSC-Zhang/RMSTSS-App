# --- 1. Load All Necessary Libraries ---
remotes::install_github("https://github.com/biomedical-data-science-internship/UTHSC-package.git")
library(shiny)
library(shinythemes)
library(ggplot2)
library(plotly)
library(kableExtra)
library(DT)
library(survival)
library(survminer)
library(dplyr)
library(magrittr)
library(future)
library(future.apply)
library(RMSTdesign) # Load the package containing the backend functions

# --- 3. UI Definition ---
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  
  titlePanel("RMSTdesign: Power and Sample Size Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Analysis Setup"),
      
      # Step 1: Upload Data
      fileInput("pilot_data_upload", "1. Upload Pilot Data (.csv)", accept = ".csv", width = "100%"),
      
      # Step 2: Main Configuration (appears after upload)
      conditionalPanel(
        condition = "output.fileUploaded",
        h4("2. Model & Column Mapping"),
        selectInput("model_selection", "Select RMST Model",
                    choices = c("Linear IPCW Model", 
                                "Additive Stratified Model", 
                                "Multiplicative Stratified Model",
                                "Semiparametric (GAM) Model",
                                "Dependent Censoring Model")),
        uiOutput("col_mapping_ui")
      ),
      
      # Step 3: Analysis Parameters (appears after model is chosen)
      conditionalPanel(
        condition = "output.fileUploaded && input.model_selection != ''",
        hr(),
        h4("3. Analysis Parameters"),
        numericInput("tau", "RMST Truncation Time (τ)", value = 365, min = 1),
        radioButtons("analysis_type", "Analysis Type", choices = c("Power Calculation", "Sample Size Search"), inline = TRUE),
        uiOutput("method_selection_ui"),
        sliderInput("alpha", "Significance Level (α)", min = 0.01, max = 0.1, value = 0.05, step = 0.01),
        uiOutput("calculation_params_ui"),
        checkboxInput("show_km", "Show Kaplan-Meier Curve", value = FALSE),
        hr(),
        actionButton("run_analysis", "Run Analysis", icon = icon("play"), class = "btn-primary", width = "100%")
      )
    ),
    
    mainPanel(
      width = 8,
      # Static tabsetPanel that will be dynamically modified by the server
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Instructions", value = "instructions_tab",
                 h3("Welcome to RMSTdesign!"),
                 p("Please follow the steps in the left panel to proceed.")
        ),
        tabPanel("Data Preview", value = "preview_tab",
                 DT::dataTableOutput("data_preview_table")
        )
      )
    )
  )
)

# --- 4. Server Definition ---
server <- function(input, output, session) {
  
  # --- Reactive Values & Triggers ---
  output$fileUploaded <- reactive({ return(!is.null(input$pilot_data_upload)) })
  outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)
  
  pilot_data_reactive <- reactive({ req(input$pilot_data_upload); read.csv(input$pilot_data_upload$datapath, stringsAsFactors = TRUE) })
  
  dynamic_tabs_ids <- reactiveVal(list())
  
  # --- Dynamic UI Rendering ---
  output$col_mapping_ui <- renderUI({ df <- pilot_data_reactive(); req(df); column_names <- names(df); tagList(selectInput("time_var", "Time-to-Event", choices = column_names, selected = intersect(column_names, c("time", "futime"))[1]), selectInput("status_var", "Status (1=event)", choices = column_names, selected = intersect(column_names, c("status", "event_primary", "pstat"))[1]), selectInput("arm_var", "Treatment Arm (1=treat)", choices = column_names, selected = intersect(column_names, "arm")[1]), if (input$model_selection %in% c("Additive Stratified Model", "Multiplicative Stratified Model")) { selectInput("strata_var", "Stratification Variable", choices = column_names, selected = intersect(column_names, c("region", "strata"))[1]) }, if (input$model_selection == "Dependent Censoring Model") { selectInput("dep_cens_var", "Dependent Censoring Status", choices = column_names, selected = intersect(column_names, c("dep_cens_status", "dstat"))[1]) }, selectizeInput("linear_terms", "Linear Covariates (Optional)", choices = column_names, multiple = TRUE), if (input$model_selection == "Semiparametric (GAM) Model") { selectizeInput("smooth_terms", "Non-Linear Covariates (GAM)", choices = column_names, multiple = TRUE) }) })
  output$method_selection_ui <- renderUI({ if (input$model_selection %in% c("Linear IPCW Model", "Multiplicative Stratified Model")) { radioButtons("calc_method", "Method", choices = c("Analytical", "Bootstrap"), selected = "Analytical", inline = TRUE) } })
  output$calculation_params_ui <- renderUI({ is_bootstrap <- (input$model_selection %in% c("Linear IPCW Model", "Multiplicative Stratified Model") && !is.null(input$calc_method) && input$calc_method == "Bootstrap") || (input$model_selection == "Semiparametric (GAM) Model"); power_inputs <- if(input$analysis_type == "Power Calculation") { textInput("sample_sizes", "Sample Sizes (per arm/stratum, comma-separated)", value = "100, 150, 200") } else { NULL }; ss_inputs <- if(input$analysis_type == "Sample Size Search") { tagList(sliderInput("target_power", "Target Power", min = 0.5, max = 0.99, value = 0.8, step = 0.01), numericInput("n_start", "Start N", value = 50, min = 10), numericInput("n_step", "Step N", value = 25, min = 5)) } else { NULL }; bootstrap_inputs <- if(is_bootstrap) { tagList(hr(), h4("Bootstrap Options"), numericInput("n_sim", "Simulations", value = 500, min = 100, step = 100), numericInput("n_cores", "Parallel Cores", value = 1, min = 1)) } else { NULL }; tagList(power_inputs, ss_inputs, bootstrap_inputs) })
  
  # --- Kaplan-Meier Analysis ---
  km_results_reactive <- reactive({ req(pilot_data_reactive(), input$time_var, input$status_var, input$arm_var); df <- pilot_data_reactive(); df$time_col <- df[[input$time_var]]; df$status_col <- df[[input$status_var]]; df$arm_col <- as.factor(df[[input$arm_var]]); surv_formula <- as.formula("Surv(time_col, status_col) ~ arm_col"); km_fit <- survfit(surv_formula, data = df); log_rank_test <- survdiff(surv_formula, data = df); log_rank_p <- 1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1); km_plot <- ggsurvplot(km_fit, data = df, pval = FALSE, conf.int = TRUE, risk.table = TRUE, ggtheme = theme_minimal(), palette = c("#E69F00", "#56B4E9"), legend.title = "Arm", legend.labs = levels(df$arm_col))$plot + labs(title = "Kaplan-Meier Survival Curve for Pilot Data", subtitle = paste0("Log-Rank Test P-Value: ", format.pval(log_rank_p, digits = 3))); return(list(plot = km_plot, p_value = log_rank_p)) })
  
  # --- Core Analysis Logic ---
  analysis_results <- eventReactive(input$run_analysis, {
    validate(
      need(pilot_data_reactive(), "Please upload pilot data."),
      need(input$tau > 0, "RMST Truncation Time (τ) must be positive.")
    )
    if(input$analysis_type == "Power Calculation"){
      validate(need(input$sample_sizes != "", "Please enter at least one sample size."))
    }
    
    withProgress(message = 'Running Analysis', value = 0, {
      
      setProgress(0.1, detail = "Preparing arguments...")
      
      method_suffix <- if (!is.null(input$calc_method) && input$calc_method == "Bootstrap") "boot" else if (input$model_selection == "Semiparametric (GAM) Model") "boot" else "analytical"
      model_prefix <- switch(input$model_selection, "Linear IPCW Model"="linear", "Additive Stratified Model"="additive", "Multiplicative Stratified Model"="MS", "Dependent Censoring Model"="DC", "Semiparametric (GAM) Model"="GAM")
      func_type <- if(input$analysis_type == "Power Calculation") "power" else "ss"
      function_to_call_name <- paste(model_prefix, func_type, method_suffix, sep = ".")
      
      args <- list(
        pilot_data = pilot_data_reactive(), time_var = input$time_var, status_var = input$status_var,
        arm_var = input$arm_var, tau = input$tau, alpha = input$alpha
      )
      
      if (!is.null(input$linear_terms) && length(input$linear_terms) > 0) args$linear_terms <- input$linear_terms
      if (func_type == "power") {
        args$sample_sizes <- as.numeric(trimws(strsplit(input$sample_sizes, ",")[[1]]))
      } else {
        args$target_power <- input$target_power
        args$n_start <- input$n_start
        args$n_step <- input$n_step
      }
      if (model_prefix %in% c("additive", "MS")) { req(input$strata_var); args$strata_var <- input$strata_var }
      if (model_prefix == "DC") { req(input$dep_cens_var); args$dep_cens_status_var <- input$dep_cens_var }
      if (model_prefix == "GAM") { if (!is.null(input$smooth_terms) && length(input$smooth_terms) > 0) args$smooth_terms <- input$smooth_terms }
      if (method_suffix == "boot") { req(input$n_sim, input$n_cores); args$n_sim <- input$n_sim; args$parallel.cores <- input$n_cores }
      
      setProgress(0.4, detail = paste("Calling", function_to_call_name, "..."))
      
      results <- tryCatch(
        do.call(function_to_call_name, args),
        error = function(e) {
          showNotification(paste("Error:", e$message), type = "error", duration = NULL)
          return(NULL)
        }
      )
      
      if (!is.null(results)) {
        log_rank_p <- km_results_reactive()$p_value
        log_rank_row <- data.frame(Statistic = "Log-Rank Test P-Value (from pilot)", Value = log_rank_p)
        if (!is.null(results$results_summary)) {
          results$results_summary <- rbind(results$results_summary, log_rank_row)
        } else {
          results$results_summary <- log_rank_row
        }
      }
      
      setProgress(1, detail = "Done.")
      req(results)
      return(results)
    })
  })
  
  # --- Observer to Manage Dynamic Tabs ---
  observeEvent(analysis_results(), {
    res <- analysis_results(); req(res)
    
    lapply(dynamic_tabs_ids(), function(id) removeTab(inputId = "main_tabs", target = id))
    
    new_ids <- c()
    
    if (input$show_km) {
      km_tab_id <- "km_tab"
      insertTab(inputId = "main_tabs",
                tabPanel("Kaplan-Meier Curve", value = km_tab_id, plotlyOutput("km_plot", height = "700px")),
                target = "preview_tab", position = "after")
      new_ids <- c(new_ids, km_tab_id)
    }
    
    plot_tab_id <- "plot_tab"
    insertTab(inputId = "main_tabs",
              tabPanel("Plot Output", value = plot_tab_id, plotlyOutput("results_plot", height = "600px")),
              target = if(input$show_km) "km_tab" else "preview_tab", position = "after")
    new_ids <- c(new_ids, plot_tab_id)
    
    table_tab_id <- "table_tab"
    insertTab(inputId = "main_tabs",
              tabPanel("Results Table", value = table_tab_id, DT::dataTableOutput("results_table")),
              target = plot_tab_id, position = "after")
    new_ids <- c(new_ids, table_tab_id)
    
    summary_tab_id <- "summary_tab"
    insertTab(inputId = "main_tabs",
              tabPanel("Summary", value = summary_tab_id, uiOutput("summary_table_ui")),
              target = table_tab_id, position = "after")
    new_ids <- c(new_ids, summary_tab_id)
    
    report_tab_id <- "report_tab"
    insertTab(inputId = "main_tabs",
              tabPanel("Download Report", value = report_tab_id, downloadButton("download_report", "Download Report")),
              target = summary_tab_id, position = "after")
    new_ids <- c(new_ids, report_tab_id)
    
    dynamic_tabs_ids(new_ids)
    updateTabsetPanel(session, "main_tabs", selected = plot_tab_id)
  })
  
  # --- Render Outputs ---
  output$data_preview_table <- DT::renderDataTable({ req(pilot_data_reactive()); DT::datatable(pilot_data_reactive(), options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE, caption = "Preview of uploaded pilot data.") })
  output$km_plot <- renderPlotly({ req(input$show_km, km_results_reactive()); plotly::ggplotly(km_results_reactive()$plot, height = 700) })
  output$results_plot <- renderPlotly({ req(analysis_results()); plotly::ggplotly(analysis_results()$results_plot, tooltip = c("x", "y")) })
  output$results_table <- DT::renderDataTable({ df <- analysis_results()$results_data; req(df); DT::datatable(df, options=list(pageLength=10, scrollX=TRUE), rownames=FALSE, caption="Calculation Results") })
  output$summary_table_ui <- renderUI({ df <- analysis_results()$results_summary; req(df); df$Value <- round(df$Value, 4); df %>% kable("html", caption = "Summary") %>% kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% HTML() })
  
  # --- Report Generation ---
  output$download_report <- downloadHandler(
    filename = function() { paste0("RMSTdesign_report_", Sys.Date(), ".html") },
    content = function(file) {
      req(analysis_results())
      temp_report <- file.path(tempdir(), "report_template.Rmd")
      
      km_section <- if (input$show_km) {
        "## Kaplan-Meier Curve of Pilot Data\n```{r km_plot, fig.height=7, fig.width=9}\nprint(params$km_res$plot)\n```"
      } else { "" }
      
      rmd_string <- c("---", "title: 'RMSTdesign Analysis Report'", "output: html_document", "params:", "  analysis_res: NA", "  km_res: NA", "---",
                      "```{r setup, include=FALSE}\nknitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)\nlibrary(ggplot2); library(kableExtra)\n```",
                      "# Analysis Results", "## Summary", "```{r summary_table}\nsummary_df <- params$analysis_res$results_summary\nsummary_df$Value <- round(summary_df$Value, 4)\nsummary_df %>% kable('html') %>% kable_styling(full_width=F)\n```",
                      km_section,
                      "## Power / Sample Size Results", "### Results Table", "```{r results_table}\nres_df <- params$analysis_res$results_data\nres_df[[names(res_df)[2]]] <- round(res_df[[names(res_df)[2]]], 3)\nres_df %>% kable('html') %>% kable_styling(full_width=F)\n```",
                      "### Results Plot", "```{r results_plot}\nprint(params$analysis_res$results_plot)\n```")
      
      writeLines(rmd_string, temp_report)
      
      rmarkdown::render(temp_report, output_file = file, 
                        params = list(
                          analysis_res = analysis_results(), 
                          km_res = if(input$show_km) km_results_reactive() else NULL
                        ), 
                        envir = new.env(parent = globalenv()))
    }
  )
}

# --- 5. Run the Application ---
shinyApp(ui = ui, server = server)