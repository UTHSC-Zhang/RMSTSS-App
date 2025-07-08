# Load necessary libraries
library(shiny)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(plotly)
library(kableExtra)
library(DT)
library(bslib)

# Source all R files in the R/ directory
r_files <- list.files(path = "R/", pattern = "\\.R$", full.names = TRUE)
sapply(r_files, source)
cat("All R scripts in the 'R/' directory have been sourced.\n")

# --- UI Definition ---
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  useShinyjs(), 
  
  titlePanel("RMSTdesign: Power and Sample Size Calculator"),
  
  sidebarLayout(
    # REVISED: Reorganized sidebar panel
    sidebarPanel(
      width = 4,
      
      wellPanel(
        h4("1. Setup & Model"),
        fileInput("pilot_data_upload", "Upload Pilot Data (.csv)", accept = ".csv"),
        selectInput("model_selection", "Select RMST Model",
                    choices = c("Linear IPCW Model", 
                                "Additive Stratified Model", 
                                "Multiplicative Stratified Model",
                                "Semiparametric (GAM) Model",
                                "Dependent Censoring Model"))
      ),
      
      wellPanel(
        h4("2. Column Mapping"),
        uiOutput("col_mapping_ui")
      ),

      wellPanel(
        h4("3. Analysis Parameters"),
        shinyjs::hidden(
          div(id = "analysis_params_panel",       
        # This fluidRow will contain the three inputs side-by-side
        fluidRow(
          column(4,
                 radioButtons("analysis_type", "Target Quantity",
                              choices = c("Power", "Sample Size"), 
                              selected = "Power")
          ),
          column(4,
                 # This will render the Analytical/Bootstrap radio buttons
                 uiOutput("method_selection_ui") 
          ),
          column(4,
                 numericInput("tau", "RMST Tau (τ)", value = 365, min = 1)
          )
        ), # End of fluidRow
        
        # The other inputs remain below
        uiOutput("analysis_inputs_ui"),
        sliderInput("alpha", "Significance Level (α)", min = 0.01, max = 0.1, value = 0.05, step = 0.01)
      )))
      ,
      uiOutput("bootstrap_options_ui"),
      
      hr(),
      actionButton("run_analysis", "Run Analysis", icon = icon("play"), class = "btn-primary btn-lg"),
      actionButton("reset_inputs", "Reset All", icon = icon("refresh"))
    ),
    
    mainPanel(
      width = 8,
      # REVISED: Reorganized tabset panel
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Instructions", 
                 h3("Welcome to RMSTdesign!"),
                 p("This application allows you to perform power and sample size calculations for clinical trials using Restricted Mean Survival Time (RMST)."),
                 tags$ol(
                   tags$li("Upload a CSV file containing your pilot study data."),
                   tags$li("Map the columns from your data to the required variables."),
                   tags$li("Select the desired statistical model and configure the analysis parameters."),
                   tags$li("Click the 'Run Analysis' button to see the results.")
                 )
        ),
        tabPanel("Data Preview", DT::dataTableOutput("data_preview_table")),
        tabPanel("Plot Output", plotlyOutput("results_plot", height = "600px")),
        
        tabPanel("Summary",
                 h4("Analysis Results"),
                 uiOutput("results_table_ui"),
                 hr(),
                 h4("Effect Size Summary"),
                 p("This summary is derived from the pilot data and forms the basis for the power/sample size calculation."),
                 uiOutput("summary_table_ui")
        ),
        
        tabPanel("Console Log", verbatimTextOutput("console_log_output")),
        
        tabPanel("Model Diagnostics",
                 h4("Diagnostic Plots for Pilot Data Models"),
                 p("These plots help assess the assumptions of the models fitted to your pilot data. ",
                   "If assumptions are violated, the power calculations may be less reliable."),
                 uiOutput("diagnostics_ui"),
                 hr(),
                 h4("Generate Analysis Report"),
                 p("Click the button below to download a complete HTML report of your analysis."),
                 downloadButton("download_report", "Download Report")
        )
      )
    )
  )
)

# --- Server Definition ---
server <- function(input, output, session) {
  # Enable live theming
  bslib::bs_themer()
  # --- A. Reactive Data and Dynamic UI ---
  pilot_data_reactive <- reactive({
    req(input$pilot_data_upload)
    tryCatch(read.csv(input$pilot_data_upload$datapath), error = function(e) {
      showNotification(paste("Error reading CSV:", e$message), type = "error"); NULL
    })
  })
  
  output$col_mapping_ui <- renderUI({
    df <- pilot_data_reactive(); req(df)
    column_names <- names(df)
    
    tagList(
      # --- ROW 1: Core Variables ---
      fluidRow(
        column(4, selectInput("time_var", "Time-to-Event", choices = column_names ,selected = column_names[2]   )),
        column(4, selectInput("status_var", "Status (1=event)", choices = column_names, selected = column_names[3])),
        column(4, selectInput("arm_var", "Treatment Arm (1=treat)", choices = column_names, selected = column_names[4]))
      ),
      
      # --- ROW 2: Conditional Model Variables ---
      fluidRow(
        # This column shows only for stratified models
        if (input$model_selection %in% c("Additive Stratified Model", "Multiplicative Stratified Model", "Semiparametric (GAM) Model")) {
          column(6, selectInput("strata_var", "Stratification Variable", choices = column_names, selected = column_names[5]))
        },
        
        # This column shows only for the dependent censoring model
        if (input$model_selection == "Dependent Censoring Model") {
          column(6, selectInput("dep_cens_var", "Dependent Censoring Status", choices = column_names, selected = column_names[5]))
        }
      ),
      
      # --- ROW 3: Covariates ---
      fluidRow(
        column(6, 
               selectizeInput("linear_terms", "Linear Covariates", choices = column_names, multiple = TRUE)
        ),
        
        # This column shows only for the GAM model
        column(6,
               if (input$model_selection == "Semiparametric (GAM) Model") {
                 selectizeInput("smooth_terms", "Non-Linear (Smooth) Covariates", choices = column_names, multiple = TRUE)
               }
        )
      )
    )
  })
  
  output$method_selection_ui <- renderUI({
    if (input$model_selection %in% c("Linear IPCW Model", "Multiplicative Stratified Model")) {
      radioButtons("calc_method", "Calculation Method", choices = c("Analytical", "Bootstrap"), selected = "Analytical")
    }
  })
  
  output$analysis_inputs_ui <- renderUI({
    if (input$analysis_type == "Power") {
      textInput("sample_sizes", "Sample Sizes (per arm/stratum, comma-separated)", value = "100, 150, 200")
    } else {
      sliderInput("target_power", "Target Power", min = 0.1, max = 1, value = 0.8, step = 0.01)
    }
  })
  
  # In app.R -> server function
  
  output$bootstrap_options_ui <- renderUI({
    
    is_bootstrap_choice <- !is.null(input$calc_method) &&
      input$calc_method == "Bootstrap" &&
      input$model_selection %in% c("Linear IPCW Model", "Multiplicative Stratified Model")
    
    is_always_bootstrap <- input$model_selection == "Semiparametric (GAM) Model"
    
    if (is_bootstrap_choice || is_always_bootstrap) {
      wellPanel(
        h4("Bootstrap Options"),
        
        # Row 1: Sample size search parameters (conditional)
        if (input$analysis_type == "Sample Size") {
          fluidRow(
            column(4, numericInput("n_start", "Start N", value = 50, min = 10)),
            column(4, numericInput("max_n_per_arm", "End N", value = 1000, min = 50)),
            column(4, numericInput("n_step", "Step", value = 25, min = 5))
          )
        },
        
        # Row 2: General simulation parameters
        fluidRow(
          column(6, numericInput("n_sim", "Simulations", value = 500, min = 100, step = 100)),
          column(6, numericInput("n_cores", "Parallel Cores", value = 1, min = 1))
        )
      )
    }
  })
  
  # --- B. UI Interactivity and Event Handlers ---
  
  # This single observer manages the visibility and enabled state of UI elements
  # based on whether data has been uploaded.
  observe({
    # This is the main condition: is data uploaded and valid?
    data_is_present <- !is.null(pilot_data_reactive()) && nrow(pilot_data_reactive()) > 0
    
    # Use shinyjs::toggle to show or hide the entire analysis parameters panel.
    shinyjs::toggle("analysis_params_panel", condition = data_is_present)
    
    # The "Run Analysis" button is only enabled if data is present AND a model is selected.
    run_button_enabled <- data_is_present && !is.null(input$model_selection) && input$model_selection != ""
    shinyjs::toggleState("run_analysis", condition = run_button_enabled)
  })
  
  # This observer handles the "Reset All" button click.
  observeEvent(input$reset_inputs, {
    # 1. Reset the file input control. This will trigger the observer below.
    shinyjs::reset("pilot_data_upload")
    
    # 2. Explicitly reset the value of the model selection dropdown to the default.
    updateSelectInput(session, "model_selection", selected = "Linear IPCW Model")
  })
  
  # This observer clears the output panels whenever the data or model changes.
  observeEvent(list(input$pilot_data_upload, input$model_selection), {
    run_output(list(results = NULL, log = "Analysis has not been run yet."))
  }, ignoreInit = TRUE) # ignoreInit prevents this from running when the app first starts.
  
  # --- C. Core Analysis Logic ---
  
  run_analysis_results <- reactive({
    
    validate(need(pilot_data_reactive(), "Please upload pilot data."))
    validate(need(input$time_var, "Please map Time-to-Event column."))
    validate(need(input$status_var, "Please map Status column."))
    validate(need(input$arm_var, "Please map Treatment Arm column."))
    
    analysis_results <- NULL
    
    log_text <- capture.output({
      analysis_results <- withProgress(message = 'Running Analysis', value = 0, {
        
        setProgress(0.1, detail = "Preparing arguments...")
        method_suffix <- if (!is.null(input$calc_method) && input$calc_method == "Bootstrap") "boot" else if (input$model_selection %in% c("Semiparametric (GAM) Model", "Additive Stratified Model")) "boot" else "analytical"
        model_prefix <- switch(input$model_selection, "Linear IPCW Model"="linear", "Additive Stratified Model"="additive", "Multiplicative Stratified Model"="MS", "Dependent Censoring Model"="DC", "Semiparametric (GAM) Model"="GAM")
        func_type <- if(input$analysis_type == "Power") "power" else "ss"
        function_to_call_name <- paste(model_prefix, func_type, method_suffix, sep = ".")
        
        args <- list(pilot_data = pilot_data_reactive(), time_var = input$time_var, status_var = input$status_var, arm_var = input$arm_var, tau = input$tau, alpha = input$alpha)
        if (!is.null(input$linear_terms) && length(input$linear_terms) > 0) args$linear_terms <- input$linear_terms
        if (input$analysis_type == "Power") {
          args$sample_sizes <- as.numeric(trimws(strsplit(input$sample_sizes, ",")[[1]]))
        } else {
          args$target_power <- input$target_power
        }
        if (model_prefix %in% c("additive", "MS", "GAM")) { req(input$strata_var); args$strata_var <- input$strata_var }
        if (model_prefix == "DC") { req(input$dep_cens_var); args$dep_cens_status_var <- input$dep_cens_var }
        if (model_prefix == "GAM" && !is.null(input$smooth_terms) && length(input$smooth_terms) > 0) { args$smooth_terms <- input$smooth_terms }
        if (method_suffix == "boot") { 
          req(input$n_sim, input$n_cores)
          args$n_sim <- input$n_sim
          args$parallel.cores <- input$n_cores
          
          # Add the sample size search parameters ONLY if it's a bootstrap SS search
          if(func_type == "ss"){
            req(input$n_start, input$max_n_per_arm, input$n_step)
            args$n_start <- input$n_start
            args$max_n_per_arm <- input$max_n_per_arm
            args$n_step <- input$n_step
          }
        }
        
        setProgress(0.4, detail = paste("Calling function:", function_to_call_name))
        
        tryCatch({
          do.call(function_to_call_name, args)
        }, error = function(e) {
          showNotification(paste("Error:", e$message), type = "error", duration = NULL); NULL
        })
      }) 
    }, type = c("output", "message"))
    
    return(list(results = analysis_results, log = paste(log_text, collapse = "\n")))
    
  }) %>% 
    bindCache(
      pilot_data_reactive(), input$model_selection, input$analysis_type,
      input$time_var, input$status_var, input$arm_var,
      input$linear_terms, input$strata_var, input$dep_cens_var, input$smooth_terms,
      input$tau, input$alpha, input$sample_sizes, input$target_power,
      input$n_sim, input$n_cores, input$calc_method
    ) %>% 
    bindEvent(input$run_analysis)
  
  run_output <- reactiveVal(list(results = NULL, log = "Analysis has not been run yet."))
  
  observe({
    run_output(run_analysis_results())
  })
  
  # --- D. Render Outputs ---
  output$data_preview_table <- DT::renderDataTable({
    req(pilot_data_reactive())
    DT::datatable(pilot_data_reactive(), options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE)
  })
  
  output$results_plot <- renderPlotly({
    req(run_output()$results$results_plot)
    plotly::ggplotly(run_output()$results$results_plot, tooltip = c("x", "y"))
  })
  
  output$results_table_ui <- renderUI({
    req(run_output()$results$results_data)
    run_output()$results$results_data %>% 
      kbl("html", caption = "Power/Sample Size Results") %>% 
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>% HTML()
  })
  
  output$summary_table_ui <- renderUI({
    req(run_output()$results$results_summary)
    run_output()$results$results_summary %>% 
      kbl("html") %>% 
      kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% HTML()
  })
  
  output$console_log_output <- renderText({
    req(run_output()$log)
    run_output()$log
  })
  
  output$diagnostics_ui <- renderUI({
    results <- run_output()$results; req(results)
    if (!is.null(results$fit_lm)) {
      tagList(h5("Linear Model (IPCW) Residuals"), plotOutput("diagnostic_plot_lm"))
    } else if (!is.null(results$fit_cens)) {
      tagList(
        h5("Censoring Model: Proportional Hazards Check"),
        p("This plot tests the proportional hazards assumption for the stratified Cox model fitted to the censoring distribution. A non-zero slope (i.e., a non-horizontal line) for a covariate suggests the assumption may be violated for that variable."),
        plotOutput("diagnostic_plot_coxzph")
      )
    } else {
      p("Diagnostic plots are not available for this model or the analysis has not been run.")
    }
  })
  
  output$diagnostic_plot_lm <- renderPlot({
    fit <- run_output()$results$fit_lm; req(fit)
    par(mfrow = c(2, 2)); plot(fit)
  })
  
  output$diagnostic_plot_coxzph <- renderPlot({
    fit <- run_output()$results$fit_cens; req(fit)
    par(mfrow = c(1, 1)); plot(survival::cox.zph(fit))
  })
  
  # --- E. Report Generation Logic ---
  output$download_report <- downloadHandler(
    filename = function() { paste0("RMSTdesign_report_", Sys.Date(), ".html") },
    content = function(file) {
      req(run_output()$results)
      withProgress(message = 'Generating Report...', value = 0, {
        
        incProgress(0.2, detail = "Gathering inputs...")
        method <- if (!is.null(input$calc_method) && input$calc_method == "Bootstrap") {"Bootstrap"} else if (input$model_selection %in% c("Semiparametric (GAM) Model", "Additive Stratified Model")) {"Bootstrap"} else {"Analytical"}
        
        report_inputs <- list(
          model_selection = input$model_selection, analysis_type = input$analysis_type,
          method = method, tau = input$tau, alpha = input$alpha,
          sample_sizes = input$sample_sizes, target_power = input$target_power,
          n_sim = input$n_sim, n_cores = input$n_cores)
        
        report_params <- list(inputs = report_inputs, results = run_output()$results, log = run_output()$log)
        
        incProgress(0.6, detail = "Rendering document...")
        rmarkdown::render("report_template.Rmd", output_file = file,
                          params = report_params, envir = new.env(parent = globalenv()))
      })
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)