# Load necessary libraries
library(shiny)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(plotly)
library(kableExtra)
library(DT)
library(RMSTdesign) # Your package

# --- UI Definition ---
ui <- fluidPage(
  theme = shinytheme("yeti"), # Apply a visual theme
  useShinyjs(), # Initialize shinyjs for the reset button
  
  titlePanel("RMSTdesign: Power and Sample Size Calculator"),
  
  sidebarLayout(
    # --- 1. Sidebar Panel for Inputs ---
    sidebarPanel(
      width = 4,
      h4("1. Setup"),
      fileInput("pilot_data_upload", "Upload Pilot Data (.csv)", accept = ".csv"),
      uiOutput("col_mapping_ui"),
      
      hr(),
      h4("2. Model Configuration"),
      selectInput("model_selection", "Select RMST Model",
                  choices = c("Linear IPCW Model", 
                              "Additive Stratified Model", 
                              "Multiplicative Stratified Model",
                              "Semiparametric (GAM) Model",
                              "Dependent Censoring Model")),
      
      uiOutput("method_selection_ui"),
      
      numericInput("tau", "RMST Truncation Time (τ)", value = 365, min = 1),
      tags$small(class = "form-text text-muted", "The time point at which the survival curve is truncated."),
      
      hr(),
      h4("3. Analysis Parameters"),
      radioButtons("analysis_type", "Select Analysis Type",
                   choices = c("Power Calculation", "Sample Size Search"), selected = "Power Calculation"),
      
      uiOutput("analysis_inputs_ui"),
      
      sliderInput("alpha", "Significance Level (α)", min = 0.01, max = 0.1, value = 0.05, step = 0.01),
      
      uiOutput("bootstrap_options_ui"),
      
      hr(),
      actionButton("run_analysis", "Run Analysis", icon = icon("play"), class = "btn-primary"),
      actionButton("reset_inputs", "Reset All Inputs", icon = icon("refresh"), class = "btn-secondary")
    ),
    
    # --- 2. Main Panel for Outputs ---
    mainPanel(
      width = 8,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Instructions", 
                 h3("Welcome to RMSTdesign!"),
                 p("This application allows you to perform power and sample size calculations for clinical trials using Restricted Mean Survival Time (RMST)."),
                 p("To get started:"),
                 tags$ol(
                   tags$li("Upload a CSV file containing your pilot study data."),
                   tags$li("Map the columns from your data to the required variables."),
                   tags$li("Select the desired statistical model and configure the analysis parameters."),
                   tags$li("Click the 'Run Analysis' button to see the results.")
                 )
        ),
        tabPanel("Data Preview", DT::dataTableOutput("data_preview_table")),
        tabPanel("Plot Output", plotlyOutput("results_plot", height = "600px")),
        tabPanel("Results Table", uiOutput("results_table_ui")),
        tabPanel("Summary", uiOutput("summary_table_ui")),
        tabPanel("Download Report",
                 h4("Generate Analysis Report"),
                 p("Click the button below to download a complete HTML report of your analysis, including inputs, results, and plots."),
                 downloadButton("download_report", "Download Report")
        )
      )
    )
  )
)

# --- Server Definition ---
server <- function(input, output, session) {
  
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
      h4("Column Mapping"),
      selectInput("time_var", "Time-to-Event", choices = column_names),
      selectInput("status_var", "Status (1=event)", choices = column_names),
      selectInput("arm_var", "Treatment Arm (1=treat)", choices = column_names),
      selectizeInput("linear_terms", "Linear Covariates (Optional)", choices = column_names, multiple = TRUE),
      if (input$model_selection %in% c("Additive Stratified Model", "Multiplicative Stratified Model")) {
        selectInput("strata_var", "Stratification Variable", choices = column_names)
      },
      if (input$model_selection == "Dependent Censoring Model") {
        selectInput("dep_cens_var", "Dependent Censoring Status", choices = column_names)
      },
      if (input$model_selection == "Semiparametric (GAM) Model") {
        selectizeInput("smooth_terms", "Non-Linear Covariates (GAM)", choices = column_names, multiple = TRUE)
      }
    )
  })
  
  output$method_selection_ui <- renderUI({
    if (input$model_selection %in% c("Linear IPCW Model", "Multiplicative Stratified Model")) {
      radioButtons("calc_method", "Calculation Method", choices = c("Analytical", "Bootstrap"), selected = "Analytical")
    }
  })
  
  output$analysis_inputs_ui <- renderUI({
    if (input$analysis_type == "Power Calculation") {
      textInput("sample_sizes", "Sample Sizes (per arm/stratum, comma-separated)", value = "100, 150, 200")
    } else {
      sliderInput("target_power", "Target Power", min = 0.5, max = 0.99, value = 0.8, step = 0.01)
    }
  })
  
  output$bootstrap_options_ui <- renderUI({
    is_bootstrap <- (input$model_selection %in% c("Linear IPCW Model", "Multiplicative Stratified Model") && !is.null(input$calc_method) && input$calc_method == "Bootstrap") ||
      (input$model_selection == "Semiparametric (GAM) Model")
    if (is_bootstrap) {
      tagList(hr(), h4("Bootstrap Options"),
              numericInput("n_sim", "Number of Simulations", value = 500, min = 100, step = 100),
              numericInput("n_cores", "Parallel Cores", value = 1, min = 1))
    }
  })
  
  # --- B. Reset Button Logic ---
  observeEvent(input$reset_inputs, {
    shinyjs::reset("pilot_data_upload")
    # You can add more resets for other specific inputs if needed
  })
  
  # --- C. Core Analysis Logic ---
  analysis_results <- eventReactive(input$run_analysis, {
    # Input Validation
    validate(need(pilot_data_reactive(), "Please upload pilot data."))
    validate(need(input$tau > 0, "RMST Truncation Time (τ) must be positive."))
    if(input$analysis_type == "Power Calculation"){
      validate(need(input$sample_sizes != "", "Please enter at least one sample size."))
    }
    
    # Use withProgress for a detailed progress bar
    withProgress(message = 'Running Analysis', value = 0, {
      
      setProgress(0.1, detail = "Preparing arguments...")
      
      # Determine method
      method_suffix <- if (!is.null(input$calc_method) && input$calc_method == "Bootstrap") "boot" else if (input$model_selection == "Semiparametric (GAM) Model") "boot" else "analytical"
      model_prefix <- switch(input$model_selection, "Linear IPCW Model"="linear", "Additive Stratified Model"="additive", "Multiplicative Stratified Model"="MS", "Dependent Censoring Model"="DC", "Semiparametric (GAM) Model"="GAM")
      func_type <- if(input$analysis_type == "Power Calculation") "power" else "ss"
      function_to_call_name <- paste(model_prefix, func_type, method_suffix, sep = ".")
      
      # Prepare arguments
      args <- list(pilot_data = pilot_data_reactive(), time_var = input$time_var, status_var = input$status_var, arm_var = input$arm_var, tau = input$tau, alpha = input$alpha)
      if (!is.null(input$linear_terms)) args$linear_terms <- input$linear_terms
      if (func_type == "power") { args$sample_sizes <- as.numeric(trimws(strsplit(input$sample_sizes, ",")[[1]])) } else { args$target_power <- input$target_power }
      if (model_prefix %in% c("additive", "MS")) { req(input$strata_var); args$strata_var <- input$strata_var }
      if (model_prefix == "DC") { req(input$dep_cens_var); args$dep_cens_status_var <- input$dep_cens_var }
      if (model_prefix == "GAM") { req(input$smooth_terms); args$smooth_terms <- input$smooth_terms }
      if (method_suffix == "boot") { req(input$n_sim); args$n_sim <- input$n_sim; args$parallel.cores <- input$n_cores }
      
      setProgress(0.4, detail = "Calling calculation function...")
      
      results <- tryCatch(do.call(function_to_call_name, args), error = function(e) {
        showNotification(paste("Error:", e$message), type = "error", duration = NULL); NULL
      })
      
      setProgress(1, detail = "Done.")
      req(results)
      return(results)
    })
  })
  
  # --- D. Render Outputs ---
  output$data_preview_table <- DT::renderDataTable({
    req(pilot_data_reactive())
    DT::datatable(pilot_data_reactive(), options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE, caption = "Preview of uploaded pilot data.")
  })
  
  output$results_plot <- renderPlotly({
    req(analysis_results())
    plotly::ggplotly(analysis_results()$results_plot, tooltip = c("x", "y"))
  })
  
  output$results_table_ui <- renderUI({
    df <- analysis_results()$results_data; req(df)
    df %>% kable("html") %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>% HTML()
  })
  
  output$summary_table_ui <- renderUI({
    df <- analysis_results()$results_summary; req(df)
    df %>% kable("html", caption = "Summary of Estimated Effect from Pilot Data") %>% kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% HTML()
  })
  
  # --- E. Report Generation Logic ---
  output$download_report <- downloadHandler(
    filename = function() { paste0("RMSTdesign_report_", Sys.Date(), ".html") },
    content = function(file) {
      # Ensure results are available before rendering
      req(analysis_results())
      
      # Create a temporary directory to work in
      temp_dir <- tempdir()
      temp_report <- file.path(temp_dir, "report_template.Rmd")
      
      # This is a placeholder for your report template file
      # In a real app, 'report_template.Rmd' would exist in your app's directory
      writeLines(c(
        "---",
        "title: 'RMSTdesign Analysis Report'",
        "output: html_document",
        "params:",
        "  results: NA",
        "---",
        "```{r setup, include=FALSE}",
        "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
        "library(ggplot2)",
        "library(kableExtra)",
        "```",
        "# Analysis Results",
        "## Summary of Estimated Effect",
        "```{r}",
        "params$results$results_summary %>% kable('html') %>% kable_styling(full_width=F)",
        "```",
        "## Results Table",
        "```{r}",
        "params$results$results_data %>% kable('html') %>% kable_styling(full_width=F)",
        "```",
        "## Results Plot",
        "```{r}",
        "print(params$results$results_plot)",
        "```"
      ), temp_report)
      
      # Render the report
      rmarkdown::render(temp_report, output_file = file,
                        params = list(results = analysis_results()),
                        envir = new.env(parent = globalenv())
      )
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)