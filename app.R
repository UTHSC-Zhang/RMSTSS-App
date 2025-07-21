# Load necessary libraries
packages <- c("shiny", "shinyjs", "bslib", "DT", "ggplot2", "plotly", "survival", "survminer", "kableExtra", "magrittr", "rmarkdown")
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
})

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
              fluidRow(
                column(4,
                       radioButtons("analysis_type", "Target Quantity",
                                    choices = c("Power", "Sample Size"),
                                    selected = "Power")
                ),
                column(4,
                       uiOutput("method_selection_ui")
                ),
                column(4,
                       numericInput("L", "RMST L (τ)", value = 365, min = 1)
                )
              ),
              
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
        
        tabPanel("Plot Output",
                 h4("Kaplan-Meier Survival Plot (Pilot Data)"),
                 p("This plot shows the survival probability over time for each treatment arm. You can hover over the lines for details."),
                 plotlyOutput("survival_plotly_output", height = "500px"),
                 hr(),
                 h4("Power vs. Sample Size Curve"),
                 p("This plot shows the calculated power for different sample sizes, or the required sample size for the target power."),
                 plotlyOutput("results_plot", height = "500px")
        ),
        
        tabPanel("Summary",
                 h4("Analysis Results"),
                 uiOutput("results_table_ui"),
                 hr(),
                 h4("Effect Size Summary"),
                 p("This summary is derived from the pilot data and forms the basis for the power/sample size calculation."),
                 uiOutput("summary_table_ui"),
                 
                 uiOutput("logrank_summary_ui"),
                 
                 uiOutput("download_button_ui")
        ),
        
        tabPanel("Console Log", verbatimTextOutput("console_log_output"))
      )
    )
  )
)

# --- Server Definition ---
#' @importFrom bslib bs_themer
#' @importFrom DT datatable renderDataTable
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_light ylim
#' @importFrom kableExtra kbl kable_styling
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly renderPlotly
#' @importFrom rmarkdown render
#' @importFrom shiny reactive req showNotification renderUI tagList fluidRow
#' @importFrom shiny column selectInput selectizeInput radioButtons textInput
#' @importFrom shiny sliderInput wellPanel h4 numericInput observe observeEvent
#' @importFrom shiny updateSelectInput reactiveVal validate need withProgress
#' @importFrom shiny setProgress renderText downloadButton downloadHandler
#' @importFrom shiny removeNotification hr p bindCache bindEvent
#' @importFrom shinyjs toggle toggleState reset
#' @importFrom stats as.formula complete.cases pchisq sd na.omit
#' @importFrom survival survfit Surv survdiff
#' @importFrom survminer ggsurvplot
server <- function(input, output, session) {
  bslib::bs_themer()
  
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
      fluidRow(
        column(4, selectInput("time_var", "Time-to-Event", choices = column_names ,selected = column_names[2]   )),
        column(4, selectInput("status_var", "Status (1=event)", choices = column_names, selected = column_names[3])),
        column(4, selectInput("arm_var", "Treatment Arm (1=treat)", choices = column_names, selected = column_names[4]))
      ),
      fluidRow(
        if (input$model_selection %in% c("Additive Stratified Model", "Multiplicative Stratified Model", "Semiparametric (GAM) Model")) {
          column(6, selectInput("strata_var", "Stratification Variable", choices = column_names, selected = column_names[5]))
        },
        if (input$model_selection == "Dependent Censoring Model") {
          column(6, selectInput("dep_cens_var", "Dependent Censoring Status", choices = column_names, selected = column_names[5]))
        }
      ),
      fluidRow(
        column(6,
               selectizeInput("linear_terms", "Linear Covariates", choices = column_names, multiple = TRUE)
        ),
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
  
  output$bootstrap_options_ui <- renderUI({
    is_bootstrap_choice <- !is.null(input$calc_method) &&
      input$calc_method == "Bootstrap" &&
      input$model_selection %in% c("Linear IPCW Model", "Multiplicative Stratified Model")
    is_always_bootstrap <- input$model_selection %in% c("Semiparametric (GAM) Model", "Additive Stratified Model")
    
    if (is_bootstrap_choice || is_always_bootstrap) {
      wellPanel(
        h4("Bootstrap Options"),
        if (input$analysis_type == "Sample Size") {
          fluidRow(
            column(4, numericInput("n_start", "Start N", value = 50, min = 10)),
            column(4, numericInput("max_n_per_arm", "End N", value = 1000, min = 50)),
            column(4, numericInput("n_step", "Step", value = 25, min = 5))
          )
        },
        fluidRow(
          column(6, numericInput("n_sim", "Simulations", value = 500, min = 100, step = 100)),
          column(6, numericInput("n_cores", "Parallel Cores", value = 1, min = 1))
        )
      )
    }
  })
  
  observe({
    data_is_present <- !is.null(pilot_data_reactive()) && nrow(pilot_data_reactive()) > 0
    shinyjs::toggle("analysis_params_panel", condition = data_is_present)
    run_button_enabled <- data_is_present && !is.null(input$model_selection) && input$model_selection != ""
    shinyjs::toggleState("run_analysis", condition = run_button_enabled)
  })
  
  observeEvent(input$reset_inputs, {
    shinyjs::reset("pilot_data_upload")
    updateSelectInput(session, "model_selection", selected = "Linear IPCW Model")
  })
  
  observeEvent(list(input$pilot_data_upload, input$model_selection), {
    run_output(list(results = NULL, log = "Analysis has not been run yet."))
  }, ignoreInit = TRUE)
  
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
        function_to_call_name <- paste(model_prefix, func_type, method_suffix,"app", sep = ".")
        
        args <- list(pilot_data = pilot_data_reactive(), time_var = input$time_var, status_var = input$status_var, arm_var = input$arm_var, L = input$L, alpha = input$alpha)
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
          
          if(func_type == "ss"){
            req(input$n_start, input$max_n_per_arm, input$n_step)
            args$n_start <- input$n_start
            args$max_n_per_arm <- input$max_n_per_arm
            args$n_step <- input$n_step
          }
        }
        
        setProgress(0.4, detail = paste("Calling function:", function_to_call_name))
        
        main_calc_results <- tryCatch({
          do.call(function_to_call_name, args)
        }, error = function(e) {
          showNotification(paste("Error:", e$message), type = "error", duration = NULL); NULL
        })
        
        setProgress(0.9, detail = "Performing survival analysis...")
        
        logrank_summary_df <- NULL
        analysis_data_for_plot <- NULL
        tryCatch({
          cat("\n\n--- Survival Analysis on Pilot Data ---\n")
          
          analysis_data <- data.frame(
            time = pilot_data_reactive()[[input$time_var]],
            status = as.numeric(pilot_data_reactive()[[input$status_var]]),
            arm = as.factor(pilot_data_reactive()[[input$arm_var]])
          )
          analysis_data_for_plot <- analysis_data
          
          fixed_formula <- as.formula("Surv(time, status) ~ arm")
          
          logrank_test <- survdiff(fixed_formula, data = analysis_data)
          cat("Log-Rank Test:\n")
          print(logrank_test)
          
          p_value <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
          logrank_summary_df <- data.frame(
            Statistic = "Chi-Square",
            Value = round(logrank_test$chisq, 3),
            DF = length(logrank_test$n) - 1,
            `P-Value` = format.pval(p_value, eps = .001, digits = 3)
          )
          
        }, error = function(e) {
          cat("\nError during survival analysis:", e$message, "\n")
        })
        
        main_calc_results$logrank_summary <- logrank_summary_df
        main_calc_results$analysis_data_for_plot <- analysis_data_for_plot
        
        return(main_calc_results)
      })
    }, type = c("output", "message"))
    
    return(list(results = analysis_results, log = paste(log_text, collapse = "\n")))
    
  }) %>%
    bindCache(
      pilot_data_reactive(), input$model_selection, input$analysis_type,
      input$time_var, input$status_var, input$arm_var,
      input$linear_terms, input$strata_var, input$dep_cens_var, input$smooth_terms,
      input$L, input$alpha, input$sample_sizes, input$target_power,
      input$n_sim, input$n_cores, input$calc_method
    ) %>%
    bindEvent(input$run_analysis)
  
  run_output <- reactiveVal(list(results = NULL, log = "Analysis has not been run yet."))
  
  observe({
    run_output(run_analysis_results())
  })
  
  output$data_preview_table <- DT::renderDataTable({
    req(pilot_data_reactive())
    DT::datatable(pilot_data_reactive(), options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE)
  })
  
  output$survival_plotly_output <- renderPlotly({
    req(run_output()$results$analysis_data_for_plot, input$alpha)
    plot_data <- run_output()$results$analysis_data_for_plot
    
    # Create the survival fit object without the conf.level argument
    fit <- survfit(Surv(time, status) ~ arm, data = plot_data)
    
    # Pass the confidence level directly to ggsurvplot
    p <- ggsurvplot(
      fit,
      data = plot_data,
      conf.int = TRUE,
      conf.int.alpha = 0.3,
      conf.int.style = "ribbon",
      conf.level = 1 - input$alpha, # Pass conf.level here
      palette = c("#007BFF", "#D9534F"),
      legend.title = input$arm_var,
      xlab = paste("Time (in units of '", input$time_var, "')"),
      ylab = "Survival Probability",
      ggtheme = theme_light()
    )
    
    ggplotly(p$plot)
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
  
  output$logrank_summary_ui <- renderUI({
    req(run_output()$results$logrank_summary)
    tagList(
      hr(),
      h4("Log-Rank Test Summary (Pilot Data)"),
      p("This is a non-parametric test comparing the survival distributions of the treatment arms from the pilot data."),
      run_output()$results$logrank_summary %>%
        kbl("html", booktabs = TRUE) %>%
        kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>%
        HTML()
    )
  })
  
  output$console_log_output <- renderText({
    req(run_output()$log)
    run_output()$log
  })
  
  output$download_button_ui <- renderUI({
    req(run_output()$results)
    tagList(
      hr(),
      h4("Generate Analysis Report"),
      p("Click the button below to download a complete PDF report of your analysis."),
      downloadButton("download_report", "Download PDF Report")
    )
  })
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("RMSTdesign_report_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(run_output()$results)
      
      # Show a notification that the report is being generated
      id <- showNotification("Generating PDF report... This may take a moment.",
                             duration = NULL, closeButton = FALSE, type = "message")
      on.exit(removeNotification(id), add = TRUE)
      
      # Gather all inputs for the report parameters
      report_inputs <- list(
        model_selection = input$model_selection,
        analysis_type = input$analysis_type,
        time_var = input$time_var,
        status_var = input$status_var,
        arm_var = input$arm_var,
        L = input$L,
        alpha = input$alpha,
        sample_sizes = input$sample_sizes,
        target_power = input$target_power
      )
      
      # Create the list of parameters to pass to the Rmd file
      report_params <- list(
        inputs = report_inputs,
        results = run_output()$results,
        log = run_output()$log
      )
      
      # Render the Rmd file
      rmarkdown::render(
        "report_template.Rmd",
        output_file = file,
        params = report_params,
        envir = new.env(parent = globalenv()) # Use a clean environment
      )
    }
  )
}

shinyApp(ui = ui, server = server)