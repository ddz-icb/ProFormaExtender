library(shiny)
library(shinyjs)
library(shinyBS)
library(readxl)
library(stringr)
library(httr)
library(jsonlite)

# zentrale Berechnung einbinden
source("core_calculation.R")

options(scipen = 999)
options(shiny.maxRequestSize = 100*1024^2)


ui <- fluidPage(
  useShinyjs(),
  
  tags$div(
    style = "display: flex; align-items: center; margin-bottom: 20px;",
    tags$img(
      src = "logo.png",
      height = "60px",
      style = "margin-right: 15px;"
    ),
    tags$div(
      tags$h2("ProForma Generation from Phospho Data", style = "margin: 0;"),
      tags$span(
        "Version 0.5.0",
        id = "version_label",
        style = "color: #555; font-size: 14px; cursor: pointer;",
        title = "This is the newest version"
      )
      
    )
  ),
  

  tags$head(
    tags$style(HTML("
      .shiny-file-input {
        border: 2px dashed #aaa;
        border-radius: 10px;
        padding: 30px;
        text-align: center;
        background-color: #fafafa;
        transition: 0.3s;
        font-size: 15px;
      }
      .shiny-file-input:hover {
        border-color: #007BFF;
        background-color: #f0f8ff;
      }
      .progress-message {
        font-style: italic;
        color: #555;
      }
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      fileInput(
        "excel",
        label = "📂 Drag & drop or select your Proteome Discoverer output file (Excel)",
        accept = c(".xlsx"),
        buttonLabel = "Browse...",
        placeholder = "No file selected"
      ),
      
      actionButton("start", "Start Calculation"),
      br(), br(),
      
      uiOutput("status_text"),
      br(),
      
      downloadButton("download", "Download Result (.tsv)")
    ),
    
    mainPanel(
      h4("Preview of calculated result"),
      fluidRow(
        column(
          width = 12,
          tags$div(
            style = "
              border: 1px solid #ccc;
              padding: 10px;
              border-radius: 8px;
              max-height: 600px; 
              overflow-y: auto; 
              overflow-x: auto;
              background: #fafafa;
              font-size: 12px;
              white-space: nowrap;
            ",
            tableOutput("preview")
          )
          
        )
      )
    )
  )
)


server <- function(input, output, session) {
  
  addTooltip(
    session,
    id = "version_label",
    title = "This is the newest version",
    placement = "right",
    trigger = "hover"
  )
  
  shinyjs::disable("download")
  result <- reactiveVal(NULL)
  
  output$preview <- renderTable({
    req(result())
    head(result(), 20)
  }, rownames = TRUE)
  
  observeEvent(input$start, {
    shinyjs::disable("download")
    req(input$excel)
    
    withProgress(message = "Starting calculation:", value = 0, {
      
      incProgress(0.1, detail = "Reading Excel...")
      
      # ---- Zentrale Berechnung ----
      incProgress(0.3, detail = "Running ProForma pipeline...")
      df <- run_proforma_pipeline(input$excel$datapath)
      
      incProgress(1, detail = "Finished")
      result(df)
    })
    
    output$status_text <- renderUI({
      tags$p("Calculation completed – result ready for download")
    })
    
    shinyjs::enable("download")
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste0(gsub("-", "", as.character(Sys.Date())), "_proforma-ext.output.tsv")
    },
    content = function(file) {
      write.table(result(), file, sep="\t", quote=FALSE, row.names=FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
