library(shiny)
library(shinyjs)
library(shinyBS)
library(readxl)
library(stringr)
library(httr)
library(jsonlite)

source("core_calculation.R")

# ── Python-Konverter aufrufen (MQ → PD-Format) ───────────────────────────────
run_mq_converter <- function(modpep_path, phospho_path,
                             peptides_path = NULL,
                             keep_reverse = FALSE,
                             keep_contaminants = FALSE) {
  
  script <- "mq_to_pd_converter.py"
  outfile <- tempfile(fileext = ".xlsx")
  
  args <- c(
    script,
    "--modpep",  shQuote(modpep_path),
    "--phospho", shQuote(phospho_path),
    "--output",  shQuote(outfile)
  )
  if (!is.null(peptides_path))
    args <- c(args, "--peptides", shQuote(peptides_path))
  if (keep_reverse)
    args <- c(args, "--keep-reverse")
  if (keep_contaminants)
    args <- c(args, "--keep-contaminants")
  
  ret <- system2("python3", args, stdout = TRUE, stderr = TRUE)
  attr(ret, "status") <- attr(ret, "status") %||% 0L
  
  list(log = paste(ret, collapse = "\n"), outfile = outfile)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

options(scipen = 999)
options(shiny.maxRequestSize = 200 * 1024^2)

# ═══════════════════════════════════════════════════════════════════════════════
# UI
# ═══════════════════════════════════════════════════════════════════════════════

ui <- fluidPage(
  useShinyjs(),
  
  # ── Styles ────────────────────────────────────────────────────────────────
  tags$head(tags$style(HTML("

    /* ---------- Basis ---------- */
    body { font-family: 'Segoe UI', Arial, sans-serif; background: #f7f8fa; }

    /* ---------- Header ---------- */
    .app-header {
      display: flex; align-items: center; gap: 14px;
      padding: 18px 0 14px;
      border-bottom: 2px solid #e0e4ec;
      margin-bottom: 24px;
    }
    .app-header img  { height: 54px; }
    .app-header h2   { margin: 0; font-size: 1.4rem; color: #1a2640; font-weight: 700; }
    .app-header .ver { font-size: 12px; color: #8896aa; margin-top: 3px; cursor: pointer; }

    /* ---------- Input-Mode-Selector ---------- */
    .mode-selector {
      display: flex; gap: 12px; margin-bottom: 24px;
    }
    .mode-card {
      flex: 1;
      border: 2px solid #d0d8e8;
      border-radius: 10px;
      padding: 16px 18px;
      cursor: pointer;
      background: #fff;
      transition: border-color .2s, box-shadow .2s, background .2s;
      user-select: none;
    }
    .mode-card:hover { border-color: #4a7fd4; box-shadow: 0 2px 8px #4a7fd420; }
    .mode-card.active {
      border-color: #2e5fa3;
      background: #eef4ff;
      box-shadow: 0 2px 10px #2e5fa325;
    }
    .mode-card .icon { font-size: 26px; margin-bottom: 6px; }
    .mode-card .title { font-weight: 700; font-size: 15px; color: #1a2640; margin-bottom: 3px; }
    .mode-card .desc  { font-size: 12px; color: #6a7a90; line-height: 1.4; }

    /* ---------- Section boxes ---------- */
    .section-box {
      background: #fff;
      border: 1px solid #dce3ef;
      border-radius: 10px;
      padding: 20px 22px;
      margin-bottom: 18px;
    }
    .section-title {
      font-weight: 700; font-size: 13px; color: #4a5a78;
      text-transform: uppercase; letter-spacing: .05em;
      margin-bottom: 14px; padding-bottom: 8px;
      border-bottom: 1px solid #eef0f5;
    }

    /* ---------- File inputs ---------- */
    .shiny-input-container label { font-weight: 600; font-size: 13px; color: #2e3d55; }

    /* ---------- Converter log ---------- */
    #converter_log {
      font-family: 'Consolas', monospace; font-size: 11px;
      background: #1a1f2e; color: #a8d8a8;
      border-radius: 6px; padding: 12px;
      max-height: 180px; overflow-y: auto;
      white-space: pre-wrap; margin-top: 8px;
    }

    /* ---------- Status text ---------- */
    .status-ok  { color: #2e7d32; font-weight: 600; }
    .status-err { color: #b00020; font-weight: 600; }

    /* ---------- Preview table ---------- */
    .preview-wrap {
      border: 1px solid #dce3ef; border-radius: 8px;
      padding: 12px; background: #fafbfd;
      max-height: 560px; overflow: auto;
      font-size: 11.5px; white-space: nowrap;
    }

    /* ---------- Buttons ---------- */
    .btn-primary { background: #2e5fa3; border-color: #2e5fa3; font-weight: 600; }
    .btn-primary:hover { background: #1e4a8a; border-color: #1e4a8a; }
    .btn-success { font-weight: 600; }

    /* ---------- Divider ---------- */
    .step-arrow {
      text-align: center; font-size: 22px; color: #b0bdd0;
      margin: -4px 0 6px; line-height: 1;
    }

  "))),
  
  # ── Header ────────────────────────────────────────────────────────────────
  tags$div(
    class = "app-header",
    tags$img(src = "logo.jpg"),
    tags$div(
      tags$h2("ProForma Generation from Phospho Data"),
      tags$div("Version 0.6.0", class = "ver", id = "version_label",
               title = "Supports direct MaxQuant input (v0.6)")
    )
  ),
  
  # ── Layout ────────────────────────────────────────────────────────────────
  sidebarLayout(
    
    # ── Sidebar ─────────────────────────────────────────────────────────────
    sidebarPanel(
      width = 4,
      
      # ── Mode selector ────────────────────────────────────────────────────
      tags$div(class = "section-box",
               tags$div("Input type", class = "section-title"),
               
               tags$div(class = "mode-selector",
                        
                        tags$div(
                          id = "card_pd", class = "mode-card active",
                          onclick = "Shiny.setInputValue('input_mode', 'pd', {priority: 'event'})",
                          
                          tags$div("Proteome Discoverer", class = "title"),
                          tags$div("Upload a PD Excel file", class = "desc")
                        ),
                        
                        tags$div(
                          id = "card_mq", class = "mode-card",
                          onclick = "Shiny.setInputValue('input_mode', 'mq', {priority: 'event'})",
                          
                          tags$div("MaxQuant", class = "title"),
                          tags$div("Convert MQ output files first, then run ProForma", class = "desc")
                        )
               )
      ),
      
      # ── PD panel ─────────────────────────────────────────────────────────
      conditionalPanel("input.input_mode === 'pd' || !input.input_mode",
                       
                       tags$div(class = "section-box",
                                tags$div("Proteome Discoverer input", class = "section-title"),
                                fileInput("excel", "Excel file (.xlsx)",
                                          accept      = ".xlsx",
                                          multiple    = FALSE,
                                          buttonLabel = "Browse…"
                                )
                       )
      ),
      
      # ── MQ panel ─────────────────────────────────────────────────────────
      conditionalPanel("input.input_mode === 'mq'",
                       
                       tags$div(class = "section-box",
                                tags$div("MaxQuant files", class = "section-title"),
                                
                                fileInput("mq_modpep",
                                          label       = "modificationSpecificPeptides.txt",
                                          accept      = ".txt",
                                          multiple    = FALSE,
                                          buttonLabel = "Browse…"
                                ),
                                fileInput("mq_phospho",
                                          label       = "Phospho (STY) Sites.txt",
                                          accept      = ".txt",
                                          multiple    = FALSE,
                                          buttonLabel = "Browse…"
                                ),
                                fileInput("mq_peptides",
                                          label       = "peptides.txt ",
                                          accept      = ".txt",
                                          multiple    = FALSE,
                                          buttonLabel = "Browse…"
                                ),
                                
                                tags$hr(style = "margin: 12px 0;"),
                                tags$div(style = "font-size:12px; color:#5a6a80; margin-bottom:8px;",
                                         "Filter options"),
                                checkboxInput("mq_keep_reverse",
                                              "Keep reverse-database hits", value = FALSE),
                                checkboxInput("mq_keep_contaminants",
                                              "Keep potential contaminants", value = FALSE),
                                
                                tags$div(class = "step-arrow", "▼"),
                                
                                actionButton("run_converter", "Convert MQ → Excel",
                                             class = "btn btn-warning btn-block",
                                             style = "width:100%; font-weight:700;"),
                                
                                uiOutput("converter_status_ui"),
                                uiOutput("converter_log_ui"),
                                uiOutput("download_converted_ui")
                       )
      ),
      
      # ── ProForma section (always visible once file ready) ─────────────────
      tags$div(class = "step-arrow", "▼"),
      
      tags$div(class = "section-box",
               tags$div("ProForma calculation", class = "section-title"),
               
               uiOutput("ready_file_hint"),
               
               actionButton("start", "Start Calculation",
                            class = "btn btn-primary btn-block",
                            style = "width:100%;"),
               br(),
               uiOutput("status_text"),
               br(),
               uiOutput("download_ui")
      )
    ),
    
    # ── Main panel ────────────────────────────────────────────────────────────
    mainPanel(
      width = 8,
      h4("Preview of calculated result",
         style = "color:#2e3d55; font-weight:700; margin-bottom:14px;"),
      tags$div(class = "preview-wrap",
               tableOutput("preview")
      )
    )
  ),
  
  # ── JS: card highlight logic ──────────────────────────────────────────────
  tags$script(HTML("
    $(document).on('shiny:inputchanged', function(e) {
      if (e.name === 'input_mode') {
        $('.mode-card').removeClass('active');
        $('#card_' + e.value).addClass('active');
      }
    });
  "))
)


# ═══════════════════════════════════════════════════════════════════════════════
# Server
# ═══════════════════════════════════════════════════════════════════════════════

server <- function(input, output, session) {
  
  addTooltip(session, "version_label", "v0.6: MaxQuant input supported",
             placement = "right", trigger = "hover")
  
  # ── Reactive state ──────────────────────────────────────────────────────────
  result            <- reactiveVal(NULL)
  converted_xlsx    <- reactiveVal(NULL)   # temp path of MQ-converted file
  converter_success <- reactiveVal(FALSE)
  
  # ── Helper: which Excel path to use for the pipeline? ─────────────────────
  active_excel_path <- reactive({
    mode <- input$input_mode %||% "pd"
    if (mode == "mq") {
      converted_xlsx()          # NULL until conversion done
    } else {
      req(input$excel)
      input$excel$datapath
    }
  })
  
  # ── Hint above the "Start Calculation" button ─────────────────────────────
  output$ready_file_hint <- renderUI({
    mode <- input$input_mode %||% "pd"
    if (mode == "mq") {
      if (isTRUE(converter_success())) {
        tags$p(class = "status-ok", style = "font-size:12px; margin-bottom:8px;",
               "✔ Converted file ready — click Start Calculation")
      } else {
        tags$p(style = "font-size:12px; color:#8896aa; margin-bottom:8px;",
               "Convert MaxQuant files first (above), then run the calculation.")
      }
    } else {
      if (!is.null(input$excel)) {
        tags$p(class = "status-ok", style = "font-size:12px; margin-bottom:8px;",
               paste0("✔ ", input$excel$name))
      } else {
        tags$p(style = "font-size:12px; color:#8896aa; margin-bottom:8px;",
               "Upload a PD Excel file above.")
      }
    }
  })
  
  # ══════════════════════════════════════════════════════════════════════════
  # MQ Converter
  # ══════════════════════════════════════════════════════════════════════════
  
  converter_log_rv <- reactiveVal("")
  
  observeEvent(input$run_converter, {
    req(input$mq_modpep, input$mq_phospho)
    
    converter_success(FALSE)
    converted_xlsx(NULL)
    result(NULL)
    converter_log_rv("")
    
    withProgress(message = "Converting MaxQuant files…", value = 0, {
      incProgress(0.2, detail = "Reading modificationSpecificPeptides.txt")
      
      res <- tryCatch({
        run_mq_converter(
          modpep_path       = input$mq_modpep$datapath,
          phospho_path      = input$mq_phospho$datapath,
          peptides_path     = if (!is.null(input$mq_peptides)) input$mq_peptides$datapath else NULL,
          keep_reverse      = isTRUE(input$mq_keep_reverse),
          keep_contaminants = isTRUE(input$mq_keep_contaminants)
        )
      }, error = function(e) {
        list(log = e$message, outfile = NULL)
      })
      
      incProgress(1, detail = "Done")
      
      converter_log_rv(res$log)
      
      if (!is.null(res$outfile) && file.exists(res$outfile)) {
        converted_xlsx(res$outfile)
        converter_success(TRUE)
      }
    })
  })
  
  output$converter_status_ui <- renderUI({
    req(input$run_converter > 0)
    if (isTRUE(converter_success())) {
      tags$p(class = "status-ok", style = "margin-top:10px;",
             "✔ Conversion successful")
    } else if (input$run_converter > 0) {
      tags$p(class = "status-err", style = "margin-top:10px;",
             "✖ Conversion failed — see log below")
    }
  })
  
  output$converter_log_ui <- renderUI({
    req(input$run_converter > 0)
    tagList(
      tags$div("Conversion log", style = "font-size:11px; color:#7a8a9a; margin-top:10px;"),
      tags$div(id = "converter_log",
               textOutput("converter_log_text", inline = TRUE)
      )
    )
  })
  
  output$converter_log_text <- renderText({
    converter_log_rv()
  })
  
  # Download converted Excel
  output$download_converted_ui <- renderUI({
    req(isTRUE(converter_success()))
    tagList(
      br(),
      downloadButton("download_converted",
                     "Download converted Excel",
                     class = "btn btn-default btn-sm btn-block",
                     style = "width:100%;")
    )
  })
  
  output$download_converted <- downloadHandler(
    filename = function() {
      paste0(gsub("-", "", Sys.Date()), "_mq_converted.xlsx")
    },
    content = function(file) {
      req(converted_xlsx())
      file.copy(converted_xlsx(), file)
    }
  )
  
  # ══════════════════════════════════════════════════════════════════════════
  # ProForma Pipeline
  # ══════════════════════════════════════════════════════════════════════════
  
  observeEvent(input$start, {
    shinyjs::disable("download")
    result(NULL)
    output$status_text <- renderUI(NULL)
    
    excel_path <- active_excel_path()
    if (is.null(excel_path)) {
      showModal(modalDialog(
        title = "No input file",
        tags$p("Please upload a Proteome Discoverer file, or convert a MaxQuant dataset first."),
        footer = tagList(modalButton("Close")), easyClose = TRUE
      ))
      return()
    }
    
    withProgress(message = "Running ProForma pipeline:", value = 0, {
      out <- tryCatch({
        incProgress(0.1, detail = "Parsing input file…")
        incProgress(0.3, detail = "Querying UniProt API…")
        df <- run_proforma_pipeline(excel_path)
        incProgress(1,   detail = "Done")
        df
      }, error = function(e) e)
    })
    
    if (inherits(out, "error")) {
      showModal(modalDialog(
        title = "Calculation failed",
        tags$p(style = "color:#b00020; font-weight:bold;",
               "An error occurred during the calculation."),
        tags$p(out$message),
        footer = tagList(modalButton("Close")), easyClose = TRUE
      ))
      return()
    }
    
    result(out)
    output$status_text <- renderUI(
      tags$p(class = "status-ok", "✔ Calculation complete — ready for download")
    )
    shinyjs::enable("download")
  })
  
  # Preview
  output$preview <- renderTable({
    req(result())
    head(result(), 20)
  }, rownames = TRUE)
  
  # Download result
  output$download <- downloadHandler(
    filename = function() {
      paste0(gsub("-", "", Sys.Date()), "_proforma_output.xlsx")
    },
    content = function(file) {
      req(result())
      writexl::write_xlsx(result(), file)
    }
  )
  
  output$download_ui <- renderUI({
    req(result())
    downloadButton("download", "Download Result (.xlsx)", class = "btn btn-success btn-block",
                   style = "width:100%;")
  })
}


shinyApp(ui = ui, server = server)
