library(shiny)
library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
library(bslib) # Using unified bslib for a robust UI

# --- STABLE ONLINE VISUALIZATION APP ---
# Run this file after running the updated Seg_data.Rmd.
# It allows dynamic gene exploration and labeling comparison.

data_file <- "../Random data/VisiumHD_Drosophila_Resolution_List.rds"

# --- UI Definition ---
ui <- page_sidebar(
  title = "Drosophila Midgut: Spatial FOV Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE, # Allow vertical scrolling so plots aren't squashed
  
  sidebar = sidebar(
    title = "Controls",
    # GENE SEARCHBAR
    selectizeInput("gene_choice", "Search Gene (FlyBase Symbol):",
                   choices = NULL, # Populated on server
                   options = list(
                     placeholder = 'Type molecular marker (e.g. Dl, nub, odd)...',
                     maxOptions = 1000,
                     openOnFocus = FALSE
                   )),
    hr(),
    # LABEL TOGGLE
    selectInput("label_choice", "Labeling System (Violin Plot):",
                choices = c("Single-Cell Transfer" = "Predicted_CellType")),
    hr(),
    helpText("How to use:"),
    tags$small(
      tags$ul(
        tags$li("Search for a gene to see its spatial 'hotspots'."),
        tags$li("Check the violin plots to see expression across cell types."),
        tags$li("The map shows high-res cell boundaries (polygons).")
      )
    ),
    verbatimTextOutput("status_info")
  ),
  
  # MAIN DISPLAY
  # Increased heights and removed 'full_screen' dependence
  layout_columns(
    col_widths = 12,
    card(
      height = "800px", # Forced height for FOV
      card_header("Spatial FOV Expression (Cell Boundaries)"),
      plotOutput("fov_plot", height = "750px")
    ),
    card(
      height = "500px", # Forced height for Violine
      card_header("Cell-Type Expression Analysis"),
      plotOutput("vln_plot", height = "450px")
    )
  )
)

# --- Server Logic ---
server <- function(input, output, session) {
  
  # Reactive object storage
  obj_reactive <- reactiveVal(NULL)
  
  # Data Initialization
  observe({
    if (file.exists(data_file)) {
      withProgress(message = 'Loading Annotated Data...', value = 0.5, {
        vlist <- readRDS(data_file)
        obj <- vlist$seg
        
        # Fallback for missing labels
        if(!"Predicted_CellType" %in% colnames(obj@meta.data)) obj$Predicted_CellType <- "Unlabeled"
        
        obj_reactive(obj)
        
        # Populate the searchbar
        updateSelectizeInput(session, "gene_choice", 
                            choices = sort(rownames(obj)), 
                            server = TRUE)
      })
    } else {
      showNotification("Data file not found! Please run Seg_data.Rmd first.", type = "error")
    }
  })
  
  # 1. FOV Spatial Plot
  output$fov_plot <- renderPlot({
    req(obj_reactive(), input$gene_choice)
    ImageFeaturePlot(obj_reactive(), features = input$gene_choice, 
                     fov = "fov", max.cutoff = "q95", size = 0.8) +
      scale_fill_viridis_c(option = "magma", name = "Log1p") +
      theme_void() +
      labs(title = paste("Spatial Distribution of", input$gene_choice))
  })
  
  # 2. Violin Plot (Comparison Mode)
  output$vln_plot <- renderPlot({
    req(obj_reactive(), input$gene_choice, input$label_choice)
    VlnPlot(obj_reactive(), features = input$gene_choice, 
            group.by = input$label_choice, pt.size = 0.1) +
      scale_fill_viridis_d(option = "viridis") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            axis.title.x = element_blank()) +
      labs(subtitle = paste("Grouping by:", input$label_choice))
  })
  
  output$status_info <- renderText({
    req(obj_reactive())
    paste("Dataset Ready:", ncol(obj_reactive()), "cells.")
  })
}

# Run
shinyApp(ui, server)
