library(shiny)
library(Seurat)
library(ggplot2)
library(viridis)
library(plotly)
library(dplyr)
library(bslib) # Using unified bslib for a robust UI

# --- STABLE ONLINE VISUALIZATION APP ---
# Run this file after running the updated Seg_data.Rmd.
# It allows dynamic gene exploration and labeling comparison.

data_file <- "/home/mattia/Thesis/03Slide4/Random data/VisiumHD_Drosophila_Resolution_List.rds"

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
    ),
    card(
      height = "800px",
      card_header("3D Midgut Expression Model"),
      plotlyOutput("mesh3d_plot", height = "750px")
    )
  )
)

# --- Server Logic ---
server <- function(input, output, session) {
  
  # Reactive object storage
  obj_reactive <- reactiveVal(NULL)
  obj3d_reactive <- reactiveVal(NULL)
  mesh_reactive <- reactiveVal(NULL)
  
  # Data Initialization
  observe({
    if (file.exists(data_file)) {
      withProgress(message = 'Loading Annotated Data...', value = 0.5, {
        vlist <- readRDS(data_file)
        obj <- vlist$seg
        
        # Fallback for missing labels
        if(!"Predicted_CellType" %in% colnames(obj@meta.data)) obj$Predicted_CellType <- "Unlabeled"
        
        obj_reactive(obj)
        
        # Load 3D Object
        obj3d_path <- "/home/mattia/Thesis/03Slide4/Random data/merged_guts.rds"
        if(file.exists(obj3d_path)) {
          obj3d_reactive(readRDS(obj3d_path))
        }
        
        # Load 3D Mesh Cache
        mesh_cache <- "/home/mattia/Thesis/03Slide4/figures/3D_mesh_plots_starvation/tube_mesh_cache.rds"
        if (file.exists(mesh_cache)) {
          mesh_reactive(readRDS(mesh_cache))
        }
        
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
  
  # 3. 3D Mesh Plotly
  output$mesh3d_plot <- renderPlotly({
    req(obj3d_reactive(), input$gene_choice, mesh_reactive())
    
    mesh_data <- mesh_reactive()
    obj3d <- obj3d_reactive()
    
    # Extract expression and pseudotime from the 3D object
    gene <- input$gene_choice
    if (!gene %in% rownames(obj3d)) return(NULL)
    
    data_df <- FetchData(obj3d, vars = c(gene, "Spatial_Pseudotime"))
    colnames(data_df) <- c("Expression", "Pseudotime")
    data_df <- data_df[!is.na(data_df$Pseudotime), ]
    
    # Compute smoothed expression for the mesh vertices
    loess_model <- tryCatch({
      loess(Expression ~ Pseudotime, data = data_df, span = 0.3)
    }, error = function(e) NULL)
    
    num_vertices <- nrow(mesh_data$v_matrix)
    vertex_smooth_expr <- rep(0, num_vertices)
    
    if (!is.null(loess_model)) {
      predicted_expr <- predict(loess_model, newdata = data.frame(Pseudotime = mesh_data$vertex_pseudotime[mesh_data$final_mask]))
      predicted_expr[predicted_expr < 0] <- 0
      vertex_smooth_expr[mesh_data$final_mask] <- predicted_expr
    }
    
    # Color mapping
    pal <- colorRamp(c("#440154", "#3b528b", "#21908c", "#5dc863", "#fde725"))
    min_expr <- min(vertex_smooth_expr[mesh_data$final_mask], na.rm = TRUE)
    max_expr <- max(vertex_smooth_expr[mesh_data$final_mask], na.rm = TRUE)
    
    vertex_colors <- rep("rgb(220, 220, 220)", num_vertices)
    
    if (max_expr > min_expr) {
      norm_expr <- (vertex_smooth_expr[mesh_data$final_mask] - min_expr) / (max_expr - min_expr)
      cols <- pal(norm_expr)
      cols[is.na(cols)] <- 220
      vertex_colors[mesh_data$final_mask] <- apply(cols, 1, function(x) paste0("rgb(", x[1], ",", x[2], ",", x[3], ")"))
    }
    
    fig <- plot_ly(
      x = mesh_data$v_matrix[, 1],
      y = mesh_data$v_matrix[, 2],
      z = mesh_data$v_matrix[, 3],
      i = mesh_data$f_matrix[, 1],
      j = mesh_data$f_matrix[, 2],
      k = mesh_data$f_matrix[, 3],
      type = "mesh3d",
      vertexcolor = vertex_colors,
      name = gene
    )
    
    fig %>% layout(
      title = paste("3D Spatial Expression:", gene),
      scene = list(
        aspectmode = "data",
        xaxis = list(visible = FALSE),
        yaxis = list(visible = FALSE),
        zaxis = list(visible = FALSE)
      )
    )
  })
  
  output$status_info <- renderText({
    req(obj_reactive())
    paste("Dataset Ready:", ncol(obj_reactive()), "cells.")
  })
}

# Run
shinyApp(ui, server)
