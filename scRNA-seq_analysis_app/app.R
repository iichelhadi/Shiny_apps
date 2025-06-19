# app.R - Main application file for scRNA-seq Shiny app

# Source global.R to load packages and global functions
source("global.R")
source("data-io.R")
source("preprocessing.R")
source("utilities.R")
source("visualization.R")
source("analysis.R")

# UI Definition --------------------------------------------------------------
ui <- fluidPage(
  theme = shinytheme("superhero"),
  
  # App title
  titlePanel("scRNA-seq Analysis App"),
  
  # Navbar layout
  navbarPage(
    "scRNA-seq Explorer",
    
    # Tab 1: Data Input -----------------------------------------------------
    tabPanel(
      "Data Input",
      sidebarLayout(
        sidebarPanel(
          # Input options
          shinyDirButton("directory_input", "Choose 10X Directory", 
                         "Select directory",class = "btn-primary"),
          tags$hr(),
          actionButton("load_data_btn", "Load Data", class = "btn-primary"),          
          tags$hr(),
          
          # QC filtering parameters
          h4("Quality Control Parameters"),
          numericInput(
            "percent_mt", "Max % mitochondrial:", 
            value = 10, min = 0, max = 100, step = 1
          ),
          numericInput(
            "min_features", "Min genes per cell:", 
            value = 500, min = 0, step = 100
          ),
          numericInput(
            "max_features", "Max genes per cell:", 
            value = 3000, min = 0, step = 100
          ),
          numericInput(
            "min_counts", "Min UMI counts per cell:", 
            value = 500, min = 0, step = 100
          ),
          numericInput(
            "max_counts", "Max UMI counts per cell:", 
            value = 10000, min = 0, step = 100
          ),
          numericInput(
            "min_genes_per_umi", "Min complexity (log10GenesPerUMI):", 
            value = 0.8, min = 0, max = 1, step = 0.05
          ),
          
          # Doublet detection
          radioButtons(
            "doublets", "Detect doublets:",
            choices = c("Yes" = "yes", "No" = "no"),
            selected = "yes"
          ),
          
          # Only show remove doublets option if detection is enabled
          conditionalPanel(
            condition = "input.doublets == 'yes'",
            checkboxInput("remove_doublets", "Remove doublets", value = TRUE)
          ),
          
          actionButton("filter_btn", "Apply Filters", class = "btn-primary"),
        ),
        
        mainPanel(
          # Data summary
          h3("Selected Directory"),
          verbatimTextOutput("directory_output"),
          
          tags$hr(),
          
          h3("Quality Control Plots"),
          plotOutput("qc_plot", height = "500px", width = "700px"),
          
          tags$hr(),
          
          h3("After Filtering"),
          plotOutput("filtered_qc_plot", height = "500px", width = "700px")
        )
      )
    ),
    
    # Tab 2: Preprocessing ------------------------------------------------------
    tabPanel(
      "Preprocessing",
      sidebarLayout(
        sidebarPanel(
          # Variable features
          numericInput(
            "n_features", "Number of variable features:", 
            value = 2000, min = 500, step = 100
          ),
          
          # Point size for plots
          numericInput(
            "pt_size", "Point size for plots:", 
            value = 1, min = 0.1, max = 5, step = 0.1
          ),
          
          # Plot dimensions
          h4("Plot Dimensions"),
          textInput(
            "dim_plot_width", "Width (e.g., '500px'):", 
            value = "1200px"
          ),
          textInput(
            "dim_plot_height", "Height (e.g., '500px'):", 
            value = "500px"
          ),
          
          actionButton("run_dim_reduction_btn", "Run Dimension Reduction", class = "btn-primary"),
          
          tags$hr(),
          
          # Clustering options
          numericInput(
            "resolution", "Clustering resolution:", 
            value = 0.5, min = 0.1, max = 2, step = 0.1
          ),
          selectInput(
            "algorithm", "Clustering algorithm:",
            choices = list(
              "Louvain" = 1, 
              "Louvain with multilevel refinement" = 2,
              "SLM" = 3, 
              "Leiden" = 4
            ),
            selected = 4
          ),
          
          actionButton("run_clustering_btn", "Run Clustering", class = "btn-primary")
        ),
        
        mainPanel(
          tabsetPanel(
            tabPanel(
              "Dimension Reduction", 
              plotOutput("dim_plot", height = "600px")
            ),
            tabPanel(
              "Clustering Results", 
              plotOutput("cluster_plot", height = "600px")
            ),
            tabPanel(
              "Quality Metrics",
              plotOutput("qc_metrics_plot", height = "600px")
            )
          )
        )
      )
    ),
    
    # Tab 3: Marker Genes -------------------------------------------------------
    tabPanel(
      "Marker Genes",
      sidebarLayout(
        sidebarPanel(
          # Marker type
          radioButtons(
            "marker_type", "Which markers to compute:",
            choices = c("All clusters" = "all", "Specific clusters" = "specific"),
            selected = "all"
          ),
          
          # Statistical options
          numericInput(
            "logfc_threshold", "Log FC threshold:", 
            value = 0.25, min = 0, max = 5, step = 0.05
          ),
          
          radioButtons(
            "only_pos", "Return only positive markers:",
            choices = c("Yes" = "TRUE", "No" = "FALSE"),
            selected = "TRUE"
          ),
          
          # Specific clusters (conditional)
          conditionalPanel(
            condition = "input.marker_type == 'specific'",
            selectInput(
              "test_cluster", "Test cluster:",
              choices = NULL
            ),
            selectInput(
              "reference_cluster", "Reference cluster:",
              choices = NULL
            )
          ),
          
          # Number of genes for heatmap
          h4("Heatmap Options"),
          numericInput(
            "n_genes_heatmap", "Number of genes per cluster:", 
            value = 5, min = 1, max = 20, step = 1
          ),
          numericInput(
            "heatmap_height", "Height (pixels):", 
            value = 800, min = 600, max = 2000, step = 50
          ),
          numericInput(
            "heatmap_width", "Width (pixels):", 
            value = 1000, min = 600, max = 2000, step = 50
          ),
          
          actionButton("find_markers_btn", "Find Markers", class = "btn-primary"),
          
          tags$hr(),
          
          downloadButton("download_markers_btn", "Download Markers")
        ),
        
        mainPanel(
          tags$style(HTML("
            .dataTables_wrapper .dataTables_length, 
            .dataTables_wrapper .dataTables_filter, 
            .dataTables_wrapper .dataTables_info, 
            .dataTables_wrapper .dataTables_processing, 
            .dataTables_wrapper .dataTables_paginate {
              color: #ffffff;
            }
            thead {
              color: #ffffff;
            }
            tbody {
              color: #ffffff;
            }
            .bootstrap4 {
              color: #ffffff !important;
            }
          ")),
          
          h3("Differentially Expressed Genes"),
          DT::dataTableOutput("markers_table"),
          
          tags$hr(),
          
          h3("Marker Gene Heatmap"),
          plotOutput("marker_heatmap", height = "800px", width = "1000px")
        )
      )
    ),
    
    # Tab 4: Cell Type Annotation -----------------------------------------------
    tabPanel(
      "Cell Annotation",
      sidebarLayout(
        sidebarPanel(
          h4("Manual Annotation"),
          helpText("Assign cell types to clusters manually:"),
          
          # Dynamic UI will be inserted here for cluster annotations
          uiOutput("cluster_annotation_ui"),
          
          actionButton("apply_manual_annotation", "Apply Annotations", class = "btn-primary"),
          
          tags$hr(),
          
          downloadButton("download_annotated_obj_btn", "Download Annotated Object")
        ),
        
        mainPanel(
          plotOutput("annotation_plot", height = "600px"),
          
          tags$hr(),
          
          h4("Cell Type Summary"),
          DT::dataTableOutput("celltype_summary_table")
        )
      )
    ),
    
    # Tab 5: Visualization ------------------------------------------------------
    tabPanel(
      "Visualization",
      sidebarLayout(
        sidebarPanel(
          # Grouping variable selection
          h4("Plot Settings"),
          selectInput(
            "group_by_viz", "Color/Group by:",
            choices = c("Clusters" = "seurat_clusters"),
            selected = "seurat_clusters"
          ),
          
          # Gene selection
          h4("Gene Expression"),
          selectizeInput(
            "gene_select", "Select Gene(s):",
            choices = NULL,
            multiple = TRUE,
            options = list(maxItems = 6)
          ),
          
          # Plot type
          radioButtons(
            "plot_type", "Plot type:",
            choices = c("Feature plot" = "feature", "Violin plot" = "violin", "UMAP by groups" = "dimplot"),
            selected = "feature"
          ),
          
          # Plot dimensions
          h4("Plot Dimensions"),
          textInput(
            "viz_plot_width", "Width (e.g., '500px'):", 
            value = "700px"
          ),
          textInput(
            "viz_plot_height", "Height (e.g., '500px'):", 
            value = "600px"
          )
        ),
        
        mainPanel(
          plotOutput("gene_plot", height = "600px", width = "700px")
        )
      )
    ),
    
    # Tab 6: About --------------------------------------------------------------
    tabPanel(
      "About",
      fluidRow(
        column(
          8, offset = 2,
          h2("About scRNA-seq Explorer"),
          
          h3("Overview"),
          p("scRNA-seq Explorer is a user-friendly Shiny application for analyzing single-cell RNA sequencing data using the Seurat framework. This application provides a graphical interface for common analysis tasks, making it accessible to researchers without extensive programming experience."),
          
          h3("Features"),
          tags$ul(
            tags$li("Data import from 10X Genomics format"),
            tags$li("Quality control and filtering"),
            tags$li("Doublet detection"),
            tags$li("Normalization and scaling"),
            tags$li("Dimensionality reduction (PCA, UMAP)"),
            tags$li("Clustering"),
            tags$li("Gene visualization"),
            tags$li("Marker gene identification"),
            tags$li("Cell type annotation")
          ),
          
          h3("Citation"),
          p("If you use this application in your research, please cite:"),
          tags$ul(
            tags$li("Seurat: Stuart et al. (2019). Comprehensive Integration of Single-Cell Data. Cell, 177(7), 1888-1902."),
            tags$li("scDblFinder: Germain et al. (2021). scDblFinder: a pipeline for filtering cells with artifactual doublet-like expression profiles. Bioinformatics, 37(19), 3333-3335."),
            tags$li("MAST: Finak et al. (2015). MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome Biology, 16, 278.")
          ),
          
          h3("Source Code"),
          p("This application is open source. The source code is available at:"),
          p(a("GitHub Repository", href = "https://github.com/iichelhadi//Shiny_apps/scRNA-seq_analysis_app", target = "_blank")),
          
          tags$hr(),
          
          p("Developed by Elhadi Iich, 2025"),
          p("License: GPL3")
        )
      )
    )
  )
)

# Server Definition ----------------------------------------------------------
server <- function(input, output, session) {
  
  # Reactive values to store Seurat objects and other data
  values <- reactiveValues(
    seurat_obj = NULL,
    seurat_filtered = NULL,
    seurat_processed = NULL,
    marker_genes = NULL,
    current_reduction = NULL,
    cluster_annotations = NULL
  )
  
  # Directory selection
  shinyDirChoose(input, "directory_input", roots = c(home = '~'))
  
  # Parse selected directory path
  selected_dir <- reactive({
    req(input$directory_input)
    parseDirPath(roots = c(home = '~'), input$directory_input)
  })
  
  # Display selected directory
  output$directory_output <- renderText({
    req(selected_dir())
    paste("Selected directory:", selected_dir())
  })
  
  # Load 10X data
  observeEvent(input$load_data_btn, {
    req(selected_dir())
    withProgress(message = 'Loading data...', value = 0, {
      tryCatch({
        # Read 10X data
        values$seurat_obj <- read_10x(selected_dir())
        
        # Update gene choices for feature visualization
        updateSelectizeInput(
          session, 
          "gene_select",
          choices = rownames(values$seurat_obj),
          server = TRUE
        )
        
        # Plot QC metrics
        output$qc_plot <- renderPlot({
          req(values$seurat_obj)
          VlnPlot(
            values$seurat_obj,
            features = c(
              'nFeature_RNA',
              'nCount_RNA',
              'percent.mt',
              'percent.ribo',
              "percent.hsp",
              'log10GenesPerUMI'
            )
          )
        })
      }, error = function(e) {
        showNotification(
          paste("Error loading data:", e$message),
          type = "error",
          duration = NULL
        )
      })
    })
  })
  
  # Filter data
  observeEvent(input$filter_btn, {
    req(values$seurat_obj)
    withProgress(message = 'Filtering data...', value = 0, {
      tryCatch({
        # Apply basic QC filtering
        setProgress(value = 0.3, message = "Applying QC filters...")
        values$seurat_filtered <- subset(
          values$seurat_obj,
          subset = percent.mt < input$percent_mt &
            log10GenesPerUMI > input$min_genes_per_umi &
            nFeature_RNA > input$min_features &
            nFeature_RNA < input$max_features &
            nCount_RNA > input$min_counts &
            nCount_RNA < input$max_counts
        )
        
        # Detect doublets if requested
        if (input$doublets == "yes") {
          setProgress(value = 0.6, message = "Detecting doublets...")
          # Run doublet detection
          values$seurat_filtered <- doublet_detection(values$seurat_filtered, col = "orig.ident")
          
          # Count cells before and after doublet filtering
          total_cells <- ncol(values$seurat_filtered)
          doublet_cells <- sum(values$seurat_filtered$scDblFinder.class == "doublet")
          singlet_cells <- sum(values$seurat_filtered$scDblFinder.class == "singlet")
          doublet_percent <- round(doublet_cells / total_cells * 100, 1)
          
          # Show notification about doublets
          showNotification(
            paste("Detected", doublet_cells, "doublets (", doublet_percent, "%) out of", total_cells, "cells"),
            type = "message",
            duration = 10
          )
          
          # Add QC plot showing doublets
          output$filtered_qc_plot <- renderPlot({
            req(values$seurat_filtered)
            # Create layout with both standard QC metrics and doublet-specific plots
            p1 <- VlnPlot(
              values$seurat_filtered,
              features = c(
                'nFeature_RNA',
                'nCount_RNA',
                'percent.mt',
                'percent.ribo',
                "percent.hsp",
                'log10GenesPerUMI'
              ),
              ncol = 3
            )
            
            # Create violin plots grouped by doublet status
            p2 <- VlnPlot(
              values$seurat_filtered,
              features = c('nFeature_RNA', 'nCount_RNA', 'log10GenesPerUMI'),
              group.by = "scDblFinder.class",
              ncol = 3
            )
            
            # Combine plots 
            p1 / p2
          })
        } else {
          # Regular QC plot without doublets
          output$filtered_qc_plot <- renderPlot({
            req(values$seurat_filtered)
            VlnPlot(
              values$seurat_filtered,
              features = c(
                'nFeature_RNA',
                'nCount_RNA',
                'percent.mt',
                'percent.ribo',
                "percent.hsp",
                'log10GenesPerUMI'
              )
            )
          })
        }
        
        # Count cells before and after filtering
        original_cells <- ncol(values$seurat_obj)
        remaining_cells <- ncol(values$seurat_filtered)
        percent_kept <- round(remaining_cells / original_cells * 100, 1)
        
        # Show notification about filtering
        showNotification(
          paste("Kept", remaining_cells, "cells (", percent_kept, "%) out of", original_cells, "cells after filtering"),
          type = "message",
          duration = 10
        )
        
      }, error = function(e) {
        showNotification(
          paste("Error filtering data:", e$message),
          type = "error",
          duration = NULL
        )
      })
    })
  })
  
  # Process data (dimensionality reduction)
  observeEvent(input$run_dim_reduction_btn, {
    req(values$seurat_filtered)
    withProgress(message = 'Running dimensionality reduction...', value = 0, {
      tryCatch({
        # Process data - use filtered data directly if doublets were already handled
        if (input$doublets == "yes" && "scDblFinder.class" %in% colnames(values$seurat_filtered@meta.data)) {
          # Apply doublet filtering if needed
          if (input$remove_doublets) {
            setProgress(value = 0.3, message = "Removing doublets...")
            data_subset <- subset(values$seurat_filtered, subset = scDblFinder.class == "singlet")
            # Show notification about filtering
            showNotification(
              paste("Removed", 
                    sum(values$seurat_filtered$scDblFinder.class == "doublet"), 
                    "doublets. Proceeding with", 
                    ncol(data_subset), 
                    "cells."),
              type = "message",
              duration = 8
            )
            # Run dimensionality reduction
            setProgress(value = 0.6, message = "Running dimension reduction...")
            values$seurat_processed <- dimred(data_subset, nfeatures = input$n_features)
          } else {
            # Keep doublets and just run dimension reduction
            setProgress(value = 0.6, message = "Running dimension reduction (including doublets)...")
            values$seurat_processed <- dimred(values$seurat_filtered, nfeatures = input$n_features)
          }
        } else {
          # If doublets weren't detected or won't be excluded, just run dimension reduction
          setProgress(value = 0.6, message = "Running dimension reduction...")
          values$seurat_processed <- dimred(values$seurat_filtered, nfeatures = input$n_features)
        }
        
        # Update current reduction
        values$current_reduction <- "umap"
        
        # Extract the number of PCs used
        pct <- values$seurat_processed[["pca"]]@stdev / sum(values$seurat_processed[["pca"]]@stdev) * 100
        cumu <- cumsum(pct)
        pcs_used <- which(cumu > 70 & pct < 5)[1]
        
        # Store the number of PCs for display
        values$pcs_used <- pcs_used
        
        # Show notification about PCs used
        showNotification(
          paste("Using", pcs_used, "principal components based on elbow method (70% variance, <5% contribution)"),
          type = "message",
          duration = 10
        )
        
        # Plot dimension reduction results
        output$dim_plot <- renderPlot({
          req(values$seurat_processed)
          
          # If we have doublet information and did not remove doublets, visualize doublets on UMAP
          if (input$doublets == "yes" && !input$remove_doublets && 
              "scDblFinder.class" %in% colnames(values$seurat_processed@meta.data)) {
            # Using standard Seurat functions
            p1 <- Seurat::DimPlot(
              values$seurat_processed,
              reduction = 'umap',
              group.by = "scDblFinder.class",
              pt.size = input$pt_size
            )
          } else {
            # Plain UMAP without grouping
            p1 <- Seurat::DimPlot(
              values$seurat_processed,
              reduction = 'umap',
              pt.size = input$pt_size
            )
          }
          
          # QC metrics plot - use separately to avoid palette issues
          p2 <- Seurat::FeaturePlot(
            values$seurat_processed,
            reduction = 'umap',
            pt.size = input$pt_size,
            features = c('nFeature_RNA')
          )
          
          p3 <- Seurat::FeaturePlot(
            values$seurat_processed,
            reduction = 'umap',
            pt.size = input$pt_size,
            features = c('nCount_RNA')
          )
          
          p4 <- Seurat::FeaturePlot(
            values$seurat_processed,
            reduction = 'umap',
            pt.size = input$pt_size,
            features = c('percent.mt')
          )
          
          p5 <- Seurat::FeaturePlot(
            values$seurat_processed,
            reduction = 'umap',
            pt.size = input$pt_size,
            features = c('log10GenesPerUMI')
          )
          
          # Combine plots separately to avoid palette conflicts
          p1 | ((p2 / p3) | (p4 / p5))
        }, height = function() {
          # Dynamic height based on input
          height_val <- as.numeric(gsub("px", "", input$dim_plot_height))
          if (is.na(height_val) || height_val < 400) height_val <- 600
          return(height_val)
        }, width = function() {
          # Dynamic width based on input
          width_val <- as.numeric(gsub("px", "", input$dim_plot_width))
          if (is.na(width_val) || width_val < 400) width_val <- 900
          return(width_val)
        })
        
        # Add an elbow plot to explain PC selection
        output$qc_metrics_plot <- renderPlot({
          req(values$seurat_processed)
          
          # Extract PCA stats for elbow plot
          pct <- values$seurat_processed[["pca"]]@stdev / sum(values$seurat_processed[["pca"]]@stdev) * 100
          
          # Create a layout with plots
          p1 <- ElbowPlot(values$seurat_processed, ndims = 50) +
            geom_vline(xintercept = values$pcs_used, linetype="dashed", color = "red") +
            annotate("text", x = values$pcs_used + 2, y = max(pct), 
                     label = paste("PC cutoff =", values$pcs_used), color = "red")
          
          # If we have doublet information, show that on the UMAP
          if (input$doublets == "yes" && "scDblFinder.class" %in% colnames(values$seurat_processed@meta.data)) {
            # Use standard Seurat functions
            p2 <- Seurat::DimPlot(
              values$seurat_processed,
              reduction = 'umap',
              group.by = "scDblFinder.class",
              pt.size = input$pt_size
            )
            
            # Split feature plots to avoid palette issues
            p3 <- Seurat::FeaturePlot(
              values$seurat_processed,
              reduction = 'umap',
              pt.size = input$pt_size,
              features = c('nFeature_RNA')
            )
            
            p4 <- Seurat::FeaturePlot(
              values$seurat_processed,
              reduction = 'umap',
              pt.size = input$pt_size,
              features = c('nCount_RNA')
            )
            
            # Combine all plots
            (p1 / p2) | (p3 / p4)
            
          } else {
            # Default QC metrics if no doublet info - split feature plots
            p2 <- Seurat::FeaturePlot(
              values$seurat_processed,
              reduction = 'umap',
              pt.size = input$pt_size,
              features = c('nFeature_RNA')
            )
            
            p3 <- Seurat::FeaturePlot(
              values$seurat_processed,
              reduction = 'umap',
              pt.size = input$pt_size,
              features = c('nCount_RNA')
            )
            
            p4 <- Seurat::FeaturePlot(
              values$seurat_processed,
              reduction = 'umap',
              pt.size = input$pt_size,
              features = c('percent.mt')
            )
            
            p5 <- Seurat::FeaturePlot(
              values$seurat_processed,
              reduction = 'umap',
              pt.size = input$pt_size,
              features = c('log10GenesPerUMI')
            )
            
            # Combine plots in a 2x3 layout
            (p1 | p2) / (p3 | p4) / p5
          }
        }, height = function() {
          # Dynamic height based on input
          height_val <- as.numeric(gsub("px", "", input$dim_plot_height))
          if (is.na(height_val) || height_val < 400) height_val <- 800
          return(height_val)
        }, width = function() {
          # Dynamic width based on input
          width_val <- as.numeric(gsub("px", "", input$dim_plot_width))
          if (is.na(width_val) || width_val < 400) width_val <- 900
          return(width_val)
        })
        
      }, error = function(e) {
        showNotification(
          paste("Error in dimensionality reduction:", e$message),
          type = "error",
          duration = NULL
        )
      })
    })
  })
  
  # Run clustering
  observeEvent(input$run_clustering_btn, {
    req(values$seurat_processed)
    withProgress(message = 'Running clustering...', value = 0, {
      tryCatch({
        # Run clustering
        values$seurat_processed <- clustering(
          values$seurat_processed,
          resolution = input$resolution,
          algorithm = input$algorithm
        )
        
        # Store number of clusters
        num_clusters <- length(levels(values$seurat_processed$seurat_clusters))
        
        # Show notification with the number of clusters found
        showNotification(
          paste("Found", num_clusters, "clusters at resolution", input$resolution),
          type = "message",
          duration = 10
        )
        
        # Plot clustering results
        output$cluster_plot <- renderPlot({
          req(values$seurat_processed)
          Seurat::DimPlot(
            values$seurat_processed,
            reduction = 'umap',
            group.by = "seurat_clusters",
            pt.size = input$pt_size,
            label = TRUE
          )
        }, height = function() {
          # Dynamic height based on input
          height_val <- as.numeric(gsub("px", "", input$dim_plot_height))
          if (is.na(height_val) || height_val < 400) height_val <- 600
          return(height_val)
        }, width = function() {
          # Dynamic width based on input
          width_val <- as.numeric(gsub("px", "", input$dim_plot_width))
          if (is.na(width_val) || width_val < 400) width_val <- 800
          return(width_val)
        })
        
        # Update cluster choices for marker gene analysis
        updateSelectInput(
          session,
          "test_cluster",
          choices = levels(values$seurat_processed$seurat_clusters),
          selected = levels(values$seurat_processed$seurat_clusters)[1]
        )
        
        updateSelectInput(
          session,
          "reference_cluster",
          choices = levels(values$seurat_processed$seurat_clusters),
          selected = levels(values$seurat_processed$seurat_clusters)[2]
        )
        
        # Update visualization grouping choices
        current_choices <- c("Clusters" = "seurat_clusters")
        # Check if cell types already exist
        if ("cell_type" %in% colnames(values$seurat_processed@meta.data)) {
          current_choices <- c(current_choices, "Cell Types" = "cell_type")
        }
        updateSelectInput(
          session,
          "group_by_viz",
          choices = current_choices,
          selected = "seurat_clusters"
        )
        
      }, error = function(e) {
        showNotification(
          paste("Error in clustering:", e$message),
          type = "error",
          duration = NULL
        )
      })
    })
  })
  
  # Find marker genes
  observeEvent(input$find_markers_btn, {
    req(values$seurat_processed)
    withProgress(message = 'Finding marker genes...', value = 0, {
      tryCatch({
        if (input$marker_type == "all") {
          # Find markers for all clusters
          values$marker_genes <- AllMarkers(
            values$seurat_processed,
            logfc.threshold = as.numeric(input$logfc_threshold),
            only.pos = as.logical(input$only_pos)
          )
        } else {
          # Find markers for specific clusters
          values$marker_genes <- Spcfc_Markers(
            values$seurat_processed,
            logfc.threshold = as.numeric(input$logfc_threshold),
            only.pos = as.logical(input$only_pos),
            cluster1 = input$test_cluster,
            cluster2 = input$reference_cluster
          )
        }
        
        # Display marker genes table with styling for dark background
        output$markers_table <- DT::renderDataTable({
          req(values$marker_genes)
          DT::datatable(
            values$marker_genes,
            options = list(
              pageLength = 15,
              scrollX = TRUE,
              scrollY = '500px',
              dom = 'Bfrtlip',
              autoWidth = TRUE,
              columnDefs = list(
                list(width = '150px', targets = c(0, 1)),
                list(className = 'dt-center', targets = "_all")
              )
            ),
            class = 'cell-border stripe hover',
            rownames = FALSE,
            selection = 'single',
            style = 'bootstrap4'
          ) %>% 
            DT::formatStyle(
              columns = names(values$marker_genes),
              backgroundColor = '#343a40',
              color = 'white',
              fontSize = '14px'
            ) %>%
            DT::formatRound(
              columns = c('p_val', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2', 'pct_diff'),
              digits = 3
            )
        })
        
        # FIXED: Plot top markers heatmap with proper gene names and cluster labels
        output$marker_heatmap <- renderPlot({
          req(values$marker_genes, values$seurat_processed)
          
          # Get top markers per cluster based on user input
          if (input$marker_type == "all") {
            # Use user-defined number of genes per cluster
            top_markers <- values$marker_genes %>%
              group_by(cluster) %>%
              top_n(input$n_genes_heatmap, wt = avg_log2FC) %>%
              pull(gene)
            
            # Ensure we have actual gene names that exist in the object
            available_genes <- top_markers[top_markers %in% rownames(values$seurat_processed)]
            
            if (length(available_genes) == 0) {
              showNotification("No marker genes found in the dataset", type = "warning")
              return(NULL)
            }
            
            # Create a much simpler and more visible heatmap
            p <- DoHeatmap(
              values$seurat_processed,
              features = available_genes,
              group.by = "seurat_clusters",
              size = 6,        # Larger gene label text
              angle = 0,       # Horizontal gene labels
              hjust = 0.5,     # Center gene labels
              draw.lines = TRUE,
              lines.width = 2,
              group.bar.height = 0.02
            ) + 
              scale_fill_gradient2(
                low = "blue", 
                mid = "white", 
                high = "red", 
                midpoint = 0,
                name = "Expression",
                guide = guide_colorbar(
                  title.position = "top",
                  title.hjust = 0.5,
                  barwidth = 1,
                  barheight = 15,
                  frame.colour = "black",
                  ticks.colour = "black"
                )
              ) +
              theme_minimal() +
              theme(
                # Y-axis (gene names) - make them very visible
                axis.text.y = element_text(
                  size = 12,
                  color = "black",  # Changed to black for better contrast
                  face = "bold"
                ),
                
                # X-axis (remove cell names, keep cluster labels)
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                
                # Legend styling - make it very visible
                legend.text = element_text(size = 14, color = "black"),
                legend.title = element_text(size = 16, color = "black", face = "bold"),
                legend.background = element_rect(fill = "white", color = "black", linewidth = 1),
                
                # Plot background - white for maximum contrast
                plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA),
                
                # Plot title
                plot.title = element_text(
                  size = 18, 
                  color = "black", 
                  face = "bold",
                  hjust = 0.5
                ),
                
                # Panel and grid
                panel.grid = element_blank(),
                axis.title = element_blank()
              ) +
              ggtitle(paste("Top", input$n_genes_heatmap, "markers per cluster"))
            
            return(p)
            
          } else {
            # For specific clusters
            top_markers <- values$marker_genes %>%
              top_n(input$n_genes_heatmap, wt = avg_log2FC) %>%
              pull(gene)
            
            # Ensure we have actual gene names that exist in the object
            available_genes <- top_markers[top_markers %in% rownames(values$seurat_processed)]
            
            if (length(available_genes) == 0) {
              showNotification("No marker genes found in the dataset", type = "warning")
              return(NULL)
            }
            
            # Create a subset with just the two clusters of interest
            Idents(values$seurat_processed) <- "seurat_clusters"
            cells_use <- WhichCells(values$seurat_processed, 
                                    idents = c(input$test_cluster, input$reference_cluster))
            
            # Subset data for just these two clusters
            sub_obj <- subset(values$seurat_processed, cells = cells_use)
            
            # Create heatmap
            p <- DoHeatmap(
              sub_obj,
              features = available_genes,
              group.by = "seurat_clusters",
              size = 6,        # Larger gene label text
              angle = 0,       # Horizontal gene labels
              hjust = 0.5,     # Center gene labels
              draw.lines = TRUE,
              lines.width = 2,
              group.bar.height = 0.02
            ) + 
              scale_fill_gradient2(
                low = "blue", 
                mid = "white", 
                high = "red", 
                midpoint = 0,
                name = "Expression",
                guide = guide_colorbar(
                  title.position = "top",
                  title.hjust = 0.5,
                  barwidth = 1,
                  barheight = 15,
                  frame.colour = "black",
                  ticks.colour = "black"
                )
              ) +
              theme_minimal() +
              theme(
                # Y-axis (gene names) - make them very visible
                axis.text.y = element_text(
                  size = 12,
                  color = "black",  # Changed to black for better contrast
                  face = "bold"
                ),
                
                # X-axis (remove cell names, keep cluster labels)
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                
                # Legend styling - make it very visible
                legend.text = element_text(size = 14, color = "black"),
                legend.title = element_text(size = 16, color = "black", face = "bold"),
                legend.background = element_rect(fill = "white", color = "black", linewidth = 1),
                
                # Plot background - white for maximum contrast
                plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA),
                
                # Plot title
                plot.title = element_text(
                  size = 18, 
                  color = "black", 
                  face = "bold",
                  hjust = 0.5
                ),
                
                # Panel and grid
                panel.grid = element_blank(),
                axis.title = element_blank()
              ) +
              ggtitle(paste("Top", input$n_genes_heatmap, "markers:", input$test_cluster, "vs", input$reference_cluster))
            
            return(p)
          }
        }, height = function() {
          max(800, input$heatmap_height)
        }, width = function() {
          max(1000, input$heatmap_width)
        })
        
      }, error = function(e) {
        showNotification(
          paste("Error finding markers:", e$message),
          type = "error",
          duration = NULL
        )
      })
    })
  })
  
  # Download marker genes
  output$download_markers_btn <- downloadHandler(
    filename = function() {
      "marker_genes.csv"
    },
    content = function(file) {
      req(values$marker_genes)
      write.csv(values$marker_genes, file)
    }
  )
  
  # Dynamic UI for manual annotation
  output$cluster_annotation_ui <- renderUI({
    req(values$seurat_processed)
    
    if (!"seurat_clusters" %in% colnames(values$seurat_processed@meta.data)) {
      return(p("Please run clustering first."))
    }
    
    clusters <- levels(as.factor(values$seurat_processed$seurat_clusters))
    
    # Create text inputs for each cluster
    input_list <- lapply(clusters, function(cluster) {
      textInput(
        inputId = paste0("cluster_", cluster, "_annotation"),
        label = paste("Cluster", cluster, ":"),
        value = paste0("Cell_type_", cluster),
        placeholder = "Enter cell type..."
      )
    })
    
    do.call(tagList, input_list)
  })
  
  # FIXED: Apply manual annotations with proper cell mapping
  observeEvent(input$apply_manual_annotation, {
    req(values$seurat_processed)
    
    tryCatch({
      # Check if clustering has been done
      if (!"seurat_clusters" %in% colnames(values$seurat_processed@meta.data)) {
        showNotification(
          "Please run clustering first before annotating cell types.",
          type = "warning",
          duration = 5
        )
        return()
      }
      
      clusters <- levels(as.factor(values$seurat_processed$seurat_clusters))
      
      # Collect annotations from text inputs
      annotations <- sapply(clusters, function(cluster) {
        input_id <- paste0("cluster_", cluster, "_annotation")
        annotation <- input[[input_id]]
        if (is.null(annotation) || annotation == "") {
          return(paste0("Cell_type_", cluster))  # Default if empty
        }
        return(annotation)
      })
      
      # Create cell type mapping
      # Get the current cluster assignments for each cell
      current_clusters <- as.character(values$seurat_processed$seurat_clusters)
      
      # Map cluster numbers to cell type names
      cell_types <- annotations[current_clusters]
      
      # Add cell types to metadata using proper cell names
      values$seurat_processed@meta.data$cell_type <- cell_types
      
      # Store the annotation mapping
      values$cluster_annotations <- data.frame(
        cluster = clusters,
        cell_type = annotations,
        stringsAsFactors = FALSE
      )
      
      # Update visualization choices to include cell types
      current_choices <- c("Clusters" = "seurat_clusters", "Cell Types" = "cell_type")
      updateSelectInput(
        session,
        "group_by_viz",
        choices = current_choices,
        selected = "cell_type"  # Auto-select the new annotations
      )
      
      # Show success notification
      showNotification(
        "Cell type annotations applied successfully!",
        type = "message",
        duration = 5
      )
      
      # Update annotation plot
      output$annotation_plot <- renderPlot({
        req(values$seurat_processed)
        Seurat::DimPlot(
          values$seurat_processed,
          reduction = 'umap',
          group.by = 'cell_type',
          pt.size = input$pt_size,
          label = TRUE,
          repel = TRUE
        ) + 
          ggtitle("Manual Cell Type Annotations") +
          theme(plot.title = element_text(hjust = 0.5))
      })
      
      # Create summary table
      output$celltype_summary_table <- DT::renderDataTable({
        req(values$seurat_processed)
        
        # Count cells per cell type
        cell_counts <- table(values$seurat_processed$cell_type)
        summary_df <- data.frame(
          Cell_Type = names(cell_counts),
          Cell_Count = as.numeric(cell_counts),
          Percentage = round(as.numeric(cell_counts) / sum(cell_counts) * 100, 2),
          stringsAsFactors = FALSE
        )
        
        DT::datatable(
          summary_df,
          options = list(
            pageLength = 15,
            dom = 't',
            autoWidth = TRUE
          ),
          style = 'bootstrap4',
          rownames = FALSE
        ) %>%
          DT::formatStyle(
            columns = names(summary_df),
            backgroundColor = '#343a40',
            color = 'white',
            fontSize = '14px'
          )
      })
      
    }, error = function(e) {
      showNotification(
        paste("Error applying annotations:", e$message),
        type = "error",
        duration = 10
      )
      print(paste("Annotation error details:", e$message))
    })
  })
  
  # UPDATED: Gene visualization with user-selected grouping
  output$gene_plot <- renderPlot({
    req(values$seurat_processed)
    
    # Check what type of plot to create
    if (input$plot_type == "dimplot") {
      # UMAP plot colored by selected grouping variable
      if (input$group_by_viz %in% colnames(values$seurat_processed@meta.data)) {
        Seurat::DimPlot(
          values$seurat_processed,
          reduction = 'umap',
          group.by = input$group_by_viz,
          pt.size = input$pt_size,
          label = TRUE,
          repel = TRUE
        ) + 
          ggtitle(paste("UMAP colored by", 
                        ifelse(input$group_by_viz == "seurat_clusters", "Clusters", "Cell Types"))) +
          theme(plot.title = element_text(hjust = 0.5))
      } else {
        # Fallback to clusters if selected grouping doesn't exist
        Seurat::DimPlot(
          values$seurat_processed,
          reduction = 'umap',
          group.by = "seurat_clusters",
          pt.size = input$pt_size,
          label = TRUE,
          repel = TRUE
        ) + 
          ggtitle("UMAP colored by Clusters") +
          theme(plot.title = element_text(hjust = 0.5))
      }
    } else if (length(input$gene_select) > 0) {
      # Gene expression plots
      if (input$plot_type == "feature") {
        # Feature plot
        Seurat::FeaturePlot(
          values$seurat_processed,
          reduction = 'umap',
          features = input$gene_select,
          pt.size = input$pt_size
        )
      } else if (input$plot_type == "violin") {
        # Violin plot with selected grouping
        group_var <- input$group_by_viz
        if (!group_var %in% colnames(values$seurat_processed@meta.data)) {
          group_var <- "seurat_clusters"  # Fallback
        }
        
        if (length(input$gene_select) <= 3) {
          Seurat::VlnPlot(
            values$seurat_processed,
            features = input$gene_select,
            group.by = group_var,
            pt.size = 0
          )
        } else {
          # For more than 3 genes, use stacked violin plots
          Seurat::VlnPlot(
            values$seurat_processed,
            features = input$gene_select,
            group.by = group_var,
            stack = TRUE,
            flip = TRUE
          )
        }
      }
    } else {
      # Default plot when no genes selected
      if (input$group_by_viz %in% colnames(values$seurat_processed@meta.data)) {
        Seurat::DimPlot(
          values$seurat_processed,
          reduction = 'umap',
          group.by = input$group_by_viz,
          pt.size = input$pt_size,
          label = TRUE,
          repel = TRUE
        ) + 
          ggtitle(paste("UMAP colored by", 
                        ifelse(input$group_by_viz == "seurat_clusters", "Clusters", "Cell Types"))) +
          theme(plot.title = element_text(hjust = 0.5))
      } else {
        # Fallback plot
        Seurat::DimPlot(
          values$seurat_processed,
          reduction = 'umap',
          group.by = "seurat_clusters",
          pt.size = input$pt_size,
          label = TRUE,
          repel = TRUE
        ) + 
          ggtitle("UMAP colored by Clusters") +
          theme(plot.title = element_text(hjust = 0.5))
      }
    }
  }, height = function() {
    # Dynamic height based on input
    height_val <- as.numeric(gsub("px", "", input$viz_plot_height))
    if (is.na(height_val) || height_val < 400) height_val <- 600
    return(height_val)
  }, width = function() {
    # Dynamic width based on input
    width_val <- as.numeric(gsub("px", "", input$viz_plot_width))
    if (is.na(width_val) || width_val < 400) width_val <- 800
    return(width_val)
  })
  
  # Download annotated Seurat object
  output$download_annotated_obj_btn <- downloadHandler(
    filename = function() {
      "annotated_seurat_object.rds"
    },
    content = function(file) {
      req(values$seurat_processed)
      saveRDS(values$seurat_processed, file)
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
