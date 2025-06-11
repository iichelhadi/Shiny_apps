# R/visualization.R - Functions for visualization

#' Generate QC plots for a Seurat object
#'
#' @param seurat_obj Seurat object with QC metrics
#' @param group_by Optional variable to group cells by
#' @param ncol Number of columns for multi-panel plots
#' @return A patchwork object with QC plots
plot_qc <- function(seurat_obj, group_by = NULL, ncol = 2) {
  # Determine QC features to plot
  qc_features <- c(
    "nFeature_RNA", 
    "nCount_RNA", 
    "percent.mt"
  )
  
  # Add optional features if they exist
  if ("percent.ribo" %in% colnames(seurat_obj@meta.data)) {
    qc_features <- c(qc_features, "percent.ribo")
  }
  if ("percent.hsp" %in% colnames(seurat_obj@meta.data)) {
    qc_features <- c(qc_features, "percent.hsp")
  }
  if ("log10GenesPerUMI" %in% colnames(seurat_obj@meta.data)) {
    qc_features <- c(qc_features, "log10GenesPerUMI")
  }
  
  # Create violin plots
  if (!is.null(group_by) && group_by %in% colnames(seurat_obj@meta.data)) {
    vln <- VlnPlot(
      seurat_obj, 
      features = qc_features, 
      group.by = group_by,
      pt.size = 0,
      ncol = ncol
    )
  } else {
    vln <- VlnPlot(
      seurat_obj, 
      features = qc_features,
      pt.size = 0,
      ncol = ncol
    )
  }
  
  # Create scatter plots for paired metrics
  scatter1 <- FeatureScatter(
    seurat_obj, 
    feature1 = "nCount_RNA", 
    feature2 = "nFeature_RNA"
  )
  
  scatter2 <- FeatureScatter(
    seurat_obj, 
    feature1 = "nCount_RNA", 
    feature2 = "percent.mt"
  )
  
  # Combine plots
  combined <- vln / (scatter1 | scatter2)
  
  return(combined)
}

#' Generate a dimensionality reduction plot for a Seurat object
#'
#' @param seurat_obj Seurat object with dimensionality reduction computed
#' @param reduction Dimensionality reduction to plot (e.g., "umap", "tsne", "pca")
#' @param group_by Variable to group cells by (e.g., "seurat_clusters", "cell_type")
#' @param split_by Variable to split the plot by
#' @param pt_size Point size
#' @param label Whether to add labels to clusters
#' @param label_size Size of labels
#' @param repel Whether to repel labels
#' @param colors Custom color palette
#' @return A ggplot object
plot_dimred <- function(seurat_obj,
                        reduction = "umap",
                        group_by = "seurat_clusters",
                        split_by = NULL,
                        pt_size = 0.5,
                        label = TRUE,
                        label_size = 5,
                        repel = TRUE,
                        colors = NULL) {
  
  # Check if reduction exists
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in Seurat object"))
  }
  
  # Check if group_by exists
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Group variable", group_by, "not found in metadata"))
  }
  
  # Check if split_by exists (if provided)
  if (!is.null(split_by) && !split_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Split variable", split_by, "not found in metadata"))
  }
  
  # Use scCustomize for more flexibility
  if (is.null(colors)) {
    # Use default custom palette
    plot <- DimPlot_scCustom(
      seurat_obj,
      reduction = reduction,
      group.by = group_by,
      split.by = split_by,
      pt.size = pt_size,
      label = label,
      label.size = label_size,
      repel = repel
    )
  } else {
    # Use custom colors
    plot <- DimPlot_scCustom(
      seurat_obj,
      reduction = reduction,
      group.by = group_by,
      split.by = split_by,
      pt.size = pt_size,
      label = label,
      label.size = label_size,
      repel = repel,
      colors_use = colors
    )
  }
  
  return(plot)
}

#' Generate a feature plot for a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param features Features to plot
#' @param reduction Dimensionality reduction to use
#' @param slot Slot to pull data from
#' @param pt_size Point size
#' @param order Whether to order points by expression
#' @param min_cutoff Minimum cutoff for color scale
#' @param max_cutoff Maximum cutoff for color scale
#' @param ncol Number of columns
#' @param colors Custom color palette
#' @return A ggplot object
plot_feature <- function(seurat_obj,
                        features,
                        reduction = "umap",
                        slot = "data",
                        pt_size = 0.5,
                        order = TRUE,
                        min_cutoff = "q1",
                        max_cutoff = "q99",
                        ncol = 2,
                        colors = NULL) {
  
  # Check if reduction exists
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in Seurat object"))
  }
  
  # Use scCustomize for more flexibility
  if (is.null(colors)) {
    # Use default custom palette
    plot <- FeaturePlot_scCustom(
      seurat_obj,
      features = features,
      reduction = reduction,
      slot = slot,
      pt.size = pt_size,
      order = order,
      min.cutoff = min_cutoff,
      max.cutoff = max_cutoff,
      ncol = ncol
    )
  } else {
    # Use custom colors
    plot <- FeaturePlot_scCustom(
      seurat_obj,
      features = features,
      reduction = reduction,
      slot = slot,
      pt.size = pt_size,
      order = order,
      min.cutoff = min_cutoff,
      max.cutoff = max_cutoff,
      ncol = ncol,
      colors_use = colors
    )
  }
  
  return(plot)
}

#' Generate a violin plot for a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param features Features to plot
#' @param group_by Variable to group cells by
#' @param split_by Variable to split the plot by
#' @param pt_size Point size (0 for no points)
#' @param ncol Number of columns
#' @param colors Custom color palette
#' @param stack Whether to stack violin plots
#' @return A ggplot object
plot_violin <- function(seurat_obj,
                       features,
                       group_by = "seurat_clusters",
                       split_by = NULL,
                       pt_size = 0,
                       ncol = 2,
                       colors = NULL,
                       stack = FALSE) {
  
  # Check if group_by exists
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Group variable", group_by, "not found in metadata"))
  }
  
  # Check if split_by exists (if provided)
  if (!is.null(split_by) && !split_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Split variable", split_by, "not found in metadata"))
  }
  
  if (stack) {
    # Use stacked violin plots from scCustomize
    plot <- Stacked_VlnPlot(
      seurat_obj,
      features = features,
      group.by = group_by,
      split.by = split_by,
      colors_use = colors
    )
  } else {
    # Use regular violin plots
    if (is.null(colors)) {
      # Use default custom palette
      plot <- VlnPlot_scCustom(
        seurat_obj,
        features = features,
        group.by = group_by,
        split.by = split_by,
        pt.size = pt_size,
        ncol = ncol
      )
    } else {
      # Use custom colors
      plot <- VlnPlot_scCustom(
        seurat_obj,
        features = features,
        group.by = group_by,
        split.by = split_by,
        pt.size = pt_size,
        ncol = ncol,
        colors_use = colors
      )
    }
  }
  
  return(plot)
}

#' Generate a dot plot for a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param features Features to plot
#' @param group_by Variable to group cells by
#' @param split_by Variable to split the plot by
#' @param scale Whether to scale expression values
#' @param colors Custom color palette
#' @return A ggplot object
plot_dot <- function(seurat_obj,
                    features,
                    group_by = "seurat_clusters",
                    split_by = NULL,
                    scale = TRUE,
                    colors = NULL) {
  
  # Check if group_by exists
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Group variable", group_by, "not found in metadata"))
  }
  
  # Check if split_by exists (if provided)
  if (!is.null(split_by) && !split_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Split variable", split_by, "not found in metadata"))
  }
  
  # Use scCustomize for more flexibility
  if (is.null(colors)) {
    # Use default custom palette
    plot <- DotPlot_scCustom(
      seurat_obj,
      features = features,
      group.by = group_by,
      split.by = split_by,
      scale = scale
    )
  } else {
    # Use custom colors
    plot <- DotPlot_scCustom(
      seurat_obj,
      features = features,
      group.by = group_by,
      split.by = split_by,
      scale = scale,
      colors_use = colors
    )
  }
  
  return(plot)
}

#' Generate a heatmap for a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param features Features to plot
#' @param group_by Variable to group cells by
#' @param assay Assay to use
#' @param slot Slot to pull data from
#' @param scale Whether to scale data
#' @param colors Custom color palette
#' @return A ggplot object
plot_heatmap <- function(seurat_obj,
                        features,
                        group_by = "seurat_clusters",
                        assay = NULL,
                        slot = "scale.data",
                        scale = TRUE,
                        colors = NULL) {
  
  # Determine which assay to use
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  
  # Check if group_by exists
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Group variable", group_by, "not found in metadata"))
  }
  
  # Use DoHeatmap from Seurat with custom parameters
  if (is.null(colors)) {
    plot <- DoHeatmap(
      seurat_obj,
      features = features,
      group.by = group_by,
      assay = assay,
      slot = slot,
      scale = scale
    )
  } else {
    plot <- DoHeatmap(
      seurat_obj,
      features = features,
      group.by = group_by,
      assay = assay,
      slot = slot,
      scale = scale,
      cols = colors
    )
  }
  
  return(plot)
}

#' Generate a ridge plot for a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param features Features to plot
#' @param group_by Variable to group cells by
#' @param ncol Number of columns
#' @param colors Custom color palette
#' @return A ggplot object
plot_ridge <- function(seurat_obj,
                      features,
                      group_by = "seurat_clusters",
                      ncol = 2,
                      colors = NULL) {
  
  # Check if group_by exists
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Group variable", group_by, "not found in metadata"))
  }
  
  # Use RidgePlot from Seurat with custom parameters
  if (is.null(colors)) {
    plot <- RidgePlot(
      seurat_obj,
      features = features,
      group.by = group_by,
      ncol = ncol
    )
  } else {
    plot <- RidgePlot(
      seurat_obj,
      features = features,
      group.by = group_by,
      ncol = ncol,
      cols = colors
    )
  }
  
  return(plot)
}

#' Generate an interactive dimensionality reduction plot with plotly
#'
#' @param seurat_obj Seurat object with dimensionality reduction computed
#' @param reduction Dimensionality reduction to plot
#' @param group_by Variable to group cells by
#' @param label_text Variables to show in tooltips
#' @param pt_size Point size
#' @param colors Custom color palette
#' @return A plotly object
plot_interactive_dimred <- function(seurat_obj,
                                  reduction = "umap",
                                  group_by = "seurat_clusters",
                                  label_text = c("orig.ident", "nFeature_RNA", "nCount_RNA", "percent.mt"),
                                  pt_size = 3,
                                  colors = NULL) {
  
  # Check if reduction exists
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in Seurat object"))
  }
  
  # Check if group_by exists
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Group variable", group_by, "not found in metadata"))
  }
  
  # Extract reduction coordinates
  reduction_coords <- Embeddings(seurat_obj, reduction = reduction)
  
  # Prepare data frame for plotting
  plot_data <- data.frame(
    x = reduction_coords[, 1],
    y = reduction_coords[, 2],
    group = seurat_obj[[group_by]]
  )
  
  # Add additional metadata for tooltips
  for (col in label_text) {
    if (col %in% colnames(seurat_obj@meta.data)) {
      plot_data[[col]] <- seurat_obj[[col]]
    }
  }
  
  # Create tooltip text
  tooltip_text <- apply(plot_data[, c("group", label_text[label_text %in% colnames(seurat_obj@meta.data)])], 
                        1, function(x) {
                          paste0(names(x), ": ", x, collapse = "<br>")
                        })
  
  # Create color mapping
  if (is.null(colors)) {
    colors <- colorRampPalette(c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725"))(length(unique(plot_data$group)))
  }
  
  color_mapping <- setNames(colors[1:length(unique(plot_data$group))], sort(unique(plot_data$group)))
  
  # Create plotly plot
  p <- plotly::plot_ly(
    data = plot_data,
    x = ~x,
    y = ~y,
    color = ~group,
    colors = color_mapping,
    type = "scatter",
    mode = "markers",
    marker = list(size = pt_size),
    text = tooltip_text,
    hoverinfo = "text"
  )
  
  # Set layout
  p <- plotly::layout(
    p,
    title = paste0(reduction, " plot colored by ", group_by),
    xaxis = list(title = paste0(reduction, "_1")),
    yaxis = list(title = paste0(reduction, "_2")),
    legend = list(title = list(text = group_by))
  )
  
  return(p)
}

#' Generate an interactive feature plot with plotly
#'
#' @param seurat_obj Seurat object
#' @param feature Feature to plot
#' @param reduction Dimensionality reduction to use
#' @param assay Assay to use
#' @param slot Slot to pull data from
#' @param pt_size Point size
#' @param label_text Variables to show in tooltips
#' @param colors Custom color palette
#' @return A plotly object
plot_interactive_feature <- function(seurat_obj,
                                   feature,
                                   reduction = "umap",
                                   assay = NULL,
                                   slot = "data",
                                   pt_size = 3,
                                   label_text = c("orig.ident", "nFeature_RNA", "nCount_RNA"),
                                   colors = NULL) {
  
  # Check if reduction exists
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in Seurat object"))
  }
  
  # Determine which assay to use
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  
  # Check if feature exists
  if (!feature %in% rownames(seurat_obj)) {
    stop(paste("Feature", feature, "not found in Seurat object"))
  }
  
  # Extract reduction coordinates
  reduction_coords <- Embeddings(seurat_obj, reduction = reduction)
  
  # Get feature expression
  feature_expr <- FetchData(seurat_obj, vars = feature, slot = slot, assay = assay)
  
  # Prepare data frame for plotting
  plot_data <- data.frame(
    x = reduction_coords[, 1],
    y = reduction_coords[, 2],
    expr = feature_expr[, 1]
  )
  
  # Add additional metadata for tooltips
  for (col in label_text) {
    if (col %in% colnames(seurat_obj@meta.data)) {
      plot_data[[col]] <- seurat_obj[[col]]
    }
  }
  
  # Add expression to tooltip
  plot_data[["feature"]] <- feature
  
  # Create tooltip text
  tooltip_text <- apply(plot_data[, c("feature", "expr", label_text[label_text %in% colnames(seurat_obj@meta.data)])], 
                        1, function(x) {
                          paste0(names(x), ": ", round(as.numeric(x), 3), collapse = "<br>")
                        })
  
  # Create color mapping
  if (is.null(colors)) {
    colors <- colorRampPalette(c("lightgrey", "#FDE725FF"))(100)
  }
  
  # Create plotly plot
  p <- plotly::plot_ly(
    data = plot_data,
    x = ~x,
    y = ~y,
    color = ~expr,
    colors = colors,
    type = "scatter",
    mode = "markers",
    marker = list(size = pt_size),
    text = tooltip_text,
    hoverinfo = "text"
  )
  
  # Set layout
  p <- plotly::layout(
    p,
    title = paste0(feature, " expression"),
    xaxis = list(title = paste0(reduction, "_1")),
    yaxis = list(title = paste0(reduction, "_2")),
    colorbar = list(title = feature)
  )
  
  return(p)
}

#' Plot elbow plot to determine optimal number of PCs
#'
#' @param seurat_obj Seurat object with PCA computed
#' @param ndims Number of dimensions to show
#' @return A ggplot object
plot_elbow <- function(seurat_obj, ndims = 50) {
  ElbowPlot(seurat_obj, ndims = ndims)
}

#' Plot PC loadings for genes contributing to principal components
#'
#' @param seurat_obj Seurat object with PCA computed
#' @param dims Principal components to plot
#' @param nfeatures Number of features to show
#' @return A ggplot object
plot_pc_loadings <- function(seurat_obj, dims = 1:2, nfeatures = 10) {
  VizDimLoadings(seurat_obj, dims = dims, nfeatures = nfeatures)
}

#' Plot heatmap of PC loadings
#'
#' @param seurat_obj Seurat object with PCA computed
#' @param dims Principal components to plot
#' @param balanced Whether to show both positive and negative loadings
#' @param ncol Number of columns
#' @return A ggplot object
plot_pc_heatmap <- function(seurat_obj, dims = 1:9, balanced = TRUE, ncol = 3) {
  DimHeatmap(seurat_obj, dims = dims, balanced = balanced, ncol = ncol)
}

#' Plot expression of multiple genes as a bubble chart grouped by cluster identity
#'
#' @param seurat_obj Seurat object
#' @param features Features to plot
#' @param group_by Variable to group cells by
#' @param dot_scale Scale for dot size
#' @param cluster_idents Whether to cluster identities
#' @param scale Whether to scale expression values
#' @return A ggplot object
plot_marker_bubble <- function(seurat_obj, 
                             features, 
                             group_by = "seurat_clusters",
                             dot_scale = 8, 
                             cluster_idents = FALSE,
                             scale = TRUE) {
  
  DotPlot(seurat_obj, 
          features = features, 
          group.by = group_by,
          dot.scale = dot_scale, 
          cluster.idents = cluster_idents,
          scale = scale) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
