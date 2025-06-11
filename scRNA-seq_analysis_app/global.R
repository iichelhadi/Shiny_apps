# global.R - Load packages and global functions for scRNA-seq app

# Suppress startup messages for cleaner console
suppressPackageStartupMessages({
  # Core Shiny packages
  library(shiny)
  library(shinythemes)
  library(shinyFiles)
  library(DT)
  
  # Bioconductor packages
  library(scDblFinder)
  library(BiocParallel)
  library(MAST)
  library(SingleR)
  
  # Visualization packages
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  
  # Analysis packages
  library(Seurat)
  library(tidyverse)
  library(reticulate)
  library(scCustomize)
  
  # For cell type reference datasets
  library(celldex)
  
  # Additional utilities
  library(parallelly)
})

# Check if leidenalg is available and import it
if(reticulate::py_module_available(module = 'leidenalg')) {
  reticulate::import('leidenalg')
} else {
  warning("Python module 'leidenalg' not available. Some clustering methods will be disabled.")
}

# Set random seed for reproducibility
set.seed(123)

# Global options
options(shiny.maxRequestSize = 1000 * 1024^2) # Increase file size limit to 1GB

############################ functions #########################################


# Detect doublets in a Seurat object using scDblFinder
doublet_detection <- function(obj, col) {
  seurat_obj <- obj
  print('##################### convert to sce ##########################')
  sce <- as.SingleCellExperiment(seurat_obj)
  
  print('##################### scDblFinder #############################')
  sce <- scDblFinder(
    sce,
    samples = col,
    BPPARAM = MulticoreParam(as.integer(parallelly::availableCores(methods = 'nproc')) - 1),
    clusters = TRUE
  )
  
  # Add doublet classification and score to metadata
  seurat_obj@meta.data <- cbind(
    seurat_obj@meta.data, 
    scDblFinder.class = sce@colData[, 'scDblFinder.class'],
    scDblFinder.score = sce@colData[, 'scDblFinder.score']
  )
  
  # Print summary stats
  n_cells <- ncol(seurat_obj)
  n_doublets <- sum(seurat_obj@meta.data$scDblFinder.class == "doublet")
  n_singlets <- sum(seurat_obj@meta.data$scDblFinder.class == "singlet")
  doublet_rate <- round(n_doublets / n_cells * 100, 2)
  
  print(paste0("Detected ", n_doublets, " doublets (", doublet_rate, "%) out of ", n_cells, " cells"))
  
  return(seurat_obj)
}

# Run dimensionality reduction on a Seurat object
dimred <- function(obj, nfeatures) {
  print('##################### normalize data ##########################')
  seurat_obj <- NormalizeData(obj)
  
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = nfeatures
  )
  
  all.genes <- rownames(seurat_obj)
  
  seurat_obj <- ScaleData(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt')
  )
  
  seurat_obj <- RunPCA(
    seurat_obj, 
    features = VariableFeatures(object = seurat_obj)
  )
  
  # Determine number of PCs to use
  pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  pcs <- which(cumu > 70 & pct < 5)[1]
  print('###### pcs ########')
  print(pcs)
  
  # Run UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:pcs)
  
  return(seurat_obj)
}

# Run clustering on a Seurat object
clustering <- function(obj, resolution, algorithm) {
  # Determine number of PCs to use
  pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  pcs <- which(cumu > 70 & pct < 5)[1]
  print('###### pcs ########')
  print(pcs)
  
  # Find neighbors and clusters
  seurat_obj <- FindNeighbors(obj, dims = 1:pcs)
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = as.numeric(resolution),
    algorithm = as.integer(algorithm)
  )
  
  return(seurat_obj)
}

# Find marker genes for all clusters
AllMarkers <- function(obj, logfc.threshold = 0.25, only.pos = TRUE) {
  de_markers <- FindAllMarkers(
    obj,
    test.use = 'MAST',
    verbose = TRUE,
    logfc.threshold = logfc.threshold,
    only.pos = only.pos
  ) %>%
    scCustomize::Add_Pct_Diff() %>% 
    arrange(desc(pct_diff))
  
  return(de_markers)
}

# Find marker genes between specific clusters
Spcfc_Markers <- function(obj, logfc.threshold = 0.25, only.pos = TRUE, cluster1, cluster2) {
  de_markers <- FindMarkers(
    obj,
    ident.1 = cluster1,
    ident.2 = cluster2,
    test.use = 'MAST',
    verbose = TRUE,
    logfc.threshold = logfc.threshold,
    only.pos = only.pos
  ) %>%
    rownames_to_column("gene") %>%
    scCustomize::Add_Pct_Diff() %>% 
    arrange(desc(pct_diff))
  
  return(de_markers)
}

# Generate a feature plot
FeaturePlot_scCustom <- function(obj, features, reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = "q1", max.cutoff = "q99", ncol = NULL, colors_use = NULL) {
  if(is.null(colors_use)) {
    # Default viridis-based color palette
    colors_use <- viridisLite::viridis(n = 100, direction = 1)
  }
  
  FeaturePlot(
    obj,
    features = features,
    reduction = reduction,
    pt.size = pt.size,
    order = order,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff,
    ncol = ncol,
    cols = colors_use
  )
}

# Generate a dimension reduction plot
DimPlot_scCustom <- function(obj, reduction = "umap", group.by = "seurat_clusters", pt.size = 1, label = TRUE, label.size = 5, repel = TRUE, colors_use = NULL) {
  if(!is.null(colors_use) && colors_use == "black") {
    DimPlot(
      obj,
      reduction = reduction,
      group.by = group.by,
      pt.size = pt.size,
      label = label,
      label.size = label.size,
      repel = repel,
      cols = "black"
    )
  } else if(!is.null(colors_use)) {
    DimPlot(
      obj,
      reduction = reduction,
      group.by = group.by,
      pt.size = pt.size,
      label = label,
      label.size = label.size,
      repel = repel,
      cols = colors_use
    )
  } else {
    # Default color palette
    DimPlot(
      obj,
      reduction = reduction,
      group.by = group.by,
      pt.size = pt.size,
      label = label,
      label.size = label.size,
      repel = repel
    )
  }
}

# Generate stacked violin plots
Stacked_VlnPlot <- function(obj, features, group.by = "seurat_clusters", pt.size = 0, colors_use = NULL) {
  VlnPlot(
    obj,
    features = features,
    group.by = group.by,
    pt.size = pt.size,
    stack = TRUE,
    flip = TRUE,
    cols = colors_use
  )
}

# Generate a violin plot with custom colors
VlnPlot_scCustom <- function(obj, features, group.by = "seurat_clusters", split.by = NULL, pt.size = 0, ncol = NULL, colors_use = NULL) {
  VlnPlot(
    obj,
    features = features,
    group.by = group.by,
    split.by = split.by,
    pt.size = pt.size,
    ncol = ncol,
    cols = colors_use
  )
}
