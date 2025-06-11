# R/preprocessing.R - Functions for data preprocessing

#' Detect doublets in a Seurat object using scDblFinder
#'
#' @param seurat_obj Seurat object
#' @param sample_col Column name in metadata that identifies samples
#' @param BPPARAM BiocParallel parameters for parallel processing
#' @return Seurat object with doublet information added
detect_doublets <- function(seurat_obj, sample_col = "orig.ident", BPPARAM = NULL) {
  # Create default BPPARAM if not provided
  if (is.null(BPPARAM)) {
    BPPARAM <- create_bioc_processor(prop = 0.75)
  }
  
  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_obj)
  
  # Run scDblFinder
  sce <- scDblFinder(
    sce,
    samples = sample_col,
    BPPARAM = BPPARAM,
    clusters = TRUE
  )
  
  # Add doublet information back to Seurat object
  seurat_obj@meta.data$scDblFinder.class <- sce@colData[, "scDblFinder.class"]
  seurat_obj@meta.data$scDblFinder.score <- sce@colData[, "scDblFinder.score"]
  
  return(seurat_obj)
}

#' Detect doublets using DoubletFinder
#'
#' @param seurat_obj Seurat object with PCA computed
#' @param pcs Number of PCs to use
#' @param sct Whether SCTransform was used
#' @return Seurat object with doublet information added
detect_doublets_df <- function(seurat_obj, pcs = NULL, sct = FALSE) {
  if (is.null(pcs)) {
    # If pcs not provided, calculate based on elbow point
    pcs <- get_suggested_pcs(seurat_obj)
  }
  
  # Get optimal pK parameter
  # This is the "proportion of artificial nearest neighbors"
  sweep_res <- paramSweep_v3(seurat_obj, PCs = 1:pcs, sct = sct)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  
  # Get the pK value at the maximum BCmetric
  optimal_pk <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) %>%
    pull(pK) %>%
    as.character() %>%
    as.numeric()
  
  # Compute doublet rate based on number of cells
  n_cells <- ncol(seurat_obj)
  expected_doublets <- 0.008 * n_cells / 1000  # Approx. doublet rate for 10X data
  
  # Run DoubletFinder
  seurat_obj <- doubletFinder_v3(
    seurat_obj,
    PCs = 1:pcs,
    pN = 0.25,  # Proportion of artificial doublets
    pK = optimal_pk,
    nExp = round(expected_doublets * n_cells),
    reuse.pANN = FALSE,
    sct = sct
  )
  
  # Standardize column names for compatibility with other functions
  # Find the DoubletFinder classification column (name varies with parameters)
  df_cols <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(df_cols) > 0) {
    # Create standardized column
    seurat_obj@meta.data$DoubletFinder.class <- seurat_obj@meta.data[[df_cols[1]]]
    # Rename to match scDblFinder
    seurat_obj@meta.data$DoubletFinder.class <- plyr::mapvalues(
      seurat_obj@meta.data$DoubletFinder.class,
      from = c("Singlet", "Doublet"),
      to = c("singlet", "doublet")
    )
    
    # Also standardize the score column
    df_score_cols <- grep("pANN", colnames(seurat_obj@meta.data), value = TRUE)
    if (length(df_score_cols) > 0) {
      seurat_obj@meta.data$DoubletFinder.score <- seurat_obj@meta.data[[df_score_cols[1]]]
    }
  }
  
  return(seurat_obj)
}

#' Filter doublets from a Seurat object
#'
#' @param seurat_obj Seurat object with doublet information
#' @param doublet_col Column name with doublet classification
#' @return Filtered Seurat object with only singlets
filter_doublets <- function(seurat_obj, doublet_col = "scDblFinder.class") {
  if (!doublet_col %in% colnames(seurat_obj@meta.data)) {
    stop("Doublet column not found in metadata. Run detect_doublets first.")
  }
  
  # Keep only singlets
  seurat_filtered <- subset(seurat_obj, subset = get(doublet_col) == "singlet")
  
  return(seurat_filtered)
}

#' Normalize and scale data in a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param normalization_method Normalization method ("LogNormalize" or "SCT")
#' @param scale_factor Scale factor for normalization
#' @param vars_to_regress Variables to regress out during scaling
#' @param verbose Print progress messages
#' @return Normalized and scaled Seurat object
normalize_and_scale <- function(seurat_obj, 
                               normalization_method = "LogNormalize", 
                               scale_factor = 10000, 
                               vars_to_regress = c("nCount_RNA", "percent.mt"),
                               verbose = TRUE) {
  
  if (normalization_method == "LogNormalize") {
    # Standard log-normalization
    seurat_obj <- NormalizeData(
      seurat_obj,
      normalization.method = normalization_method,
      scale.factor = scale_factor,
      verbose = verbose
    )
    
    # Find variable features
    seurat_obj <- FindVariableFeatures(
      seurat_obj,
      selection.method = "vst",
      nfeatures = 2000,
      verbose = verbose
    )
    
    # Scale data
    seurat_obj <- ScaleData(
      seurat_obj,
      features = VariableFeatures(seurat_obj),
      vars.to.regress = vars_to_regress,
      verbose = verbose
    )
    
  } else if (normalization_method == "SCT") {
    # SCTransform normalization (more robust to technical variation)
    seurat_obj <- SCTransform(
      seurat_obj,
      vars.to.regress = vars_to_regress,
      verbose = verbose
    )
  } else {
    stop("Unsupported normalization method. Use 'LogNormalize' or 'SCT'.")
  }
  
  return(seurat_obj)
}

#' Run dimensionality reduction on a Seurat object
#'
#' @param seurat_obj Normalized and scaled Seurat object
#' @param reduction_method Dimensionality reduction method ("pca", "umap", "tsne", or "all")
#' @param dims_to_use Number of dimensions to use for UMAP and t-SNE
#' @param assay Assay to use for dimensionality reduction
#' @param seed Random seed for reproducibility
#' @param verbose Print progress messages
#' @return Seurat object with dimensionality reduction results
run_dim_reduction <- function(seurat_obj, 
                             reduction_method = "all", 
                             dims_to_use = NULL,
                             assay = NULL,
                             seed = 42,
                             verbose = TRUE) {
  
  # Determine which assay to use
  if (is.null(assay)) {
    if ("SCT" %in% names(seurat_obj@assays)) {
      assay <- "SCT"
    } else {
      assay <- DefaultAssay(seurat_obj)
    }
  }
  
  # Set default assay
  DefaultAssay(seurat_obj) <- assay
  
  # Run PCA if requested or needed for other reductions
  if (reduction_method %in% c("pca", "all") || reduction_method %in% c("umap", "tsne")) {
    if (verbose) message("Running PCA...")
    seurat_obj <- RunPCA(
      seurat_obj,
      features = VariableFeatures(object = seurat_obj),
      npcs = 50,  # Calculate more PCs than we'll use for flexibility
      verbose = verbose
    )
  }
  
  # Determine number of dimensions to use
  if (is.null(dims_to_use)) {
    dims_to_use <- get_suggested_pcs(seurat_obj)
    if (verbose) message(paste("Using", dims_to_use, "PCs based on elbow method"))
  }
  
  # Run UMAP if requested
  if (reduction_method %in% c("umap", "all")) {
    if (verbose) message("Running UMAP...")
    set.seed(seed)
    seurat_obj <- RunUMAP(
      seurat_obj,
      dims = 1:dims_to_use,
      verbose = verbose
    )
  }
  
  # Run t-SNE if requested
  if (reduction_method %in% c("tsne", "all")) {
    if (verbose) message("Running t-SNE...")
    set.seed(seed)
    seurat_obj <- RunTSNE(
      seurat_obj,
      dims = 1:dims_to_use,
      verbose = verbose
    )
  }
  
  return(seurat_obj)
}

#' Run clustering on a Seurat object
#'
#' @param seurat_obj Seurat object with dimensionality reduction
#' @param resolution Resolution parameter for clustering
#' @param dims_to_use Number of dimensions to use
#' @param algorithm Algorithm to use (1=Louvain, 2=Louvain with multilevel refinement, 3=SLM, 4=Leiden)
#' @param verbose Print progress messages
#' @return Seurat object with clustering results
run_clustering <- function(seurat_obj, 
                          resolution = 0.8, 
                          dims_to_use = NULL,
                          algorithm = 4,
                          verbose = TRUE) {
  
  # Determine number of dimensions to use
  if (is.null(dims_to_use)) {
    dims_to_use <- get_suggested_pcs(seurat_obj)
    if (verbose) message(paste("Using", dims_to_use, "PCs based on elbow method"))
  }
  
  # Find nearest neighbors
  seurat_obj <- FindNeighbors(
    seurat_obj,
    dims = 1:dims_to_use,
    verbose = verbose
  )
  
  # Find clusters
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = resolution,
    algorithm = algorithm,
    verbose = verbose
  )
  
  return(seurat_obj)
}

#' Run batch correction using Harmony
#'
#' @param seurat_obj Seurat object with PCA computed
#' @param batch_var Batch variable in metadata to correct for
#' @param dims_to_use Number of dimensions to use
#' @param verbose Print progress messages
#' @return Seurat object with batch-corrected reduction
run_harmony <- function(seurat_obj, 
                        batch_var, 
                        dims_to_use = NULL,
                        verbose = TRUE) {
  
  # Check if batch variable exists
  if (!batch_var %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Batch variable", batch_var, "not found in metadata"))
  }
  
  # Determine number of dimensions to use
  if (is.null(dims_to_use)) {
    dims_to_use <- get_suggested_pcs(seurat_obj)
    if (verbose) message(paste("Using", dims_to_use, "PCs based on elbow method"))
  }
  
  # Run Harmony
  seurat_obj <- harmony::RunHarmony(
    seurat_obj,
    group.by.vars = batch_var,
    dims.use = 1:dims_to_use,
    verbose = verbose
  )
  
  # Run UMAP on harmony embeddings
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "harmony",
    dims = 1:dims_to_use,
    verbose = verbose
  )
  
  return(seurat_obj)
}

#' Annotate cell types using SingleR
#'
#' @param seurat_obj Seurat object
#' @param ref_dataset Reference dataset name
#' @param labels_to_use Labels to use from reference (main or fine)
#' @param assay Assay to use for annotation
#' @param clusters Whether to predict labels for clusters instead of cells
#' @param cluster_col Column name with cluster IDs if clusters=TRUE
#' @param BPPARAM BiocParallel parameters for parallel processing
#' @return Seurat object with cell type annotations
annotate_cell_types <- function(seurat_obj,
                               ref_dataset = "BlueprintEncodeData",
                               labels_to_use = "main",
                               assay = "RNA",
                               clusters = FALSE,
                               cluster_col = "seurat_clusters",
                               BPPARAM = NULL) {
  
  # Create default BPPARAM if not provided
  if (is.null(BPPARAM)) {
    BPPARAM <- create_bioc_processor(prop = 0.75)
  }
  
  # Load reference dataset
  ref_data <- switch(
    ref_dataset,
    "HumanPrimaryCellAtlasData" = celldex::HumanPrimaryCellAtlasData(),
    "BlueprintEncodeData" = celldex::BlueprintEncodeData(),
    "MouseRNAseqData" = celldex::MouseRNAseqData(),
    "ImmuneCellExpressionData" = celldex::ImmuneCellExpressionData(),
    "DatabaseImmuneCellExpressionData" = celldex::DatabaseImmuneCellExpressionData(),
    stop("Unsupported reference dataset")
  )
  
  # Get expression data
  expr_mat <- GetAssayData(seurat_obj, slot = "data", assay = assay)
  
  # Run SingleR
  if (clusters) {
    # Get cluster IDs
    clusters <- seurat_obj@meta.data[[cluster_col]]
    
    # Aggregate expression by cluster
    cluster_expr <- matrix(
      0,
      nrow = nrow(expr_mat),
      ncol = length(unique(clusters))
    )
    
    rownames(cluster_expr) <- rownames(expr_mat)
    colnames(cluster_expr) <- sort(unique(clusters))
    
    for (cl in colnames(cluster_expr)) {
      cluster_expr[, cl] <- rowMeans(expr_mat[, clusters == cl, drop = FALSE])
    }
    
    # Run SingleR on clusters
    singler_results <- SingleR(
      test = cluster_expr,
      ref = ref_data,
      labels = ref_data[[labels_to_use]],
      BPPARAM = BPPARAM
    )
    
    # Add results to metadata
    cluster_annotations <- data.frame(
      cluster = colnames(cluster_expr),
      cell_type = singler_results$labels,
      score = singler_results$scores[cbind(1:nrow(singler_results), max.col(singler_results$scores))],
      stringsAsFactors = FALSE
    )
    
    # Map cluster annotations to cells
    cell_types <- cluster_annotations$cell_type[match(clusters, cluster_annotations$cluster)]
    scores <- cluster_annotations$score[match(clusters, cluster_annotations$cluster)]
    
    seurat_obj$SingleR.clusters <- clusters
    seurat_obj$SingleR.cell_type <- cell_types
    seurat_obj$SingleR.score <- scores
    
  } else {
    # Run SingleR on individual cells
    singler_results <- SingleR(
      test = expr_mat,
      ref = ref_data,
      labels = ref_data[[labels_to_use]],
      BPPARAM = BPPARAM
    )
    
    # Add results to metadata
    seurat_obj$SingleR.cell_type <- singler_results$labels
    seurat_obj$SingleR.score <- singler_results$scores[cbind(1:nrow(singler_results), max.col(singler_results$scores))]
  }
  
  return(seurat_obj)
}

#' Add manual cell type annotations to a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param annotations Data frame with cluster to cell type mappings
#' @param cluster_col Column in seurat_obj metadata with cluster IDs
#' @param annotation_col Column in annotations data frame with cell type labels
#' @param name Name to use for the new metadata column
#' @return Seurat object with manual annotations added
add_manual_annotations <- function(seurat_obj,
                                 annotations,
                                 cluster_col = "seurat_clusters",
                                 annotation_col = "cell_type",
                                 name = "cell_type") {
  
  # Check if cluster column exists
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Cluster column", cluster_col, "not found in metadata"))
  }
  
  # Check if annotation column exists
  if (!annotation_col %in% colnames(annotations)) {
    stop(paste("Annotation column", annotation_col, "not found in annotations data frame"))
  }
  
  # Get cluster ID column from annotations
  cluster_id_col <- setdiff(colnames(annotations), annotation_col)[1]
  
  # Create mapping from cluster IDs to annotations
  cluster_to_annotation <- setNames(
    annotations[[annotation_col]],
    annotations[[cluster_id_col]]
  )
  
  # Map cluster IDs to annotations
  seurat_obj[[name]] <- cluster_to_annotation[as.character(seurat_obj[[cluster_col]])]
  
  return(seurat_obj)
}
