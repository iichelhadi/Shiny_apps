# R/analysis.R - Functions for data analysis

#' Find marker genes for all clusters
#'
#' @param seurat_obj Seurat object
#' @param assay Assay to use
#' @param test.use Statistical test to use
#' @param min.pct Minimum percentage of cells expressing the gene
#' @param logfc.threshold Log fold-change threshold
#' @param only.pos Only return positive markers
#' @param add.pct.diff Add percentage difference between clusters
#' @return A data frame with marker genes
find_all_markers <- function(seurat_obj,
                            assay = NULL,
                            test.use = "MAST",
                            min.pct = 0.1,
                            logfc.threshold = 0.25,
                            only.pos = TRUE,
                            add.pct.diff = TRUE) {
  
  # Determine which assay to use
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  
  # Set active assay
  DefaultAssay(seurat_obj) <- assay
  
  # Find markers
  markers <- FindAllMarkers(
    seurat_obj,
    assay = assay,
    test.use = test.use,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    only.pos = only.pos,
    verbose = TRUE
  )
  
  # Add percentage difference if requested
  if (add.pct.diff && "scCustomize" %in% .packages(TRUE)) {
    markers <- scCustomize::Add_Pct_Diff(markers)
    markers <- markers %>% arrange(cluster, desc(pct_diff))
  } else {
    markers <- markers %>% arrange(cluster, desc(avg_log2FC))
  }
  
  return(markers)
}

#' Find marker genes for specific clusters
#'
#' @param seurat_obj Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 Identity class to compare to
#' @param assay Assay to use
#' @param test.use Statistical test to use
#' @param min.pct Minimum percentage of cells expressing the gene
#' @param logfc.threshold Log fold-change threshold
#' @param only.pos Only return positive markers
#' @param add.pct.diff Add percentage difference between clusters
#' @return A data frame with marker genes
find_markers <- function(seurat_obj,
                        ident.1,
                        ident.2 = NULL,
                        assay = NULL,
                        test.use = "MAST",
                        min.pct = 0.1,
                        logfc.threshold = 0.25,
                        only.pos = TRUE,
                        add.pct.diff = TRUE) {
  
  # Determine which assay to use
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  
  # Set active assay
  DefaultAssay(seurat_obj) <- assay
  
  # Find markers
  markers <- FindMarkers(
    seurat_obj,
    ident.1 = ident.1,
    ident.2 = ident.2,
    assay = assay,
    test.use = test.use,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    only.pos = only.pos,
    verbose = TRUE
  )
  
  # Add percentage difference if requested
  if (add.pct.diff && "scCustomize" %in% .packages(TRUE)) {
    markers <- markers %>% 
      rownames_to_column("gene") %>%
      scCustomize::Add_Pct_Diff() %>%
      arrange(desc(pct_diff))
  } else {
    markers <- markers %>% 
      rownames_to_column("gene") %>%
      arrange(desc(avg_log2FC))
  }
  
  return(markers)
}

#' Find conserved markers across conditions
#'
#' @param seurat_obj Seurat object
#' @param cluster.col Column containing cluster IDs
#' @param group.col Column containing group/condition IDs
#' @param assay Assay to use
#' @param test.use Statistical test to use
#' @param min.pct Minimum percentage of cells expressing the gene
#' @param logfc.threshold Log fold-change threshold
#' @param only.pos Only return positive markers
#' @return A data frame with conserved marker genes
find_conserved_markers <- function(seurat_obj,
                                 cluster.col = "seurat_clusters",
                                 group.col,
                                 assay = NULL,
                                 test.use = "MAST",
                                 min.pct = 0.1,
                                 logfc.threshold = 0.25,
                                 only.pos = TRUE) {
  
  # Determine which assay to use
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  
  # Check if columns exist
  if (!cluster.col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Cluster column", cluster.col, "not found in metadata"))
  }
  if (!group.col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Group column", group.col, "not found in metadata"))
  }
  
  # Set cluster identity
  Idents(seurat_obj) <- cluster.col
  
  # Get all clusters
  clusters <- unique(seurat_obj@meta.data[[cluster.col]])
  
  # Initialize list to store markers
  conserved_markers_list <- list()
  
  # Find conserved markers for each cluster
  for (cluster in clusters) {
    conserved_markers_list[[as.character(cluster)]] <- FindConservedMarkers(
      seurat_obj,
      ident.1 = cluster,
      grouping.var = group.col,
      assay = assay,
      test.use = test.use,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      only.pos = only.pos,
      verbose = TRUE
    ) %>%
      rownames_to_column("gene") %>%
      mutate(cluster = cluster)
  }
  
  # Combine all markers
  conserved_markers <- bind_rows(conserved_markers_list)
  
  return(conserved_markers)
}

#' Run gene set enrichment analysis (GSEA) on marker genes
#'
#' @param markers Data frame with marker genes
#' @param gene_col Column name containing gene names/IDs
#' @param score_col Column name containing gene scores (log2FC)
#' @param organism Organism for MSigDB gene sets ("human" or "mouse")
#' @param category MSigDB gene set category (default: "H" for hallmark)
#' @param subcategory MSigDB gene set subcategory (default: NULL)
#' @param min_size Minimum size of gene set
#' @param max_size Maximum size of gene set
#' @param nperm Number of permutations
#' @return A list with GSEA results
run_gsea <- function(markers,
                    gene_col = "gene",
                    score_col = "avg_log2FC",
                    organism = "human",
                    category = "H",
                    subcategory = NULL,
                    min_size = 15,
                    max_size = 500,
                    nperm = 1000) {
  
  # Check if required packages are available
  if (!requireNamespace("fgsea", quietly = TRUE) || 
      !requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Packages 'fgsea' and 'msigdbr' are required for GSEA analysis")
  }
  
  # Get MSigDB gene sets
  if (is.null(subcategory)) {
    gene_sets <- msigdbr::msigdbr(species = organism, category = category) %>%
      dplyr::select(gs_name, gene_symbol) %>%
      dplyr::group_by(gs_name) %>%
      dplyr::summarize(genes = list(gene_symbol))
  } else {
    gene_sets <- msigdbr::msigdbr(species = organism, 
                                 category = category, 
                                 subcategory = subcategory) %>%
      dplyr::select(gs_name, gene_symbol) %>%
      dplyr::group_by(gs_name) %>%
      dplyr::summarize(genes = list(gene_symbol))
  }
  
  # Convert gene sets to list format
  gene_sets_list <- setNames(gene_sets$genes, gene_sets$gs_name)
  
  # Create ranked gene list
  ranked_genes <- markers %>%
    dplyr::select(!!sym(gene_col), !!sym(score_col)) %>%
    dplyr::arrange(desc(!!sym(score_col))) %>%
    deframe()
  
  # Run GSEA
  gsea_results <- fgsea::fgsea(
    pathways = gene_sets_list,
    stats = ranked_genes,
    minSize = min_size,
    maxSize = max_size,
    nperm = nperm
  )
  
  # Sort results by significance
  gsea_results <- gsea_results %>%
    dplyr::arrange(pval)
  
  return(gsea_results)
}

#' Run over-representation analysis (ORA) on marker genes
#'
#' @param markers Data frame with marker genes
#' @param gene_col Column name containing gene names/IDs
#' @param universe Vector of all genes in the dataset (background genes)
#' @param organism Organism for MSigDB gene sets ("human" or "mouse")
#' @param category MSigDB gene set category (default: "H" for hallmark)
#' @param subcategory MSigDB gene set subcategory (default: NULL)
#' @param p_cutoff P-value cutoff for significance
#' @param q_cutoff Q-value cutoff for significance
#' @return A data frame with ORA results
run_ora <- function(markers,
                   gene_col = "gene",
                   universe = NULL,
                   organism = "human",
                   category = "H",
                   subcategory = NULL,
                   p_cutoff = 0.05,
                   q_cutoff = 0.2) {
  
  # Check if required packages are available
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || 
      !requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Packages 'clusterProfiler' and 'msigdbr' are required for ORA analysis")
  }
  
  # Get marker genes
  gene_list <- markers[[gene_col]]
  
  # Get MSigDB gene sets
  if (is.null(subcategory)) {
    msigdb_df <- msigdbr::msigdbr(species = organism, category = category)
  } else {
    msigdb_df <- msigdbr::msigdbr(species = organism, 
                                 category = category, 
                                 subcategory = subcategory)
  }
  
  # Create term2gene and term2name data frames
  term2gene <- msigdb_df %>% dplyr::select(gs_name, gene_symbol)
  term2name <- msigdb_df %>% dplyr::select(gs_name, gs_description) %>% distinct()
  
  # Run ORA
  ora_results <- clusterProfiler::enricher(
    gene = gene_list,
    universe = universe,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff
  )
  
  # Convert to data frame
  if (!is.null(ora_results) && nrow(ora_results@result) > 0) {
    ora_df <- ora_results@result
  } else {
    ora_df <- data.frame()
  }
  
  return(ora_df)
}

#' Calculate cell cycle scores
#'
#' @param seurat_obj Seurat object
#' @param assay Assay to use
#' @param s.features S phase genes
#' @param g2m.features G2/M phase genes
#' @return Seurat object with cell cycle scores
calculate_cell_cycle_scores <- function(seurat_obj,
                                       assay = NULL,
                                       s.features = NULL,
                                       g2m.features = NULL) {
  
  # Determine which assay to use
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  
  # Set default gene lists if not provided
  if (is.null(s.features) || is.null(g2m.features)) {
    data("cc.genes.updated.2019", package = "Seurat", envir = environment())
    s.features <- cc.genes.updated.2019$s.genes
    g2m.features <- cc.genes.updated.2019$g2m.genes
  }
  
  # Calculate cell cycle scores
  seurat_obj <- CellCycleScoring(
    seurat_obj,
    s.features = s.features,
    g2m.features = g2m.features,
    assay = assay
  )
  
  return(seurat_obj)
}

#' Perform differential abundance analysis
#'
#' @param seurat_obj Seurat object
#' @param cluster_col Column containing cluster IDs
#' @param condition_col Column containing condition/group IDs
#' @param reference_level Reference level for condition
#' @param min_cells Minimum number of cells per cluster-condition group
#' @return A data frame with differential abundance results
calculate_differential_abundance <- function(seurat_obj,
                                            cluster_col = "seurat_clusters",
                                            condition_col,
                                            reference_level = NULL,
                                            min_cells = 10) {
  
  # Check if required packages are available
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for differential abundance analysis")
  }
  
  # Check if columns exist
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Cluster column", cluster_col, "not found in metadata"))
  }
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Condition column", condition_col, "not found in metadata"))
  }
  
  # Create contingency table
  cell_counts <- table(seurat_obj@meta.data[[cluster_col]], 
                      seurat_obj@meta.data[[condition_col]])
  
  # Filter out clusters with too few cells
  keep_clusters <- rowSums(cell_counts >= min_cells) == ncol(cell_counts)
  cell_counts <- cell_counts[keep_clusters, ]
  
  # Create design matrix
  conditions <- colnames(cell_counts)
  if (is.null(reference_level)) {
    reference_level <- conditions[1]
  }
  
  # Reorder conditions to put reference level first
  conditions <- c(reference_level, setdiff(conditions, reference_level))
  cell_counts <- cell_counts[, conditions]
  
  # Create DGEList object
  y <- edgeR::DGEList(counts = cell_counts)
  
  # Create design matrix
  design <- model.matrix(~ factor(conditions))
  colnames(design) <- c("Intercept", paste0(conditions[-1], "_vs_", reference_level))
  
  # Estimate dispersion
  y <- edgeR::estimateDisp(y, design)
  
  # Fit model
  fit <- edgeR::glmQLFit(y, design)
  
  # Test for differential abundance
  results_list <- list()
  
  for (i in 2:ncol(design)) {
    contrast_name <- colnames(design)[i]
    qlf <- edgeR::glmQLFTest(fit, coef = i)
    results <- edgeR::topTags(qlf, n = nrow(y$counts))$table
    results$contrast <- contrast_name
    results$cluster <- rownames(results)
    results_list[[contrast_name]] <- results
  }
  
  # Combine results
  all_results <- bind_rows(results_list)
  
  # Calculate fold changes in proportions
  prop_table <- prop.table(cell_counts, margin = 2)
  
  for (contrast in unique(all_results$contrast)) {
    condition <- strsplit(contrast, "_vs_")[[1]][1]
    ref <- strsplit(contrast, "_vs_")[[1]][2]
    
    clusters <- all_results$cluster[all_results$contrast == contrast]
    for (cluster in clusters) {
      idx <- which(all_results$cluster == cluster & all_results$contrast == contrast)
      
      if (length(idx) == 1) {
        # Calculate proportional fold change
        prop_fc <- prop_table[cluster, condition] / prop_table[cluster, ref]
        all_results$prop_FC[idx] <- prop_fc
      }
    }
  }
  
  # Sort results
  all_results <- all_results %>%
    dplyr::arrange(contrast, PValue)
  
  return(all_results)
}

#' Perform integration of multiple Seurat objects
#'
#' @param seurat_list List of Seurat objects
#' @param integration_method Integration method ("CCA", "RPCA", or "harmony")
#' @param features.to.use Features to use for integration (default: 2000)
#' @param dims Number of dimensions to use
#' @param normalization_method Normalization method ("LogNormalize" or "SCT")
#' @param batch_var Variable name for harmony integration
#' @param reference Optional reference dataset index for integration
#' @param verbose Print progress messages
#' @return Integrated Seurat object
integrate_seurat_objects <- function(seurat_list,
                                    integration_method = "CCA",
                                    features.to.use = 2000,
                                    dims = 30,
                                    normalization_method = "LogNormalize",
                                    batch_var = "orig.ident",
                                    reference = NULL,
                                    verbose = TRUE) {
  
  # Check input
  if (length(seurat_list) < 2) {
    stop("At least two Seurat objects required for integration")
  }
  
  # Check if the integration method is valid
  if (!integration_method %in% c("CCA", "RPCA", "harmony")) {
    stop("Invalid integration method. Use 'CCA', 'RPCA', or 'harmony'")
  }
  
  # Use SCTransform if requested
  if (normalization_method == "SCT") {
    if (verbose) message("Applying SCTransform to all objects...")
    
    # Process each object with SCTransform
    seurat_list <- lapply(seurat_list, function(x) {
      SCTransform(x, verbose = verbose)
    })
    
    if (integration_method %in% c("CCA", "RPCA")) {
      # Select integration features
      features <- SelectIntegrationFeatures(
        object.list = seurat_list,
        nfeatures = features.to.use,
        verbose = verbose
      )
      
      # Prepare for integration
      seurat_list <- PrepSCTIntegration(
        object.list = seurat_list,
        anchor.features = features,
        verbose = verbose
      )
      
      # Find integration anchors
      if (is.null(reference)) {
        anchors <- FindIntegrationAnchors(
          object.list = seurat_list,
          normalization.method = "SCT",
          anchor.features = features,
          dims = 1:dims,
          reduction = integration_method,
          verbose = verbose
        )
      } else {
        anchors <- FindIntegrationAnchors(
          object.list = seurat_list,
          normalization.method = "SCT",
          anchor.features = features,
          dims = 1:dims,
          reduction = integration_method,
          reference = reference,
          verbose = verbose
        )
      }
      
      # Integrate data
      integrated <- IntegrateData(
        anchorset = anchors,
        normalization.method = "SCT",
        dims = 1:dims,
        verbose = verbose
      )
      
      # Set default assay
      DefaultAssay(integrated) <- "integrated"
      
    } else if (integration_method == "harmony") {
      # Merge objects
      merged <- merge(
        x = seurat_list[[1]],
        y = seurat_list[-1],
        add.cell.ids = names(seurat_list),
        project = "integrated"
      )
      
      # Run PCA
      merged <- RunPCA(
        merged,
        npcs = dims,
        verbose = verbose
      )
      
      # Run Harmony
      integrated <- RunHarmony(
        merged,
        group.by.vars = batch_var,
        dims.use = 1:dims,
        verbose = verbose
      )
      
      # Update reduction
      integrated@reductions$pca <- integrated@reductions$harmony
    }
    
  } else { # Use standard normalization
    # Normalize and find variable features for each object
    if (verbose) message("Normalizing and finding variable features...")
    
    seurat_list <- lapply(seurat_list, function(x) {
      x <- NormalizeData(x, verbose = verbose)
      x <- FindVariableFeatures(
        x,
        nfeatures = features.to.use,
        verbose = verbose
      )
      return(x)
    })
    
    if (integration_method %in% c("CCA", "RPCA")) {
      # Select integration features
      features <- SelectIntegrationFeatures(
        object.list = seurat_list,
        nfeatures = features.to.use,
        verbose = verbose
      )
      
      # Find integration anchors
      if (is.null(reference)) {
        anchors <- FindIntegrationAnchors(
          object.list = seurat_list,
          dims = 1:dims,
          anchor.features = features,
          reduction = integration_method,
          verbose = verbose
        )
      } else {
        anchors <- FindIntegrationAnchors(
          object.list = seurat_list,
          dims = 1:dims,
          anchor.features = features,
          reduction = integration_method,
          reference = reference,
          verbose = verbose
        )
      }
      
      # Integrate data
      integrated <- IntegrateData(
        anchorset = anchors,
        dims = 1:dims,
        verbose = verbose
      )
      
      # Set default assay
      DefaultAssay(integrated) <- "integrated"
      
    } else if (integration_method == "harmony") {
      # Merge objects
      merged <- merge(
        x = seurat_list[[1]],
        y = seurat_list[-1],
        add.cell.ids = names(seurat_list),
        project = "integrated"
      )
      
      # Normalize and scale the merged object
      merged <- NormalizeData(merged, verbose = verbose)
      merged <- FindVariableFeatures(
        merged,
        nfeatures = features.to.use,
        verbose = verbose
      )
      merged <- ScaleData(merged, verbose = verbose)
      
      # Run PCA
      merged <- RunPCA(
        merged,
        npcs = dims,
        verbose = verbose
      )
      
      # Run Harmony
      integrated <- RunHarmony(
        merged,
        group.by.vars = batch_var,
        dims.use = 1:dims,
        verbose = verbose
      )
    }
  }
  
  # Run dimensional reduction
  if (verbose) message("Running dimensional reduction...")
  
  # Scale data if not using SCTransform
  if (normalization_method != "SCT") {
    integrated <- ScaleData(integrated, verbose = verbose)
  }
  
  # Run UMAP
  if (integration_method == "harmony") {
    # Use harmony embeddings for UMAP
    integrated <- RunUMAP(
      integrated,
      reduction = "harmony",
      dims = 1:dims,
      verbose = verbose
    )
  } else {
    # Use PCA for UMAP
    integrated <- RunUMAP(
      integrated,
      reduction = "pca",
      dims = 1:dims,
      verbose = verbose
    )
  }
  
  return(integrated)
}
