# R/data_io.R

#’ Read 10X Genomics data into a Seurat object
#’
#’ @param path Path to the 10X directory containing matrix, barcodes, and features
#’ @param min.cells Include features detected in at least this many cells
#’ @param min.features Include cells where at least this many features are detected
#’ @param project Project name for the Seurat object
#’ @return A Seurat object with QC metrics calculated
read_10x <- function(path, min.cells = 3, min.features = 200, project = "scRNA") {
  raw_counts <- Seurat::Read10X(data.dir = path)
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts     = raw_counts,
    min.cells  = min.cells,
    min.features = min.features,
    project    = project
  )
  rm(raw_counts); gc()
  
  # compute QC metrics
  seurat_obj[["percent.mt"]]   <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  seurat_obj[["percent.hsp"]]  <- PercentageFeatureSet(seurat_obj, pattern = "^HSP")
  seurat_obj$log10GenesPerUMI  <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  
  return(seurat_obj)
}

#' Filter Seurat object based on QC metrics
#'
#' @param seurat_obj Seurat object with QC metrics
#' @param min.features Minimum number of features per cell
#' @param max.features Maximum number of features per cell
#' @param min.counts Minimum number of counts per cell
#' @param max.counts Maximum number of counts per cell
#' @param max.mt Maximum percentage of mitochondrial genes
#' @param min.genes_per_umi Minimum complexity (log10GenesPerUMI)
#' @return Filtered Seurat object
filter_seurat <- function(seurat_obj, 
                         min.features = 500, 
                         max.features = 5000,
                         min.counts = 1000,
                         max.counts = 40000,
                         max.mt = 20,
                         min.genes_per_umi = 0.8) {
  
  # Apply filters
  seurat_filtered <- subset(
    seurat_obj,
    subset = nFeature_RNA >= min.features &
      nFeature_RNA <= max.features &
      nCount_RNA >= min.counts &
      nCount_RNA <= max.counts &
      percent.mt <= max.mt &
      log10GenesPerUMI >= min.genes_per_umi
  )
  
  return(seurat_filtered)
}

#' Save a Seurat object to an RDS file
#'
#' @param seurat_obj Seurat object to save
#' @param file File path to save to
#' @param verbose Print message about saving
#' @return Invisibly returns TRUE on success
save_seurat_object <- function(seurat_obj, file, verbose = TRUE) {
  saveRDS(seurat_obj, file = file, compress = TRUE)
  if (verbose) {
    message(paste("Saved Seurat object to", file))
  }
  invisible(TRUE)
}

#' Load a Seurat object from an RDS file
#'
#' @param file File path to load from
#' @param verbose Print message about loading
#' @return Loaded Seurat object
load_seurat_object <- function(file, verbose = TRUE) {
  if (verbose) {
    message(paste("Loading Seurat object from", file))
  }
  readRDS(file)
}

#' Export Seurat data to various formats
#'
#' @param seurat_obj Seurat object
#' @param export_type Type of data to export ("counts", "metadata", "markers", "embeddings")
#' @param file File path to save to
#' @param slots Slots to export if exporting counts (defaults to "counts")
#' @param markers Markers data frame if exporting markers
#' @return Invisibly returns TRUE on success
export_seurat_data <- function(seurat_obj, export_type, file, slots = "counts", markers = NULL) {
  switch(
    export_type,
    "counts" = {
      counts <- GetAssayData(seurat_obj, slot = slots)
      write.csv(as.matrix(counts), file = file)
    },
    "metadata" = {
      metadata <- seurat_obj@meta.data
      write.csv(metadata, file = file)
    },
    "markers" = {
      if (!is.null(markers)) {
        write.csv(markers, file = file)
      } else {
        stop("Markers data frame must be provided for export_type='markers'")
      }
    },
    "embeddings" = {
      # Get all dimension reduction embeddings
      reductions <- names(seurat_obj@reductions)
      embeddings_list <- lapply(reductions, function(red) {
        as.data.frame(Embeddings(seurat_obj, reduction = red))
      })
      names(embeddings_list) <- reductions
      
      # Combine with metadata
      result <- cbind(
        seurat_obj@meta.data,
        do.call(cbind, embeddings_list)
      )
      
      write.csv(result, file = file)
    },
    {
      stop("Invalid export_type. Must be one of 'counts', 'metadata', 'markers', or 'embeddings'")
    }
  )
  
  invisible(TRUE)
}

#' Load a gene set from a GMT file
#'
#' @param file File path to GMT file
#' @return List of gene sets
load_gene_set <- function(file) {
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }
  
  lines <- readLines(file)
  gene_sets <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    list(
      name = parts[1],
      description = parts[2],
      genes = parts[3:length(parts)]
    )
  })
  
  names(gene_sets) <- sapply(gene_sets, function(x) x$name)
  
  return(gene_sets)
}

#' Create sample metadata from file names
#'
#' @param seurat_obj Seurat object
#' @param sample_pattern Regular expression pattern to extract sample info
#' @param sample_column Name of the column to store sample info
#' @return Seurat object with sample metadata added
add_sample_metadata <- function(seurat_obj, sample_pattern = NULL, sample_column = "sample") {
  if (is.null(sample_pattern)) {
    # If no pattern provided, use cell names to guess samples
    cell_names <- colnames(seurat_obj)
    
    # Try to extract sample names from cell barcodes based on common formats
    # This assumes format like SAMPLENAME_BARCODE or SAMPLENAME-BARCODE
    sample_names <- sub("_.*", "", cell_names)
    sample_names <- sub("-.*", "", sample_names)
    
    # If all sample names are the same, they probably don't have sample identifiers
    if (length(unique(sample_names)) == 1) {
      warning("Could not automatically determine samples. Using 'sample1' for all cells.")
      sample_names <- rep("sample1", length(cell_names))
    }
  } else {
    # Extract sample names using the provided pattern
    cell_names <- colnames(seurat_obj)
    sample_names <- rep(NA_character_, length(cell_names))
    
    # Apply the pattern to each cell name
    for (i in seq_along(cell_names)) {
      match <- regexpr(sample_pattern, cell_names[i], perl = TRUE)
      if (match > 0) {
        sample_names[i] <- regmatches(cell_names[i], match)
      }
    }
    
    # Check if pattern worked
    if (all(is.na(sample_names))) {
      warning("Pattern did not match any cell names. Using 'sample1' for all cells.")
      sample_names <- rep("sample1", length(cell_names))
    } else {
      # Replace NAs with default
      sample_names[is.na(sample_names)] <- "unknown"
    }
  }
  
  # Add to metadata
  seurat_obj[[sample_column]] <- sample_names
  
  return(seurat_obj)
}
