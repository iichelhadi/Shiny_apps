# R/utilities.R - Helper functions for the scRNA-seq app

#' Create a reusable tooltip with help text
#'
#' @param id The ID of the element to attach the tooltip to
#' @param title The tooltip text
#' @param placement Where to place the tooltip (top, bottom, left, right)
#' @return A bsTooltip object
create_tooltip <- function(id, title, placement = "right") {
  bsTooltip(id = id, title = title, placement = placement, trigger = "hover")
}

#' Format numeric values with a specified number of decimal places
#'
#' @param x Numeric value to format
#' @param digits Number of decimal places
#' @return Formatted string
format_numeric <- function(x, digits = 2) {
  format(round(x, digits), nsmall = digits)
}

#' Get the number of processors for parallel computation
#'
#' @param prop Proportion of available cores to use (0-1)
#' @return Number of cores to use
get_n_processors <- function(prop = 0.75) {
  max(1, floor(parallelly::availableCores(methods = 'nproc') * prop))
}

#' Create a multi-core BiocParallel processor
#'
#' @param prop Proportion of available cores to use (0-1)
#' @return BiocParallelParam object
create_bioc_processor <- function(prop = 0.75) {
  n_cores <- get_n_processors(prop)
  BiocParallel::MulticoreParam(workers = n_cores)
}

#' Create a progress object for long-running tasks
#'
#' @param session The Shiny session
#' @param min Minimum progress value
#' @param max Maximum progress value
#' @param message Progress message to display
#' @return A Shiny Progress object
create_progress <- function(session, min = 0, max = 1, message = "Processing...") {
  progress <- shiny::Progress$new(session = session, min = min, max = max)
  progress$set(message = message, value = min)
  return(progress)
}

#' Convert a CSS size string (e.g., "400px") to numeric
#'
#' @param size_str Size string with units
#' @return Numeric value without units
parse_size <- function(size_str) {
  as.numeric(gsub("[^0-9.]", "", size_str))
}

#' Check if a string matches a pattern
#'
#' @param string String to check
#' @param pattern Pattern to match
#' @return TRUE if the string matches the pattern, FALSE otherwise
string_matches <- function(string, pattern) {
  grepl(pattern, string)
}

#' Generate a time-stamped filename
#'
#' @param prefix Prefix for the filename
#' @param extension File extension
#' @return Filename with timestamp
timestamped_filename <- function(prefix, extension = "csv") {
  paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", extension)
}

#' Format a data frame for display in DT::datatable
#'
#' @param df Data frame to format
#' @param digits Number of decimal places for numeric columns
#' @return Formatted data frame
format_dt_table <- function(df, digits = 3) {
  df_out <- df
  numeric_cols <- sapply(df, is.numeric)
  df_out[, numeric_cols] <- lapply(df_out[, numeric_cols, drop = FALSE], function(x) {
    format(round(x, digits), nsmall = digits)
  })
  return(df_out)
}

#' Get suggested number of PCs based on elbow method
#'
#' @param seurat_obj Seurat object with PCA results
#' @param variance_cutoff Cumulative variance cutoff
#' @param pct_cutoff Individual PC variance contribution cutoff
#' @return Suggested number of PCs
get_suggested_pcs <- function(seurat_obj, variance_cutoff = 70, pct_cutoff = 5) {
  pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  pcs <- which(cumu > variance_cutoff & pct < pct_cutoff)[1]
  if(is.na(pcs)) {
    # If no PCs meet the criteria, use a minimum of 10 or 1/3 of total PCs
    pcs <- min(10, floor(length(pct)/3))
  }
  return(pcs)
}

#' Save a Seurat object to a file with compression
#'
#' @param seurat_obj Seurat object to save
#' @param file File path to save to
#' @return Invisibly returns the file path
save_seurat <- function(seurat_obj, file) {
  saveRDS(seurat_obj, file = file, compress = TRUE)
  invisible(file)
}

#' Load a Seurat object from a file
#'
#' @param file File path to load from
#' @return Loaded Seurat object
load_seurat <- function(file) {
  readRDS(file)
}

#' Check if required packages are installed
#'
#' @param packages Character vector of package names
#' @return TRUE if all packages are installed, FALSE otherwise
check_packages <- function(packages) {
  installed <- vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  return(all(installed))
}
