
#' Inverse Age Transformation for Developmental Clocks
#'
#' Reverts a transformed age prediction (e.g., from a log-scale) back to
#' the chronological age scale, based on the Horvath (2013) method.
#'
#' @param x A numeric vector of transformed age predictions.
#' @param adult.age The numeric age defining the threshold for adulthood (default: 20).
#' @export
#' @return A numeric vector of age predictions reverted to the chronological age scale.
anti.trafo <- function(x, adult.age = 20) {
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
}




#' Impute and Append Missing CpGs for Clock Analysis
#'
#' Pre-processes a DNA methylation matrix for clock analysis by handling
#' both missing values (NA) and entirely missing CpGs.
#'
#' This function performs two key operations:
#' 1. Imputes missing values (NA) within existing CpGs using the specified method.
#' 2. Appends any CpGs required by the clock that are absent from the
#'    input matrix, filling them with the provided default values from \code{CpGs}.
#'
#' @param DNAm A numeric matrix of DNA methylation (beta) values.
#'   Rows must be samples and columns must be CpG IDs.
#' @param method The imputation method for missing (NA) values.
#'   Currently, only "mean" is supported.
#' @param CpGs A named numeric vector where names are the required CpG IDs
#'   and values are their default (e.g., population mean) beta values.
#' @param subset A logical value. If TRUE (default), NA imputation is
#'   only applied to the intersection of CpGs in \code{DNAm} and \code{CpGs}.
#'   If FALSE, it is applied to all CpGs present in \code{DNAm}.
#' @export
#' @return A numeric matrix, subsetted to the CpGs specified in \code{CpGs},
#'   with all missing values and missing columns imputed.

PCClocks_impute_DNAm <- function (DNAm, method = c("mean"), CpGs = NULL, subset = TRUE){
  #suppressWarnings(check_DNAm(DNAm))
  method <- match.arg(method)
  checkmate::assert_numeric(CpGs, min.len = 1, names = "unique",
                            null.ok = TRUE)
  checkmate::assert_logical(subset, len = 1, any.missing = FALSE,
                            null.ok = FALSE)
  if (is.null(CpGs)) {
    CpGs <- numeric(length = ncol(DNAm))
    names(CpGs) <- colnames(DNAm)
  }
  needed_cpgs <- setdiff(names(CpGs), colnames(DNAm))
  needed_matrix <- matrix(CpGs[needed_cpgs], ncol = length(needed_cpgs),
                          nrow = nrow(DNAm), byrow = TRUE, dimnames = list(row.names(DNAm),
                                                                           needed_cpgs))
  to_be_imputed <- if (subset) {
    intersect(names(CpGs), colnames(DNAm))
  }
  else {
    colnames(DNAm)
  }
  n_miss <- colSums(is.na(DNAm[, to_be_imputed, drop = F]))
  imputed_matrix <- if (sum(n_miss) == 0) {
    DNAm[, to_be_imputed, drop = F]
  }
  else if (method == "mean") {
    cbind(DNAm[, names(which(n_miss == 0)), drop = F], mean_impute(DNAm[,
                                                                        names(which(n_miss > 0)), drop = F]))
  }
  else {
    stop("Unsupported Method")
  }
  return(cbind(imputed_matrix, needed_matrix)[, c(to_be_imputed,
                                                  needed_cpgs), drop = F])
}



#' Robustly Load Serialized Package Data (.qs2)
#'
#' An internal helper function to safely load .qs2 data files from a local
#' path, validate their existence, and provide user-friendly error messages.
#'
#' @param object_name The base name of the data object (e.g., "PCClocks_data")
#'   to be loaded, excluding the ".qs2" extension.
#' @param path An optional character string specifying either the full path
#'   to the .qs2 file or the path to the directory containing it. If NULL,
#'   defaults to the path from \code{get_OmniAgeR_path()}.
#'
#' @return The deserialized R object (e.g., a list or data.frame)
#'   from the .qs2 file.
#' @export
#'
#'

#load_PCClocks_data
load_OmniAgeR_data <- function (object_name,
                                path = NULL){

  checkmate::assert_string(object_name)

  full_object_name <- paste0(object_name, ".qs2")
  checkmate::assert_character(path, len = 1, null.ok = TRUE)

  file_path <- if (!is.null(path)) {
    stopifnot(`Provided file/dir doesn't exists` = dir.exists(path) ||
                file.exists(path))
    info <- file.info(path)
    if (info$isdir) {
      file.path(path, full_object_name)
    }
    else {
      path
    }
  }
  else {
    file.path(get_OmniAgeR_path(), full_object_name)
  }


  tryCatch({
    obj <- qs2::qs_read(file_path, validate_checksum = TRUE)
  }, error = function(e) {

    stop(sprintf(paste("Failed to read '%s'.",
                       "Did you download '%s' and passed the path or loaded object to this function?",
                       "See function `download_OmniAgeR_data()` to download the required data.",
                       "Function failed: %s", sep = "\n"),
                 full_object_name, full_object_name, e))
  })

  return(obj)
}





#' Get Default Data Storage Directory
#'
#' Retrieves the OS-specific, cross-platform data directory for the
#' 'OmniAgeR' package using \code{tools::R_user_dir}.
#' Creates the directory if it does not already exist.
#'
#' @return A character string representing the file path to the
#'   package's data directory.
#'
#' @export
#'
get_OmniAgeR_path <- function() {

  path <- tools::R_user_dir("OmniAgeR", which = "data")

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  return(path)
}





#' Download OmniAgeR Data
#'
#' Downloads large data files (.qs2) required for clock calculations.
#' Data is downloaded from the Zenodo repository.
#'
#' @param clocks A character vector of clock datasets to download.
#'   Valid options (based on the function's default) are "SystemsAge" and "PCClocks".
#' @param path A character string specifying the directory to save the data.
#'   If NULL (default), uses the path from \code{get_OmniAgeR_path()}.
#' @param force A logical value. If TRUE, data will be re-downloaded
#'   and will overwrite any existing files. Default is FALSE.
#' @param ... Additional arguments passed to \code{zen4R::ZenodoRecord$downloadFiles()}.
#'
#' @return Invisibly, a logical vector where TRUE indicates
#'   a successful download for each file.
#' @export
download_OmniAgeR_data <- function (clocks = c("SystemsAge", "PCClocks"),
                                          path = NULL,
                                          force = FALSE,
                                          ...) {


  BIOMARKER_MANIFEST <- data.frame(
    clock_name = c("SystemsAge", "PCClocks"),
    file_name = c("SystemsAge_data.qs2", "PCClocks_data.qs2")
  )

  ZENODO_DOI <- "10.5281/zenodo.17162604" #


  clocks <- match.arg(clocks, several.ok = TRUE)

  if (is.null(path)) {

    path <- get_OmniAgeR_path()
    if (!dir.exists(path)) {
      message("Creating a new folder at ", path)
      if (!dir.create(path, showWarnings = TRUE, recursive = TRUE)) {
        stop("Failed to create directory at ", path)
      }
    }
  }
  checkmate::assert_directory_exists(path, access = "rw", .var.name = "path")


  if ("all" %in% clocks) {
    clocks <- BIOMARKER_MANIFEST$clock_name
  }
  clocks <- unique(clocks)


  manifest_subset <- BIOMARKER_MANIFEST[BIOMARKER_MANIFEST$clock_name %in% clocks, ]

  if(nrow(manifest_subset) == 0 && length(clocks) > 0) {
    warning("No valid clock names provided. Check spelling.")
    return(invisible(FALSE))
  }

  download_name <- manifest_subset$file_name
  download_to <- file.path(path, download_name)

  exists <- if (force) {
    rep(FALSE, times = length(download_to))
  } else {
    file.exists(download_to)
  }

  if (any(exists)) {
    message(paste("Skipping", sum(exists), "file(s) that already exist. Use force = TRUE to re-download."))

    clocks <- clocks[!exists]
    download_name <- download_name[!exists]
    download_to <- download_to[!exists]
  }

  if (length(clocks) == 0) {
    message("No new files to download.")
    return(invisible(TRUE))
  }


  if (!requireNamespace("zen4R", quietly = TRUE)) {
    stop("Please install the 'zen4R' package to download these files (`install.packages('zen4R')`).")
  }

  message("Connecting to Zenodo...")
  zenodo <- zen4R::ZenodoManager$new(logger = "INFO")

  rec <- zenodo$getRecordByDOI(ZENODO_DOI)

  if(is.null(rec)) {
    stop(paste("Could not find Zenodo record for DOI:", YOUR_ZENODO_DOI))
  }


  message(paste("Attempting to download", length(clocks), "file(s) to:", path))
  download_success <- logical(length(clocks))
  names(download_success) <- download_name

  for (i in seq_along(clocks)) {
    tryCatch({
      current_file_name <- download_name[i]
      message("Downloading ", current_file_name, "...")


      rec$downloadFiles(path = path,
                        files = list(current_file_name),
                        ...)

      if (!file.exists(download_to[i])) {
        stop(paste("File does not exist after download:", download_to[i]))
      }

      message("Successfully downloaded ", current_file_name)
      download_success[i] <- TRUE

    }, error = function(e) {
      message("Error: ", e$message)
      warning(paste("Failed to download:", download_name[i]), call. = FALSE)
      download_success[i] <- FALSE
    })
  }

  return(invisible(download_success))
}






