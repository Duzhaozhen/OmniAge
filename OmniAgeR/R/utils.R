#' Load Pre-trained Models and Example Data for OmniAgeR
#'
#' @description
#' This function seamlessly loads specific aging omic clock models. It checks 
#' locally bundled models first, fetches large models from Zenodo (caching 
#' them in the local user data directory), and falls back to OmniAgeRData if needed.
#'
#' @param title A character string specifying the exact name of the model.
#' @param verbose A logical flag. If `TRUE` (default), prints status messages.
#' @param force A logical flag. If `TRUE`, forces re-downloading the large model even if it exists locally.
#'
#' @return An R object containing the requested model parameters or reference data.
#' @export
#' @importFrom utils download.file
#' @examples
#' # Load the Horvath2013 model weights
#' horvath2013Model <- loadOmniAgeRdata("omniager_horvath2013_coef")

loadOmniAgeRdata <- function(title, verbose = FALSE, force = FALSE) {
  
  # -------------------------------------------------------------------------
  # STEP 1: Check for Small Models bundled directly within OmniAgeR
  # -------------------------------------------------------------------------
  local_path <- system.file("extdata", paste0(title, ".rds"), package = "OmniAgeR")
  
  if (local_path != "" && file.exists(local_path)) {
    if (verbose) message("[OmniAgeR] Loading file: ", title)
    return(readRDS(local_path))
  }
  
  # -------------------------------------------------------------------------
  # STEP 2: Check for Large Models hosted on Zenodo (Base R with Shadow File Caching)
  # -------------------------------------------------------------------------
  zenodo_urls <- list(
    "omniager_hannum_example" = "https://zenodo.org/api/records/18832408/files/omniager_hannum_example.rds/content",
    "PCClocks_data"           = "https://zenodo.org/api/records/18832408/files/PCClocks_data.qs2/content",
    "SystemsAge_data"         = "https://zenodo.org/api/records/18832408/files/SystemsAge_data.qs2/content"
  )
  
  if (title %in% names(zenodo_urls)) {
    targetUrl <- zenodo_urls[[title]]
    file_ext <- ifelse(grepl("\\.qs2", targetUrl), ".qs2", ".rds")
    filename <- paste0(title, file_ext)
    
    cache_dir <- tools::R_user_dir("OmniAgeR", which = "data")
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    dest_file <- file.path(cache_dir, filename)
    
    if (force || !file.exists(dest_file)) {
      if (verbose) {
        message("[OmniAgeR] Downloading large data: ", filename)
        message("          Saving to: ", cache_dir)
      }
      
      # 1. Extend the download timeout of the basic R
      old_timeout <- getOption("timeout")
      options(timeout = max(60, old_timeout)) 
      
      # 2. Shadow File Path
      tmp_file <- paste0(dest_file, ".tmp")
      
      tryCatch({
        # 3. Check the download status and file validity
        status <- download.file(url = targetUrl, destfile = tmp_file, mode = "wb", quiet = !verbose)
        
        # Check the download status and file validity
        if (status != 0) {
          stop("download.file returned a non-zero status.")
        }
        if (!file.exists(tmp_file) || file.info(tmp_file)$size == 0) {
          stop("Downloaded file is missing or 0 bytes.")
        }
        
        # 4. The download was a perfect success! Rename.tmp to the official file
        file.rename(from = tmp_file, to = dest_file)
        
      }, error = function(e) {
        # 5. If an error occurs, delete this unfinished.tmp file
        if (file.exists(tmp_file)) {
          if (verbose) message("[OmniAgeR] Cleaning up partial shadow file...")
          unlink(tmp_file)
        }
        error_msg <- paste0(
          "\n",
          "====================================================================\n",
          "[OmniAgeR] AUTOMATED DOWNLOAD FAILED \n",
          "====================================================================\n",
          "Your internet connection to Zenodo was interrupted. \n",
          "Please manually download the data to continue:\n\n",
          "STEP 1: Open your web browser and download this file:\n",
          "        ", targetUrl, "\n\n",
          "STEP 2: Make sure the downloaded file is named EXACTLY:\n",
          "        ", filename, "\n\n",
          "STEP 3: Move the downloaded file into this specific directory:\n",
          "        ", cache_dir, "\n\n",
          "STEP 4: Run your R code again. OmniAgeR will automatically detect it.\n",
          "====================================================================\n"
        )
        
        stop(error_msg, call. = FALSE)
        
      }, finally = {
        # 6. Restore the user's default timeout setting
        options(timeout = old_timeout)
      })
      
    } else {
      if (verbose) message("[OmniAgeR] Loading large data from local cache: ", dest_file)
    }
    

    if (file_ext == ".qs2") {
      if (!requireNamespace("qs2", quietly = TRUE)) {
        stop("[OmniAgeR] The 'qs2' package is required to load .qs2 data. Please install it.")
      }
      return(qs2::qs_read(dest_file))
    } else {
      return(readRDS(dest_file))
    }
  }
  
  # -------------------------------------------------------------------------
  # STEP 3: If all fail
  # -------------------------------------------------------------------------
  stop(
    "[OmniAgeR] Data '", title, 
    "' could not be loaded. It is not built-in, not found on Zenodo."
  )
}
#' Developmental Age Transformation
#'
#' @param x A vector of sample ages
#' @param adultAge The age considered to be the cutoff for adulthood
#'
#' @return transformed age prediction
#' @noRd
.antiTrafo <- function(x, adultAge = 20) {
    ifelse(x < 0, (1 + adultAge) * exp(x) - 1, (1 + adultAge) * x + adultAge)
}
