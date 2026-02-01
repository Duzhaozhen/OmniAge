#' Developmental Age Transformation
#'
#' @param x A vector of sample ages
#' @param adult.age The age considered to be the cutoff for adulthood
#'
#' @return transformed age prediction
anti.trafo <- function(x, adult.age = 20) {
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
}



#' Download OmniAgeR Example Data
#'
#' Downloads large example datasets (.rda) required for vignettes and examples.
#' Data is stored in the package's user data directory.
#'
#' @param datasets A character vector of datasets to download (e.g., "seu_gabitto_2024_filtered").
#' @param path Path to save data. Defaults to \code{get_OmniAgeR_path()}.
#' @param force Logical. If TRUE, overwrites existing files.
#' @export
download_OmniAgeR_example <- function(datasets = c("Hannum_example.rda","GA_example.rda","LungInv.rda","Tursiops_example.rda","TZH_example_CTF.rda","brain_frohlich_control_example_15donors.rda",
                                                   "seu_gabitto_2024_filtered","Yazar_CD4T_CD8T_example.rda","CTS_ExampleData_Liver.rda","CTS_MurphyGSE88890.rda",
                                                   "CTS_PaiGSE112179.rda"),
                                      path = NULL,
                                      force = FALSE, ...) {
  
  # List defining the sample data
  EXAMPLE_MANIFEST <- data.frame(
    name = c("Hannum_example", "GA_example","LungInv", "Tursiops_example", 
             "TZH_example_CTF", "brain_frohlich_control_example_15donors",
             "seu_gabitto_2024_filtered", "Yazar_CD4T_CD8T_example","CTS_ExampleData_Liver","CTS_MurphyGSE88890",
             "CTS_PaiGSE112179"),
    file = c("Hannum_example.rda", "GA_example.rda","LungInv.rda", "Tursiops_example.rda", 
             "TZH_example_CTF.rda", "brain_frohlich_control_example_15donors.rda",
             "seu_gabitto_2024_filtered.rda", "Yazar_CD4T_CD8T_example.rda","CTS_ExampleData_Liver.rda","CTS_MurphyGSE88890.rda",
             "CTS_PaiGSE112179.rda")
  )
  
  ZENODO_DOI <- "10.5281/zenodo.18287372" 
  
  if (is.null(path)) path <- get_OmniAgeR_path()
  
  files_to_download <- EXAMPLE_MANIFEST$file[EXAMPLE_MANIFEST$name %in% datasets]
  
  if(length(files_to_download) == 0) {
    files_to_download <- EXAMPLE_MANIFEST$file[EXAMPLE_MANIFEST$file %in% datasets]
  }
  
  if (length(files_to_download) == 0) {
    stop(paste("The dataset", datasets, "is not in the manifest. Please check the spelling."))
  }
  
  
  if (!force) {
    files_to_download <- files_to_download[!file.exists(file.path(path, files_to_download))]
  }
  
  if (length(files_to_download) == 0) {
    message("All requested example datasets already exist.")
    return(invisible(TRUE))
  }
  
  message("Connecting to Zenodo to download example datasets...")
  zenodo <- zen4R::ZenodoManager$new(logger = "INFO")
  rec <- zenodo$getRecordByDOI(ZENODO_DOI)
  
  for (f in files_to_download) {
    message("Downloading example: ", f, "...")
    rec$downloadFiles(path = path, files = list(f), ...)
  }
  
  return(invisible(TRUE))
}



#' Load OmniAgeR Example Data
#'
#' Safely load .rda example datasets from the local storage path.
#'
#' @param dataset_name Base name of the dataset (e.g., "seu_gabitto_2024_filtered").
#' @param envir The environment where the data should be loaded. Defaults to parent.frame().
#' @export
load_OmniAgeR_example <- function(dataset_name, envir = parent.frame()) {
  
  full_name <- paste0(dataset_name, ".rda")
  file_path <- file.path(get_OmniAgeR_path(), full_name)
  
  if (!file.exists(file_path)) {
    stop(sprintf("Example data '%s' not found. Please run `download_OmniAgeR_example('%s')` first.", 
                 full_name, dataset_name))
  }
  
  # load .rda file
  load(file_path, envir = envir)
  message(sprintf("Dataset '%s' loaded into environment.", dataset_name))
}
