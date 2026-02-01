#' @title Predict DNA methylation age using CTS clocks
#' @description
#' This is the function for computing DNA methylation age using CTS (cell type
#' specific) clocks. The inputs including DNAm matrix, the CTS clock you want to
#' use, is your DNAm data from bulk tissue sample or sorted cell sample, cell
#' type fraction matrix if you want to use Neu-In/Glia-In/Brain clock, tissue of
#' your DNAm data samples and the number of cores if you want to do parallel computing.
#'
#' @details
#' This function supports a variety of Cell-Type-Specific (CTS) clocks.
#'
#' **Available `CTSclocks` include:**
#' * `'Neu-In'`
#' * `'Glia-In'`
#' * `'Brain'`
#' * `'Neu-Sin'`
#' * `'Glia-Sin'`
#' * `'Hep'`
#' * `'Liver'`
#''
#'  The clocks are grouped below based on their biological target:
#'
#' **1. Cell-Type Specific Clocks**
#'
#' These clocks are trained to measure aging in specific cell populations.
#'
#' * `'Neu-In'` (Intrinsic): Measures cell-intrinsic aging of **Neurons**.
#'     (Uses processed data: residuals/Z-scores).
#' * `'Glia-In'` (Intrinsic): Measures cell-intrinsic aging of **Glial cells**.
#'     (Uses processed data: residuals/Z-scores).
#' * `'Neu-Sin'` (Semi-intrinsic): Measures aging of **Neurons** using raw data.
#'     (Reflects both intrinsic aging and composition changes).
#' * `'Glia-Sin'` (Semi-intrinsic): Measures aging of **Glial cells** using raw data.
#'     (Reflects both intrinsic aging and composition changes).
#' * `'Hep'` (Semi-intrinsic): Measures aging of **Hepatocytes** using raw data.
#'
#' **2. Non Cell-Type Specific Clocks**
#'
#' * `'Brain'` (Intrinsic): A intrinsic clock for **whole brain tissue**.
#'     (Uses processed data: residuals/Z-scores).
#' * `'Liver'` (Semi-intrinsic): A intrinsic clock for **whole liver tissue**.
#'     (Uses raw data).
#'
#' @param data.m  A DNAm matrix (row: CpGs, column: samples) of the samples you want to get  DNAm age predicted by a CTS clock.
#' @param CTSclocks A character vector of one or more clocks to apply.
#'                  (e.g., 'Neu-In', 'Hep', c('Neu-In', 'Neu-Sin', 'Brain')).
#' @param dataType Type of the samples ('bulk' or 'sorted').
#' @param CTF.m Optional cell type fraction matrix (rows: samples, columns: cell types).
#'              Required for 'Intrinsic' bulk clocks if tissue is not 'brain'.
#' @param tissue What tissue are your samples from ('brain' or 'otherTissue').
#' @param coreNum Number of cores for parallel computation in preprocessing.
#'
#' @param verbose Logical. If `TRUE` (default), the function will print
#'   progress messages to the console.
#'
#' @return A data.frame of predicted DNAm ages (samples x clocks).
#'
#' @import glmnet
#' @import parallel
#'
#' @export
#'
#' @references
#' Tong H, Guo X, Jacques M, Luo Q, Eynon N, Teschendorff AE.
#' Cell-type specific epigenetic clocks to quantify biological age at cell-type resolution.
#' \emph{Aging} 2024
#'
#' @examples
#' download_OmniAgeR_example("CTS_MurphyGSE88890")
#' load_OmniAgeR_example("CTS_MurphyGSE88890")
#' agePred_df <- CTS_Clocks(Murphy_beta_m, CTSclock = c('Neu-In','Neu-Sin'), dataType = 'bulk', CTF.m = NULL, tissue = 'brain',coreNum = 1)
#'
#' download_OmniAgeR_example("CTS_PaiGSE112179")
#' load_OmniAgeR_example("CTS_PaiGSE112179")
#' agePred_df <- CTS_Clocks(Pai_beta_m, CTSclock = c('Neu-In','Neu-Sin'), dataType = 'sorted', CTF.m = NULL, tissue = 'brain',coreNum = 1)
#'
#' download_OmniAgeR_example("CTS_ExampleData_Liver")
#' load_OmniAgeR_example("CTS_ExampleData_Liver")
#' agePred_df <- CTS_Clocks(Liver_beta_m, CTSclock = c('Hep','Liver'), dataType = 'bulk', CTF.m = NULL, tissue = 'otherTissue',coreNum = 1)
#'
#'



CTS_Clocks <- function(data.m,
                      CTSclocks = c('Neu-In'),
                      dataType = c('bulk', 'sorted'),
                      CTF.m = NULL,
                      tissue = c('brain', 'otherTissue'),
                      coreNum = NULL,
                      verbose = TRUE){

  if (verbose) {
    print(paste0("[CTS_Clocks] Starting CTS_Clocks calculation..."))
  }

  # --- 1. Initialization and Setup ---

  data("CTS_Clocks_Coef")
  # Set default core number for parallel computation in ProcessData
  if (is.null(coreNum)) {
    coreNum <- ceiling(parallel::detectCores() / 2)
  }

  # Initialize a list to store prediction results for each clock
  results.ls <- list()

  # Define clock types based on their preprocessing requirements
  # 'Intrinsic' clocks require data to be processed (residuals or Z-score)
  intrinsic_clock_models <- c('Neu-In', 'Glia-In', 'Brain')
  # 'Semi-intrinsic' linear clocks use raw data
  semi_intrinsic_linear_models <- c('Neu-Sin', 'Glia-Sin')
  # 'Semi-intrinsic' glmnet clocks use raw data and have a different prediction logic
  semi_intrinsic_glmnet_models <- c('Hep', 'Liver')


  # --- 2. Optimized Preprocessing (Run-Once Logic) ---

  # Check if any requested clock requires the 'Intrinsic' preprocessing step
  needs_processing <- any(intrinsic_clock_models %in% CTSclocks)

  # Lazily initialize the processed data object.
  # It will only be computed if 'needs_processing' is TRUE.
  processed_data.m <- NULL

  if (needs_processing) {
    if (dataType == 'bulk') {
      # This is the most computationally expensive step.
      # We calculate residuals only ONCE for all 'Intrinsic' bulk clocks.

      processed_data.m <- CTS_Clocks_ProcessData(data.m,
                                      dataType = 'bulk',
                                      tissue = tissue,
                                      CTF.m = CTF.m,
                                      coreNum = coreNum)
    } else if (dataType == 'sorted') {
      # For 'sorted' data, 'Intrinsic' clocks require Z-score standardization.

      processed_data.m <- CTS_Clocks_ProcessData(data.m,
                                      dataType = 'sorted',
                                      coreNum = coreNum)
    }
  }

  # --- 3. Iterate and Predict for Each Requested Clock ---

  for (clock_name in CTSclocks) {

    if (clock_name %in% intrinsic_clock_models) {
      # --- Path A: 'Intrinsic' Clocks (e.g., Neu-In) ---

      # Load the clock coefficients data file (e.g., Neu-InCoef.rda)
      #data(list = paste0("CTS_",clock_name, 'Coef'), envir = environment())
      clockCoef.df <- CTS_Clocks_Coef[[clock_name]]


      # Predict using the pre-computed processed data (residuals or Z-scored)
      intercept <- clockCoef.df$coef[1]
      weights_df <- clockCoef.df[-1, ]
      weights_vec <- weights_df$coef
      names(weights_vec) <- weights_df$probe
      coef_list <- list(intercept, weights_vec)

      results.ls[[clock_name]] <- calculateLinearPredictor(
        beta.m = processed_data.m,
        coef.lv = coef_list,
        clock.name = clock_name,
        verbose
      )


    } else if (clock_name %in% semi_intrinsic_linear_models) {
      # --- Path B: 'Semi-intrinsic' Linear Clocks (e.g., Neu-Sin) ---

      # Load the clock coefficients (e.g., Neu-SinCoef.rda)
      #data(list = paste0("CTS_",clock_name, 'Coef'), envir = environment())
      clockCoef.df <- CTS_Clocks_Coef[[clock_name]]

      intercept <- clockCoef.df$coef[1]
      weights_df <- clockCoef.df[-1, ]
      weights_vec <- weights_df$coef
      names(weights_vec) <- weights_df$probe
      coef_list <- list(intercept, weights_vec)

      # Predict using the new function and the raw, unprocessed data.m
      results.ls[[clock_name]] <- calculateLinearPredictor(
        beta.m = data.m,
        coef.lv = coef_list,
        clock.name = clock_name,
        verbose
      )

    } else if (clock_name %in% semi_intrinsic_glmnet_models) {
      # --- Path C: 'Semi-intrinsic' glmnet Clocks (e.g., Hep) ---

      # Load the .glm model object (e.g., HepClock.rda)
      #data(list = paste0("CTS_",clock_name, 'Clock'), envir = environment())

      # Evaluate the loaded object name (e.g., HepClock.glm)
      # We create a copy so that trimming does not affect the loaded object
      # in subsequent loops.
      current_clock.glm <- CTS_Clocks_Coef[[clock_name]]

      # --- Begin glmnet-specific probe matching ---
      # This logic finds the intersection of CpGs between the clock and input data

      # Get CpGs from the trained clock model
      ClockCpGs.v <- rownames(current_clock.glm$beta)
      # Get CpGs from the user's input data
      TestSetCpGs.v <- rownames(data.m)

      ClockCpGs_df <- as.data.frame(coef(current_clock.glm))
      non_zero_cpg <- rownames(ClockCpGs_df)[which(ClockCpGs_df[[1]] !=0 )]
      non_zero_cpg <- non_zero_cpg[-1]

      if (verbose) {
        print(paste0("[",clock_name,
                     "] Number of represented ", clock_name, " CpGs (max=",
                     length(non_zero_cpg), ")=", sum(TestSetCpGs.v %in% non_zero_cpg)))
      }

      # Find the indices of clock CpGs in the user's data
      idx <- match(ClockCpGs.v, TestSetCpGs.v)

      # Trim the clock's beta matrix to only include probes present in the user's data
      # We re-evaluate the original object to get a fresh beta matrix to trim
      #original_beta <- eval(parse(text = paste0(clock_name, 'Clock.glm')))$beta
      original_beta <- current_clock.glm$beta
      current_clock.glm$beta <- original_beta[!is.na(idx), , drop = FALSE]

      # Update the clock's dimensions
      current_clock.glm$dim <- c(sum(!is.na(idx)), 1)

      # Trim the user's data matrix to match the clock's probes and order
      data_trimmed.m <- data.m[na.omit(idx), , drop = FALSE]

      # Perform prediction using the trimmed clock and trimmed data
      results.ls[[clock_name]] <- as.vector(predict(current_clock.glm,
                                                           newx = t(data_trimmed.m)))
    } else {
      warning(paste("Clock name '", clock_name, "' not recognized. Skipping."))
    }
  }

  # --- 4. Format and Return Output ---

  if (length(CTSclocks) == 0) {
    warning("No clocks were specified or recognized.")
    return(NULL)
  }


  # return a data.frame (samples x clocks).
  return(as.data.frame(results.ls, row.names = colnames(data.m)))

}

#' @title Internal Preprocessing for CTS Clocks
#' @description
#' This is an internal helper function for \code{CTS_Clocks}. It preprocesses
#' the DNAm matrix for 'Intrinsic' clock models. It performs one of two
#' operations based on \code{dataType}:
#' \enumerate{
#'   \item \strong{'sorted'}: Applies Z-score standardization.
#'   \item \strong{'bulk'}: Calculates residuals by regressing out cell type
#'     fractions (CTFs), then applies Z-score standardization to the residuals.
#' }
#'
#' @param data.m A numeric matrix of DNAm (beta) values (CpGs x Samples).
#' @param dataType Character string: 'bulk' or 'sorted'.
#' @param tissue Character string: 'brain' or 'otherTissue'.
#' @param CTF.m An optional numeric matrix of cell type fractions
#'   (Samples x CellTypes).
#' @param coreNum Number of cores for parallel \code{mclapply}.
#' @param verbose Logical. If `TRUE` (default), the function will print
#'   progress messages to the console.
#'
#' @return A processed numeric matrix (CpGs x Samples), either Z-scored data
#'   or Z-scored residuals.
#'
#' @import HiBED
#' @importFrom parallel mclapply
#'
#' @noRd


CTS_Clocks_ProcessData <- function(data.m, dataType = c('bulk', 'sorted'),
                                   tissue = c('brain', 'otherTissue'), CTF.m = NULL,
                                   coreNum = coreNum,verbose = TRUE){

  if(dataType == 'sorted'){
    ## Normalize the data
    dataSD.v = unlist(parallel::mclapply(1:nrow(data.m),function(i) sd(data.m[i,]), mc.cores = coreNum))
    dataZ.m = (data.m - rowMeans(data.m))/dataSD.v
    return(dataZ.m)
  }else if(dataType == 'bulk'){
    if (is.null(CTF.m)){
      if(tissue == 'brain'){
        ## Estimate the cell type fractions
        estF.m = HiBED::HiBED_deconvolution(data.m, h=1)
        estF.m = estF.m/100
        colnames(estF.m) = c('EndoStrom', 'Glia', 'Neu')
        estF.m = estF.m[, c('Neu', 'Glia', 'EndoStrom')]
        CTF.m = as.matrix(estF.m)
      }else if(tissue == 'otherTissue'){
        stop("Cell type fraction matrix (CTF.m) is missing. If you don't have it, you are recommended to use R package EpiSCORE, EpiDISH or some other deconvolution algorithms to estimate the cell type fractions of your samples.")
      }
    }
    ## Adjust for cell type fractions and normalize the data
    if (verbose) {
      print("[CTS_Clocks] Processing the data may take some time. Don't worry.")
      }
    lm.o = lm(t(data.m) ~ CTF.m)
    res.m = t(lm.o$res)
    resSD.v = unlist(parallel::mclapply(1:nrow(res.m),function(i) sd(res.m[i,]), mc.cores = coreNum))
    dataZ.m = (res.m - rowMeans(res.m))/resSD.v
    return(dataZ.m)
  }
}

