#' @title Calculate Universal Pan-Mammalian Epigenetic Clocks
#' @description
#' Applies the three universal pan-mammalian epigenetic clocks (Clock 1, 2, and 3)
#' from Lu et al. (2023) to a given dataset of DNA methylation values.
#'
#' @details
#' This function encapsulates the complete logic from the original publication's
#' script. It merges sample metadata with AnAge data, processes methylation
#' values, calculates the three clocks by matrix multiplication with pre-trained
#' models, and performs the necessary inverse mathematical transformations
#' to report age in years.
#'
#' @param betaM A numeric matrix (Rows: CpGs, Cols: Samples).
#' @param speciesName A character string or vector specifying the Latin species
#' name(s) (e.g., "Homo sapiens").
#' @param anageData A data.frame containing the AnAge database information.
#'   Must include 'SpeciesLatinName', 'GestationTimeInYears',
#'   'averagedMaturity.yrs', and 'maxAge'.
#' @param minCoverage Numeric (0-1). Minimum required proportion of CpGs present.
#' Default is 0.
#' @param verbose Logical. Whether to print status messages.
#'
#' @return A data.frame containing the 'Sample', 'SpeciesLatinName',
#'   and the calculated ages: 'DNAmAgePanMammalianClock1','DNAmRelativeAge',
#'    'DNAmAgePanMammalianClock2', 'DNAmRelativeAdultAge' and
#'   'DNAmAgePanMammalianClock3'. Any additional columns
#'   from the input `sampleInfo` (like 'Age', 'Tissue') will also be returned.
#'
#' @export
#'
#' @references
#' Lu, A.T., Fei, Z., Haghani, A. et al.
#' Universal DNA methylation age across mammalian tissues.
#' \emph{Nat Aging.} 2023
#'
#'
#'
#' @examples
#' tursiopsExample <- loadOmniAgeRdata(
#'     "omniager_tursiops_example",
#'     verbose = FALSE
#' )
#'
#' ## This anage_data is from
#' ## https://github.com/shorvath/MammalianMethylationConsortium
#' anageData <- loadOmniAgeRdata(
#'     "omniager_anage_data",
#'     verbose = FALSE
#' )
#' # Run the calculation with progress messages
#' clockResults <- universalPanMammalianClocks(
#'     betaM = tursiopsExample$beta_m,
#'     speciesName = tursiopsExample$PhenoTypes$SpeciesLatinName,
#'     anageData = anageData
#' )
#'
universalPanMammalianClocks <- function(betaM,
                                        speciesName,
                                        anageData = NULL,
                                        minCoverage = 0,
                                        verbose = TRUE) {
    if (verbose) message("[UniversalPanMammalianClocks] Initializing calculation...")

    # --- 1. Load Internal Data ---
    panMammalianClockCoef <- loadOmniAgeRdata(
        "omniager_pan_mammalian_clock_coef",
        verbose = verbose
    )
    # --- 2. Data Preparation & Merging ---
    sampleInfo <- data.frame(
        Sample = colnames(betaM),
        SpeciesLatinName = speciesName,
        stringsAsFactors = FALSE
    )

    anageSubset <- anageData[, c(
        "SpeciesLatinName", "GestationTimeInYears",
        "averagedMaturity.yrs", "maxAge"
    )]
    info <- merge(sampleInfo, anageSubset, by = "SpeciesLatinName", all.x = TRUE)

    if (any(is.na(info$maxAge))) {
        missingSp <- unique(info$SpeciesLatinName[is.na(info$maxAge)])
        warning("[UniversalPanMammalianClocks] Missing AnAge data for species: ", paste(missingSp, collapse = ", "))
    }

    # The HighmaxAge required for processing Clock 2
    info$HighmaxAge <- info$maxAge * 1.3
    info$HighmaxAge[info$SpeciesLatinName %in% c("Homo sapiens", "Mus musculus")] <-
        info$maxAge[info$SpeciesLatinName %in% c("Homo sapiens", "Mus musculus")]

    # --- 3. Prediction (Step 1: Linear Scores) ---
    yNames <- c("Y.pred1", "Y.pred2", "Y.pred3")

    for (k in seq_along(yNames)) {
        modelDf <- panMammalianClockCoef[[k]]

        info[[yNames[k]]] <- .calLinearClock(
            betaM = betaM[, info$Sample, drop = FALSE],
            coefData = modelDf,
            clockLabel = paste0("PanMammalian_Clock", k),
            minCoverage = minCoverage,
            verbose = verbose
        )
    }

    # --- 4. Post-Processing (Step 2: Inverse Transformations) ---
    if (verbose) message("[UniversalPanMammalianClocks] Applying age transformations...")

    # Clock 1: Simple log-linear
    info$DNAmAgePanMammalianClock1 <- exp(info$Y.pred1) - 2

    # Clock 2: Relative Age
    info$DNAmRelativeAge <- exp(-exp(-info$Y.pred2))

    info$DNAmAgePanMammalianClock2 <- info$DNAmRelativeAge * (info$HighmaxAge + info$GestationTimeInYears) - info$GestationTimeInYears

    # Clock 3: Relative Adult Age

    a2 <- info$GestationTimeInYears / info$averagedMaturity.yrs
    m1 <- 5 * (a2^0.38)

    # Reversal conversion Clock 3
    y3 <- info$Y.pred3
    relAdultAge <- ifelse(y3 < 0, (exp(y3) - 1) * m1 + m1, y3 * m1 + m1)

    info$DNAmRelativeAdultAge <- relAdultAge
    info$DNAmAgePanMammalianClock3 <- relAdultAge * (info$averagedMaturity.yrs + info$GestationTimeInYears) - info$GestationTimeInYears

    # --- 5. Return Output ---
    finalResults <- info[, setdiff(names(info), c("Y.pred1", "Y.pred2", "Y.pred3", "HighmaxAge"))]

    return(finalResults)
}
