#' @title Get CpG Probes from Various DNAm Biomarkers
#' @description
#' Retrieves the lists of CpG probes used by the clocks specified
#' in the `clock_names` argument.
#'
#' @param clock_names A character vector specifying which clocks to retrieve probes for. Accepts all keywords from `EpiAge()` (e.g., "all","chronological", "biological") or individual clock names.
#'
#' @param PCClocks_RData
#' **Required if any PC Clock is requested.** This argument must be the
#' **pre-loaded data object** (a list) from
#' `load_OmniAgeR_data(object_name = "PCClocks_data")`.
#'
#' @param SystemsAge_RData
#' **Required if 'SystemsAge' is requested.** This argument must be the
#' **pre-loaded data object** (a list) from
#' `load_OmniAgeR_data(object_name = "SystemsAge_data")`.
#'
#' @return A named list where each element is a character vector
#' of CpG probes for a requested clock.
#'
#'
#' @export
#'
#' @examples
#' # Get probes for all clocks
#'PCClocks_RData <- load_OmniAgeR_data(object_name = "PCClocks_data")
#'SystemsAge_RData <- load_OmniAgeR_data(object_name = "SystemsAge_data")
#' probes <- getClockProbes(PCClocks_RData = PCClocks_RData,SystemsAge_RData=SystemsAge_RData)



#' @title Retrieve CpG Probes for Various DNAm Biomarkers
#'
#' @description
#' This function extracts the specific CpG probes used to construct various DNA methylation-based models.
#' These models include:
#' \itemize{
#'   \item **Epigenetic Clocks**: Chronological, biological, and mitotic age estimators (e.g., Horvath, GrimAge).
#'   \item **Surrogate Biomarkers**: DNAm predictors for plasma proteins (e.g., CRP, IL6), mitotic rates, telomere length and CHIP score.
#'   \item **Phenotypic Traits**: Estimators for lifestyle factors (e.g., BMI, Smoking).
#' }
#' Use this function to filter your methylation matrix (beta values) to only the necessary probes before calculating scores.
#'
#' @param clock_names A character vector specifying which models to retrieve probes for.
#'   Options include specific model names (e.g., `"GrimAge2"`, `"McCartney_Trait"`) or keywords for groups:
#'   \itemize{
#'     \item `"all"`: (Default) Retrieve probes for all available models.
#'     \item `"chronological"`: chronological age predictors.
#'     \item `"biological"`: Biological age predictors (e.g., PhenoAge, DunedinPACE).
#'     \item `"mitotic"`: Estimates of cell division rates (e.g., epiTOC).
#'     \item `"surrogate_biomarkers"`: Proxies for proteins, CHIP, and traits.
#'     \item `"celltype_specific"`: Neuron, Glia, or tissue-specific clocks.
#'   }
#'
#' @param clock_names A character vector specifying which models to retrieve probes for.
#'   Options include specific model names (e.g., `"GrimAge2"`, `"McCartney_Trait"`) or keywords for groups:
#'   \itemize{
#'     \item `"all"`: (Default) Retrieve probes for all available models.
#'     \item `"chronological"`: Chronological age predictors (First-generation clocks, e.g., Horvath, Hannum).
#'     \item `"biological"`: Biological age predictors focused on mortality/morbidity (e.g., PhenoAge, GrimAge, DunedinPACE).
#'     \item `"cellular_aging"`: A bundle including both mitotic clocks and DNAm telomere length estimators.
#'     \item `"mitotic"`: Estimates of cumulative cell division rates (e.g., epiTOC, HypoClock).
#'     \item `"dnamtl"`: DNA methylation-based estimators of telomere length (e.g., DNAmTL).
#'     \item `"causal"`: Clocks trained using causal inference methods to distinguish damage from adaptation (e.g., DamAge, AdaptAge).
#'     \item `"stochastic"`: Stochastic analogues of the Horvath, PhenoAge, and Zhang clocks (e.g., StocH, StocZ).
#'     \item `"cross_species"`: Clocks applicable across species or specific to mouse (e.g., UniversalPanMammalian).
#'     \item `"gestational"`: Predictors of gestational age (GA) (e.g., Bohlin_GA, Knight_GA, EPIC_GA).
#'     \item `"celltype_specific"`: Cell-type specific clocks (e.g., Neuron, Glia, Hep).
#'     \item `"surrogate_biomarkers"`: Proxies for plasma proteins, lifestyle traits, and clonal hematopoiesis (CHIP).
#'   }
#'
#'
#' @param PCClocks_RData A list object containing PC-Clock data (loaded via \code{load_OmniAgeR_data}).
#'   Required if requesting PC-based clocks (e.g., PCHorvath, PCGrimAge).
#' @param SystemsAge_RData A list object containing SystemsAge data (loaded via \code{load_OmniAgeR_data}).
#'   Required if requesting SystemsAge.
#'
#' @return A named \code{list} where each element is a \code{character vector} containing the CpG probe IDs required for a specific model.
#'
#' @export
#'
#' @examples
#' # 1. Get probes for a specific chronological aging predictor
#' chronological_probes <- getClockProbes(clock_names = "Horvath2013")
#'
#' # 2. Get probes for a specific biological aging clocks
#' bio_probes <- getClockProbes(clock_names = "GrimAge1")
#'
#' # 3. Get probes for all biomaker
#' PCClocks_RData <- load_OmniAgeR_data(object_name = "PCClocks_data")
#' SystemsAge_RData <- load_OmniAgeR_data(object_name = "SystemsAge_data")
#' all_probes <- getClockProbes(clock_names = "all",PCClocks_RData,SystemsAge_RData)
#'





getClockProbes <- function(clock_names = "all",
                           PCClocks_RData = NULL,
                           SystemsAge_RData = NULL) {

  # ===================================================================
  # 1. REPLICATE EpiAge Clock Definitions
  # ===================================================================
  mitotic_clocks <- c("epiTOC1", "epiTOC2", "epiTOC3", "stemTOCvitro", "stemTOC",
                      "RepliTali", "HypoClock", "EpiCMIT_Hyper", "EpiCMIT_Hypo")
  dnamtl_clocks <- c("DNAmTL","PCDNAmTL")
  cellular_aging_clocks <- c(mitotic_clocks, dnamtl_clocks)

  chronological_clocks <- c("Horvath2013", "Hannum","Lin","VidalBralo","ZhangClock",
                            "Horvath2018","Bernabeu_cAge","CorticalClock","PedBE",
                            "CentenarianClock","Retro_age","PCHorvath2013","PCHorvath2018","PCHannum")

  biological_clocks <- c("Zhang10","PhenoAge", "DunedinPACE", "GrimAge1","GrimAge2","PCPhenoAge","PCGrimAge1",
                         "DNAmFitAge",
                         "IC_Clock","SystemsAge")

  causal_clocks <- c("CausalAge","DamAge","AdaptAge")
  stochastic_clocks <- c("StocH", "StocZ", "StocP")
  #ensemble_clocks <- c("EnsembleAgeHumanMouse", "EnsembleAgeMouse_Static", "EnsembleAgeMouse_Dynamic")
  #cross_species_clocks <- c("UniversalPanMammalianClocks","PanMammalianBlood","PanMammalianSkin")
  cross_species_supprot <- c("EnsembleAgeHumanMouse", "EnsembleAgeMouse_Static", "EnsembleAgeMouse_Dynamic",
                             "UniversalPanMammalianClocks","PanMammalianBlood","PanMammalianSkin")

  cts_clocks <- c("Neu-In", "Glia-In", "Neu-Sin", "Glia-Sin", "Hep") # 'Brain' and 'Liver' are often bundles

  GA_clocks <- c("Bohlin_GA", "EPIC_GA", "Knight_GA", "Lee_GA", "Mayne_GA")

  EpiAgeClocks <-  unique(c(chronological_clocks, biological_clocks, causal_clocks, cellular_aging_clocks,
                            stochastic_clocks, cross_species_supprot, cts_clocks, GA_clocks))

  Surrogate_Biomarkers <- c("CompCRP","CompCHIP","CompIL6","CompEpiScores","McCartney_Trait")

  all_individual_clocks <- unique(c(EpiAgeClocks, Surrogate_Biomarkers))
  all_keywords <- c("all","cellular_aging", "mitotic", "dnamtl", "chronological", "biological", "causal",
                    "stochastic","cross_species","celltype_specific","gestational","surrogate_biomarkers")

  # ===================================================================
  # 2. REPLICATE EpiAge Clock Name Processing
  # ===================================================================
  # (This logic MUST be kept in sync with EpiAge.R)

  clocks_to_process <- c()
  if ("all" %in% clock_names) { clocks_to_process <- c(clocks_to_process, all_individual_clocks) }
  if ("cellular_aging" %in% clock_names) { clocks_to_process <- c(clocks_to_process, cellular_aging_clocks) }
  if ("mitotic" %in% clock_names) { clocks_to_process <- c(clocks_to_process, mitotic_clocks) }
  if ("dnamtl" %in% clock_names) { clocks_to_process <- c(clocks_to_process, dnamtl_clocks) }
  if ("chronological" %in% clock_names) { clocks_to_process <- c(clocks_to_process, chronological_clocks) }
  if ("biological" %in% clock_names) { clocks_to_process <- c(clocks_to_process, biological_clocks) }
  if ("causal" %in% clock_names) { clocks_to_process <- c(clocks_to_process, causal_clocks) }
  if ("stochastic" %in% clock_names) { clocks_to_process <- c(clocks_to_process, stochastic_clocks) }
  if ("celltype_specific" %in% clock_names) { clocks_to_process <- c(clocks_to_process, cts_clocks) }
  if ("gestational" %in% clock_names) { clocks_to_process <- c(clocks_to_process, GA_clocks) }
  if ("surrogate_biomarkers" %in% clock_names) { clocks_to_process <- c(clocks_to_process, Surrogate_Biomarkers) }

  specific_clocks <- setdiff(clock_names, all_keywords)
  clocks_to_process <- unique(c(clocks_to_process, specific_clocks))

  invalid_names <- setdiff(clocks_to_process, all_individual_clocks)
  if (length(invalid_names) > 0) {
    warning("The following clock names are not valid and will be ignored: ", paste(invalid_names, collapse = ", "))
  }
  clocks_to_run <- intersect(clocks_to_process, all_individual_clocks)

  if (length(clocks_to_run) == 0) {
    stop("No valid clocks were selected.")
  }

  # ===================================================================
  # 3. The "Probe Manifest"
  # ===================================================================

  # ===================================================================
  # Helper Function: Load Internal Package Data safely
  # ===================================================================
  probes.ls <- list()

  load_pkg_data <- function(data_name) {
    e <- new.env()
    pkg_name <- "OmniAgeR"

    # Step 1: Try to load the data file.
    try(
      data(list = data_name, package = pkg_name, envir = e),
      silent = TRUE
    )
    # Step 2: Get the names of all objects that were *actually* loaded into 'e'.
    loaded_object_names <- ls(e)

    # Step 3: Check if the loading failed.
    if (length(loaded_object_names) == 0) {
      warning(paste("Data file", shQuote(data_name),
                    "not found in package", shQuote(pkg_name)))
      return(NULL)
    }
    # Step 4: Retrieve all loaded objects from 'e'.
    object_list <- mget(loaded_object_names, envir = e)
    # Step 5: Return the complete list of objects.
    return(object_list)
  }


  # ===================================================================
  # PART A: Complex/Bundled Clocks (PC Clocks & SystemsAge)
  # ===================================================================

  # PCClocks Bundle
  pc_clock_group <- c("PCHorvath2013", "PCHorvath2018", "PCHannum", "PCPhenoAge", "PCGrimAge1", "PCDNAmTL")
  requested_pc_clocks <- intersect(clocks_to_run, pc_clock_group)

  if (length(requested_pc_clocks) > 0) {
    if (is.null(PCClocks_RData)) {
      warning("Cannot get PC Clock probes. `PCClocks_RData` must be provided.")
    } else if (is.character(PCClocks_RData)) {
      PCClocks_RData <- load_OmniAgeR_data(object_name = "PCClocks_data", path = PCClocks_RData)
    }
    if (rlang::hash(PCClocks_RData) != "46386ec4be2b2a5239cf67b242d7dc24") {
      stop(" The downloaded PCClocks data is corrupted or the wrong data (e.g., SystemsAge) was passed. See `?download_methylCIPHER()`.")
    }
    # All PC clocks use one combined set of probes
    probes.ls$PCClocks_ALL <- PCClocks_RData$CpGs

    # Remove from individual processing
    clocks_to_run <- setdiff(clocks_to_run, pc_clock_group)
  }

  # SystemsAge Bundle
  if ("SystemsAge" %in% clocks_to_run) {
    if (is.null(SystemsAge_RData)) {
      warning("Cannot get SystemsAge probes. `SystemsAge_RData` must be provided.")
    } else if (is.character(SystemsAge_RData)) {
      SystemsAge_RData <- load_OmniAgeR_data(object_name = "SystemsAge_data", path = SystemsAge_RData)
    }
    if (rlang::hash(SystemsAge_RData) != "d984914ff6aa17d8a6047fed5f9f6e4d") {
      stop(" The downloaded SystemsAge_RData data is corrupted or the wrong data (e.g., SystemsAge) was passed. See `?download_methylCIPHER()`.")
    }
    # All PC clocks use one combined set of probes
    probes.ls$SystemsAge <- SystemsAge_RData$CpGs

    # Remove from individual processing
    clocks_to_run <- setdiff(clocks_to_run, "SystemsAge")
  }

  # ===================================================================
  # PART B: Individual Model Retrieval
  # ===================================================================

  for (clock in clocks_to_run) {

    # Use 'try' to prevent one error from stopping the whole loop
    tryCatch({

      if (clock == "Horvath2013") {
        d <- load_pkg_data("Horvath2013CpG")
        d <- as.data.frame(d)
        probes.ls$Horvath2013 <- as.vector(d[2:nrow(d), 1])
      } else if (clock == "Hannum") {
        d <- load_pkg_data("Hannum")
        probes.ls$Hannum <- d$Hannum_coef$CpGmarker
      } else if (clock == "Lin") {
        d <- load_pkg_data("Lin")
        probes.ls$Lin <- d$Lin_coef$CpGmarker[-length(d[[1]][[1]])]
      } else if (clock == "VidalBralo") {
        d <- load_pkg_data("VidalBraloCoef")
        probes.ls$VidalBralo <-  d$VidalBralo_Coef$probe[-1]
      } else if (clock == "ZhangClock") {
        d <- load_pkg_data("ZhangCpG")
        probes.ls$ZhangClock <- d$Zhang_coef_df$CpG[-1]
      } else if (clock == "Horvath2018") {
        d <- load_pkg_data("Horvath_SkinBlood_CpG")
        probes.ls$Horvath2018 <- d$Horvath2_coef$CpGmarker[-1]
      } else if (clock == "Bernabeu_cAge") {
        d <- load_pkg_data("Bernabeu_cAge_Coef")
        probes.ls$Bernabeu_cAge <- rownames(d$Bernabeu_cAge_means)
      } else if (clock == "CorticalClock") {
        d <- load_pkg_data("CorticalClockCoef")
        probes.ls$CorticalClock <-  d$CorticalClock_coef$probe[-1]
      } else if (clock == "PedBE") {
        d <- load_pkg_data("PedBECoef")
        probes.ls$PedBE <- d$PedBECoef$Probe[-1]
      } else if (clock == "CentenarianClock") {
        d <- load_pkg_data("CentenarianENCoef")
        probes.ls$CentenarianClock_40 <- d$CentenarianENCoef$ENCen40$prob[-1]
        probes.ls$CentenarianClock_100 <- d$CentenarianENCoef$ENCen100$prob[-1]
      } else if (clock == "Retro_age") {
        d <- load_pkg_data("Retro_age_Coef")
        probes.ls$Retro_age_V1 <- d$Retro_age_Coef$V1$prob[-1]
        probes.ls$Retro_age_V2 <- d$Retro_age_Coef$V2$prob[-1]
      } else if (clock == "Zhang10") {
        d <- load_pkg_data("Zhang10Coef")
        probes.ls$Zhang10 <- d$Zhang10Coef$Probe
      } else if (clock == "PhenoAge") {
        d <- load_pkg_data("PhenoAgeCpG");
        probes.ls$PhenoAge <- d$PhenoAge_coef_df$CpGmarker[-1]
      } else if (clock == "DunedinPACE") {
        d <- load_pkg_data("DunedinPACECpG")
        probes.ls$DunedinPACE <- d$mPACE_Models_list$model_probes$DunedinPACE
      } else if (clock == "GrimAge1") {
        d <- load_pkg_data("GrimAge1CpG")
        probes.ls$GrimAge1 <- unique(setdiff(d$grimage1$dnam.all$var,c("Intercept", "Age")))
      } else if (clock == "GrimAge2") {
        d <- load_pkg_data("GrimAge2CpG")
        probes.ls$GrimAge2 <- unique(setdiff(d$grimage2$dnam.all$var,c("Intercept", "Age")))
      } else if (clock == "DNAmFitAge") { ################# 需要再判断一下
        d <- load_pkg_data("DNAmFitAgeCoef");
        d_temp <- d$DNAmFitnessModels$AllCpGs
        probes.ls$DNAmFitAge <- d$DNAmFitnessModels$AllCpGs ## 需不需要grimage1呢
      } else if (clock == "IC_Clock") {
        d <- load_pkg_data("IC_Clock_Coef")
        probes.ls$IC_Clock <- d$IC_Clock_Coef$probe[-1]
      } else if (clock %in% c("CausalAge", "DamAge", "AdaptAge")) {
        d <- load_pkg_data("CausalCpG")
        probes.ls$CausalAge <- d$causalClock.l$causal$term[-1]
        probes.ls$DamAge <- d$causalClock.l$dam$term[-1]
        probes.ls$AdaptAge <- d$causalClock.l$adapt$term[-1]
      } else if (clock %in% c("StocH", "StocZ", "StocP")) {
        d <- load_pkg_data("glmStocALL")
        ### StocH
        temp_coef <- coef(d$glmStocALL.lo$H)
        probes.ls$StocH <- names(temp_coef[-1,ncol(temp_coef)])
        ### StocZ
        temp_coef <- coef(d$glmStocALL.lo$Z)
        probes.ls$StocZ <- names(temp_coef[-1,ncol(temp_coef)])
        ### StocP
        temp_coef <- coef(d$glmStocALL.lo$P)
        probes.ls$StocP <- names(temp_coef[-1,ncol(temp_coef)])

      } else if (clock %in% c("EnsembleAge")) {
        d <- load_pkg_data("EnsembleAgeCoef")
        d_tmep <- d$EnsembleAgeCoef$Dynamic
        probes.ls$EnsembleAgeMouse_Dynamic <- setdiff(unique(unlist(lapply(d_tmep, function(x) x$CGid))),"Intercept")
        probes.ls$EnsembleAgeMouse_Static <- d$EnsembleAgeCoef$Static$Static$CGid[-1]
        probes.ls$EnsembleAgeMouse_Static_Top <- d$EnsembleAgeCoef$Static$Static_Top$CGid[-1]
        probes.ls$EnsembleAgeHumanMouse <- d$EnsembleAgeCoef$HumanMouse$HumanMouse$CGid[-1]
        rm(d_tmep)
      } else if (clock == "epiTOC1") {
        d <- load_pkg_data("dataETOC3");
        probes.ls$epiTOC1 <- d$dataETOC3.l$epiTOC
      } else if (clock == "epiTOC2") {
        d <- load_pkg_data("dataETOC3")
        probes.ls$epiTOC2 <- rownames(d$dataETOC3.l$epiTOC2)
      } else if (clock == "epiTOC3") {
        d <- load_pkg_data("dataETOC3");
        probes.ls$epiTOC3 <- rownames(d$dataETOC3.l$epiTOC3)
      } else if (clock == "RepliTali") {
        d <- load_pkg_data("Replitali")
        probes.ls$RepliTali <- d$replitali.cpg.v[-1]
      } else if (clock == "HypoClock") {
        d <- load_pkg_data("dataETOC3")
        probes.ls$HypoClock <- d$dataETOC3.l$`cm-solo-WCGW`
      } else if (clock == "EpiCMIT_Hyper") {
        d <- load_pkg_data("EpiCMITcpgs")
        probes.ls$EpiCMIT_Hyper <- d$epiCMIT.df$Name[d$epiCMIT.df$epiCMIT.class == "epiCMIT.hyper"]
      } else if (clock == "EpiCMIT_Hypo") {
        d <- load_pkg_data("EpiCMITcpgs")
        probes.ls$EpiCMIT_Hypo <- d$epiCMIT.df$Name[d$epiCMIT.df$epiCMIT.class == "epiCMIT.hypo"]
      } else if (clock == "DNAmTL") {
        d <- load_pkg_data("DNAmTLCoef")
        probes.ls$DNAmTL <- d$DNAmTLCoef$probe[-1]
      } else if (clock %in% c("UniversalPanMammalianClocks")) {
        d <- load_pkg_data("PanMammalianClockCoef")
        probes.ls$UniversalPanMammalianClock1 <- d$PanMammalianClockCoef$PanMammalianClock1$prob[-1]
        probes.ls$UniversalPanMammalianClock2 <- d$PanMammalianClockCoef$PanMammalianClock2$prob[-1]
        probes.ls$UniversalPanMammalianClock3 <- d$PanMammalianClockCoef$PanMammalianClock3$prob[-1]
      }else if (clock %in% c("PanMammalianBlood")) { ###稍后再用
        d <- load_pkg_data("PanMammalianBloodCoef")
        probes.ls$PanMammalianBlood2 <- d$PanMammalianBloodCoef$PanMammalianBlood2$probe[-1]
        probes.ls$PanMammalianBlood3 <-  d$PanMammalianBloodCoef$PanMammalianBlood3$probe[-1]
      }else if (clock %in% c("PanMammalianSkin")) { ###稍后再用
        d <- load_pkg_data("PanMammalianSkinCoef")
        probes.ls$PanMammalianSkin2 <- d$PanMammalianSkinCoef$PanMammalianSkin2$probe[-1]
        probes.ls$PanMammalianSkin3 <- d$PanMammalianSkinCoef$PanMammalianSkin3$probe[-1]
      } else if (clock == "Neu-In") {
        d <- load_pkg_data("CTS_Clocks_Coef");  # d <- load_pkg_data("CTS_Neu-InCoef");
        probes.ls$Neu_In <- d$CTS_Clocks_Coef$`Neu-In`$probe[-1]
      } else if (clock == "Glia-In") {
        d <- load_pkg_data("CTS_Clocks_Coef");
        probes.ls$Glia_In <- d$CTS_Clocks_Coef$`Glia-In`$probe[-1]
      } else if (clock == "Neu-Sin") {
        d <- load_pkg_data("CTS_Clocks_Coef");
        probes.ls$Neu_Sin <- d$CTS_Clocks_Coef$`Neu-Sin`$probe[-1]
      } else if (clock == "Glia-Sin") {
        d <- load_pkg_data("CTS_Clocks_Coef");
        probes.ls$Glia_Sin <- d$CTS_Clocks_Coef$`Glia-Sin`$probe[-1]
      } else if (clock == "Hep") {
        d <- load_pkg_data("CTS_Clocks_Coef");
        d_tmep <- as.matrix(d$CTS_Clocks_Coef$Hep$beta)
        probes.ls$Hep <-  rownames(d_tmep)[d_tmep[,1] !=0]  # Coefs from glm
        rm(d_tmep)
      }else if (clock == "Bohlin_GA") {    ## add the Gestational Age
        d <- load_pkg_data("BohlinGACoef");
        probes.ls$Bohlin_GA <- d$BohlinGACoef$probe[-1]
      }else if (clock == "EPIC_GA") {
        d <- load_pkg_data("EPICGACoef");
        probes.ls$EPIC_GA <- d$EPICGACoef$probe[-1]
      }else if (clock == "Lee_GA") {
        d <- load_pkg_data("LeeGACoef");
        probes.ls$LeeControl <- d$LeeGACoef$LeeControl$probe[-1]
        probes.ls$LeeRobust <-  d$LeeGACoef$LeeRobust$probe[-1]
        probes.ls$LeeRefinedRobust <-  d$LeeGACoef$LeeRefinedRobust$probe[-1]
      }else if (clock == "Knight_GA") {
        d <- load_pkg_data("KnightCoef");
        probes.ls$Knight_GA <- d$KnightCoef$probe[-1]
      }else if (clock == "Mayne_GA") {
        d <- load_pkg_data("MayneGACoef");
        probes.ls$Mayne_GA <- d$MayneGACoef$probe[-1]
      }else if (clock == "CompCRP") {    ## add Surrogate Biomarkers & Scores
        d <- load_pkg_data("crpCpG");
        probes.ls$CRP <- names(d$crpCpG.lv$crp)
        probes.ls$intCRP <- names(d$crpCpG.lv$intCRP)
      }else if (clock == "CompCHIP") {
        d <- load_pkg_data("chipCpG")
        probes.ls$AnyCHIP <- unique(names(d$chipCpG.lv$AnyCHIP))
        probes.ls$DNMT3A <- unique(names(d$chipCpG.lv$DNMT3A))
        probes.ls$TET2 <- unique(names(d$chipCpG.lv$TET2))
        probes.ls$ASXL1 <- unique(names(d$chipCpG.lv$ASXL1))
      }else if (clock == "CompIL6") {
        d <- load_pkg_data("IL6Coef")
        probes.ls$IL6 <- d$IL6Coef$prob
      }else if (clock == "CompEpiScores") {
        d <- load_pkg_data("EpiScoresCoef")
        # split
        d_temp <- split(d$EpiScoresCoef$CpG_Site, d$EpiScoresCoef$Predictor)
        names(d_temp) <- paste0("EpiScores_", names(d_temp))
        probes.ls <- c(probes.ls, d_temp)
        rm(d_temp)
      }else if (clock == "McCartney_Trait") {
        d <- load_pkg_data("McCartneyTraitCoef");
        d_temp <- lapply(d$McCartneyTraitCoef, function(df) {
          return(df$CpG)
        })
        names(d_temp) <- paste0("McCartney_", names(d_temp))
        probes.ls <- c(probes.ls, d_temp)
      }

    }, error = function(e) {
      warning(paste("Failed to get probes for clock:", clock, ". Error:", e$message))
    })
  }

  return(probes.ls)
}
