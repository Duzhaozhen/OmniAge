#' @title List all available epigenetic clocks by category
#'
#' @description
#' Returns a named list containing the names of all individual clocks
#' implemented in the `EpiMarker` function, grouped by their functional category.
#' This is useful for knowing which strings are valid for the `clockNames`
#' argument.
#'
#' @return A named list where each name is a clock category (e.g., "mitotic",
#'   "chronological") and each value is a character vector of clock names in
#'   that category.
#'
#' @examples
#' clockCategories <- listEpiMarker()
#' @export

listEpiMarker <- function() {
    mitotic_clocks <- c(
        "epiTOC1", "epiTOC2", "epiTOC3", "stemTOCvitro", "stemTOC",
        "RepliTali", "HypoClock", "EpiCMIT_Hyper", "EpiCMIT_Hypo"
    )
    dnamtl_clocks <- c("DNAmTL", "PCDNAmTL")
    # cellular_aging_clocks <- c(mitotic_clocks, dnamtl_clocks)
    chronological_clocks <- c(
        "Horvath2013", "Hannum", "Lin", "VidalBralo", "ZhangClock",
        "Horvath2018", "Bernabeu_cAge", "CorticalClock", "PedBE",
        "CentenarianClock", "Retro_age", "ABEC", "eABEC", "cABEC",
        "PipekElasticNet", "PipekFilteredh", "PipekRetrainedh", "WuClock",
        "Weidner", "IntrinClock", "Garagnani",
        "PCHorvath2013", "PCHorvath2018", "PCHannum"
    )
    biological_clocks <- c(
        "Zhang10", "PhenoAge", "DunedinPACE", "GrimAge1", "GrimAge2", "PCPhenoAge",
        "PCGrimAge1", "DNAmFitAge", "IC_Clock", "SystemsAge"
    )
    causal_clocks <- c("CausalAge", "DamAge", "AdaptAge")
    stochastic_clocks <- c("StocH", "StocZ", "StocP")
    cts_clock <- c("Neu-In", "Neu-Sin", "Glia-In", "Glia-Sin", "Hep")

    ga_clocks <- c("BohlinGA", "EPICGA", "KnightGA", "LeeGA", "MayneGA")
    
    ensembleAge <- c("EnsembleAge_HumanMouse", 
                     "EnsembleAge_Static", 
                     "EnsembleAge_Dynamic")
    
    crossSpeciesSupport <- c(
        ensembleAge,
        "UniversalPanMammalianClocks",
        "PanMammalianBlood", "PanMammalianSkin",
        "PetkovichMouse", "ThompsonMouse", "MeerMouse"
    )
    
    cts_clock <- c("Neu-In", "Neu-Sin", "Glia-In", "Glia-Sin", "Hep")
    surrogateBiomarkers <- c("CRP", "CHIP", "IL6", "EpiScores")
    traitPred <- c("McCartneyTrait")
    diseaseRisk <- c("SmokeIndex", "HepatoXuRisk")

    return(list(
        mitotic = mitotic_clocks,
        dnamtl = dnamtl_clocks,
        chronological = chronological_clocks,
        biological = biological_clocks,
        causal = causal_clocks,
        cellTypeSpecific = cts_clock,
        stochastic = stochastic_clocks,
        gestationalAge = ga_clocks,
        surrogateBiomarkers = surrogateBiomarkers,
        traitPred = traitPred,
        diseaseRisk = diseaseRisk,
        crossSpeciesSupport = crossSpeciesSupport
    ))
}


#' @title Estimate epigenetic marker
#'
#' @description
#' A unified wrapper to calculate various epigenetic clocks from DNAm data.
#' This wrapper function provides a unified interface to run multiple clocks,
#' such as "Horvath2013" or "PhenoAge", either individually or as predefined
#' groups (e.g., "cellular_aging").
#'
#' @param betaM DNAm beta value matrix with rows labeling Illumina 450k/EPIC
#' CpGs and columns labeling samples.
#' @param clockNames A character vector specifying which clocks to calculate.
#'   This can include individual clock names (e.g., "Horvath2013") or special
#'   keywords:
#'   - `"all"`: (Default) Calculates all available clocks:
#'     `cellular_aging` + `chronological` + `biological` + `causal`+
#'     `stochastic` + `ensemble` + `cross_species` + `cellTypeSpecific`.
#'   - `"cellular_aging"`: Calculates all cellular aging clocks
#'   (both `mitotic` and `dnamtl`).
#'   - `"mitotic"`: Calculates only the mitotic clocks (e.g., `epiTOC1`).
#'   - `"dnamtl"`: Calculates only the DNAm Telomere Length clock (`DNAmTL`,
#'   `PCDNAmTL`).
#'   - `"mitotic"`: Calculates all defined mitotic clocks (e.g., `epiTOC1`).
#'   - `"chronological"`: Calculates chronological age clocks (e.g.,
#'   `Horvath2013`, `Hannum`, `ZhangClock`).
#'   - `"biological"`: Calculates biological age clocks (e.g., `PhenoAge`,
#'   `DunedinPACE`, `GrimAge1`, `GrimAge2`).
#'   - `"causal"`: Calculates causality-related clocks (`CausalAge`).
#'   - `"stochastic"`: Calculates the stochastic clocks (StocH, StocP, StocZ)
#'   as a single bundle.
#'   - `"cellTypeSpecific"`: Calculates all cell-type-specific clocks
#'   (`Neu-In`, `Glia-Sin`, etc.) as a single bundle.
#'
#'   - `"gestationalAge"`: Calculates all gestational age clocks
#'   (`BohlinGA`, `KnightGA`, etc.) as a single bundle.
#'
#'   - `"surrogateBiomarkers"`: Calculates all surrogate biomarkers
#'   (`CRP`, `CHIP`, etc.) as a single bundle.
#'   - `"traitPred"`: Calculates trait predictor (`McCartneyTrait`, etc.)
#'   as a single bundle.
#'   - `"diseaseRisk"`: Calculates disease risk predictions (`SmokeIndex`, etc.)
#'    as a single bundle.
#'   - `"crossSpeciesSupport"`: Calculates all cross species clocks
#'   (`UniversalPanMammalianClocks`, etc.) as a single bundle.
#'
#'   Multiple keywords and individual names can be combined (e.g.,
#'   `c("Horvath2013", "GrimAge2")`).
#'   Use `listEpiMarker()` to see a categorized list of all available clocks.
#'
#' @param chronAge Optional numeric vector of chronological ages.
#'   \bold{Note:} This is mandatory for \code{GrimAge}, \code{DNAmFitAge}, and \code{epiTOC} versions 2 & 3.
#' @param sexVec Optional character vector ('Male', 'Female').
#'   \bold{Note:} Required for \code{GrimAge} and \code{DNAmFitAge}.
#' 
#' @param minCoverage Numeric (0-1). Minimum proportion of required CpGs
#' present. Default 0.
#' @param verbose Logical. Whether to print progress messages.
#' @param ... Additional context-specific arguments:
#'   \describe{
#'     \item{\code{pcClockData}}{External data object for PC clocks.}
#'     \item{\code{systemsAgeData}}{External data object for SystemsAge calculations.}
#'     \item{\code{ctsDataType}}{Character. Required for \code{"cellTypeSpecific"} clocks. 
#'     Either \code{"bulk"} (default) or \code{"sorted"}.}
#'     \item{\code{ctsTissue}}{Character. Required for \code{"cellTypeSpecific"} clocks. 
#'     Either \code{"brain"} (allows auto-deconvolution) or \code{"otherTissue"}.}
#'     \item{\code{ctfM}}{Numeric matrix. Optional cell type fraction matrix (samples x cell types). 
#'     Required for \code{"Intrinsic"} CTS clocks (e.g., \code{"Neu-In"}) when \code{ctsDataType = "bulk"} and 
#'     \code{ctsTissue = "otherTissue"}.}     
#'     \item{\code{speciesName}}{Character string (e.g., "Mus musculus") for cross-species clocks.}
#'     \item{\code{anageData}}{A data.frame containing mammalian lifespan information. 
#'       Required for panMammalian clocks. It is recommended to use the internal database: 
#'       \code{anageData = loadOmniAgeRdata("omniager_anage_data")}. 
#'       If providing custom data, it must include columns: \code{SpeciesLatinName}, 
#'       \code{GestationTimeInYears}, \code{averagedMaturity.yrs}, and \code{maxAge}.}
#'   }

#' 
#' @details
#' The EpiMarker function implements a wide range of epigenetic clocks.
#' This includes cellular aging clocks (e.g., `epiTOC2`, `DNAmTL`),
#' chronological age clocks (e.g., `Horvath2013`, `Hannum`, `ZhangClock`),
#' biological age clocks (e.g., `PhenoAge`, `DunedinPACE`, `GrimAge2`),
#' causality-related clocks ( `CausalAge`, `DamAge`, `AdaptAge`),
#' stochastic clocks(`StocH`, `StocP`, `StocZ`), cell-type specific clocks
#' (e.g., `Neu-In`, `Glia-In`), gestational age clocks(e.g., `BohlinGA`,
#' `KnightGA`), surrogate biomarkers(e.g., `CRP`, `IL6`), trait prediction
#' (`McCartneyTrait`), disease risk prediction(e.g., `HepatoXuRisk`),
#' and cross species clocks(e.g., `UniversalPanMammalianClocks`).
#'
#'
#' === Cellular Aging Clocks ===
#' This category includes clocks measuring cell proliferation (mitotic clocks)
#' and telomere length.
#' **Mitotic Clocks**
#' * `epiTOC1`: (Yang et al. 2016) Average beta value of 385 promoter CpGs.
#' * `epiTOC2`: (Teschendorff et al. 2020) Estimates stem cell divisions using
#' a dynamic model.
#' * `epiTOC3`: Estimates stem cell divisions based on unmethylated population
#' doubling associated CpGs.
#' * `stemTOCvitro`: (Zhu et al. 2024 The 0.95 upper quantile of the 629
#' stemTOCvitro CpGs. The stemTOCvitro CpGs are promoters CpGs that are
#' unmethylated in fetal tissue-types and undergo DNA hypermethylation with
#' increased population-doublings.
#' * `stemTOC`: (Zhu et al. 2024) The 0.95 upper quantile of the 371 stemTOC
#' CpGs. Compared to stemTOCvitro CpGs, the stemTOC CpGs are filtered for
#' significant DNA hypermethylation with chronological age in large in-vivo
#' datasets.
#' * `RepliTali`: (Endicott et al. 2022) Based on 87 population doubling
#' associated hypomethylated CpGs.
#' * `HypoClock`: (Teschendorff et al. 2020) Based on hypomethylation at 678
#' solo-WCGW sites.
#' * `EpiCMIT_hyper`: (Duran-Ferrer et al. 2020) Average beta value of 184
#' age-associated hypermethylated CpGs.
#' * `EpiCMIT_hypo`: (Duran-Ferrer et al. 2020) Average beta value of 1164
#' age-associated hypomethylated CpGs.
#'
#' **DNAm Telomere Length (TL) Clocks**
#' * `DNAmTL`: (Lu AT et al. 2019) Calculating the Leukocyte telomere length.
#' * `PCDNAmTL`: (Lu AT et al. 2019) (Higgins-Chen et al. 2022)
#' Computationally-bolstered versions of the original clocks. These "PC clocks"
#' are retrained using principal components (PCs) derived from a large
#' CpG set to minimize technical noise and improve reliability
#'
#' === Chronological Clocks ===
#' * `Horvath2013`: (Horvath 2013) The original pan-tissue clock from 353 CpGs.
#' * `Hannum`: (Hannum et al. 2013) Blood-specific clock from 71 CpGs.
#' * `Lin`: (Lin et al. 2016) Clock based on 99 CpGs.
#' * `VidalBralo`: (Vidal-Bralo et al. 2016) Clock based on 8 CpGs.
#' * `ZhangClock`: (Zhang et al. 2019) Improved-precision clock from 514 CpGs.
#' * `Horvath2018`: (Horvath et al. 2018) Skin & blood clock.
#' * `Bernabeu_cAge`: (Bernabeu et al. 2022) Clock based on 3225 CpGs.
#' * `CorticalClock`: (Shireby et al. 2020) Crain cortical clock based on 347
#' CpGs.
#' * `PedBE`: (McEwen et al. 2019) Pediatric buccal epithelial clock in
#' Children.
#' * `CentenarianClock`: (Eric Dec et al. 2023) Centenarian Epigenetic Clocks.
#' * `Retro_age`: (Ndhlovu et al. 2024) Retroelement-based clock developed on
#' the EPIC v1/v2 arrays.
#' * `ABEC`: (Lee et al. 2020) Adult Blood-based EPIC Clock.
#' * `eABEC`: (Lee et al. 2020) Extended Adult Blood-based EPIC Clock.
#' * `cABEC`: (Lee et al. 2020) Common Adult Blood-based EPIC Clock.
#' * `PipekElasticNet`: (Pipek et al. 2023) Pipek's Multi-tissue Elastic Net
#' Epigenetic Clock (239 CpGs).
#' * `PipekFilteredh`: (Pipek et al. 2023) Pipek's Filtered Horvath Epigenetic
#' Clock (272 CpGs).
#' * `PipekRetrainedh`: (Pipek et al. 2023) Pipek's Retrained Horvath Epigenetic
#'  Clock (308 CpGs).
#' * `WuClock`: (Wu et al. 2019) Epigenetic Clock for Pediatric Age Estimation.
#' * `Weidner`: (Weidner et al. 2014) Calculate Weidner Epigenetic Age (3 CpGs).
#' * `IntrinClock`: (Tomusiak et al. 2024) Calculates the intrinsic cellular age
#' (380 CpGs).
#' * `Garagnani`: (Garagnani et al. 2012) The Garagnani ELOVL2-based Epigenetic
#' Age Score(1 CpG).
#' * `PCHorvath2013`, `PCHorvath2018`, `PCHannum`: (Higgins-Chen et al. 2022)
#' Computationally-bolstered versions of the original clocks. These "PC clocks"
#' are retrained using principal components (PCs) derived from a large
#' CpG set to minimize technical noise and improve reliability
#'
#'
#' === Biological Clocks ===
#' * `Zhang10`: (Zhang et al. 2017) 10-CpG clock associated with mortality.
#' * `PhenoAge`: (Levine et al. 2018) Predicts phenotypic age from 513 CpGs.
#' * `DunedinPACE`: (Belsky et al. 2022) Quantifies the pace of
#' biological aging.
#' * `GrimAge1`: The original GrimAge clock (Lu et al. 2019).
#' * `GrimAge2`: (Lu et al. 2022) Updated composite biomarker of mortality risk.
#' * `PCPhenoAge`, `PCGrimAge1`: (Higgins-Chen et al. 2022)
#'   Computationally-bolstered PC versions of the PhenoAge and GrimAge1
#'   clocks, retrained on principal components for enhanced reliability.
#' * `DNAmFitAge`: (McGreevy et al. 2023) Biological age indicator incorporating
#' DNAmGrimAge and 3 DNAm-based physical fitness markers.
#' * `IC_Clock`: (Fuentealba et al. 2025) The Intrinsic Capacity (IC) Clock
#' based on 91 CpGs.
#' * `SystemsAge`: (Sehgal et al. 2025) The Systems Age and 11
#' System-Specific Scores.
#'
#' === Causal Clocks ===
#' * `CausalAge`, `DamAge`, `AdaptAge`: (Ying et al. 2024) Three clocks
#' derived from a causality-enriched model. They dissect aging into distinct
#' components: 'Causal' (CausalAge), 'Damage' (DamAge),
#' and 'Adaptation' (AdaptAge).
#'
#' === Stochastic Clocks ===
#' * `StocH`, `StocP`, `StocZ`: (Tong et al. 2024) Stochastic analogues of the
#'   Horvath, PhenoAge, and Zhang clocks, trained on artificial cohorts to
#'   quantify the stochastic component of epigenetic aging.
#'
#' === Cell-Type Specific Clocks ===
#' * `Neu-In`, `Glia-In`, `Brain`: (Tong et al. 2024) Intrinsic clocks for
#' neurons, glia, and whole brain.
#' * `Neu-Sin`, `Glia-Sin`, `Hep`, `Liver`: (Tong et al. 2024) Semi-intrinsic
#' clocks for neurons, glia, hepatocytes, and whole liver.
#'
#' === Gestational Age ===
#' * `BohlinGA`: (Bohlin et al. 2016) The Bohlin Gestational Age (Cord Blood)
#' * `EPICGA`: (Haftorn et al. 2021) The Gestational Age clock based on 176
#' Illumina EPIC CpGs.
#' * `LeeGA`: (Lee et al. 2019) The Lee gestational age.
#' * `KnightGA`: (Knight et al. 2016) The Knight gestational age.
#' * `MayneGA`: (Mayne et al. 2017) The Mayne gestational age.
#'
#' === Surrogate Biomarkers ===
#' * `CRP`: (Wielscher et al. 2022) The DNAm-based C-Reactive Protein Score.
#' * `CHIP`: (Kirmani et al. 2025) The CHIP-related Methylation Scores.
#' * `IL6`: (Stevenson et al. 2021) The DNAm-Based Proxy for IL-6.
#' * `EpiScores`: (Gadd et al. 2022) The Epigenetic Scores for the
#' Circulating Proteome
#'
#' === Trait Prediction ===
#' * `McCartneyTrait`: (McCartney et al. 2018) The Epigenetic Surrogates for
#' Complex Traits.
#'
#' === Disease Risk Prediction ===
#' * `SmokeIndex`: (Teschendorff et al. 2015) The DNAm-based Smoking Index.
#' * `HepatoXuRisk`: (Xu et al. 2017) The HepatoXu ctDNA Methylation
#' Scores for Hepatocellular Carcinoma.
#'
#' === Cross Species Support===
#' * `EnsembleAge`: (Haghani et al. 2025) The EnsembleAge Epigenetic Clocks.
#' 
#' * `UniversalPanMammalianClocks`: (Lu et al. 2023) The Universal Pan-Mammalian
#'  Epigenetic Clocks.
#' * `PanMammalianBlood`: (Lu et al. 2023) The Universal Pan-Mammalian Blood
#' Epigenetic Clocks.
#' * `PanMammalianSkin`: (Lu et al. 2023) The Universal Pan-Mammalian Skin
#' Epigenetic Clocks.
#'
#' * `PetkovichMouse`: (Petkovich et al. 2017)  Petkovich Blood Mouse DNA 
#' Methylation Age Predictor (RRSB clock).
#' * `ThompsonMouse`: (Thompson et al. 2018)  Thompson Mouse Multi-Tissue DNA 
#' Methylation Age Predictor (RRSB clock).
#' * `MeerMouse`: (Meer et al. 2018)  Meer Whole Lifespan Multi-Tissue Mouse DNA
#'  Methylation Age Predictor (RRSB clock).
#'
#' Each vector is named with the sample IDs from the `rownames` of `betaM`.
#'
#' @return A list containing the entries for the calculated clocks.
#'
#' * `epiTOC1`: The epiTOC1 score.
#' * `epiTOC2`: A list (tnsc, tnsc2, irS, etc.).
#' * `epiTOC3`: A list (tnsc, tnsc2, irS, etc.).
#' * `stemTOCvitro`: The stemTOCvitro score.
#' * `stemTOC`: The stemTOC score.
#' * `RepliTali`: The RepliTali score.
#' * `HypoClock`: The HypoClock score.
#' * `EpiCMIT_hyper`: The EpiCMIT_hyper score.
#' * `EpiCMIT_hypo`: The EpiCMIT_hypo score.
#' * `DNAmTL`: The Leukocyte telomere length.
#' * `Horvath2013`: Horvath epigenetic clock age.
#' * `Hannum`: Hannum epigenetic clock age.
#' * `Lin`: Lin epigenetic clock age.
#' * `VidalBralo`: VidalBralo epigenetic clock age.
#' * `ZhangClock`: Zhang epigenetic clock age.
#' * `Horvath2018`: Horvath skin & blood clock age.
#' * `Bernabeu_cAge`: Bernabeu cAge.
#' * `CorticalClock`: CorticalClock age.
#' * `PedBE`: PedBE clock age.
#' * `Zhang10`: Zhang 10-CpG score.
#' * `PhenoAge`: PhenoAge score.
#' * `DunedinPACE`: DunedinPACE score.
#' * `GrimAge1`: A data.frame (Sample, Age, Female, DNAm..., DNAmGrimAge1)
#' * `GrimAge2`: A data.frame (Sample, Age, Female, DNAm..., DNAmGrimAge2).
#' * `CausalClock`: A list (Causal, Damage, Adaptation scores).
#' * `IC_Clock`: The Intrinsic Capacity (IC) score
#' * `ABEC`: ABEC age.
#' * `eABEC`: eABEC age.
#' * `cABEC`: cABEC age
#' * `PipekElasticNet`: PipekElasticNet age
#' * `PipekFilteredh`:  PipekFilteredh age
#' * `PipekRetrainedh`:  PipekRetrainedh age
#' * `WuClock`:  WuClock age
#' * `Weidner`: Weidner age
#' * `IntrinClock`: IntrinClock Age Prediction.
#' * `Garagnani`: Garagnani Age Score.
#' * `BohlinGA`: Bohlin Gestational Age
#' * `EPICGA`: EPIC Gestational Age
#' * `KnightGA`:  Knight Gestational Age
#' * `LeeGA`:  A list (LeeControl, LeeRobust, LeeRefinedRobust).
#' * `MayneGA`: Mayne Gestational Age
#' * `CRP`: A list (CRP, intCRP).
#' * `CHIP`: A list (AnyCHIP, DNMT3A, TET2, ASXL1).
#' * `IL6`:  IL-6 proxy score
#' * `EpiScores`:  A list including 109 protein score.
#' * `McCartneyTrait`: A list containing 10 named elements, one for each trait.
#' * `SmokeIndex`: A DNAm-based Smoking Index.
#' * `HepatoXuRisk`: The HepatoXuRisk scores for hepatocellular carcinoma.
#' * `EnsembleAge_HumanMouse`: A list where each element is a named numeric vector of
#' predicted values.
#' * `EnsembleAge_Static`: A list where each element is a named numeric vector of
#' predicted values.
#' * `EnsembleAge_Dynamic`: A list where each element is a named numeric vector of
#' predicted values.
#' * `UniversalPanMammalianClocks`: A data.frame (Sample, SpeciesLatinName, ...,
#' DNAmAgePanMammalianClock3).
#' * `PanMammalianBlood`: A data.frame (Sample, SpeciesLatinName, ...,
#' DNAmAgePanMammalianBlood3).
#' * `PanMammalianSkin`: A data.frame (Sample, SpeciesLatinName, ...,
#' DNAmAgePanMammalianSkin3).
#' * `PetkovichMouse`: A numeric vector of predicted mouses chronological ages in months.
#' * `ThompsonMouse`: A numeric vector of predicted mouses chronological ages in months.
#' * `MeerMouse`: A numeric vector of predicted mouses chronological ages in months.
#' * `PCClocks`: (If requested) A data.frame containing the original
#' `SampleID` columns, appended with 14 new columns for
#' the calculated PC clock values (including `PCHorvath2013`, `PCPhenoAge`,
#' `PCDNAmTL`, etc.).
#' * `SystemsAge`: (If requested) A data.frame where the first column is
#' `SampleID` and the subsequent 13 columns contain the calculated scores
#' * `DNAmFitAge`: (If requested) A single data.frame containing DNAmFitAge, 
#' and all 6 related fitness biomarkers.
#' * `StochClocks`: A list containing three numeric vectors: StocH, StocP, and
#' StocZ, representing predicted DNAm ages.#'
#'
#' @references
#' Yang Z, Wong A, Kuh D, et al.
#' Correlation of an epigenetic mitotic clock with cancer risk.
#' \emph{Genome Biol.} 2016
#'
#' Teschendorff AE.
#' A comparison of epigenetic mitotic-like clocks for cancer risk prediction.
#' \emph{Genome Med.} 2020
#'
#' Endicott JL, Nolte PA, Shen H, Laird PW.
#' Cell division drives DNA methylation loss in late-replicating domains in
#' primary human cells.
#' \emph{Nat Commun.} 2022
#'
#' Duran-Ferrer M, Clot G, Nadeu F, et al.
#' The proliferative history shapes the DNA methylome of B-cell tumors and
#' predicts clinical outcome.
#' \emph{Nat Cancer} 2020
#' Horvath S.
#' DNA methylation age of human tissues and cell types.
#' \emph{Genome Biol.} 2013
#'
#' Zhang Q, Vallerga CL, Walker RM, et al.
#' Improved precision of epigenetic clock estimates across tissues and its
#' implication for biological ageing.
#' \emph{Genome Med.} 2019
#'
#' Levine ME, Lu AT, Quach A, et al.
#' An epigenetic biomarker of aging for lifespan and healthspan.
#' \emph{Aging} 2018
#'
#' Belsky DW, Caspi A, Corcoran DL, et al.
#' DunedinPACE, a DNA methylation biomarker of the pace of aging.
#' \emph{eLife} 2022
#'
#' Lu AT, Quach A, Wilson JG, et al.
#' DNA methylation GrimAge strongly predicts lifespan and healthspan
#' \emph{Aging} 2019
#'
#' Lu AT, Binder AM, Zhang J, et al.
#' DNA methylation GrimAge version 2.
#' \emph{Aging} 2022
#'
#' Ying K, Liu H, Tarkhov AE, et al.
#' Causality-enriched epigenetic age uncouples damage and adaptation.
#' \emph{Nat Aging} 2024
#'
#' Hannum G, Guinney J, Zhao L, et al.
#' Genome-wide methylation profiles reveal quantitative views of human aging rates.
#' \emph{Mol Cell.} 2013
#'
#' Lin Q, Weidner CI, Costa IG, et al.
#' DNA methylation levels at individual age-associated CpG sites can be
#' indicative for life expectancy.
#' \emph{Aging} 2016
#'
#' Vidal-Bralo L, Lopez-Golan Y, Gonzalez A.
#' Simplified Assay for Epigenetic Age Estimation in Whole Blood of Adults.
#' \emph{Front Genet.} 2016
#'
#' Zhang Y, Wilson R, Heiss J, et al.
#' DNA methylation signatures in peripheral blood strongly predict all-cause mortality.
#' \emph{Nat Commun.} 2017
#'
#' Horvath S, Oshima J, Martin GM, et al.
#' Epigenetic clock for skin and blood cells applied to Hutchinson Gilford
#' Progeria Syndrome and ex vivo studies.
#' \emph{Aging} 2018
#'
#' McEwen LM, O'Donnell KJ, McGill MG, et al.
#' The PedBE clock accurately estimates DNA methylation age in pediatric buccal cells.
#' \emph{Proc Natl Acad Sci U S A.} 2020
#'
#' Dec, E., Clement, J., Cheng, K. et al.
#' Centenarian clocks: epigenetic clocks for validating claims of exceptional longevity.
#' \emph{GeroScience} 2023
#'
#' Ndhlovu LC, Bendall ML, Dwaraka V, et al.
#' Retro-age: A unique epigenetic biomarker of aging captured by DNA methylation
#' states of retroelements.
#' \emph{Aging Cell.} 2024
#'
#' Higgins-Chen AT, Thrush KL, Wang Y, et al.
#' A computational solution for bolstering reliability of epigenetic clocks:
#' Implications for clinical trials and longitudinal tracking.
#' \emph{Nat Aging.} (2022).
#'
#' Sehgal, R., Markov, Y., Qin, C. et al.
#' Systems Age: a single blood methylation test to quantify aging heterogeneity
#' across 11 physiological systems.
#' \emph{Nat Aging} (2025).
#'
#' Tong, H., Dwaraka, V.B., Chen, Q. et al.
#' Quantifying the stochastic component of epigenetic aging.
#' \emph{Nat Aging} (2024). \doi{10.1038/s43587-024-00636-6}
#'
#' Fuentealba M, Rouch L, Guyonnet S, et al.
#' A blood-based epigenetic clock for intrinsic capacity predicts mortality and
#' is associated with clinical, immunological and lifestyle factors.
#' \emph{Nature Aging.} 2025
#'
#' McGreevy KM, Radak Z, Torma F, et al.
#' DNAmFitAge: biological age indicator incorporating physical fitness.
#' \emph{Aging} 2023
#'
#' Lu AT, Seeboth A, Tsai PC, et al.
#' DNA methylation-based estimator of telomere length
#' \emph{Aging} 2019
#'
#' Tong H, Guo X, Jacques M, Luo Q, Eynon N, Teschendorff AE.
#' Cell-type specific epigenetic clocks to quantify biological age at cell-type
#' resolution.
#' \emph{Aging} 2024
#'
#' Lee, Y., Haftorn, K.L., Denault, W.R.P. et al.
#' Blood-based epigenetic estimators of chronological age in human adults using
#' DNA methylation data from the Illumina MethylationEPIC array.
#' \emph{BMC Genomics} 2020
#'
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for
#' methylation array data.
#' \emph{J Math Chem} 2023
#'
#' Wu, Xiaohui et al.
#' DNA methylation profile is a quantitative measure of biological aging in
#' children
#' \emph{Aging} 2019
#'
#' Weidner, C.I., Lin, Q., Koch, C.M. et al.
#' Aging of blood can be tracked by DNA methylation changes at just three
#' CpG sites.
#' \emph{Genome Biol} 2014
#'
#' Tomusiak, A., Floro, A., Tiwari, R. et al.
#' Development of an epigenetic clock resistant to changes in immune cell
#' composition.
#' \emph{Commun Biol} 2024
#'
#' Garagnani, P. et al.
#' Methylation of ELOVL2 gene as a new epigenetic marker of age.
#' \emph{Aging Cell} 2012
#'
#' Bohlin J, Håberg SE, Magnus P, et al.
#' Prediction of gestational age based on genome-wide differentially
#' methylated regions. \emph{Genome Biol.} 2016
#'
#' Haftorn KL, Lee Y, Denault WRP, et al.
#' An EPIC predictor of gestational age and its application to newborns
#' conceived by assisted reproductive technologies.
#' \emph{Clin Epigenetics.} 2021
#'
#' Lee Y, Choufani S, Weksberg R, et al.
#' Placental epigenetic clocks: estimating gestational age using placental
#' DNA methylation levels.
#' \emph{Aging} 2019
#'
#' Knight AK, Craig JM, Theda C, et al.
#' An epigenetic clock for gestational age at birth
#' based on blood methylation data.
#' \emph{Genome Biol.} 2016
#'
#' Mayne BT, Leemaqz SY, Smith AK, Breen J, Roberts CT, Bianco-Miotto T.
#' Accelerated placental aging in early onset preeclampsia pregnancies
#' identified by DNA methylation.
#' \emph{Epigenomics} 2017
#'
#' Wielscher, M., Mandaviya, P.R., Kuehnel, B. et al.
#' DNA methylation signature of chronic low-grade inflammation and its role in
#' cardio-respiratory diseases.
#' \emph{Nat Commun} 2022
#'
#' Kirmani, S., Huan, T., Van Amburg, J.C. et al.
#' Epigenome-wide DNA methylation association study of CHIP provides insight
#' into perturbed gene regulation.
#' \emph{Nat Commun} 2025
#'
#' Stevenson AJ et al.
#' Creating and Validating a DNA Methylation-Based Proxy for Interleukin-6.
#' \emph{J Gerontol A Biol Sci Med Sci.} 2021
#'
#' Gadd DA, Hillary RF, McCartney DL, et al.
#' Epigenetic scores for the circulating proteome as tools for
#' disease prediction.
#' \emph{Elife.} 2022
#'
#' McCartney DL, Hillary RF, Stevenson AJ, et al.
#' Epigenetic prediction of complex traits and death.
#' \emph{Genome Biol.} 2018
#'
#' Teschendorff, Andrew E et al.
#' Correlation of Smoking-Associated DNA Methylation Changes in Buccal Cells
#' With DNA Methylation Changes in Epithelial Cancer
#' \emph{JAMA oncology} 2015
#'
#' Xu, Rh., Wei, W., Krawczyk, M. et al.
#' Circulating tumour DNA methylation markers for diagnosis and prognosis of
#' hepatocellular carcinoma.
#' \emph{Nature Mater} 2017
#'
#' Haghani, A., Lu, A.T., Yan, Q. et al.
#' EnsembleAge: enhancing epigenetic age assessment with a
#' multi-clock framework.
#' \emph{GeroScience} 2025.
#'
#' Lu, A.T., Fei, Z., Haghani, A. et al.
#' Universal DNA methylation age across mammalian tissues.
#' \emph{Nat Aging.} 2023
#'
#' @examples
#' lungInv <- loadOmniAgeRdata(
#'     "omniager_lung_inv",
#'     verbose = FALSE
#' )
#' lungInvM <- lungInv$bmiq_m
#' phenoDf <- lungInv$PhenoTypes
#' epiMarkerOut <- epiMarker(
#'     betaM = lungInvM,
#'     clockNames = "mitotic",
#'     chronAge = phenoDf$Age,
#'     minCoverage = 0
#' )
#' ## Downloading "PCClocks_data" and "SystemsAge_data" will take a very long time.
#' \donttest{
#' hannumExample <- loadOmniAgeRdata(
#'     "omniager_hannum_example",
#'     verbose = FALSE
#' )
#' hannumBmiqM <- hannumExample[[1]]
#' phenoTypesHannum <- hannumExample[[2]]
#' age <- phenoTypesHannum$Age
#' sex <- ifelse(phenoTypesHannum$Sex == "F", "Female", "Male")
#'
#' pcClockData <- loadOmniAgeRdata(
#'     "PCClocks_data",
#'     verbose = FALSE
#' )
#' systemsAgeData <- loadOmniAgeRdata(
#'     "SystemsAge_data",
#'     verbose = FALSE
#' )
#'
#' clockNames <- c(
#'     "cellularAging", "chronological", "biological", "causal",
#'     "stochastic"
#' )
#'
#' epiMarkerOut <- epiMarker(
#'     betaM = hannumBmiqM,
#'     clockNames = clockNames,
#'     chronAge = age,
#'     sexVec = sex,
#'     minCoverage = 0,
#'     pcClockData = pcClockData,
#'     systemsAgeData = systemsAgeData
#' )
#'}
#' @export
#'


epiMarker <- function(betaM,
                      clockNames = "all",
                      chronAge = NULL,
                      sexVec = NULL,
                      minCoverage = 0,
                      verbose = TRUE,
                      ...) {
    # --- 1. Obtain the classification mapping and expand the clock list ---
    clockMap <- .getClockCategories()
    clocksToRun <- .resolveClockLogic(clockNames, clockMap)

    # --- 2. Obtain the classification mapping and expand the clock list ---
    clocksToRun <- .enforceDependencies(clocksToRun)
    clocksToRun <- .filterByParameters(clocksToRun, chronAge, sexVec, verbose)

    if (length(clocksToRun) == 0) {
        stop("[EpiMarker] No valid clocks selected or missing required Age/Sex metadata.")
    }

    # --- 3. Prepare the list for storing the results ---
    resultsList <- list()
    extraArgs <- list(...)
    # --- 4. Perform prediction ---
    # 1. EpiCMIT
    epicmitRequested <- intersect(clocksToRun, clockMap$epicmitGroup)
    if (length(epicmitRequested) > 0) {
        epicmitObj <- epiCMIT(betaM, minCoverage, verbose)
        if ("EpiCMIT_Hyper" %in% epicmitRequested) resultsList$EpiCMIT_Hyper <- epicmitObj$hyperScore
        if ("EpiCMIT_Hypo" %in% epicmitRequested) resultsList$EpiCMIT_Hypo <- epicmitObj$hypoScore
        clocksToRun <- setdiff(clocksToRun, clockMap$epicmitGroup)
    }
    # 2. CausalClock
    causalRequested <- intersect(clocksToRun, clockMap$causal)
    if (length(causalRequested) > 0) {
        resultsList$CausalClock <- causalClock(betaM, minCoverage, verbose)
        clocksToRun <- setdiff(clocksToRun, clockMap$causal)
    }
    # 3. PC Clocks
    pcGroup <- intersect(clocksToRun, clockMap$pcClocks)
    if (length(pcGroup) > 0 && !is.null(extraArgs$pcClockData)) {
        if (verbose) message("[EpiMarker] Calculating PCClocks bundle...")
        resultsList$PCClocks <- pcClocks(
            betaM, chronAge, sexVec,
            extraArgs$pcClockData,
            minCoverage, verbose
        )
        clocksToRun <- setdiff(clocksToRun, clockMap$pcClocks)
    }
    # 4. SystemsAge
    if ("SystemsAge" %in% clocksToRun && !is.null(extraArgs$systemsAgeData)) {
        resultsList$SystemsAge <- systemsAge(
            betaM, extraArgs$systemsAgeData,
            minCoverage, verbose
        )
        clocksToRun <- setdiff(clocksToRun, "SystemsAge")
    }
    # 5. Stochastic Clocks
    stochRequested <- intersect(clocksToRun, clockMap$stochClocks)
    if (length(stochRequested) > 0) {
        resultsList$StochClocks <- stochClocks(betaM, minCoverage, verbose)
        clocksToRun <- setdiff(clocksToRun, clockMap$stochClocks)
    }
    # 6. EnsembleAge
    ensembleRequested <- intersect(clocksToRun, clockMap$ensembleAge)
    
    if (length(ensembleRequested) > 0) {
      if (verbose) message("[EpiMarker] Calculating requested EnsembleAge clocks...")
      
      for (ensClock in ensembleRequested) {
        versionName <- sub("EnsembleAge_", "", ensClock)
        
        resultsList[[ensClock]] <- ensembleAge(
          betaM = betaM,
          clockVersion = versionName,
          minCoverage = minCoverage,
          verbose = verbose
        )
      }
      clocksToRun <- setdiff(clocksToRun, clockMap$ensembleAge)
    }
 
    # 7. Pan-Mammalian
    panRequested <- intersect(clocksToRun, clockMap$panMammalian)
    if (length(panRequested) > 0) {
        if (is.null(extraArgs$speciesName) || is.null(extraArgs$anageData)) {
            warning("[EpiMarker] Skipping Pan-Mammalian Clocks: speciesName or anageData missing.")
        } else {
            resultsList <- .runPanMammalianLogic(
                resultsList, betaM, panRequested,
                extraArgs, minCoverage, verbose
            )
        }
        clocksToRun <- setdiff(clocksToRun, clockMap$panMammalian)
    }
    # 8. Cell type-specific clock
    ctsGroup <- intersect(clocksToRun, clockMap$cellTypeSpecific)
    if (length(ctsGroup) > 0) {
        resultsList$CTSClocks <- ctsClocks(betaM,
            compClocks = ctsGroup,
            dataType = extraArgs$ctsDataType,
            tissue = extraArgs$ctsTissue,
            ctfM = extraArgs$ctfM,
            minCoverage = minCoverage,
            verbose = verbose
        )
        clocksToRun <- setdiff(clocksToRun, clockMap$cellTypeSpecific)
    }
    # 9. Remaining independent clock calculation
    clocksToRunClean <- setdiff(clocksToRun, "DNAmFitAge")
    for (clockLabel in clocksToRunClean) {
        resultsList[[clockLabel]] <- .dispatchIndividualClock(
            clockLabel, betaM,
            chronAge, sexVec,
            minCoverage, verbose
        )
    }
    # 10. DNAmFitAge
    if ("DNAmFitAge" %in% clocksToRun) {
        grimObj <- resultsList$GrimAge1

        if (!is.null(grimObj) && "DNAmGrimAge1" %in% colnames(grimObj)) {
            if (verbose) message("[EpiMarker] Calculating DNAmFitAge using calculated DNAmGrimAge1...")

            resultsList$DNAmFitAge <- dnamFitAge(
                betaM = betaM,
                age = chronAge,
                sex = sexVec,
                grimageVector = grimObj$DNAmGrimAge1,
                minCoverage = minCoverage,
                verbose = verbose
            )
        } else {
            warning("[EpiMarker] DNAmFitAge skipped: GrimAge1 result not found. DNAmFitAge requires GrimAge1.")
        }
    }

    if (verbose) message("[EpiMarker] All requested calculations successfully completed.")
    return(resultsList)
}

#' Clock Category Mapping
#'
#' @description
#' Defines the hierarchical structure and groups of all supported epigenetic
#' clocks. This centralizes the clock registry for easier maintenance.
#'
#' @return A named list where each element is a character vector of clock names.
#' @keywords internal
.getClockCategories <- function() {
    mitotic <- c(
        "epiTOC1", "epiTOC2", "epiTOC3", "stemTOCvitro", "stemTOC",
        "RepliTali", "HypoClock", "EpiCMIT_Hyper", "EpiCMIT_Hypo"
    )
    dnamtl <- c("DNAmTL", "PCDNAmTL")
    stochastic <- c("StocH", "StocZ", "StocP")
    panMammalian <- c("UniversalPanMammalianClocks", "PanMammalianBlood", 
                      "PanMammalianSkin")
    ensembleAge <- c("EnsembleAge_HumanMouse", 
                     "EnsembleAge_Static", 
                     "EnsembleAge_Dynamic")
    cts <- c("Neu-In", "Glia-In", "Neu-Sin", "Glia-Sin", "Hep")
    causal <- c("CausalAge", "DamAge", "AdaptAge")

    list(
        mitotic = mitotic,
        dnamtl = dnamtl,
        cellularAging = c(mitotic, dnamtl),
        chronological = c(
            "Horvath2013", "Hannum", "Lin", "VidalBralo", "ZhangClock",
            "Horvath2018", "Bernabeu_cAge", "CorticalClock", "PedBE",
            "ABEC", "eABEC", "cABEC", "PipekElasticNet",
            "PipekFilteredh", "PipekRetrainedh", "WuClock",
            "Weidner", "IntrinClock", "Garagnani",
            "CentenarianClock", "Retro_age", "PCHorvath2013", "PCHorvath2018",
            "PCHannum"
        ),
        biological = c(
            "Zhang10", "PhenoAge", "DunedinPACE", "GrimAge1", "GrimAge2", "PCPhenoAge",
            "PCGrimAge1", "DNAmFitAge", "IC_Clock", "SystemsAge"
        ),
        causal = causal,
        stochastic = stochastic,
        cellTypeSpecific = cts,
        gestationalAge = c("BohlinGA", "EPICGA", "KnightGA", "LeeGA", "MayneGA"),
        surrogateBiomarkers = c("CRP", "CHIP", "IL6", "EpiScores"),
        traitPred = c("McCartneyTrait"),
        diseaseRisk = c("SmokeIndex", "HepatoXuRisk"),
        crossSpeciesSupport = c(ensembleAge, panMammalian, "PetkovichMouse", "ThompsonMouse", "MeerMouse"),
        pcClocks = c("PCHorvath2013", "PCHorvath2018", "PCHannum", "PCPhenoAge", "PCGrimAge1", "PCDNAmTL"),
        stochClocks = stochastic,
        panMammalian = panMammalian,
        ensembleAge = ensembleAge,
        epicmitGroup = c("EpiCMIT_Hyper", "EpiCMIT_Hypo")
    )
}


#' Resolve Clock Selection Logic
#'
#' @description
#' Expands high-level category keywords into individual clock labels and
#' merges them with explicitly requested clocks.
#'
#' @param clockNames A character vector of requested clocks or categories.
#' @param clockMap A named list of categories from \code{.getClockCategories}.
#'
#' @return A unique character vector of all individual clock labels to run.
#' @keywords internal
.resolveClockLogic <- function(clockNames, clockMap) {
    if ("all" %in% clockNames) {
        return(unique(unlist(clockMap)))
    }

    resolved <- c()
    for (category in names(clockMap)) {
        if (category %in% clockNames) {
            resolved <- c(resolved, clockMap[[category]])
        }
    }
    finalClocks <- unique(c(resolved, setdiff(clockNames, names(clockMap))))
    return(finalClocks)
}


#' Enforce Clock Dependencies
#'
#' @description
#' Automatically adds precursor clocks to the execution list if a downstream
#' clock (e.g., DNAmFitAge) depends on their results.
#'
#' @param clocksToRun A character vector of requested individual clocks.
#'
#' @return A character vector including necessary dependent clocks.
#' @keywords internal
.enforceDependencies <- function(clocksToRun) {
    if ("DNAmFitAge" %in% clocksToRun) {
        if (!any(c("GrimAge1", "GrimAge2") %in% clocksToRun)) {
            message("[EpiMarker] Dependency: Adding 'GrimAge1' to run list as
              it is required by 'DNAmFitAge'.")
            clocksToRun <- c(clocksToRun, "GrimAge1")
        }
    }
    return(unique(clocksToRun))
}


#' Filter Clocks by Available Metadata
#'
#' @description
#' Removes clocks from the execution list if their required metadata
#' (chronological age or sex) is missing, preventing runtime errors.
#'
#' @param clocks Character vector of clocks to check.
#' @param chronAge Numeric vector of chronological ages (optional).
#' @param sexVec Character/factor vector of biological sex (optional).
#' @param verbose Logical, whether to issue warnings for skipped clocks.
#'
#' @return A filtered character vector of compatible clocks.
#' @keywords internal
.filterByParameters <- function(clocks, chronAge, sexVec, verbose) {
    needsAge <- c(
        "epiTOC2", "epiTOC3", "GrimAge1", "GrimAge2",
        "DNAmFitAge", "PCGrimAge1"
    )
    needsSex <- c("GrimAge1", "GrimAge2", "DNAmFitAge", "PCGrimAge1")

    removeList <- c()

    if (is.null(chronAge)) {
        incompatible <- intersect(clocks, needsAge)
        if (length(incompatible) > 0) {
            if (verbose) warning("[EpiMarker] Missing 'chronAge'. Skipping: ", paste(incompatible, collapse = ", "))
            removeList <- c(removeList, incompatible)
        }
    }

    if (is.null(sexVec)) {
        incompatible <- intersect(clocks, needsSex)
        if (length(incompatible) > 0) {
            if (verbose) warning("[EpiMarker] Missing 'sexVec'. Skipping: ", paste(incompatible, collapse = ", "))
            removeList <- c(removeList, incompatible)
        }
    }

    return(setdiff(clocks, removeList))
}


#' Dispatch Individual Clock Calculations
#'
#' @description
#' The central routing unit that maps a clock label to its specific
#' implementation function.
#'
#' @param clockLabel String, the name of the clock to calculate.
#' @param betaM Numeric DNA methylation matrix (probes x samples).
#' @param chronAge Numeric vector of chronological ages.
#' @param sexVec Vector of biological sex.
#' @param minCoverage Numeric, minimum probe coverage required (0 to 1).
#' @param verbose Logical, whether to print progress messages.
#'
#' @return The result of the specific clock function (typically a vector or data.frame).
#' @keywords internal
.dispatchIndividualClock <- function(clockLabel, betaM, chronAge, sexVec, minCoverage, verbose) {
    res <- switch(clockLabel,
        "epiTOC1" = epiTOC1(betaM, minCoverage, verbose),
        "epiTOC2" = epiTOC2(betaM, chronAge, minCoverage, verbose),
        "epiTOC3" = epiTOC3(betaM, chronAge, minCoverage, verbose),
        "stemTOCvitro" = stemTOCvitro(betaM, minCoverage, verbose),
        "stemTOC" = stemTOC(betaM, minCoverage, verbose),
        "HypoClock" = hypoClock(betaM, minCoverage, verbose),
        "RepliTali" = repliTali(betaM, minCoverage, verbose),
        "DNAmTL" = dnamTL(betaM, minCoverage, verbose),
        "Horvath2013" = horvath2013Clock(betaM, minCoverage, verbose),
        "Hannum" = hannumClock(betaM, minCoverage, verbose),
        "Lin" = linClock(betaM, minCoverage, verbose),
        "VidalBralo" = vidalBraloClock(betaM, minCoverage, verbose),
        "ZhangClock" = zhangClock(betaM, minCoverage, verbose),
        "Horvath2018" = horvath2018Clock(betaM, minCoverage, verbose),
        "Bernabeu_cAge" = bernabeuCAge(betaM, minCoverage, verbose),
        "PedBE" = pedBEClock(betaM, minCoverage, verbose),
        "CorticalClock" = corticalClock(betaM, minCoverage, verbose),
        "CentenarianClock" = centenarianClock(betaM, minCoverage, verbose),
        "Retro_age" = retroAge(betaM, minCoverage, verbose),
        "ABEC" = leeABEC(betaM, minCoverage, verbose),
        "eABEC" = leeExtendedABEC(betaM, minCoverage, verbose),
        "cABEC" = leeCommonABEC(betaM, minCoverage, verbose),
        "PipekElasticNet" = pipekElasticNet(betaM, minCoverage, verbose),
        "PipekFilteredh" = pipekFilteredh(betaM, minCoverage, verbose),
        "PipekRetrainedh" = pipekRetrainedh(betaM, minCoverage, verbose),
        "WuClock" = wuClock(betaM, minCoverage, verbose),
        "Weidner" = weidnerClock(betaM, minCoverage, verbose),
        "IntrinClock" = intrinClock(betaM, minCoverage, verbose),
        "Garagnani" = garagnaniClock(betaM, minCoverage, verbose),
        "Zhang10" = zhang10(betaM, minCoverage, verbose),
        "PhenoAge" = phenoAge(betaM, minCoverage, verbose),
        "DunedinPACE" = dunedinPACE(betaM, minCoverage, verbose),
        "GrimAge1" = grimAge1(betaM, chronAge, sexVec, minCoverage, verbose),
        "GrimAge2" = grimAge2(betaM, chronAge, sexVec, minCoverage, verbose),
        "IC_Clock" = icClock(betaM, minCoverage, verbose),
        "BohlinGA" = bohlinGa(betaM, minCoverage, verbose),
        "EPICGA" = epicGa(betaM, minCoverage, verbose),
        "KnightGA" = knightGa(betaM, minCoverage, verbose),
        "LeeGA" = LeeGa(betaM, minCoverage, verbose),
        "MayneGA" = mayneGa(betaM, minCoverage, verbose),
        "CRP" = compCRP(betaM, minCoverage, verbose),
        "CHIP" = compCHIP(betaM, minCoverage, verbose),
        "IL6" = compIL6(betaM, minCoverage, verbose),
        "EpiScores" = compEpiScores(betaM, minCoverage, verbose),
        "McCartneyTrait" = mcCartneyTrait(betaM, minCoverage, verbose),
        "SmokeIndex" = compSmokeIndex(betaM, minCoverage, verbose),
        "HepatoXuRisk" = hepatoXuRisk(betaM, minCoverage, verbose),
        "PetkovichMouse"  = petkovichMouse(betaM, minCoverage, verbose),
        "ThompsonMouse"  = thompsonMouse(betaM, minCoverage, verbose),
        "MeerMouse" = meerMouse(betaM, minCoverage, verbose),
        
        {
            warning("[EpiMarker] Clock", clockLabel, "is not yet implemented in dispatcher.")
            NULL
        }
    )
    return(res)
}


#' Pan-Mammalian Logic Dispatcher
#'
#' @description
#' Handles the specific requirements for cross-species clocks, which require
#' species-specific metadata and lifespan data.
#'
#' @param resultsList Existing list of results to append to.
#' @param betaM Numeric DNA methylation matrix.
#' @param requested Character vector of requested pan-mammalian clocks.
#' @param args List of additional arguments (extraArgs) containing species data.
#' @param minCoverage Numeric, minimum probe coverage.
#' @param verbose Logical, whether to print progress messages.
#'
#' @return Updated resultsList containing pan-mammalian results.
#' @keywords internal
.runPanMammalianLogic <- function(resultsList, betaM, requested, args, minCoverage, verbose) {
    if ("UniversalPanMammalianClocks" %in% requested) {
        resultsList$UniversalPanMammalianClocks <- universalPanMammalianClocks(
            betaM, args$speciesName, args$anageData, minCoverage, verbose
        )
    }
    if ("PanMammalianBlood" %in% requested) {
        resultsList$PanMammalianBlood <- panMammalianBlood(
            betaM, args$speciesName, args$anageData, minCoverage, verbose
        )
    }
    if ("PanMammalianSkin" %in% requested) {
        resultsList$PanMammalianSkin <- panMammalianSkin(
            betaM, args$speciesName, args$anageData, minCoverage, verbose
        )
    }
    return(resultsList)
}
