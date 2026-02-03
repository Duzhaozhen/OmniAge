#' @title List all available epigenetic clocks by category
#'
#' @description
#' Returns a named list containing the names of all individual clocks
#' implemented in the `EpiAge` function, grouped by their functional category.
#' This is useful for knowing which strings are valid for the `clock_names`
#' argument.
#'
#' @return A named list where each name is a clock category (e.g., "mitotic",
#'   "chronological") and each value is a character vector of clock names in
#'   that category.
#'
#' @examples
#' clock_categories <- listEpiAge()
#' @export

listEpiAge <- function() {

  mitotic_clocks <- c("epiTOC1", "epiTOC2", "epiTOC3", "stemTOCvitro", "stemTOC",
                             "RepliTali", "HypoClock", "EpiCMIT_Hyper", "EpiCMIT_Hypo")
  dnamtl_clocks <- c("DNAmTL","PCDNAmTL")
  #cellular_aging_clocks <- c(mitotic_clocks, dnamtl_clocks)


  chronological_clocks <- c("Horvath2013", "Hannum","Lin","VidalBralo","ZhangClock",
                            "Horvath2018","Bernabeu_cAge","CorticalClock","PedBE",
                            "CentenarianClock","Retro_age","ABEC","eABEC","cABEC",
                            "PipekElasticNet","PipekFilteredh","PipekRetrainedh","WuClock",
                            "PCHorvath2013","PCHorvath2018","PCHannum")

  biological_clocks <- c("Zhang10","PhenoAge", "DunedinPACE", "GrimAge1","GrimAge2","PCPhenoAge","PCGrimAge1","DNAmFitAge"
                         ,"IC_Clock","SystemsAge")

  causal_clocks <- c("CausalAge","DamAge","AdaptAge")

  stochastic_clocks <- c("StocH", "StocZ", "StocP")

  #cross_species_supprot <-  c("EnsembleAgeMouse_Static","EnsembleAgeMouse_Dynamic","EnsembleAgeHumanMouse",
  #                            "UniversalPanMammalianClocks","PanMammalianBlood","PanMammalianSkin")

  #ensemble_clocks <- c("EnsembleAgeMouse_Static","EnsembleAgeMouse_Dynamic","EnsembleAgeHumanMouse") #"EnsembleAgeMouse_Static","EnsembleAgeMouse_Dynamic"

  #cross_species_clocks <- c("UniversalPanMammalianClocks","PanMammalianBlood","PanMammalianSkin")

  CTS_clock <- c("Neu-In","Neu-Sin", "Glia-In","Glia-Sin","Hep")

  return(list(
    #cellular_aging = cellular_aging_clocks,
    mitotic = mitotic_clocks,
    dnamtl = dnamtl_clocks,
    chronological = chronological_clocks,
    biological = biological_clocks,
    causal = causal_clocks,
    celltype_specific = CTS_clock
    #ensemble = ensemble_clocks,
    #cross_species = cross_species_supprot
  ))
}


#' @title List available Gestational Age (GA) epigenetic clocks
#'
#' @description
#' Retrieves a character vector containing the names of all implemented
#' Gestational Age (GA) epigenetic clocks available within the package.
#' This utility function is primarily used to identify the valid clock names
#' (strings) that can be supplied to the `clock_names` argument in other
#' functions, such as `EpiGA`.
#'
#' @return A character vector where each element is a string corresponding
#'   to the name of an implemented GA clock.
#'
#' @examples
#' available_ga_clocks <- listEpiGA()
#' @export
#'
listEpiGA <- function() {

  GA_clocks <- c("Bohlin_GA", "EPIC_GA", "Knight_GA", "Lee_GA", "Mayne_GA")

  return(GA_clocks)
}



#' @title Estimate epigenetic age
#'
#' @aliases EpiAge
#' @description
#' Estimates epigenetic age by applying a suite of established DNA methylation
#' (DNAm) aging clocks to a DNAm matrix. This wrapper function provides a
#' unified interface to run multiple clocks, such as "Horvath2013" or
#' "PhenoAge", either individually or as predefined groups (e.g., "cellular_aging").
#'
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @param clock_names A character vector specifying which clocks to calculate.
#'   This can include individual clock names (e.g., "Horvath2013") or special
#'   keywords:
#'   - `"all"`: (Default) Calculates all available clocks:
#'     `cellular_aging` + `chronological` + `biological` + `causal`+ `stochastic` + `ensemble` + `cross_species` + `celltype_specific`.
#'   - `"cellular_aging"`: Calculates all cellular aging clocks (both `mitotic` and `dnamtl`).
#'   - `"mitotic"`: Calculates only the mitotic clocks (e.g., `epiTOC1`, `HypoClock`).
#'   - `"dnamtl"`: Calculates only the DNAm Telomere Length clock (`DNAmTL`,`PCDNAmTL`).
#'   - `"mitotic"`: Calculates all defined mitotic clocks (e.g., `epiTOC1`, `epiTOC2`, `HypoClock`).
#'   - `"chronological"`: Calculates chronological age clocks (e.g., `Horvath2013`, `Hannum`, `ZhangClock`).
#'   - `"biological"`: Calculates biological age clocks (e.g., `PhenoAge`, `DunedinPACE`, `GrimAge1`, `GrimAge2`).
#'   - `"causal"`: Calculates causality-related clocks (`CausalClock`).
#'   - `"stochastic"`: Calculates the stochastic clocks (StocH, StocP, StocZ) as a single bundle.
#'   - `"celltype_specific"`: Calculates all cell-type-specific clocks (`Neu-In`, `Glia-Sin`, etc.) as a single bundle.
#'
#'   Multiple keywords and individual names can be combined (e.g., `c("Horvath2013", "GrimAge2")`).
#'   Use `listEpiAge()` to see a categorized list of all available clocks.
#'
#'
#' @param ages.v
#' Optional argument representing the chronological ages or surrogates thereof of the samples. Vector must be of same length as the number of samples.
#'
#' @param sex.v
#' Optional argument representing the sex of the samples. Vector must be of same length as the number of samples.
#'
#' @param PCClocks_RData
#' **Required if any PC Clock is requested.** This argument specifies the
#' PC clock model data (the `PCClocks_data.qs2` file). This can be:
#' \itemize{
#'   \item A `character` string: The path to the *directory* containing
#'     the `PCClocks_data.qs2` file.
#'   \item A `list`: The pre-loaded data object (from
#'     `load_OmniAgeR_data(object_name = "PCClocks_data")`).
#' }
#' See `?PCClocks` for details.
#'
#' @param SystemsAge_RData
#' **Required if 'SystemsAge' is requested.** This argument specifies the
#' SystemsAge model data (the `SystemsAge_data.qs` file). This can be:
#' \itemize{
#'   \item A `character` string: The path to the *directory* containing
#'     the `SystemsAge_data.qs` file.
#'   \item A `list`: The pre-loaded data object (from
#'     `load_OmniAgeR_data(object_name = "SystemsAge_data")`).
#' }
#' See `?SystemsAge` for details.
#'
#' @param Stoch_refM.m
#' **Optional (for Stochastic Clocks).** A reference matrix for cell-type
#' deconvolution (e.g., from `EpiDISH`). If provided (along with `ages.v`),
#' this will be used to calculate Intrinsic Age Acceleration (IAA). (Default: `NULL`).
#'
#' @param CTS_dataType
#' **Required if any Cell-Type-Specific Clock is requested.**
#' Character string. The type of input sample for CTS_Clocks: `'bulk'` or `'sorted'`.
#'
#' @param CTS_CTF.m
#' **Optional (for CTS_Clocks).** Cell type fraction matrix (Samples x CellTypes).
#' See `?CTS_Clocks` for details on when this is required.
#'
#' @param CTS_tissue
#' **Required if any Cell-Type-Specific Clock is requested.**
#' Character string. The source tissue for CTS_Clocks: `'brain'` or `'otherTissue'`.
#'
#' @param CTS_coreNum
#' **Optional (for CTS_Clocks).** Number of cores for parallel computation
#' in `CTS_Clocks` preprocessing.
#'
#'
#' @details
#' The EpiAge function implements a wide range of epigenetic clocks.
#' This includes cellular aging clocks (e.g., `epiTOC2`, `DNAmTL`), chronological age clocks (e.g., `Horvath2013`, `Hannum`, `ZhangClock`),
#' biological age clocks (e.g., `PhenoAge`, `DunedinPACE`, `GrimAge2`), causality-related clocks ( `CausalAge`, `DamAge`, `AdaptAge`),
#' Stochastic Clocks(`StocH`, `StocP`, `StocZ`), and Cell-Type Specific Clocks(e.g., `Neu-In`, `Glia-In`).
#'
#'
#' === Cellular Aging Clocks ===
#' This category includes clocks measuring cell proliferation (mitotic clocks) and telomere length.
#' **Mitotic Clocks**
#' * `epiTOC1`: (Yang et al. 2016) Average beta value of 385 promoter CpGs.
#' * `epiTOC2`: (Teschendorff et al. 2020) Estimates stem cell divisions using a dynamic model.
#' * `epiTOC3`: Estimates stem cell divisions based on unmethylated population doubling associated CpGs.
#' * `stemTOCvitro`: (Yang et al. 2016) The 0.95 upper quantile of the 629 stemTOCvitro CpGs. The stemTOCvitro CpGs are promoters CpGs that are unmethylated in fetal tissue-types and undergo DNA hypermethylation with increased population-doublings.
#' * `stemTOC`: (Yang et al. 2016) The 0.95 upper quantile of the 371 stemTOC CpGs. Compared to stemTOCvitro CpGs, the stemTOC CpGs are filtered for significant DNA hypermethylation with chronological age in large in-vivo datasets.
#' * `RepliTali`: (Endicott et al. 2022) Based on 87 population doubling associated hypomethylated CpGs.
#' * `HypoClock`: (Yang et al. 2016) Based on hypomethylation at 678 solo-WCGW sites.
#' * `EpiCMIT_hyper`: (Duran-Ferrer et al. 2020) Average beta value of 184 age-associated hypermethylated CpGs.
#' * `EpiCMIT_hypo`: (Duran-Ferrer et al. 2020) Average beta value of 1164 age-associated hypomethylated CpGs.
#' * `EpiCMIT_hypo`: (Duran-Ferrer et al. 2020) Average beta value of 1164 age-associated hypomethylated CpGs.
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
#' * `CorticalClock`: (Shireby et al. 2020) Crain cortical clock based on 347 CpGs.
#' * `PedBE`: (McEwen et al. 2019) Pediatric buccal epithelial clock in Children.
#' * `CentenarianClock`: (Eric Dec et al. 2023) Centenarian Epigenetic Clocks.
#' * `Retro_age`: (Ndhlovu et al. 2024) Retroelement-based clock developed on the EPIC v1/v2 arrays.
#' * `ABEC`: (Lee et al. 2020) Adult Blood-based EPIC Clock.
#' * `eABEC`: (Lee et al. 2020) Extended Adult Blood-based EPIC Clock.
#' * `cABEC`: (Lee et al. 2020) Common Adult Blood-based EPIC Clock.
#' * `PipekElasticNet`: (Pipek et al. 2023) Pipek's Multi-tissue Elastic Net Epigenetic Clock (239 CpGs).
#' * `PipekFilteredh`: (Pipek et al. 2023) Pipek's Filtered Horvath Epigenetic Clock (272 CpGs).
#' * `PipekRetrainedh`: (Pipek et al. 2023) Pipek's Retrained Horvath Epigenetic Clock (308 CpGs).
#' * `WuClock`: (Wu et al. 2019) Wu's Epigenetic Clock for Pediatric Age Estimation.
#' * `PCHorvath2013`, `PCHorvath2018`, `PCHannum`: (Higgins-Chen et al. 2022)
#' Computationally-bolstered versions of the original clocks. These "PC clocks"
#' are retrained using principal components (PCs) derived from a large
#' CpG set to minimize technical noise and improve reliability
#'
#'
#' === Biological Clocks ===
#' * `Zhang10`: (Zhang et al. 2017) 10-CpG clock associated with mortality.
#' * `PhenoAge`: (Levine et al. 2018) Predicts phenotypic age from 513 CpGs.
#' * `DunedinPACE`: (Belsky et al. 2022) Quantifies the pace of biological aging.
#' * `GrimAge1`: The original GrimAge clock (Lu et al. 2019).
#' * `GrimAge2`: (Lu et al. 2022) Updated composite biomarker of mortality risk.
#' * `PCPhenoAge`, `PCGrimAge1`: (Higgins-Chen et al. 2022)
#'   Computationally-bolstered PC versions of the PhenoAge and GrimAge1
#'   clocks, retrained on principal components for enhanced reliability.
#' * `DNAmFitAge`: (McGreevy et al. 2023) Biological age indicator incorporating DNAmGrimAge and 3 DNAm-based physical fitness markers.
#' * `IC_Clock`: (Fuentealba et al. 2025) The Intrinsic Capacity (IC) Clock based on 91 CpGs.
#' * `SystemsAge`: (Sehgal et al. 2025) The Systems Age and 11 System-Specific Scores.
#'
#' === Causal Clocks ===
#' * `CausalAge`, `DamAge`, `AdaptAge`: (Ying et al. 2024) Three clocks derived from a causality-enriched model. They dissect aging into distinct components: 'Causal' (CausalAge), 'Damage' (DamAge), and Adaptation' (AdaptAge).
#'
#' === Stochastic Clocks ===
#' * `StocH`, `StocP`, `StocZ`: (Tong et al. 2024) Stochastic analogues of the
#'   Horvath, PhenoAge, and Zhang clocks, trained on artificial cohorts to
#'   quantify the stochastic component of epigenetic aging.
#'
#' === Cell-Type Specific Clocks ===
#' * `Neu-In`, `Glia-In`, `Brain`: (Tong et al. 2024) Intrinsic clocks for neurons, glia, and whole brain.
#' * `Neu-Sin`, `Glia-Sin`, `Hep`, `Liver`: (Tong et al. 2024) Semi-intrinsic clocks for neurons, glia, hepatocytes, and whole liver.
#'
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
#' * `PCClocks`: (If requested) A data.frame containing the original `Sample_ID`, `Age`, and `Female` columns, appended with 14 new columns for the calculated PC clock values (including `PCHorvath2013`, `PCPhenoAge`, `PCDNAmTL`, etc.).
#' * `SystemsAge`: (If requested) A data.frame where the first column is `Sample_ID` and the subsequent 13 columns contain the calculated scores
#' * `DNAmFitAge`: (If requested) A single data.frame containing DNAmFitAge, FitAgeAccel, and all 6 related fitness biomarkers.
#' * `StochClocks`: (If requested) A list containing the `mage` (DNAm ages)  and `eaa` (Extrinsic Age Acceleration, if `ages.v` was provided).Note: Intrinsic Age Acceleration (IAA) is not calculated via the `EpiAge` wrapper. Please call `StochClocks()` directly if IAA is needed.
#'
#' @seealso
#' `listEpiAge()`
#'
#' `PCClocks()`,
#' `SystemsAge()`,
#' `DNAmFitAge()`,
#' `StochClocks()`
#'
#' `epiTOC1()`, `epiTOC2()`, `epiTOC3()`, `RepliTali()`, `HypoClock()`,
#' `stemTOCvitro()`, `stemTOC()`, `EpiCMIT()`,`DNAmTL()`,
#' `Horvath2013()`, `Hannum()`, `Lin()`, `VidalBralo()`, `ZhangClock()`,
#' `Horvath2018()`, `Bernabeu_cAge()`, `CorticalClock()`, `PedBE()`,
#' `CentenarianClock()`, `Retro_age()`,
#' `Zhang10()`, `PhenoAge()`, `DunedinPACE()`, `GrimAge1()`, `GrimAge2()`,
#' `IC_Clock()`, `CausalClock()`, `CausalClock()`,`ABEC()`,`eABEC()`,`cABEC()`, 
#' `PipekElasticNet()`,`PipekFilteredh()`,`PipekRetrainedh()`,`WuClock()`.
#'
#'
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
#' Cell division drives DNA methylation loss in late-replicating domains in primary human cells.
#' \emph{Nat Commun.} 2022
#'
#' Duran-Ferrer M, Clot G, Nadeu F, et al.
#' The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome.
#' \emph{Nat Cancer} 2020
#' Horvath S.
#' DNA methylation age of human tissues and cell types.
#' \emph{Genome Biol.} 2013
#'
#' Zhang Q, Vallerga CL, Walker RM, et al.
#' Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing.
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
#' DNA methylation levels at individual age-associated CpG sites can be indicative for life expectancy.
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
#' Epigenetic clock for skin and blood cells applied to Hutchinson Gilford Progeria Syndrome and ex vivo studies.
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
#' Retro-age: A unique epigenetic biomarker of aging captured by DNA methylation states of retroelements.
#' \emph{Aging Cell.} 2024
#'
#' Higgins-Chen AT, Thrush KL, Wang Y, et al.
#' A computational solution for bolstering reliability of epigenetic clocks: Implications for clinical trials and longitudinal tracking.
#' \emph{Nat Aging.} (2022).
#'
#' Sehgal, R., Markov, Y., Qin, C. et al.
#' Systems Age: a single blood methylation test to quantify aging heterogeneity across 11 physiological systems.
#' \emph{Nat Aging} (2025).
#'
#' Tong, H., Dwaraka, V.B., Chen, Q. et al.
#' Quantifying the stochastic component of epigenetic aging.
#' \emph{Nat Aging} (2024). \doi{10.1038/s43587-024-00636-6}
#'
#' Fuentealba M, Rouch L, Guyonnet S, et al.
#' A blood-based epigenetic clock for intrinsic capacity predicts mortality and is associated with clinical, immunological and lifestyle factors.
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
#' Cell-type specific epigenetic clocks to quantify biological age at cell-type resolution.
#' \emph{Aging} 2024
#'
#' Lee, Y., Haftorn, K.L., Denault, W.R.P. et al.
#' Blood-based epigenetic estimators of chronological age in human adults using DNA methylation data from the Illumina MethylationEPIC array. 
#' \emph{BMC Genomics} 2020
#'
#' Pipek, O.A., Csabai, I.
#' A revised multi-tissue, multi-platform epigenetic clock model for methylation array data. 
#' \emph{J Math Chem} 2023
#' 
#' Wu, Xiaohui et al. 
#' DNA methylation profile is a quantitative measure of biological aging in children
#' \emph{Aging} 2019
#' 
#' @examples
#' download_OmniAgeR_example("LungInv")
#' load_OmniAgeR_example("LungInv")
#' EpiAge.o<-EpiAge(data.m = bmiq.m,clock_names="mitotic",df$Age)
#'
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' age <- PhenoTypesHannum_lv$Age
#' sex <- ifelse(PhenoTypesHannum_lv$Sex=="F","Female","Male")
#' PCClocks_RData <- load_OmniAgeR_data(object_name = "PCClocks_data")
#' SystemsAge_RData <- load_OmniAgeR_data(object_name = "SystemsAge_data")
#'
#' clock_names<-  c("cellular_aging", "chronological", "biological", "causal", "stochastic")
#'
#' EpiAge.o <- EpiAge(data.m = hannum_bmiq_m,clock_names = clock_names, ages.v = age, sex.v = sex,PCClocks_RData=PCClocks_RData,SystemsAge_RData=SystemsAge_RData)
#'
#' @export
#'

EpiAge <- function(data.m, clock_names = "all", ages.v=NULL, sex.v=NULL,
                   PCClocks_RData = NULL,
                   SystemsAge_RData = NULL,
                   Stoch_refM.m = NULL,
                   CTS_dataType = NULL,
                   CTS_CTF.m = NULL,
                   CTS_tissue = NULL,
                   CTS_coreNum = NULL,
                   #sample_info = NULL,
                   #anage_data = NULL,
                   #Ensemble_version = "HumanMouse",
                   CTF_m) {


  # Define clock categories

  mitotic_clocks <- c("epiTOC1", "epiTOC2", "epiTOC3", "stemTOCvitro", "stemTOC",
                      "RepliTali", "HypoClock", "EpiCMIT_Hyper", "EpiCMIT_Hypo")
  dnamtl_clocks <- c("DNAmTL","PCDNAmTL")

  cellular_aging_clocks <- c(mitotic_clocks, dnamtl_clocks)

  chronological_clocks <- c("Horvath2013", "Hannum","Lin","VidalBralo","ZhangClock",
                            "Horvath2018","Bernabeu_cAge","CorticalClock","PedBE",
                            "ABEC","eABEC","cABEC","PipekElasticNet",
                            "PipekFilteredh","PipekRetrainedh","WuClock",
                            "CentenarianClock","Retro_age","PCHorvath2013","PCHorvath2018","PCHannum")

  biological_clocks <- c("Zhang10","PhenoAge", "DunedinPACE", "GrimAge1","GrimAge2","PCPhenoAge","PCGrimAge1","DNAmFitAge"
                         ,"IC_Clock","SystemsAge")

  causal_clocks <- c("CausalAge","DamAge","AdaptAge")

  stochastic_clocks <- c("StocH", "StocZ", "StocP")

  #ensemble_clocks <- c("EnsembleAgeHumanMouse")

  #cross_species_clocks <- c("UniversalPanMammalianClocks","PanMammalianBlood","PanMammalianSkin")

  cts_clocks <- c("Neu-In", "Glia-In", "Neu-Sin", "Glia-Sin", "Hep")

  # Define the master list of all individual, runnable clocks
  # --- Define composite groups (for keywords) ---
  all_individual_clocks <- unique(c(chronological_clocks, biological_clocks, causal_clocks, cellular_aging_clocks, stochastic_clocks, cts_clocks))

  # --- Define all valid keywords ---
  all_keywords <- c("all","cellular_aging", "mitotic", "chronological", "biological", "causal", "stochastic", "ensemble","cross_species","celltype_specific")
  # ---- Process input to determine which clocks to run ----

  # Start with an empty vector for the final list of clocks
  clocks_to_process <- c()

  # Expand special keywords into lists of clock names
  if ("all" %in% clock_names) {
    clocks_to_process <- c(clocks_to_process, all_individual_clocks)
  }
  if ("cellular_aging" %in% clock_names) {
    clocks_to_process <- c(clocks_to_process, cellular_aging_clocks)
  }
  if ("mitotic" %in% clock_names) {
    clocks_to_process <- c(clocks_to_process, mitotic_clocks)
  }
  if ("dnamtl" %in% clock_names) {
    clocks_to_process <- c(clocks_to_process, dnamtl_clocks)
  }
  if ("chronological" %in% clock_names) {
    clocks_to_process <- c(clocks_to_process, chronological_clocks)
  }
  if ("biological" %in% clock_names) {
    clocks_to_process <- c(clocks_to_process, biological_clocks)
  }
  if ("causal" %in% clock_names) {
    clocks_to_process <- c(clocks_to_process, causal_clocks)
  }
  if ("stochastic" %in% clock_names) {
    clocks_to_process <- c(clocks_to_process, stochastic_clocks)
  }
  # if ("ensemble" %in% clock_names) {
  #   clocks_to_process <- c(clocks_to_process, ensemble_clocks)
  # }
  # if ("cross_species" %in% clock_names) { #
  #   clocks_to_process <- c(clocks_to_process, cross_species_clocks)
  # }
  if ("celltype_specific" %in% clock_names) {
    clocks_to_process <- c(clocks_to_process, cts_clocks)
  }

  # Add any specific clock names provided by the user
  specific_clocks <- setdiff(clock_names, all_keywords)
  clocks_to_process <- unique(c(clocks_to_process, specific_clocks))

  # ---- Validate the final list and run the calculations ----

  # Warn about any invalid clock names and select only the valid ones
  invalid_names <- setdiff(clocks_to_process, all_individual_clocks)
  if (length(invalid_names) > 0) {
    warning("The following clock names are not valid and will be ignored: ", paste(invalid_names, collapse = ", "))
  }
  clocks_to_run <- intersect(clocks_to_process, all_individual_clocks)

  # -------------------------------------------------------------------
  # ---- Dependency Enforcement  ----
  # -------------------------------------------------------------------
  # DNAmFitAge *requires* a GrimAge result. We will enforce GrimAge2.
  if ("DNAmFitAge" %in% clocks_to_run) {
    if (!"GrimAge1" %in% clocks_to_run && !"GrimAge2" %in% clocks_to_run) {
      message("Adding 'GrimAge1' to the calculation list as it is required by 'DNAmFitAge'.")
      clocks_to_run <- c(clocks_to_run, "GrimAge1")
    }
  }


  # Stop if no valid clocks were selected
  if (length(clocks_to_run) == 0) {
    stop("No valid clocks were selected to run. Please check the `clock_names` argument.")
  }

  # ---- Validate arguments and filter clocks ----

  clocks_need_age <- c("epiTOC2", "epiTOC3", "GrimAge1", "GrimAge2",
                       "DNAmFitAge")
  clocks_need_sex <- c("GrimAge1", "GrimAge2", "DNAmFitAge")

  # Check for missing ages.v
  if (is.null(ages.v)) {
    clocks_to_skip_age <- intersect(clocks_to_run, clocks_need_age)
    if (length(clocks_to_skip_age) > 0) {
      warning("`ages.v` is NULL. Skipping clocks that require age: ",
              paste(clocks_to_skip_age, collapse = ", "))
      clocks_to_run <- setdiff(clocks_to_run, clocks_to_skip_age)
    }
  }

  # Check for missing or invalid sex.v
  # We only need to do this if a clock needing sex is still in the list
  clocks_to_check_sex <- intersect(clocks_to_run, clocks_need_sex)

  if (length(clocks_to_check_sex) > 0) {
    if (is.null(sex.v)) {
      warning("`sex.v` is NULL. Skipping clocks that require sex: ",
              paste(clocks_to_check_sex, collapse = ", "))
      clocks_to_run <- setdiff(clocks_to_run, clocks_to_check_sex)
    } else {
      # Validate sex.v format if it's not NULL
      valid_sex_values <- c("Female", "Male")
      # Check for any non-NA values that are NOT "Female" or "Male"
      invalid_sex_entries <- setdiff(unique(stats::na.omit(sex.v)), valid_sex_values)

      if (length(invalid_sex_entries) > 0) {
        warning("Invalid values found in `sex.v`: ",
                paste(shQuote(invalid_sex_entries), collapse = ", "),
                ". Expected 'Female' or 'Male'. Skipping clocks: ",
                paste(clocks_to_check_sex, collapse = ", "))
        clocks_to_run <- setdiff(clocks_to_run, clocks_to_check_sex)
      }
    }
  }

  # Stop if no valid clocks are left to run after filtering
  if (length(clocks_to_run) == 0) {
    stop("No valid clocks to run after filtering for missing arguments (age/sex).")
  }


  # Initialize an empty list to store the results
  results.ls <- list()

  # -------------------------------------------------------------------
  # ---- PCClocks Bundle Calculation ----
  # -------------------------------------------------------------------

  # Define the group of clocks calculated *together* by your PCClocks function
  pc_clock_group <- c("PCHorvath2013", "PCHorvath2018", "PCHannum",
                      "PCPhenoAge", "PCGrimAge1","PCDNAmTL")

  # Check which PC clocks the user *actually* requested
  requested_pc_clocks <- intersect(clocks_to_run, pc_clock_group)

  # Initialize a variable to hold the results
  pc_clock_results_df <- NULL

  if (length(requested_pc_clocks) > 0) {

    # Check PCClocks_RData
    if (is.null(PCClocks_RData)) {
      warning(
        "Skipping all PC Clocks (", paste(requested_pc_clocks, collapse=", "),
        ") because the `PCClocks_RData` argument was not provided. ",
        "Please provide the path or the pre-loaded data object. See `?EpiAge`."
      )
      clocks_to_run <- setdiff(clocks_to_run, pc_clock_group)

    } else if (is.null(ages.v) || is.null(sex.v)) {
      warning(
        "Skipping all PC Clocks (", paste(requested_pc_clocks, collapse=", "),
        ") because the bundle calculation requires `ages.v` and `sex.v` for PCGrimAge."
      )
      clocks_to_run <- setdiff(clocks_to_run, pc_clock_group)

    } else {

      message("Calculating PC Clock bundle (PCHorvath, PCPhenoAge, PCGrimAge, etc.)...")

      pc_clock_results_df <- PCClocks(
        DNAm = data.m,
        age_v = ages.v,
        sex_v = sex.v,
        RData = PCClocks_RData
      )

      results.ls$PCClocks <- pc_clock_results_df
      message("PCClocks results are stored in the 'PCClocks' data.frame.")
    }

    clocks_to_run <- setdiff(clocks_to_run, pc_clock_group)
  }


  # -------------------------------------------------------------------
  # ---- SystemsAge Bundle Calculation ----
  # -------------------------------------------------------------------

  systems_age_clock <- "SystemsAge"

  if (systems_age_clock %in% clocks_to_run) {

    if (is.null(SystemsAge_RData)) {
      warning(
        "Skipping 'SystemsAge' because the `SystemsAge_RData` argument was not provided. ",
        "Please provide the path or the pre-loaded data object. See `?EpiAge`."
      )
    } else {
      #message("Calculating SystemsAge bundle...")

      systems_age_results_df <- SystemsAge(
        DNAm = data.m,
        RData = SystemsAge_RData
      )

      results.ls$SystemsAge <- systems_age_results_df
      message("SystemsAge results are stored in the 'SystemsAge' data.frame.")
    }

    clocks_to_run <- setdiff(clocks_to_run, systems_age_clock)
  }

  # -------------------------------------------------------------------
  # ---- Stochastic Clocks Bundle Calculation ----
  # -------------------------------------------------------------------
  stoch_clock_group <- c("StocH", "StocZ", "StocP")

  if (any(stoch_clock_group %in% clocks_to_run)) {
    message("Calculating Stochastic Clock bundle (StocH, StocP, StocZ)...")

    results.ls$StochClocks <- StochClocks(
      data.m = data.m,
      ages.v = ages.v,
      refM.m = Stoch_refM.m,
      verbose = FALSE)

    message("StochClocks results are stored in the 'StochClocks' list (IAA not calculated).")

    clocks_to_run <- setdiff(clocks_to_run, stoch_clock_group)
  }


  # -------------------------------------------------------------------
  # ---- EnsembleAge Bundle Calculation ----
  # -------------------------------------------------------------------
  # ensemble_clock_name <- "EnsembleAgeHumanMouse"
  #
  # if (ensemble_clock_name %in% clocks_to_run) {
  #   message("Calculating EnsembleAge (HumanMouse version)...")
  #
  #   results.ls$EnsembleAgeHumanMouse <- EnsembleAge(
  #     beta.m = data.m,
  #     clock_version = Ensemble_version,
  #     verbose = TRUE
  #   )
  #   clocks_to_run <- setdiff(clocks_to_run, ensemble_clock_name)
  # }


  # -------------------------------------------------------------------
  # ---- Pan-Mammalian Clocks Calculation ----
  # -------------------------------------------------------------------

  # 1. Define Pan-Mammalian Clocks
  # pan_mammal_group <- c("UniversalPanMammalianClocks", "PanMammalianBlood", "PanMammalianSkin")
  #
  # requested_pan_mammal_clocks <- intersect(clocks_to_run, pan_mammal_group)
  #
  # if (length(requested_pan_mammal_clocks) > 0) {
  #
  #   if (is.null(sample_info) || is.null(anage_data)) {
  #     warning(
  #       "Skipping all Pan-Mammalian Clocks (", paste(requested_pan_mammal_clocks, collapse=", "),
  #       ") because the `sample_info` and `anage_data` arguments were not provided. ",
  #        "Please provide both data.frames. See `?EpiAge`."
  #     )
  #   } else {
  #
  #     if ("UniversalPanMammalianClocks" %in% requested_pan_mammal_clocks) {
  #       message("Calculating UniversalPanMammalianClocks...")
  #       results.ls$UniversalPanMammalianClocks <- UniversalPanMammalianClocks(
  #         sample_info = sample_info,
  #         beta.m = data.m,
  #         anage_data = anage_data,
  #         verbose = TRUE
  #       )
  #     }
  #
  #     if ("PanMammalianBlood" %in% requested_pan_mammal_clocks) {
  #       message("Calculating PanMammalianBlood Clocks...")
  #       results.ls$PanMammalianBlood <- PanMammalianBlood(
  #         sample_info = sample_info,
  #         beta.m = data.m,
  #         anage_data = anage_data,
  #         verbose = TRUE
  #       )
  #     }
  #
  #     if ("PanMammalianSkin" %in% requested_pan_mammal_clocks) {
  #       message("Calculating PanMammalianSkin Clocks...")
  #       results.ls$PanMammalianSkin <- PanMammalianSkin(
  #         sample_info = sample_info,
  #         beta.m = data.m,
  #         anage_data = anage_data,
  #         verbose = TRUE
  #       )
  #     }
  #   }
  #   clocks_to_run <- setdiff(clocks_to_run, pan_mammal_group)
  # }

  # -------------------------------------------------------------------
  # ---- Cell-Type-Specific (CTS) Clocks Calculation ----
  # -------------------------------------------------------------------
  cts_clock_group <- c("Neu-In", "Glia-In", "Neu-Sin", "Glia-Sin", "Hep")
  requested_cts_clocks <- intersect(clocks_to_run, cts_clock_group)

  if (length(requested_cts_clocks) > 0) {
    # Check for required arguments
    if (is.null(CTS_dataType) || is.null(CTS_tissue)) {
      warning(
        "Skipping all CTS Clocks (", paste(requested_cts_clocks, collapse=", "),
        ") because `CTS_dataType` and `CTS_tissue` arguments were not provided. ",
        "Please provide 'bulk'/'sorted' and 'brain'/'otherTissue'. See `?EpiAge`."
      )
    } else {
      message(paste("Calculating Cell-Type-Specific (CTS) Clock bundle for:", paste(requested_cts_clocks, collapse=", ")))

      # Call the CTS_Clocks function
      cts_results_df <- CTS_Clocks(
        data.m = data.m,
        CTSclocks = requested_cts_clocks, # Pass only the ones the user asked for
        dataType = CTS_dataType,
        CTF.m = CTS_CTF.m,
        tissue = CTS_tissue,
        coreNum = CTS_coreNum,
        verbose = TRUE # You could also pass the main verbose flag if you add one
      )

      results.ls$CTS_Clocks <- cts_results_df
      message("CTS_Clocks results are stored in the 'CTS_Clocks' data.frame.")
    }

    # Remove these clocks from further individual processing
    clocks_to_run <- setdiff(clocks_to_run, cts_clock_group)
  }


  # --- Cellular Aging Clock Calculations ---
  if ("epiTOC1" %in% clocks_to_run) { results.ls$epiTOC1 <- epiTOC1(data.m) }
  if ("epiTOC2" %in% clocks_to_run) { results.ls$epiTOC2 <- epiTOC2(data.m, ages.v) }
  if ("epiTOC3" %in% clocks_to_run) { results.ls$epiTOC3 <- epiTOC3(data.m, ages.v) }
  if ("stemTOCvitro" %in% clocks_to_run) { results.ls$stemTOCvitro <- stemTOCvitro(data.m) }
  if ("stemTOC" %in% clocks_to_run) { results.ls$stemTOC <- stemTOC(data.m) }
  if ("HypoClock" %in% clocks_to_run) { results.ls$HypoClock <- HypoClock(data.m) }
  if ("RepliTali" %in% clocks_to_run) { results.ls$RepliTali <- RepliTali(data.m) }
  if ("DNAmTL" %in% clocks_to_run) { results.ls$DNAmTL <- DNAmTL(data.m) }

  # EpiCMIT clocks are calculated together
  if (any(c("EpiCMIT_Hyper", "EpiCMIT_Hypo") %in% clocks_to_run)) {
    epicmit.o <- EpiCMIT(data.m)
    if ("EpiCMIT_Hyper" %in% clocks_to_run) { results.ls$EpiCMIT_Hyper <- epicmit.o$hyperSC }
    if ("EpiCMIT_Hypo" %in% clocks_to_run) { results.ls$EpiCMIT_Hypo <- epicmit.o$hypoSC }
  }

  # --- Chronological Clock Calculations ---
  if ("Horvath2013" %in% clocks_to_run) { results.ls$Horvath2013 <- Horvath2013(data.m) }
  if ("Hannum" %in% clocks_to_run) { results.ls$Hannum <- Hannum(data.m) }
  if ("Lin" %in% clocks_to_run) { results.ls$Lin <- Lin(data.m) }
  if ("VidalBralo" %in% clocks_to_run) { results.ls$VidalBralo <- VidalBralo(data.m) }
  if ("ZhangClock" %in% clocks_to_run) { results.ls$ZhangClock <- ZhangClock(data.m) }
  if ("Horvath2018" %in% clocks_to_run) { results.ls$Horvath2018 <- Horvath2018(data.m) }
  if ("Bernabeu_cAge" %in% clocks_to_run) { results.ls$Bernabeu_cAge <- Bernabeu_cAge(data.m) }
  if ("CorticalClock" %in% clocks_to_run) { results.ls$CorticalClock <- CorticalClock(data.m) }
  if ("PedBE" %in% clocks_to_run) { results.ls$PedBE <- PedBE(data.m) }
  if ("CentenarianClock" %in% clocks_to_run) { results.ls$CentenarianClock <- CentenarianClock(data.m) }
  if ("Retro_age" %in% clocks_to_run) { results.ls$Retro_age <- Retro_age(data.m,version="V2") }
  if ("ABEC" %in% clocks_to_run) { results.ls$ABEC <- ABEC(data.m) }
  if ("eABEC" %in% clocks_to_run) { results.ls$eABEC <- eABEC(data.m) }
  if ("cABEC" %in% clocks_to_run) { results.ls$cABEC <- cABEC(data.m) }
  if ("PipekElasticNet" %in% clocks_to_run) { results.ls$PipekElasticNet <- PipekElasticNet(data.m) }
  if ("PipekFilteredh" %in% clocks_to_run) { results.ls$PipekFilteredh <- PipekFilteredh(data.m) }
  if ("PipekRetrainedh" %in% clocks_to_run) { results.ls$PipekRetrainedh <- PipekRetrainedh(data.m) }
  if ("WuClock" %in% clocks_to_run) { results.ls$WuClock <- WuClock(data.m) }

  # --- Biological Clock Calculations ---
  if ("Zhang10" %in% clocks_to_run) { results.ls$Zhang10 <- Zhang10(data.m) }
  if ("PhenoAge" %in% clocks_to_run) { results.ls$PhenoAge <- PhenoAge(data.m) }
  if ("DunedinPACE" %in% clocks_to_run) { results.ls$DunedinPACE <- DunedinPACE(data.m) }
  if ("GrimAge1" %in% clocks_to_run) { results.ls$GrimAge1 <- GrimAge1(data.m, ages.v, sex.v) }
  if ("GrimAge2" %in% clocks_to_run) { results.ls$GrimAge2 <- GrimAge2(data.m, ages.v, sex.v) }

  if ("IC_Clock" %in% clocks_to_run) { results.ls$IC_Clock <- IC_Clock(data.m) }

  # --- Causal Clock Calculations ---
  if (any(causal_clocks %in% clocks_to_run)) {
    results.ls$CausalClock <- CausalClock(data.m)
  }

  if ("DNAmFitAge" %in% clocks_to_run) {

    grimage_input_vector <- NULL

    if ("GrimAge1" %in% names(results.ls)) {
      grimage_input_vector <- results.ls$GrimAge1$DNAmGrimAge1
      message("Using 'GrimAge1' output as input for 'DNAmFitAge'.")

    } else if ("GrimAge2" %in% names(results.ls)) {
      grimage_input_vector <- results.ls$GrimAge1$DNAmGrimAge2
      message("Using 'GrimAge2' output as input for 'DNAmFitAge'.")
    } else {
      warning("Skipping 'DNAmFitAge': Required 'GrimAge1' or 'GrimAge2' result is not available.")
    }

    if (!is.null(grimage_input_vector)) {
      message("Calculating DNAmFitAge bundle...")
      results.ls$DNAmFitAge <- DNAmFitAge(
        beta_matrix = data.m,
        age_vector = ages.v,
        sex_vector = sex.v,
        grimage_vector = grimage_input_vector,
        verbose = FALSE
      )
      message("DNAmFitAge results are stored in the 'DNAmFitAge' data.frame.")
     }
  }


  return(results.ls)
}



#' @title Estimate Gestational Age (GA)
#'
#' @description
#' Estimates epigenetic Gestational Age (GA) by applying one or more
#' established GA clocks to a DNA methylation (DNAm) matrix.
#' **All calculated estimates are returned in units of gestational weeks.**
#'
#' @details
#' This wrapper function provides a unified interface to run multiple
#' GA clocks. It processes the `clock_names` argument and calls the
#' appropriate underlying clock functions.
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs
#' and columns labeling samples.
#'
#' @param clock_names
#' A character vector specifying which GA clocks to calculate.
#' Defaults to `"all"`.
#'   - `"all"`: Calculates all available GA clocks.
#'   - `c("Knight_GA", "Lee_GA", ...)`: Calculates only the specified clocks.
#'
#' Use `listEpiGA()` to see a character vector of all available GA clock names.
#'
#' @seealso
#' [listEpiGA()]
#' `Bohlin_GA()`, `EPIC_GA()`, `Knight_GA()`, `Lee_GA()`, `Mayne_GA()`
#'
#'
#' @return A list containing the entries for the calculated clocks.
#'
#' * `Bohlin_GA`: The Bohlin Gestational Age.
#' * `EPIC_GA`: The EPIC Gestational Age.
#' * `Lee_GA`: A list (LeeControl, LeeRobust, LeeRefinedRobust).
#' * `Knight_GA`: he Knight Gestational Age
#' * `Mayne_GA`: The Mayne Placental Gestational Age.
#'
#' Each vector is named with the sample IDs from the `rownames` of `beta.m`.
#'
#'
#' @seealso
#' [listEpiGA()]
#' `Bohlin_GA()`, `EPIC_GA()`, `Knight_GA()`, `Lee_GA()`, `Mayne_GA()`
#'
#'
#'
#' @references
#' Bohlin J, Håberg SE, Magnus P, et al.
#' Prediction of gestational age based on genome-wide differentially methylated regions.
#' \emph{Genome Biol.} 2016
#'
#' Haftorn KL, Lee Y, Denault WRP, et al.
#' An EPIC predictor of gestational age and its application to newborns conceived by assisted reproductive technologies.
#' \emph{Clin Epigenetics.} 2021
#'
#' Knight AK, Craig JM, Theda C, et al.
#' An epigenetic clock for gestational age at birth based on blood methylation data.
#' \emph{Genome Biol.} 2016
#'
#' Lee Y, Choufani S, Weksberg R, et al.
#' Placental epigenetic clocks: estimating gestational age using placental DNA methylation levels.
#' \emph{Aging} 2019
#'
#' Mayne BT, Leemaqz SY, Smith AK, Breen J, Roberts CT, Bianco-Miotto T.
#' Accelerated placental aging in early onset preeclampsia pregnancies identified by DNA methylation.
#' \emph{Epigenomics} 2017
#'
#'
#' @examples
#' data("GA_example")
#' download_OmniAgeR_example("GA_example")
#' load_OmniAgeR_example("GA_example")
#' EpiGA.o <- EpiGA(GA_m)
#'
#'
#' @export
#'

EpiGA <- function(data.m, clock_names = "all") {

  # 1. Define the master list of all available GA clocks
  # (This is the "Source of Truth", matching listEpiGA())
  all_ga_clocks <- c("Bohlin_GA", "EPIC_GA", "Knight_GA", "Lee_GA", "Mayne_GA")

  # 2. Define valid keywords (just "all" in this case)
  all_keywords <- c("all")

  # 3. Process the clock_names argument
  clocks_to_process <- c()

  # Expand "all" keyword
  if ("all" %in% clock_names) {
    clocks_to_process <- c(clocks_to_process, all_ga_clocks)
  }

  # Add any specific, non-keyword clock names provided by the user
  specific_clocks <- setdiff(clock_names, all_keywords)
  clocks_to_process <- unique(c(clocks_to_process, specific_clocks))

  # 4. Validate the final list of clocks to run
  invalid_names <- setdiff(clocks_to_process, all_ga_clocks)
  if (length(invalid_names) > 0) {
    warning("The following GA clock names are not valid and will be ignored: ",
            paste(invalid_names, collapse = ", "))
  }

  # Final, validated list of clocks to execute
  clocks_to_run <- intersect(clocks_to_process, all_ga_clocks)

  # Stop if no valid clocks were selected
  if (length(clocks_to_run) == 0) {
    stop("No valid GA clocks were selected to run. Please check the `clock_names` argument.")
  }

  # 5. Initialize results list and run calculations
  results.ls <- list()

  if ("Bohlin_GA" %in% clocks_to_run) {
    results.ls$Bohlin_GA <- Bohlin_GA(data.m)
  }
  if ("EPIC_GA" %in% clocks_to_run) {
    results.ls$EPIC_GA <- EPIC_GA(data.m)
  }
  if ("Knight_GA" %in% clocks_to_run) {
    results.ls$Knight_GA <- Knight_GA(data.m)
  }
  if ("Lee_GA" %in% clocks_to_run) {
    results.ls$Lee_GA <- Lee_GA(data.m)
  }
  if ("Mayne_GA" %in% clocks_to_run) {
    results.ls$Mayne_GA <- Mayne_GA(data.m)
  }

  # 6. Return the list of results
  return(results.ls)
}






