test_that("epiMarker calculates clocks correctly on Hannum dataset", {
  # Skip test if the data package is not available
  skip_if_not_installed("OmniAgeRData")
  
  # 1. Load the required example data
  hannumExample <- loadOmniAgeRdata(
    "omniager_hannum_example",
    verbose = FALSE
  )
  
  hannumBmiqM <- hannumExample[[1]]
  phenoTypesHannum <- hannumExample[[2]]
  age <- phenoTypesHannum$Age
  sex <- ifelse(phenoTypesHannum$Sex == "F", "Female", "Male")
  
  # 2. Define the exact clocks to be tested
  clockCategories <- c("cellularAging", "chronological", "biological", 
                       "causal", "stochastic", "gestationalAge",
                       "surrogateBiomarkers", "traitPred", "diseaseRisk")
  
  # Get the list of available clocks and format them
  useClocks <- listEpiMarker()
  useClocks$crossSpecies <- NULL
  useClocks <- unname(unlist(useClocks))
  
  # Exclude specific PC clocks and SystemsAge for this test suite
  excludedClocks <- c("PCHorvath2013", "PCHorvath2018", "PCHannum", 
                      "PCPhenoAge", "PCGrimAge1", "PCDNAmTL", "SystemsAge")
  useClocks <- useClocks[!useClocks %in% excludedClocks]
  
  # 3. Execute the main epiMarker function
  hannumEpiAgeRes <- epiMarker(
    betaM = hannumBmiqM,
    clockNames = useClocks,
    chronAge = age,
    sexVec = sex,
    minCoverage = 0,
    verbose = FALSE
  )
  
  # 4. Extract clocks that return simple numeric vectors (non-lists)
  simpleClocks <- hannumEpiAgeRes[!sapply(hannumEpiAgeRes, is.list)]
  
  # 5. Extract clocks that return nested lists or data frames
  listClocks <- hannumEpiAgeRes[sapply(hannumEpiAgeRes, is.list)]
  
  # Reformat nested clock outputs into flat named vectors for comparison
  refNames <- names(listClocks$CentenarianClock$ENCen40)
  
  nestedClocks <- list(
    epiTOC2_irS = listClocks$epiTOC2$irS,
    epiTOC3_irS = listClocks$epiTOC3$irS,
    CentenarianClock_40 = listClocks$CentenarianClock$ENCen40,
    CentenarianClock_100 = listClocks$CentenarianClock$ENCen100,
    
    GrimAge1 = setNames(listClocks$GrimAge1$DNAmGrimAge1, 
                        listClocks$GrimAge1$Sample),
    GrimAge2 = setNames(listClocks$GrimAge2$DNAmGrimAge2, 
                        listClocks$GrimAge2$Sample),
    DNAmFitAge = setNames(listClocks$DNAmFitAge$DNAmFitAge, 
                          listClocks$DNAmFitAge$SampleID),
    
    StocH = setNames(listClocks$StochClocks$StocH, refNames),
    StocZ = setNames(listClocks$StochClocks$StocZ, refNames),
    StocP = setNames(listClocks$StochClocks$StocP, refNames),
    
    CausalAge = setNames(listClocks$CausalClock$CausalAge, refNames),
    DamAge = setNames(listClocks$CausalClock$DamAge, 
                      listClocks$PCClocks$SampleID),
    AdaptAge = setNames(listClocks$CausalClock$AdaptAge, 
                        listClocks$SystemsAge$SampleID),
    
    Retro_age_v1 = listClocks$Retro_age$retroAgeV1,
    Retro_age_v2 = listClocks$Retro_age$retroAgeV2
  )
  
  # Append the remaining list items
  nestedClocks <- c(nestedClocks, listClocks$LeeGA, 
                    listClocks$CRP, listClocks$CHIP)
  
  # Combine simple and processed nested results
  finaleHannumResList <- c(simpleClocks, nestedClocks) 
  
  # 6. Load the expected results from the fixtures directory
  expectedRes <- readRDS(test_path("fixtures", "expected_result.rds"))
  
  # 7. Assert equality for each clock result
  for (i in seq_along(finaleHannumResList)) {
    tmpName <- names(finaleHannumResList)[i]
    expect_equal(unname(finaleHannumResList[[tmpName]]), 
                 unname(expectedRes[[tmpName]]), 
                 tolerance = 1e-5)
  }
})
