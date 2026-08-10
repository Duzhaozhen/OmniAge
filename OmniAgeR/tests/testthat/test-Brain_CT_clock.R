test_that("brain cell-type pipeline resolves model and imputation siblings", {
  
  # Create a minimal 5-fold BrainCT model
  foldModel <- data.frame(
    feature_name = c("intercept", "GENE1"),
    coefficient = c(1, 2)
  )
  
  folds <- setNames(
    replicate(5, foldModel, simplify = FALSE),
    paste0("fold_", seq_len(5))
  )
  
  # Mimic the structure of the real BrainCT model object
  modelData <- list(
    brain_ct_clocks_coef = list(
      SC_Excitatory_Neurons = folds
    ),
    Brain_CT_imputation_data_list = list(
      SC_Excitatory_Neurons = data.frame(
        feature_name = "GENE1",
        imputation_value = 3
      )
    )
  )
  
  # Replace external data loading and Seurat preprocessing
  # with small deterministic mock objects
  local_mocked_bindings(
    loadOmniAgeRdata = function(...) modelData,
    getDfSeurat = function(...) {
      data.frame(
        age = c(40, 50),
        donorId = c("donor_1", "donor_2")
      )
    },
    #.env = environment(runPredictionPipelineBrainCt)
  )
  
  result <- runPredictionPipelineBrainCt(
    sampleType = "SC",
    seuratObj = NULL,
    cellTypes = "Excitatory Neurons",
    verbose = FALSE
  )
  
  # GENE1 is absent from the input and should therefore be imputed as 3.
  # Prediction = intercept + coefficient * imputed value
  #            = 1 + 2 * 3
  #            = 7
  expect_equal(nrow(result), 2)
  expect_equal(result$prediction, c(7, 7))
  expect_equal(
    result$celltype,
    rep("Excitatory Neurons", 2)
  )
})


test_that("brain cell-type pipeline supports a single sample", {
  
  # Create the same minimal 5-fold model
  foldModel <- data.frame(
    feature_name = c("intercept", "GENE1"),
    coefficient = c(1, 2)
  )
  
  folds <- setNames(
    replicate(5, foldModel, simplify = FALSE),
    paste0("fold_", seq_len(5))
  )
  
  modelData <- list(
    brain_ct_clocks_coef = list(
      SC_Oligodendrocytes = folds
    ),
    Brain_CT_imputation_data_list = list(
      SC_Oligodendrocytes = data.frame(
        feature_name = "GENE1",
        imputation_value = 3
      )
    )
  )
  
  local_mocked_bindings(
    loadOmniAgeRdata = function(...) modelData,
    getDfSeurat = function(...) {
      data.frame(
        age = 40,
        donorId = "donor_1"
      )
    },
    #.env = environment(runPredictionPipelineBrainCt)
  )
  
  result <- runPredictionPipelineBrainCt(
    sampleType = "SC",
    seuratObj = NULL,
    cellTypes = "Oligodendrocytes",
    verbose = FALSE
  )
  
  expect_equal(nrow(result), 1)
  expect_equal(result$prediction, 7)
  expect_equal(result$age, 40)
  expect_equal(result$donorId, "donor_1")
  expect_equal(result$celltype, "Oligodendrocytes")
})