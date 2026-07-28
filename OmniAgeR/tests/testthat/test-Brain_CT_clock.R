test_that("brain cell-type pipeline resolves model and imputation siblings", {
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
      SC_Excitatory_Neurons = folds
    ),
    Brain_CT_imputation_data_list = list(
      SC_Excitatory_Neurons = data.frame(
        feature_name = "GENE1",
        imputation_value = 3
      )
    )
  )

  local_mocked_bindings(
    loadOmniAgeRdata = function(...) modelData,
    getDfSeurat = function(...) {
      data.frame(
        age = c(40, 50),
        donorId = c("donor_1", "donor_2")
      )
    },
    .env = environment(runPredictionPipelineBrainCt)
  )

  result <- runPredictionPipelineBrainCt(
    sampleType = "SC",
    seuratObj = NULL,
    cellTypes = "Excitatory Neurons",
    verbose = FALSE
  )

  expect_equal(nrow(result), 2)
  expect_equal(result$prediction, c(7, 7))
  expect_equal(result$celltype, "Excitatory Neurons")
})
