#' @title Calculate GrimAge1
#'
#' @description
#' Calculates DNA methylation GrimAge1, a composite biomarker of mortality risk and biological aging.
#'
#'@details
#' This function calculates DNAm GrimAge1 in a multi-step process. First, it predicts DNAm-based surrogate biomarkers for several plasma proteins from the input beta values. These predicted biomarkers, along with chronological age and sex, are then used to calculate a composite mortality risk score. This score is calibrated to the scale of chronological age to produce the final `DNAmGrimAge2`.
#'
#' @param beta_matrix A numeric matrix of DNA methylation beta values. Rows should
#'   represent CpG sites and columns should represent individual samples.
#' @param age A numeric vector of chronological ages for the samples corresponding
#'   to the columns in `beta_matrix`.
#' @param sex A character vector of sample sexes. Must contain "Male" or "Female"
#'   for each sample.
#'
#' @return
#' A data.frame containing the following columns:
#' \itemize{
#'   \item `Sample`: Identifier for each sample.
#'   \item `Age`: The input chronological age.
#'   \item `Female`: A numeric indicator for sex (1 = Female, 0 = Male).
#'   \item `DNAm...`: Columns for each of the predicted surrogate biomarkers (e.g., `DNAmADM`, `DNAmGDF15`).
#'   \item `DNAmGrimAge1`: The final calibrated GrimAge1 score.
#' }
#'
#'
#'
#' @references
#' Lu AT, Quach A, Wilson JG, et al.
#' DNA methylation GrimAge strongly predicts lifespan and healthspan
#' \emph{Aging} 2019
#'
#' @export
#'
#' @examples
#' download_OmniAgeR_example("Hannum_example")
#' load_OmniAgeR_example("Hannum_example")
#' age <- PhenoTypesHannum_lv$Age
#' sex <- ifelse(PhenoTypesHannum_lv$Sex=="F","Female","Male")
#' GrimAge1.out <- GrimAge1(beta_matrix = hannum_bmiq_m,age,sex)



GrimAge1 <- function(beta_matrix,age,sex) {
  # Load the model file
  data("GrimAge1CpG")
  cpgs <- grimage1[[1]]          # CpG sites and coefficients
  glmnet.final1 <- grimage1[[2]] # Final model parameters
  gold <- grimage1[[3]]          # Calibration parameters

  # --- CpG Coverage Calculation ---
  # 1. Get the list of all CpG probes required by the model (excluding non-CpG variables)
  required_cpgs_list <- setdiff(cpgs$var, c("Intercept", "Age"))
  n_required_cpgs <- length(required_cpgs_list)

  available_cpgs_in_data <- intersect(required_cpgs_list, rownames(beta_matrix))
  n_available_cpgs <- length(available_cpgs_in_data)
  print(paste0("[GrimAge1] Number of represented GrimAge1 CpGs (max=",n_required_cpgs,")=",
               n_available_cpgs))

  # Transpose beta_matrix so that samples are rows and CpGs are columns
  dat_meth <- as.data.frame(t(beta_matrix))

  # Handle sex: Convert 'sex' to 'Female' (1 for female, 0 for male)
  Female <- ifelse(sex == "Female", 1, 0)
  sample_info <- data.frame(Sample=colnames(beta_matrix),Age=age,Female=Female)
  # Merge sample information
  dat_meth <- cbind(sample_info,dat_meth)
  #dat_meth <- cbind(dat_meth,sample_info)
  # Add 'Intercept' column
  dat_meth$Intercept <- 1

  # Filter for CpG sites required by the model that are present in beta_matrix
  available_cpgs <- intersect(cpgs$var, colnames(dat_meth))
  cpgs <- subset(cpgs, var %in% c(available_cpgs,"Intercept","Age"))

  # Generate DNAm protein predictors
  Ys <- unique(cpgs$Y.pred)
  for (k in 1:length(Ys)) {
    cpgs1 <- subset(cpgs, Y.pred == Ys[k])
    if (nrow(cpgs1) > 0) {
      Xs <- dat_meth[, cpgs1$var, drop = FALSE]
      Y.pred <- as.numeric(as.matrix(Xs) %*% cpgs1$beta)
      dat_meth[, Ys[k]] <- Y.pred
    }
  }

  # Filter variables in glmnet.final1 to keep only those present in dat_meth
  available_vars <- intersect(glmnet.final1$var, colnames(dat_meth))
  glmnet.final1 <- subset(glmnet.final1, var %in% available_vars)

  # Generate the raw GrimAge2 variable ('COX')
  if (length(available_vars) > 0) {
    output.all <- dat_meth[, c('Sample', 'Age', 'Female', Ys)]
    output.all$COX <- as.numeric(as.matrix(output.all[, available_vars, drop = FALSE]) %*% glmnet.final1$beta)
  } else {
    stop("[GrimAge1] No available variables to calculate 'COX'")
  }

  # Calibrate 'COX' to 'DNAmGrimAge2'
  F_scale <- function(INPUT0, Y.pred0.name, Y.pred.name, gold) {
    out.para <- subset(gold, var == 'COX')
    out.para.age <- subset(gold, var == 'Age')
    m.age <- out.para.age$mean
    sd.age <- out.para.age$sd
    Y0 <- INPUT0[, Y.pred0.name]
    Y <- (Y0 - out.para$mean) / out.para$sd
    INPUT0[, Y.pred.name] <- as.numeric((Y * sd.age) + m.age)
    return(INPUT0)
  }
  output.all <- F_scale(output.all, 'COX', 'DNAmGrimAge1', gold)

  # Calculate age acceleration 'AgeAccelGrim1'
  #output.all$DNAmtemp <- output.all[, 'DNAmGrimAge1']
  #output.all[, 'AgeAccelGrim2'] <- residuals(lm(DNAmtemp ~ Age, data = output.all, na.action = na.exclude))
  #output.all$DNAmtemp <- NULL
  output.all$COX <- NULL



  return(output.all)
}













