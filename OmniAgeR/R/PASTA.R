#' @title Calculate PASTA-based transcriptomic age scores
#'
#' @description
#' This function computes transcriptomic age acceleration scores based on the
#' models described by Salignon et al. (2025), including the primary PASTA
#' score, a standard regression (REG) score, and the CT46 score.
#'
#' @param mat An expression matrix (numeric) with **genes as rows** and
#'   **samples as columns**.
#' @param filter_genes Logical. If \code{TRUE} (default), the matrix is
#'   subsetted to retain only the genes utilized by the pre-trained models.
#' @param rank_norm Logical. If \code{TRUE} (default), applies a rank-based
#'   inverse normal transformation (rank-normalization) to the expression data.
#' @param REG Logical. If \code{TRUE} (default), computes the REG (regression)
#'   age score.
#' @param PASTA Logical. If \code{TRUE} (default), computes the PASTA age score.
#' @param CT46 Logical. If \code{TRUE} (default), computes the CT46 age score.
#'
#' @return
#' A \code{list} where each element is a numeric vector of predicted age scores
#' for the requested model(s) (e.g., \code{res_list$PASTA}, \code{res_list$REG}).
#'
#' @export
#'
#' @references
#' Salignon, J., Tsiokou, M., Marqués, P. et al.
#' Pasta, an age-shift transcriptomic clock, maps the chemical and genetic determinants of aging and rejuvenation.
#' \emph{bioRxiv.} 2025
#'
#' @examples
#' library(magrittr)
#' library(Seurat)
#' library(glmnet)
#' download_OmniAgeR_example("seu_gabitto_2024_filtered")
#' load_OmniAgeR_example("seu_gabitto_2024_filtered")
#'
#' seu$age <- seu$development_stage %>%
#'  gsub("-year.*", "", .) %>%
#'  gsub("-", " ", .) %>%
#'  gsub("80 year old and over stage", "85", .)
#'
#' set.seed(42)
#' seu_bulk <- making_pseudobulks_from_seurat(seu,pool_by = c("cell_type", "age"), chunk_size = 512, verbose = FALSE)
#'
#' lognorm_matrix <- GetAssayData(seu_bulk, assay = "RNA", layer = "data")
#' lognorm_matrix <- as.matrix(lognorm_matrix)
#'
#' seu_bulk_meta <- seu_bulk[[c('chunk_size', 'cell_type', 'age')]]
#' seu_bulk_meta$age <- as.numeric(seu_bulk_meta$age)
#' PASTA_res <- PASTA_Scores(lognorm_matrix, filter_genes = TRUE, rank_norm = TRUE)
#'


PASTA_Scores <- function(mat, filter_genes = TRUE,
                         rank_norm = TRUE, REG = TRUE,
                         PASTA = TRUE, CT46 = TRUE) {

  res_list <- list()
  data("PASTA_Gene")

  if (filter_genes) mat <- filtering_age_model_genes(mat,v_genes_model)
  if (rank_norm)    mat <- applying_rank_normalization(mat)
  mat_t<- t(mat)
  if (REG)   res_list[["REG"]] <- predicting_age_score(mat_t, model_type = 'REG')
  if (PASTA) res_list[["PASTA"]] <- predicting_age_score(mat_t, model_type = 'PASTA')
  if (CT46) res_list[["CT46"]] <- predicting_age_score(mat_t, model_type = 'CT46')
  return(res_list)
}








#' Filter Age Model Genes from Count Matrix
#'
#' Subsets the count matrix to include only genes used in the age prediction model.
#'
#' @param mat Matrix. Count matrix.
#' @return Matrix. Filtered count matrix with median imputation.
#' @export
filtering_age_model_genes <- function(mat,v_genes_model) {
  #data(v_genes_model, envir = environment())
  mat <- mat[match(v_genes_model, rownames(mat)), ]
  median_value <- stats::median(c(mat), na.rm = TRUE)
  mat[is.na(mat)] <- median_value
  rownames(mat) <- v_genes_model
  return(mat)
}


#' Apply Rank Normalization to Matrix
#'
#' Applies rank normalization across each column of the matrix.
#'
#' @param mat Matrix. Count matrix.
#' @return Matrix. Rank-normalized matrix.
#' @export
applying_rank_normalization <- function(mat) {
  mat <- apply(mat, 2, rank, ties.method = 'average')
  return(mat)
}




#' Predict Age Score from Gene Expression Matrix
#'
#' Uses pre-trained models to predict age scores based on gene expression.
#'
#' @param mat Matrix. Processed count matrix.
#' @param model_type Character. Model type ('PASTA', 'REG', or 'CT46').
#' @return Numeric vector. Predicted age scores.
#' @export

predicting_age_score <- function(mat, model_type = 'PASTA') {
  requireNamespace('glmnet', quietly = TRUE)
  #requireNamespace('dplyr', quietly = TRUE)
  # data(cvfit_REG, envir = environment())
  # data(cvfit_PASTA, envir = environment())
  # data(beta_PASTA, envir = environment())
  # data(cvfit_C46, envir = environment())
  # data(beta_C46, envir = environment())
  data("cvfit_REG")
  data("cvfit_PASTA")
  data("beta_PASTA")
  data("cvfit_C46")
  data("beta_C46")


  if (model_type == 'PASTA') cur_model <- cvfit_PASTA
  if (model_type == 'REG')   cur_model <- cvfit_REG
  if (model_type == 'CT46')  cur_model <- cvfit_C46
  if (!model_type %in% c('PASTA', 'REG', 'CT46'))
    stop('Specify a valid model; either PASTA, REG, or CT46')

  v_age_scores <- as.numeric(stats::predict(cur_model, mat, s = 'lambda.min',
                                            type = 'link')[, 1])
  if (model_type == 'PASTA') v_age_scores <- v_age_scores * beta_PASTA
  if (model_type == 'CT46')  v_age_scores <- v_age_scores * beta_C46

  return(v_age_scores)
}


#' Create Pseudobulk Samples from Seurat Object
#'
#' Aggregates single-cell expression data into pseudobulk samples based on
#' user-defined metadata variables and a specified chunk size.
#'
#' @param seu A Seurat object.
#' @param pool_by A character vector of column names from `@meta.data`. Cells are
#'   grouped by unique combinations of these variables prior to chunking.
#'   Default: `c("cell_type", "age")`.
#' @param chunk_size A numeric value. The target number of cells per pseudobulk
#'   sample. If set to 1, no aggregation is performed and the original
#'   object is returned (with metadata updated). Default: `1000`.
#' @param verbose Logical. If `TRUE`, prints a summary table detailing the
#'   number of pseudobulk samples generated per group. Default: `TRUE`.
#'
#' @return A new Seurat object where columns represent pseudobulk samples.
#'   The `@meta.data` slot includes the original `pool_by` variables and the
#'   `chunk_size` used for aggregation.
#' @export
#'
making_pseudobulks_from_seurat <- function(seu,
                                           pool_by = c("cell_type", "age"),
                                           chunk_size = 1000,
                                           verbose = TRUE) {

  # --- 1. Validate Dependencies and Inputs ---
  if (!require('Seurat', quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }

  if (!all(pool_by %in% colnames(seu@meta.data))) {
    stop("Error: Not all variables in 'pool_by' found in Seurat metadata.")
  }

  # --- 2. Handle 'chunk_size == 1' (no aggregation) ---
  if (chunk_size == 1) {
    seu_bulk <- seu
    # Ensure metadata is still present even if no aggregation is done
    seu_bulk$chunk_size <- chunk_size
    for (col in pool_by) {
      if (!col %in% colnames(seu_bulk@meta.data)) {
        seu_bulk[[col]] <- seu[[col, drop=TRUE]]
      }
    }
    return(seu_bulk)
  }

  # --- 3. Core Chunking Logic (Base R replacement for data.table) ---

  # Extract metadata as a data.frame
  pdata1 <- seu[[pool_by]]

  # Create a single grouping factor (equivalent to data.table's 'by = pool_by')
  # interaction() creates unique IDs for "T-cell.50", "T-cell.60", etc.
  group_factor <- interaction(pdata1[, pool_by], drop = TRUE)

  # Use ave() to perform chunking and random sampling independently within each group
  # This is the Base R equivalent of data.table's '[, chunk := sample(...), by = ...]'
  pdata1$chunk <- ave(
    seq_len(nrow(pdata1)), # Pass a simple numeric vector
    group_factor,         # Group by this factor
    FUN = function(idx_in_group) {
      # 'idx_in_group' is a vector of row indices for the current group
      n_in_group <- length(idx_in_group)
      num_chunks <- ceiling(n_in_group / chunk_size)

      # Create the vector of chunk IDs for this group
      chunk_ids_for_group <- rep(
        1:num_chunks,
        each = chunk_size,
        length.out = n_in_group
      )

      # Shuffle and return the chunk IDs
      sample(chunk_ids_for_group)
    }
  )

  # Create unique chunk_id (Base R replacement for .SDcols and apply)
  # Paste all 'pool_by' columns (e.g., "L5.neuron-85")
  temp_id <- do.call(paste, c(pdata1[, pool_by], sep = "-"))
  # Add the chunk number (e.g., "L5.neuron-85-1")
  temp_id_with_chunk <- paste(temp_id, pdata1$chunk, sep = "-")

  # Create R-syntactically-valid names (e.g., "L5.neuron_85-1" -> "L5.neuron_85.1")
  pdata1$chunk_id <- make.names(temp_id_with_chunk)

  # Add the new chunk_id back to the Seurat object
  seu$chunk_id <- pdata1$chunk_id

  # --- 4. Perform Aggregation ---
  seu_bulk <- Seurat::AggregateExpression(
    seu,
    group.by = "chunk_id",
    return.seurat = TRUE,
    verbose = FALSE
  )

  # --- 5. Manually Restore Metadata (Base R replacement for data.table) ---

  # 1. Create a mapping data.frame from chunk_id to all pool_by variables
  #    (Using !duplicated is the Base R equivalent of data.table's 'unique(..., by = ...)')
  meta_map <- pdata1[!duplicated(pdata1$chunk_id), c("chunk_id", pool_by)]

  # 2. Get the metadata from the newly created pseudobulk object
  meta_new <- seu_bulk@meta.data

  # 3. Re-order the meta_map to match the exact cell order in seu_bulk
  ordered_meta_map <- meta_map[match(rownames(meta_new), meta_map$chunk_id), ]

  # 4. Add the metadata columns back to the seu_bulk object
  for (col in pool_by) {
    seu_bulk[[col]] <- ordered_meta_map[[col]]
  }

  # --- 6. Final Cleanup and Return ---
  seu_bulk$chunk_size <- chunk_size

  if (verbose) {
    print("Summary of created pseudobulk samples:")
    # Use two 'aggregate' calls to replicate the data.table summary
    # Step A: Count cells per chunk_id
    cells_per_chunk_df <- aggregate(
      list(cell_count = pdata1$chunk_id),
      by = pdata1[, c(pool_by, "chunk_id")],
      FUN = length
    )
    # Step B: Count chunks per group
    summary_table <- aggregate(
      list(chunk_count = cells_per_chunk_df$chunk_id),
      by = cells_per_chunk_df[, pool_by],
      FUN = length
    )
    print(summary_table)
  }

  return(seu_bulk)
}










