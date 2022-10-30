#' Tidy RCTD results
#'
#' @param rctd an RCTD object containing the results of the RCTD algorithm.
#' documentation for more information on interpreting the content of the
#' RCTD object
#' @param rctd_mode the mode of RCTD analysis, such as in 'multi' mode,
#' 'doublet' mode, or 'full' mode
#' @param conf whether only return the confident weights when user choose '
#' multi' mode
#'
#' @return a data.frame of cell type weights for every pixel
#' @export
rctd4weight <- function(rctd,
                        rctd_mode = c("multi", "full", "doublet"),
                        conf = TRUE) {
  rctd_mode <- match.arg(rctd_mode)
  if (!is.null(rctd) & !is.null(rctd_mode)) {
    spot_bc <- colnames(rctd@spatialRNA@counts)
    celltypes2use <- colnames(rctd@cell_type_info$info[[1]])

    if (rctd_mode == "multi") {
      if (conf) {
        nspot <- length(spot_bc)
        ncelltypes <- length(celltypes2use)
        rctd_weights_mtx <- matrix(rep(0, length(nspot * ncelltypes)),
                                   nrow = nspot, ncol = ncelltypes)
        rownames(rctd_weights_mtx) <- spot_bc
        colnames(rctd_weights_mtx) <- celltypes2use
        for (spot in seq_len(nrow(rctd_weights_mtx))) {
          rctd_weights_mtx[spot, match(rctd@results[[spot]]$cell_type_list[rctd@results[[spot]]$conf_list],
                                       colnames(rctd_weights_mtx))] <- rctd@results[[spot]]$sub_weights[rctd@results[[spot]]$conf_list]
        }
        rctd_weights <- data.frame(rctd_weights_mtx[, apply(rctd_weights_mtx, 2,
                                                            sum) > 0])
      }else{
        nspot <- length(spot_bc)
        ncelltypes <- length(celltypes2use)
        rctd_weights_mtx <- matrix(rep(0, length(nspot * ncelltypes)),
                                   nrow = nspot, ncol = ncelltypes)
        rownames(rctd_weights_mtx) <- spot_bc
        colnames(rctd_weights_mtx) <- celltypes2use
        for (spot in seq_len(nrow(rctd_weights_mtx))) {
          rctd_weights_mtx[spot, match(rctd@results[[spot]]$cell_type_list,
                                       colnames(rctd_weights_mtx))] <- rctd@results[[spot]]$sub_weights
        }
        rctd_weights <- data.frame(rctd_weights_mtx[, apply(rctd_weights_mtx, 2,
                                                            sum) > 0])
      }
    } else if (rctd_mode == "doublet") {
      doublet_results <- rctd@results$results_df
      doublet_weights <- as.data.frame(rctd@results$weights_doublet)
      rctd_doublet_weights <- data.frame(
        first_type = doublet_results@first_type,
        first_weight = doublet_weights@first_type,
        second_type = doublet_results@second_type,
        second_weight = doublet_weights@second_type,
        spot_class = doublet_results@spot_class)
      rownames(rctd_doublet_weights) <- rownames(doublet_results)
      doublet_celltype <- stringr::str_sort(union(doublet_results$first_celltype,
                                                  doublet_results$second_celltype))

      rctd_weights <- data.frame(ID = rownames(doublet_results))
      rctd_weights[, doublet_celltype] <- 0
      rownames(rctd_weights) <- rctd_weights$ID
      rctd_weights <- rctd_weights[, -1]
      for (ct in doublet_celltype) {
        rctd_weights[which(rctd_doublet_weights$first_type == ct), ct] <- rctd_doublet_weights[which(RCTD_doublet_weights$first_type == ct), "first_weight"]
        rctd_weights[which(rctd_doublet_weights$second_type == ct), ct] <- rctd_doublet_weights[which(rctd_doublet_weights$second_type == ct), "second_weight"]
      }
    } else if (rctd_mode == "full") {
      rctd_weights <- data.frame(rctd@results[[1]]$all_weights)
      for (id in seq(2:length(rctd@results))) {
        rctd_weights <- cbind(rctd_weights, data.frame(rctd@results[[id]]$all_weights))
      }
    }
  }

  return(rctd_weights)
}
