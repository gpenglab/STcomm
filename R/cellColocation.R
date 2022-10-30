#' Get the co-localized cell type pairs according to the weights for each pixel
#' after RCTD pipeline
#'
#' @param rctd An RCTD object or a data.frame containing the confident results
#' of the RCTD algorithm.
#' @param rctd_mode  RCTD analysis, such as in 'multi' mode, 'doublet' mode,
#' or 'full' mode.
#' @param method The method to identify co-localized celltypes, including 'pcc'
#' (calculating Pearson Correlation Coefficient) and 'jac'
#' (calculating Jaccard similarity coefficient)
#' @param topn Threshold of the number of cell-type pairs to keep according to
#' the method of identifying co-localized cell type pairs
#' @param pcc Threshold of the Pearson Correlation Coefficient for identifying
#' significant co-localized cell type pairs
#' @param pval Threshold of the p-value for identifying significant
#' co-localized cell type pairs. Only for 'pcc' method
#' @param padj Threshold of the ajusted p-value for identifying significant
#' co-localized cell type pairs. Only for 'pcc' method
#' @param jac Threshold of the Jaccard similarity coefficient for identifying
#' significant co-localized cell type pairs
#' @importFrom Hmisc rcorr
#' @importFrom dplyr mutate group_by arrange filter top_n ungroup desc distinct
#' @importFrom tidyr gather
#'
#' @return a data.frame contains co-localized cell-type pairs
#' @export
cellColocation <- function(rctd,
                           rctd_mode = "multi",
                           method = c("pcc", "jac"),
                           topn = 10,
                           pcc = 0.05,
                           pval = 0.05,
                           padj = 0.05,
                           jac = 0.05) {
  # data frame as input  # rctd object as  input
  if(is(rctd,'RCTD') & !is.null(rctd_mode)){
    message("Input a RCTD weights from a RCTD object")
    rctd_df <- rctd4weight(rctd,
                           rctd_mode = rctd_mode,
                           conf = TRUE)
  }else{
    message("Input a RCTD weights from a data frame")
    rctd_df <- as.data.frame(rctd)
  }

  sort_repair <- function(row1) {
    cells <- sort(row1[1:2])
    repair <- paste(cells[1], cells[2], sep = "_")
    return(repair)
  }

  method <- match.arg(method)
  if (method == "pcc") {
    message("Identify co-localized celltypes by calculating Pearson Correlation Coefficient")
    pearson_cor <- Hmisc::rcorr(as.matrix(rctd_df), type = "pearson")
    pearson_cor_df <- as.data.frame(pearson_cor$r)
    pearson_cor_df$cellTYpe_1 <- rownames(pearson_cor_df)
    pearson_cor_long <- pearson_cor_df %>%
      tidyr::gather(cellTYpe_2, PCC, -cellTYpe_1) %>%
      dplyr::mutate("cellTYpe_pair" = paste(cellTYpe_1, cellTYpe_2, sep = "_"))


    pvalue_cor_df <- as.data.frame(pearson_cor$P)
    pvalue_cor_df$cellTYpe_1 <- rownames(pvalue_cor_df)
    pvalue_cor_long <- pvalue_cor_df %>%
      tidyr::gather(cellTYpe_2, Pval, -cellTYpe_1) %>%
      dplyr::mutate("cellTYpe_pair" = paste(cellTYpe_1, cellTYpe_2, sep = "_"))
    pearson_cor_long_pval <- data.frame(pearson_cor_long,
                                        Pval = pvalue_cor_long$Pval)

    pval_adj_cor_df <- as.data.frame(apply(pearson_cor$P, 1,
                                           function(Pval) p.adjust(Pval, method = "BH")))
    pval_adj_cor_df$cellTYpe_1 <- rownames(pval_adj_cor_df)
    pval_adj_cor_long <- pval_adj_cor_df %>%
      tidyr::gather(cellTYpe_2, P.adj, -cellTYpe_1) %>%
      dplyr::mutate("cellTYpe_pair" = paste(cellTYpe_1, cellTYpe_2, sep = "_"))
    pearson_cor_long_pval_adj <- data.frame(pearson_cor_long_pval,
                                            P.adj = pval_adj_cor_long$P.adj)

    # filter the pairs of celltype with pcc, pvalue and ajusted pvalue thresholds.
    topn_celltypePair_padj <- pearson_cor_long_pval_adj %>%
      dplyr::group_by(cellTYpe_1) %>%
      dplyr::arrange(desc(PCC)) %>%
      dplyr::top_n(topn, PCC) %>%
      dplyr::filter(PCC != 1, PCC > pcc, Pval < pval, P.adj < padj) %>%
      dplyr::arrange(cellTYpe_1) %>%
      dplyr::ungroup()

    # filter the repetitive pairs of celltype
    topn_celltypePair_padj$reorder_pair <- apply(topn_celltypePair_padj, 1,
                                                 FUN = sort_repair)
    topn_celltypePair_padj <- dplyr::distinct(topn_celltypePair_padj,
                                              reorder_pair, .keep_all = TRUE)
    cells_colocation <- topn_celltypePair_padj
  } else if (method == "jac") {
    message("Identify co-localized celltypes by calculating Jaccard Similarity Coefficient")
    # binary the matrix of rctd
    rctd_mtx <- as.matrix(rctd_df)
    for (id in seq_len(nrow(rctd_mtx))) {
      rctd_mtx[id, which(rctd_mtx[id, ] > 0)] <- 1
    }
    # calculate Jaccard Index
    jaccard <- function(a, b) {
      intersection <- length(intersect(a, b))
      unions <- length(a) + length(b) - intersection
      s <- intersection / unions
      return(s)
    }

    jac_dist <- data.frame(
      cellType_1 = "cellType_1",
      cellType_2 = "cellType_2",
      celltype_pairs = "celltype_pairs",
      Jac = 0
    )
    celltypes <- colnames(rctd_mtx)
    for (ct in celltypes) {
      for (id in celltypes) {
        Jac <- jaccard(which(rctd_mtx[, ct] != 0),
                       which(as.matrix(rctd_mtx)[, id] != 0))
        cellType_1 <- ct
        cellType_2 <- id
        celltype_pairs <- paste(cellType_1, cellType_2, sep = "_")
        df <- data.frame(cellType_1, cellType_2, celltype_pairs, Jac)
        jac_dist <- rbind(jac_dist, df)
      }
    }
    jac_dist <- jac_dist[-1, ]

    # filter with jac threshold
    sub_jac_dist <- jac_dist %>%
      dplyr::group_by(cellType_1) %>%
      dplyr::arrange(dplyr::desc(Jac)) %>%
      dplyr::filter(Jac != 1, Jac > jac) %>%
      dplyr::ungroup()

    sub_jac_dist$reorder_pair <- apply(sub_jac_dist, 1, FUN = sort_repair)
    sub_jac_dist <- dplyr::distinct(sub_jac_dist, reorder_pair, .keep_all = TRUE)
    cells_colocation <- sub_jac_dist
  }

  return(cells_colocation)
}
