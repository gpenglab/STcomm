#' Get the significant cell-cell communication restricted by ligand-receptor
#' co-expressed in co-localized cell type pairs
#'
#' @param object Seurat object
#' @param weights.df a data frame of cell type confident weights for every pixel
#' @param ctpairs a data frame contains co-localized cell-type pairs
#' @param cellchat CellChat object or a net data.frame with single-cell data,
#' which contained the communication probability and infered cellular
#' communication network
#' @param db the species to use for CellChatDB to choose, including 'mouse' and
#' 'human'
#' @param fisher.pavl threshold of the p-value for fisher exact test
#' @param fisher.padj threshold of the ajusted p-value for fisher exact test
#' @importFrom stringr str_to_title
#' @importFrom Seurat GetAssayData
#' @importFrom CellChat subsetCommunication
#' @importFrom dplyr bind_rows mutate
#'
#' @return a data frame contains the significant cell-cell communication
#' restricted by ligand-receptor co-expressed in co-localized cell type pairs
#' @export
st_comm <- function(object,
                    weights.df,
                    ctpairs,
                    cellchat,
                    db = c("mouse", "human"),
                    fisher.pavl = 0.05,
                    fisher.padj = 0.05) {
  if (!is.null(object)) {
    exp_st <- t(as.matrix(Seurat::GetAssayData(object, assay = "Spatial", slot = "data")))
  } else {
    print("It is necessary to input the Seurat object of spatial transcriptome data!")
  }

  if (!is.null(weights.df)) {
    rctd_conf <- weights.df

    # binary the rctd_conf
    rctd_conf_bina <- as.matrix(rctd_conf)
    for (ct in colnames(rctd_conf)) {
      rctd_conf_bina[rctd_conf_bina[, ct] > 0, ct] <- 1
    }
  } else {
    print("It is necessary to input the confident weights matrix or data frame of RCTD analysis!")
  }

  if (!is.null(ctpairs)) {
    colocal_cts <- as.data.frame(ctpairs)
  } else {
    print("It is necessary to input the data.frame of confident co-localized cell type pairs!")
  }

  # get the LR pairs from CellChatDB
  message("Now is getting the LR pairs from CellChatDB.")
  db <- match.arg(db)
  if (is(cellchat, 'CellChat') & !is.null(db)) {
    if(length(cellchat@net) != 0 & length(cellchat@netP) != 0){
      net.df <- CellChat::subsetCommunication(cellchat)
      net.df$ct_pairs <- paste0(net.df$source, "_", net.df$target)

      interaction_input <- cellchat@DB$interaction
    }else{
      print("A CellChat object must computed the communication probability and infered cellular communication network!")
    }
  }else{
    net.df <- as.data.frame(cellchat)
    net.df$ct_pairs <- paste0(net.df$source, "_", net.df$target)

    if (db == "mouse") {
      CellChatDB <- CellChat::CellChatDB.mouse
    } else if (db == "human") {
      CellChatDB <- CellChat::CellChatDB.human
    }
    interaction_input <- CellChatDB$interaction
  }

  lrs <- strsplit(interaction_input$interaction_name, "_")
  for (i in seq_len(length(lrs))) {
    lrs[[i]] <- stringr::str_to_title(unlist(lrs[[i]]))
    if (length(lrs[[i]]) == 3) {
      lrs[[i]] <- data.frame( # ligand = lrs[[i]][1],
        recepter1 = lrs[[i]][2],
        recepter2 = lrs[[i]][3]
      )
    } else {
      lrs[[i]] <- data.frame( # ligand = lrs[[i]][1],
        recepter1 = lrs[[i]][2]
      )
    }
  }
  lrs_df <- dplyr::bind_rows(lrs)
  rownames(lrs_df) <- rownames(interaction_input)
  interaction_input <- cbind(interaction_input, lrs_df)
  interaction_input2 <- interaction_input[, c("interaction_name",
                                              "pathway_name",
                                              "ligand", "recepter1",
                                              "recepter2", "annotation")]
  lr_genes <- unique(stringr::str_to_title(unlist(strsplit(interaction_input$interaction_name, "_"))))

  # fisher exact test
  message("Now is Fisher exact test between the co-localizaed cell type pairs and the co-expressed LR pairs.")
  lrst_df <- data.frame(
    interaction_name = "test",
    pathway_name = "test",
    ligand = "test",
    recepter1 = "test",
    recepter2 = "test",
    annotation = "test",
    interaction_name2 = "test",
    source = "ct1", target = "ct2",
    CT.pair.colocated = 0,
    LR.pair.coexp = 0,
    CT.LR.pair_colocated = 0,
    CT.LR.pair_co.ratio = 0,
    FT.pval = 0
  )
  exp_lrs_lst <- list()
  for (ctp in seq_len(nrow(colocal_cts))) {
    ct1 <- colocal_cts[ctp, 1]
    ct2 <- colocal_cts[ctp, 2]
    tmp_ctp <- rctd_conf_bina[, c(ct1, ct2)]
    ctp_colocal <- data.frame(rctd_conf_bina[, c(ct1, ct2)])[apply(rctd_conf_bina[, c(ct1, ct2)], 1, sum) > 1, ]
    abcd <- nrow(tmp_ctp)
    ab <- nrow(ctp_colocal)
    if (!is.null(ab)) {
      exp_lrst_df <- data.frame(spot_bc = rownames(exp_st))
      rownames(exp_lrst_df) <- rownames(exp_st)
      for (id in seq_len(nrow(interaction_input2))) {
        lr <- stringr::str_to_title(strsplit(interaction_input2$interaction_name, "_")[[id]])
        if (length(intersect(lr, colnames(exp_st))) == length(lr)) {
          lr_df <- interaction_input2[id, ]
          exp_lr <- as.data.frame(exp_st[, lr])
          if (length(lr) == 2) {
            exp_lr$LR <- exp_lr[, 1] * exp_lr[, 2]
            colnames(exp_lr)[3] <- paste(lr, collapse = "_")
          } else {
            exp_lr$LR <- exp_lr[, 1] * exp_lr[, 2] * exp_lr[, 3]
            colnames(exp_lr)[4] <- paste(lr, collapse = "_")
          }
          # LR pair in st
          ac <- length(which(exp_lr[, ncol(exp_lr)] > 0))
          # LR pair in st-colocated CT pairs
          a <- length(which(exp_lr[rownames(ctp_colocal), ncol(exp_lr)] > 0))
          # LR pair not in st-colocated CT pairs
          b <- ab - a
          # LR pair-Non_colocated CT pairs
          c <- ac - a
          # LR pair not in st-Non_colocated CT pairs
          d <- abcd - ab - c
          x <- matrix(c(a, c, b, d), ncol = 2, nrow = 2)
          # fisher test
          ft.pval <- fisher.test(x)$p.value
          lr_per_spot <- a / ab
          LR <- paste(lr, collapse = "_")
          lrs_spot <- data.frame(lr_df,
                                 interaction_name2 = LR, source = ct1, target = ct2,
                                 CT.pair.colocated = ab,
                                 LR.pair.coexp = ac,
                                 CT.LR.pair_colocated = a,
                                 CT.LR.pair_co.ratio = lr_per_spot,
                                 FT.pval = ft.pval
          )
          lrst_df <- rbind(lrst_df, lrs_spot)
          exp_lrst_df <- cbind(exp_lrst_df, exp_lr)
        } else {
          next
        }
      }
    }

    exp_lrst_df <- list(exp_lrst_df)
    names(exp_lrst_df) <- paste0(ct1, "_", ct2)
    exp_lrs_lst <- c(exp_lrs_lst, exp_lrst_df)
  }
  lrst_df <- lrst_df[-1, ]

  # caculate adjusted pvlue
  lrst_df$FT.padj <- p.adjust(lrst_df$FT.pval, "fdr", n = nrow(lrst_df))
  lrst_df$CT.pairs <- paste0(lrst_df$source, "_", lrst_df$target)

  # filter by FT.pval and FT.padj
  lrst_df <- lrst_df[lrst_df$FT.pval < fisher.pavl & lrst_df$FT.padj < fisher.padj, ]
  lrst_df <- lrst_df[lrst_df$CT.LR.pair_co.ratio > 0, ]

  # filter the net of CellChat
  net.df_filter <- net.df[1, ]
  net.df_filter <- net.df_filter %>%
    dplyr::mutate(CT.pair.colocated = 0, LR.pair.coexp = 0,
                  CT.LR.pair_colocated = 0,
                  CT.LR.pair_co.ratio = 0,
                  FT.pval = 0, FT.padj = 0)

  for (ct in intersect(net.df$reoder_pair, lrst_df$CT.pairs)) {
    lr_common <- intersect(
      net.df[net.df$reorder_pair == ct, "interaction_name"],
      lrst_df[lrst_df$CT.pairs == ct, "interaction_name"]
    )

    for (interaction in lr_common) {
      tmp.st_chat <- lrst_df[lrst_df$CT.pairs == ct & lrst_df$interaction_name == interaction, ]
      tmp <- net.df[net.df$reorder_pair == ct & net.df$interaction_name == interaction, ]
      tmp <- tmp %>% dplyr::mutate(
        CT.pair.colocated = tmp.st_chat$CT.pair.colocated,
        LR.pair.coexp = tmp.st_chat$LR.pair.coexp,
        CT.LR.pair_colocated = tmp.st_chat$CT.LR.pair_colocated,
        CT.LR.pair_co.ratio = tmp.st_chat$CT.LR.pair_co.ratio,
        FT.pval = tmp.st_chat$FT.pval,
        FT.padj = tmp.st_chat$FT.padj
      )
      net.df_filter <- rbind(net.df_filter, tmp)
    }
  }

  return(net.df_filter)
}
