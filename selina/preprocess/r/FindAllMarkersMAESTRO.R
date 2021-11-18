FindAllMarkersMAESTRO <- function(object, test.use = 'presto',
                                  features = NULL, min.pct = 0.1, logfc.threshold = 0.25,
                                  only.pos = FALSE, verbose = TRUE, return.thresh = 1e-2,
                                  min.cells.feature = 3, min.cells.group = 3,
                                  slot = "data", latent.vars = NULL) {
  ident.all <- sort(unique(Idents(object = object)))
  if (test.use != 'presto') {
    all.res <- lapply(0:(length(ident.all) - 1), function(x) {
      if (verbose) {
        message("Calculating cluster ", x)
      }
      test.res <- tryCatch(expr = FindMarkersMAESTRO(object, ident.1 = x, ident.2 = NULL, test.use = test.use,
                                                      features = features, min.pct = min.pct, logfc.threshold = logfc.threshold,
                                                      only.pos = only.pos, verbose = verbose),
                            error = function(cond) { return(cond$message) })
      if (class(test.res) == "character") {
        message(paste0(test.res, ", no differential features identified"))
        test.res <- NULL
      } else {
        test.res$cluster <- as.character(x)
        test.res$gene <- rownames(test.res)
      }
      return(test.res)
    })
    all.res.df <- do.call("rbind", all.res)
    all.res.df <- all.res.df[which(all.res.df$p_val < return.thresh),]
    rownames(all.res.df) = NULL
  } else {
    all.res.df = FindAllMarkersPresto(object, features = features, min.pct = min.pct,
                                      logfc.threshold = logfc.threshold,
                                      only.pos = only.pos, return.thresh = return.thresh,
                                      slot = slot)
  }
  return(all.res.df)
}

FindAllMarkersPresto <- function(object, features = NULL, min.pct = 0.1, logfc.threshold = 0.25,
                                 only.pos = FALSE, return.thresh = 1e-2,
                                 slot = "data") {
  features = if (is.null(features)) { rownames(object) } else { features }
  X_matrix = GetAssayData(object, slot = slot)
  y <- Idents(object) %>% unlist %>% as.character()
  test.res = wilcoxauc(X_matrix, y)

  # Calculate logFC
  y = factor(y)
  if (slot != "scale.data") {
    if (slot == "data") {
      X = expm1(X_matrix)
    }
    group_sums = sumGroups(X, y, 1)
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
      log(group_means[, g] + 1) - log(((cs - group_sums[g,]) / (length(y) - gs[g])) + 1)
    }))
  } else {
    group_sums = sumGroups(X, y, 1)
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
      group_means[, g] - ((cs - group_sums[g,]) / (length(y) - gs[g]))
    }))
  }

  test.res$avg_logFC = as.vector(lfc)
  res = test.res[, c("pval", "avg_logFC", "pct_in", "pct_out", "padj", "group", "feature")]
  res[, c("pct_in", "pct_out")] = round(res[, c("pct_in", "pct_out")] / 100, digits = 3)
  colnames(res) = c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
  res <- res %>% dplyr::filter(.data$p_val < return.thresh &
                         abs(.data$avg_logFC) > logfc.threshold &
                         (.data$pct.1 > min.pct |
                         .data$pct.2 > min.pct) &
                         .data$gene %in% features)
  if (only.pos) {
    res <- res %>% dplyr::filter(.data$avg_logFC > 0)
  }
  res <- res %>% dplyr::arrange(.data$cluster, .data$p_val, desc(.data$avg_logFC))
  return(res)
}