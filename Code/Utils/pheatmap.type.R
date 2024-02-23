# taken from Pigengene::pheatmap.type, https://doi.org/doi:10.18129/B9.bioc.Pigengene
# Fixed: the clustering was done ignoring scale, i.e. assuming scale="none", which lead to wrong results
pheatmap.type <- function (Data, annRow, type = colnames(annRow)[1], doTranspose = FALSE, 
                           conditions = "Auto", scale = scale, ...) 
{
  res <- list()
  annRow <- annRow[, type, drop = FALSE]
  if (is.null(rownames(annRow))) 
    stop("annrow must have row names!")
  if (any(!rownames(annRow) %in% rownames(Data))) 
    stop("annRow has rows that are not present in rows of Data!")
  annRow <- annRow[order(annRow[, 1]), , drop = FALSE]
  samplesOriginalOrder <- rownames(Data)
  if (conditions[1] == "Auto") 
    conditions <- unique(as.character(annRow[, 1]))
  if (any(!conditions %in% unique(as.character(annRow[, 1])))) {
    warning("Some of the conditions are not in annRow.")
    conditions <- intersect(conditions, unique(as.character(annRow[, 1])))
  }
  pheatmapS <- list()
  dataPlot <- c()
  ann1 <- c()
  for (cond in conditions) {
    condSamples <- rownames(annRow)[which(annRow == cond)]
    if (length(condSamples) > 1) {
      pa <- pheatmap(Data[condSamples, , drop = FALSE], 
                     cluster_cols = FALSE, silent = TRUE, scale = scale)
      pheatmapS[[as.character(cond)]] <- pa
      o2 <- pa$tree_row$order
    }
    else {
      o2 <- 1
    }
    dataPlot <- rbind(dataPlot, Data[condSamples[o2], , drop = FALSE])
    ann1 <- rbind(ann1, annRow[condSamples[o2], , drop = FALSE])
  }
  if (!doTranspose) {
    pAll <- pheatmap(dataPlot, annotation_row = ann1, cluster_rows = FALSE, scale = scale, 
                     ...)
  }
  else {
    pAll <- pheatmap(t(dataPlot), annotation_col = ann1, 
                     cluster_cols = FALSE, scale = scale, ...)
  }
  res[["pheatmapS"]] <- pheatmapS
  res[["pheat"]] <- pAll
  res[["ordering"]] <- match(rownames(dataPlot), samplesOriginalOrder)
  res[["annRowAll"]] <- ann1
  invisible(res)
}