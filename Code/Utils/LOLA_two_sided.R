require(reshape2)
require(data.table)
require(LOLA)

# this is code taken from the LOLA package and changed to compute two-sided fisher tests for each database entry instead of one-sided ones.
# This function is based on older code that wrongly reports the OddsRatio column to be be logOddsRatio
# This problem is currently fixed in downstream code. 
runLOLA <- function (userSets, userUniverse, regionDB, minOverlap = 1, cores = 1, redefineUserSets = FALSE) 
{
  support = d = b = userSet = pValueLog = rnkSup = rnkPV = rnkLO = NULL
  logOddsRatio = maxRnk = meanRnk = dbSet = description = NULL
  annotationDT = regionDB$regionAnno
  testSetsGRL = regionDB$regionGRL
  annotationDT[, `:=`(dbSet, seq_len(nrow(annotationDT)))]
  setkey(annotationDT, dbSet)
  userSets = LOLA:::listToGRangesList(userSets)
  testSetsGRL = LOLA:::listToGRangesList(testSetsGRL)
  LOLA:::setLapplyAlias(cores)
  if (any(is.null(names(testSetsGRL)))) {
    names(testSetsGRL) = seq_along(testSetsGRL)
  }
  if (redefineUserSets) {
    userSets = redefineUserSets(userSets, userUniverse, cores = cores)
    userSets = LOLA:::listToGRangesList(userSets)
  }
  userSetsLength = unlist(LOLA:::lapplyAlias(as.list(userSets), length))
  if (!any(isDisjoint(userSets))) {
    message("You have non-disjoint userSets.")
  }
  message("Calculating unit set overlaps...")
  geneSetDatabaseOverlap = LOLA:::lapplyAlias(as.list(userSets), LOLA:::countOverlapsRev, 
                                              testSetsGRL, minoverlap = minOverlap)
  olmat = do.call(cbind, geneSetDatabaseOverlap)
  message("Calculating universe set overlaps...")
  testSetsOverlapUniverse = countOverlaps(testSetsGRL, userUniverse, 
                                          minoverlap = minOverlap)
  universeLength = length(userUniverse)
  scoreTable = data.table(reshape2::melt(t(olmat), variable.factor = FALSE))
  setnames(scoreTable, c("Var1", "Var2", "value"), c("userSet", 
                                                     "dbSet", "support"))
  if ("factor" %in% class(scoreTable[, userSet])) {
    scoreTable$userSet = as.character(scoreTable$userSet)
  }
  message("Calculating Fisher scores...")
  scoreTable[, `:=`(c("b", "c"), list(b = testSetsOverlapUniverse[match(dbSet, 
                names(testSetsOverlapUniverse))] - support, c = userSetsLength - support))]
  scoreTable[, `:=`(d, universeLength - support - b - c)]
  if (scoreTable[, any(b < 0)]) {
    warning(cleanws("Negative b entry in table. This means either: 1) Your user sets\n\t\tcontain items outside your universe; or 2) your universe has a region that\n\t\toverlaps multiple user set regions, interfering with the universe set overlap\n\t\tcalculation."))
    return(scoreTable)
  }
  if (scoreTable[, any(c < 0)]) {
    warning("Negative c entry in table. Bug with userSetsLength; this should not happen.")
    return(scoreTable)
  }
  scoreTable[, `:=`(c("pValue", "logOddsRatio"), fisher.test(matrix(c(support, 
                    b, c, d), 2, 2), alternative = "two.sided")[c("p.value", 
                    "estimate")]), by = list(userSet, dbSet)]
  scoreTable[, pValue := ifelse(pValue > 1, 1, pValue)]
  if (requireNamespace("qvalue", quietly = TRUE)) {
    tryCatch({
      scoreTable[, `:=`(qValue, qvalue::qvalue(pValue)$qvalue)]
    }, error = function(e) {
      warning(paste0("Problem in FDR calculation with qvalue [",e$message,"]. p-Value range:", str(range(scoreTable$pValue)) ) )
    })
  }
  else {
  }
  scoreTable[, `:=`(pValueLog, -log10(pValue))]
  scoreTable[, `:=`(rnkSup, rank(-support, ties.method = "min")), 
             by = userSet]
  scoreTable[, `:=`(rnkPV, rank(-pValueLog, ties.method = "min")), 
             by = userSet]
  scoreTable[, `:=`(rnkLO, rank(-logOddsRatio, ties.method = "min")), 
             by = userSet]
  scoreTable[, `:=`(maxRnk, max(c(rnkSup, rnkPV, rnkLO))), 
             by = list(userSet, dbSet)]
  scoreTable[, `:=`(meanRnk, signif(mean(c(rnkSup, rnkPV, rnkLO)), 
                                    3)), by = list(userSet, dbSet)]
  setkeyv(scoreTable, "dbSet")
  scoreTable = scoreTable[annotationDT]
  scoreTable[, `:=`(description, substr(description, 0, 80))]
  orderedCols = c("userSet", "dbSet", "collection", "pValueLog", 
                  "logOddsRatio", "support", "rnkPV", "rnkLO", "rnkSup", 
                  "maxRnk", "meanRnk", "b", "c", "d", "description", "cellType", 
                  "tissue", "antibody", "treatment", "dataSource", "filename")
  unorderedCols = setdiff(colnames(scoreTable), orderedCols)
  setcolorder(scoreTable, c(orderedCols, unorderedCols))
  scoreTable[order(pValueLog, -meanRnk, decreasing = TRUE), 
  ]
}
