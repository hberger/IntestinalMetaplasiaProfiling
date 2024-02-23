# Code taken from https://github.com/YuLab-SMU/enrichplot/blob/master/R/dotplot.R

#' default_labeller
#'
#' default labeling function that uses the
#' internal string wrapping function `ep_str_wrap`
#' @noRd
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

#' ep_str_wrap internal string wrapping function
#' @param string the string to be wrapped
#' @param width the maximum number of characters before wrapping to a new line
#' @noRd
ep_str_wrap <- function(string, width) {
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}


##' @rdname dotplot
##' @param object compareClusterResult object
##' @param by one of "geneRatio", "Percentage" and "count"
##' @param split ONTOLOGY or NULL
##' @param includeAll logical
##' @param font.size font size
##' @param title figure title
##' @param group a logical value, whether to connect the 
##' nodes of the same group with wires.
##' @param shape a logical value, whether to use nodes of 
##' different shapes to distinguish the group it belongs to
##' @param colorBy variable that used to color enriched terms,
##' e.g. 'pvalue', 'p.adjust' or 'qvalue'
#dotplot.compareClusterResult
my_dotplot <- function(object, x= "Cluster", colorBy="p.adjust",
                       showCategory=5, by="geneRatio", size="geneRatio",
                       split=NULL, includeAll=TRUE,
                       font.size=12, title="", label_format = 30,
                       group = FALSE, shape = FALSE) {
  color <- NULL
  df <- fortify(object, showCategory=showCategory, by=by,
                includeAll=includeAll, split=split)
  df[[x]] = gsub("\n"," ",df[[x]])
  
  if (by != "geneRatio")
    df$GeneRatio <- parse_ratio(df$GeneRatio)
  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }
  if (is.null(size)) size <- by
  by2 <- switch(size, rowPercentage = "Percentage", 
                count         = "Count", 
                geneRatio     = "GeneRatio")    
  p <- ggplot(df, aes_string(x = x, y = "Description", size = by2))      
  if (group) {
    p <- p + geom_line(aes_string(color = "Cluster", group = "Cluster"), size=.3) + 
      ggnewscale::new_scale_colour()     
  }
  
  if (shape) {
    ggstar <- "ggstar"
    require(ggstar, character.only=TRUE)
    # p <- p + ggsymbol::geom_symbol(aes_string(symbol = "Cluster", fill = colorBy)) +
    p <- p + ggstar::geom_star(aes_string(starshape="Cluster", fill=colorBy)) + 
      scale_fill_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))
  }  else {
    p <- p +  geom_point(aes_string(color = colorBy)) 
  }  
  #suppressMessages(print(
  p + scale_color_continuous(low="red", high="blue",
                             guide=guide_colorbar(reverse=TRUE)) +
    ylab(NULL) + ggtitle(title) + DOSE::theme_dose(font.size) +
    scale_size_continuous(range=c(3, 8)) + 
    scale_y_discrete(labels = label_func)
  #))
}

