#' selfssGSEA function
#'
# This function allow you to perform ssGSEA using your self gene list
# There is a born gene list including ten clock genes:ARNTL,CLOCK,PER1,PER2,PER3,CRY1,CRY2,NR1D1,NPAS2,RORA
#' @param pathGenes Enter a vector,the self-help genes.Defult pathGenes is circadian clock gene list
#' @param data Enter your expression matrix,whose row is gene and column is sample

selfssGSEA <- function(pathGenes=geneList, data){
  geneList <- list("score" = c(NA, pathGenes))
  g <- gsva(data, geneList, method="ssgsea", kcdf="Gaussian", abs.ranking=T)
  s <- g-min(g)/(max(g)-min(g))
  return(t(s))
}
