#' From gprofiler2 to Cytoscape network
#'
#' The \code{gprofiler2cytoscape} package enables the visualization 
#' of enrichment results obtained by \code{gprofiler2} 
#' (\url{https://cran.r-project.org/web/packages/gprofiler2/index.html})  
#' under the form of \code{Cytoscape} network (\url{https://cytoscape.org/}). 
#' 
#' In those networks, both gene datasets (GO terms/pathways/protein complexes)
#' and genes are represented as nodes. A edge connect a gene to its datasets.
#' In the current version, only genes present in at least one gene dataset 
#' are retained.
#'
#' @docType package
#'
#' @name gprofiler2cytoscape-package
#'
#' @aliases gprofiler2cytoscape-package gprofiler2cytoscape
#'
#' @author Astrid Deschênes and
#' Pascal Belleau
#'
#' Maintainer:
#' Astrid Deschênes <adeschen@hotmail.com>
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{createNetwork}} {for transforming functional 
#'     enrichment results from gprofiler2 into a Cytoscape network}
#' }
#' 
#' @encoding UTF-8
#' @keywords package
NULL


#' The result of a functional enrichment analysis done with \code{gprofiler2} 
#' (\url{https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html}).
#' 
#' The object is a \code{list} with 2 entries. It contains the results of the
#' enrichment analysis as well as the metadata related to the analysis.
#'
#' This dataset can be
#' used to test the \code{\link{createNetwork}} function.
#'
#' @name demoGOST
#'
#' @docType data
#'
#' @aliases demoGOST
#'
#' @format The \code{list} contains two entries. The \code{result} entry 
#' contains a \code{data.frame} with the significant results obtained by
#' an enrichment analysis done with \code{gprofiler2}. The \code{meta} entry 
#' contains a named list with all the metadata for the query.
#'
#' @return  A \code{list} containing two entries. The \code{result} entry 
#' contains a \code{data.frame} with the significant results obtained by
#' an enrichment analysis done with \code{gprofiler2}. 
#' The \code{meta} entry contains a named list with all the 
#' metadata for the query.
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{createNetwork}} {for transforming functional 
#'     enrichment results from gprofiler2 into a Cytoscape network}
#' }
#'
#' @usage data(demoGOST)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(demoGOST)
#'
#' \donttest{
#' 
#' ## Create network for WikiPathways results
#' createNetwork(gostObject = demoGOST, source="WP")
#' 
#' }
#'
#'
NULL