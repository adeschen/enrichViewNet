#' gprofiler2cytoscape: TODO
#'
#' The gprofiler2cytoscape package TODO
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
#'     \item \code{\link{createNetwork}} {TODO}
#' }
#' 
#' @encoding UTF-8
#' @keywords package
NULL


#' All samples information, formated by \code{methylKit}, in a
#' \code{methylRawList} format (for demo purpose).
#'
#' The object is a \code{list} with 3 entries. Each entry corresponds to the
#' information for one generation (first entry = first generation, etc..)
#' stored in a \code{methylRawList}.
#' There are 12 samples (6 controls and 6 cases) for each generation. Each
#' sample information is stored in a \code{methylRaw} object.
#'
#' This dataset can be
#' used to test the \code{runPermutation} function.
#'
#' @name demoGOST
#'
#' @docType data
#'
#' @aliases demoGOST
#'
#' @format A \code{list} containing three \code{methylRawList} objects. Each
#' \code{methylRawList} contains the information for one generation
#' (first entry = first generation, etc..). Each sample information is
#' stored in a \code{methylRaw} object. There is \code{methylRaw} objects
#' (6 controls and 6 cases) in each generation.
#'
#' @return A \code{list} containing three \code{methylRawList} objects. Each
#' \code{methylRawList} contains the information for one generation
#' (first entry = first generation, etc..). Each sample information is
#' stored in a \code{methylRaw} object. There is \code{methylRaw} objects
#' (6 controls and 6 cases) in each generation.
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{createNetwork}} {for TODO}
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