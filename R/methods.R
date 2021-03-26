
#' @title TODO
#' 
#' @description TODO
#' 
#' @param gostObject a \code{list} created by grofiler2 that contains
#' the results from an enrichment analysis
#' 
#' @param source a \code{character} string representing the selected source 
#' that will be used to generate the network. To hand-pick the terms to be 
#' used, "TERM_ID" should be used and the list of selected term IDs should
#' be passed through the \code{termIDs} parameter. The possible sources are 
#' "GO:BP" for Gene Ontology Biological Process, "GO:CC" for Gene Ontology  
#' Cellular Component, "GO:MF" for Gene Ontology Molecular Function, 
#' "KEGG" for Kegg, "REAC" for Reactome, "TF" for TRANSFAC, "MIRNA" for 
#' miRTarBase, "CORUM" for CORUM database, "HP" for Human phenotype ontology
#' and "WP" for WikiPathways.  Default: "TERM_ID".
#' 
#' @param termIDs a \code{array} of \code{character} strings that contains the
#' term IDS retained for the creation of the network. Default: \code{NULL}.
#' 
#' @param title a \code{character} string TODO
#' 
#' @param collection a \code{character} string representing the  TODO
#' 
#' 
#' @return TODO
#' 
#' @examples
#'
#' ## TODO
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom RCy3 cytoscapePing 
#' @encoding UTF-8
#' @export
createNetwork <- function(gostObject, source=c("TERM_ID", "GO:MF", "GO:CC",
    "GO:BP", "KEGG", "REAC", "TF", "MIRNA", "HPA", "CORUM", "HP","WP"), 
    termIDs=NULL, title="gprofiler network", collection="enrichment results") {
    
    ## Test that gostObject is a gprofiler2 result 
    if (!("list" %in% class(gostObject) && "result" %in% names(gostObject) &&
          "meta" %in% names(gostObject)))   {
        stop(paste0("The gostObject object should be a list with meta ", 
                        "and result as entries corresponding to gprofiler2 ", 
                        "enrichment output."))
    } 
    
    ## Test that Cytoscape is running
    isRunning <- isCytoscapeRunning()
    
    ## Extract results
    gostResults <- gostObject$result
    
    ## Filter results
    if (source == "TERM_ID") {
        gostResults <- gostResults[gostResults$term_id %in% termIDs,]
    } else {
        gostResults <- gostResults[gostResults$source == source,]
    }
    
    
    if (isRunning) {
        createCytoscapeNetwork(gostResults=gostResults, gostObject=gostObject,
                                title=title, collection=collection)
    }
}
        