
#' @title Using functional enrichment results from gprofiler2 to create a 
#' Cytoscape network
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
#' @param termIDs a \code{vector} of \code{character} strings that contains the
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
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(demoGOST)
#'
#' \donttest{
#' 
#' ## Create network for Gene Ontology - Molecular Function related results
#' createNetwork(gostObject = demoGOST, source="GO:MF")
#' 
#' }
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom RCy3 cytoscapePing 
#' @importFrom strex match_arg
#' @encoding UTF-8
#' @export
createNetwork <- function(gostObject, source=c("TERM_ID", "GO:MF", "GO:CC",
    "GO:BP", "KEGG", "REAC", "TF", "MIRNA", "HPA", "CORUM", "HP","WP"), 
    termIDs=NULL, title="gprofiler network", collection="enrichment results") {
    
    ## Validate source is among the possible choices
    source <- match_arg(source, ignore_case=TRUE)
    
    ## Validate parameters
    validateCreateNetworkArguments(gostObject=gostObject, source=source,
                                    termIDs=termIDs)
    
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
        