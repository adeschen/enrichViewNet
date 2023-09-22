
#' @title Using functional enrichment results in  gprofiler2 format to create a 
#' enrichment map
#' 
#' @description User selected enrichment terms are used to create a enrichment 
#' map. The selection of the term can by specifying by the 
#' source of the terms (GO:MF, REAC, TF, etc...) or by listing the selected 
#' term IDs. The map is only generated when there is at least on 
#' significant term to graph.
#' 
#' @param gostObject a \code{list} corresponding to gprofiler2 enrichment 
#' output that contains and that contains 
#' the results from an enrichment analysis.
#' 
#' @param query a \code{character} string representing the name of the query 
#' that is going to be used to generate the graph. The query must exist in the 
#' \code{gostObject} object.
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
#' @param removeRoot a \code{logical} that specified if the root terms of 
#' the selected source should be removed (when present). Default: \code{TRUE}.
#' 
#' @param termIDs a \code{vector} of \code{character} strings that contains the
#' term IDS retained for the creation of the network. Default: \code{NULL}.
#' 
#' @param title a \code{character} string representing the name TODO
#' 
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{vector} of terms, 
#' the selected terms will be displayed. Default: \code{30L}.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped. Default: \code{FALSE}.
#' 
#' @param cexLabelCategory a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1). Default: \code{1}.s
#' 
#' @param cexCategory a positive \code{numeric} representing the amount by 
#' which plotting category nodes should be scaled relative to the default (1). 
#' Default: \code{1}.
#' 
#' @return a \code{ggplot} object which is the enrichment map for enrichment 
#' results.
#' 
#' @examples
#'
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#'
#' ## Extract query information (only one in this dataset)
#' query <- unique(parentalNapaVsDMSOEnrichment$result$query)
#' 
#' ## Create graph for Gene Ontology - Cellular Component related results
#' createEnrichMap(gostObject=parentalNapaVsDMSOEnrichment, 
#'     query=query, source="GO:CC", removeRoot=TRUE,
#'     title="GO Cellular Component")
#' 
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom strex match_arg
#' @encoding UTF-8
#' @export
createEnrichMap <- function(gostObject, query, source=c("TERM_ID", "GO:MF", 
        "GO:CC", "GO:BP", "KEGG", "REAC", "TF", "MIRNA", "HPA", "CORUM", 
        "HP", "WP"), 
        termIDs=NULL, removeRoot=TRUE, title="gprofiler network", 
        showCategory=30, groupCategory=FALSE, cexLabelCategory=1,
        cexCategory=1) {
    
    ## Validate source is among the possible choices
    source <- match_arg(source, ignore_case=TRUE)
    
    ## Validate parameters
    validateCreateEnrichMapArguments(gostObject=gostObject, query=query, 
        source=source, termIDs=termIDs, removeRoot=removeRoot, title=title, 
        showCategory=showCategory, cexLabelCategory=cexLabelCategory,
        groupCategory=groupCategory, cexCategory=cexCategory)
    
    ## Extract results
    gostResults <- gostObject$result
    
    ## Retain results associated to query
    gostResults <- gostResults[gostResults$query == query,]
    
    ## Filter results
    if (source == "TERM_ID") {
        gostResults <- gostResults[gostResults$term_id %in% termIDs,]
    } else {
        gostResults <- gostResults[gostResults$source == source,]
    }
    
    ## Remove root term if required
    if (removeRoot) {
        gostResults <- removeRootTerm(gostResults)
        if (nrow(gostResults) == 0) {
            stop("With removal of the root term, there is no ", 
                 "enrichment term left")
        }
    }
    
    meta <- gostObject$meta
    
    backgroundGenes <- meta$query_metadata$queries[[query]]
    
    significantMethod <- meta$query_metadata$significance_threshold_method
    
    ## Create basic emap
    emap <- createBasicEmap(gostResults=gostResults, 
                backgroundGenes=backgroundGenes, title=title, 
                showCategory=showCategory, cexLabelCategory=cexLabelCategory,
                groupCategory=groupCategory, cexCategory=cexCategory, 
                significantMethod=significantMethod)
    
    return(emap)
}
