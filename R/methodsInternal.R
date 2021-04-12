
#' @title Create network and load it into Cytoscape
#' 
#' @description Create network from gprofiler2 results and load it 
#' into Cytoscape
#' 
#' @param gostResults a \code{GRangesList} TODO
#' 
#' @param gostObject TODO
#' 
#' @param title a \code{GRangesList} TODO
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
#' @importFrom RCy3 createNetworkFromDataFrames setNodeColorMapping
#' @importFrom RCy3 setNodeLabelBypass setNodeWidthMapping 
#' @encoding UTF-8
#' @keywords internal
createCytoscapeNetwork <- function(gostResults, gostObject, title, collection) {
    
    done <- array()   
    nodes <- list()
    edges <- list()
    
    for(i in seq_len(nrow(gostResults))) {
        term <- gostResults$term_id[i]
        query <- gostResults$query[i]
        termName <- gostResults$term_name[i]
        #termShort <- str_replace(selectedTF[i, "term_name"], pattern="Factor: ", "")
        #termShort <- str_replace(termShort, pattern="; motif:.+$", "")
        nodes[[length(nodes) + 1]] <- data.frame(id=c(term),
                                            group=c("TERM"),
                                            alias=c(termName),
                                            stringsAsFactors=FALSE)
        
        res <- gconvert(query = c(term))
        genes <- gostObject$meta$query_metadata$queries[[query]]
        
        for (g in genes) {
            if (g %in% res$target) {
                geneName <- res[res$target == g, c("name")]
                if (! g %in% done) {
                    nodes[[length(nodes) + 1]] <- data.frame(id=c(g),
                                                        group=c("GENE"),
                                                        alias=c(geneName),
                                                        stringsAsFactors=FALSE)
                    done <- c(done, g)
                }
                edges[[length(edges) + 1]] <- data.frame(source=c(term),
                                                        target=c(g),
                                                        interaction=c("contains"),  # optional
                                                        stringsAsFactors=FALSE)
            }
        }
    }
    
    nodeInfo <- do.call(rbind, nodes)
    edgeInfo <- do.call(rbind, edges)
    
    createNetworkFromDataFrames(nodes=nodeInfo, edges=edgeInfo, 
                                title=title, collection=collection)
    
    column <- 'group'
    control.points <- c("TERM", "GENE")
    setNodeColorMapping(table.column=column, 
        table.column.values=control.points, colors=c('#ffba42', '#99CCFF'), 
        mapping.type="discrete", style.name="default")
    
    setNodeLabelBypass(node.names=nodeInfo$id, new.labels=nodeInfo$alias)
    
    setNodeWidthMapping(table.column=column, 
        table.column.values=control.points, widths = c(100, 75),
        mapping.type="discrete", style.name="default")
}


#' @title Verifying that Cytoscape is running
#' 
#' @description Verifying that Cytoscape is running
#' 
#' @return a \code{logical} indicating if Cytoscape is running.
#' 
#' @examples
#'
#' ## TODO
#' 
#' gprofiler2cytoscape:::isCytoscapeRunning()
#' 
#' @author Astrid Deschênes
#' @importFrom RCy3 cytoscapePing 
#' @encoding UTF-8
#' @keywords internal
isCytoscapeRunning <- function() {
    
    out <- tryCatch(
        {
            cytoscapePing()
            #message("Cytoscape is running.\n")
            return(TRUE)
        },
        error=function(cond) {
            message(paste0("Unable to connect to Cytoscape. \n", 
                                "Sif file will be created.\n"))
            return(FALSE)
        },
        warning=function(cond) {
            message(paste0("Unable to connect to Cytoscape. \n", 
                                "Sif file will be created."))
            return(FALSE)
        }
    ) 
    
    return(out)
}


#' @title TODO
#' 
#' @description TODO
#' 
#' @param gostObject a \code{list} created by grofiler2 that contains
#' the results from an enrichment analysis.
#' 
#' @param source a \code{character} string representing the selected source 
#' that will be used to generate the network. To hand-pick the terms to be 
#' used, "TERM_ID" should be used and the list of selected term IDs should
#' be passed through the \code{termIDs} parameter. The possible sources are 
#' "GO:BP" for Gene Ontology Biological Process, "GO:CC" for Gene Ontology  
#' Cellular Component, "GO:MF" for Gene Ontology Molecular Function, 
#' "KEGG" for Kegg, "REAC" for Reactome, "TF" for TRANSFAC, "MIRNA" for 
#' miRTarBase, "CORUM" for CORUM database, "HP" for Human phenotype ontology
#' and "WP" for WikiPathways. 
#' 
#' @param termIDs a \code{vector} of \code{character} strings that contains the
#' term IDS retained for the creation of the network. 
#' 
#' @return TRUE
#' 
#' @examples
#'
#' ## TODO
#' 
#' gostObject <- list()
#' gostObject[["meta"]] <- list()
#' gostObject[["result"]] <- list()
#' 
#' gprofiler2cytoscape:::validateCreateNetworkArguments(gostObject=gostObject,
#'     source="GO:MF", termIDs=NULL)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
validateCreateNetworkArguments <- function(gostObject, source, termIDs) {
    
    ## Test that gostObject is a gprofiler2 result 
    if (!("list" %in% class(gostObject) && "result" %in% names(gostObject) &&
            "meta" %in% names(gostObject)))   {
        stop(paste0("The gostObject object should be a list with meta ", 
                    "and result as entries corresponding to gprofiler2 ", 
                    "enrichment output."))
    } 
    
    return(TRUE)   
}

