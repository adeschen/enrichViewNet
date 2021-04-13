
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
        #termShort <- str_replace(selectedTF[i, "term_name"], 
        #                    pattern="Factor: ", "")
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
                                                    interaction=c("contains"), 
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
        table.column.values=control.points, widths=c(100, 75),
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
#' ## Test if Cytoscape is running
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


#' @title Validate arguments passed to creatNetwork() function
#' 
#' @description Validate the arguments passed to creatNetwork() function.
#' First, the object containing the enrichment results must correspond to a 
#' object created by  \code{gprofiler2} software. Second, the selected 
#' source must at least have one enriched term in the results. Then, if the
#' source is 'TERM_ID', the listed terms must be present in the enrichment
#' results.
#' 
#' @param gostObject a \code{list} created by \code{gprofiler2} that contains
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
#' @param termIDs a \code{vector} of \code{character} strings that contains 
#' the term IDs retained for the creation of the network. This parameter is 
#' only used when \code{source} is set to "TERM_ID".
#' 
#' @return \code{TRUE} when all arguments are valid
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(demoGOST)
#' 
#' ## Check that all arguments are valid
#' gprofiler2cytoscape:::validateCreateNetworkArguments(gostObject=demoGOST,
#'     source="GO:BP", termIDs=NULL)
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
    
    if (source != "TERM_ID") {
        if (sum(gostObject$result$source == source) < 1) {
            stop(paste0("There is no enriched term for the selected ", 
                    "source \'", source, "\'."))    
        }
    } else {
        if (is.null(termIDs)) {
            stop(paste0("A vector of terms should be given through the ",
                    "\'termIDs\' parameter when source is \'TERM_ID\'."))  
        }
        else {
            if(!all(termIDs %in% gostObject$result$term_id)) {
                stop(paste0("Not all listed terms are present in the  ",
                                "enrichment results.")) 
            }
        }
    }
    
    return(TRUE)   
}

