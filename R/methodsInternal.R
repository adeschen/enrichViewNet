
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