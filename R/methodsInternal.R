
#' @title Create network and load it into Cytoscape
#' 
#' @description Create network from gprofiler2 results and load it 
#' into Cytoscape
#' 
#' @param gostResults a \code{data.frame} containing the terms retained
#' for the creation of the network.
#' 
#' @param gostObject a \code{list} created by gprofiler2 that contains
#' the results from an enrichment analysis.
#' 
#' @param title a \code{character} string representing the name assigned to 
#' the network.
#' 
#' @param collection a \code{character} string representing the collection 
#' name assigned to the network.
#' 
#' @return \code{TRUE}
#' 
#' @examples
#'
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(demoGOST)
#' 
#' ## Only retained the GO - Molecular Function results
#' results <- demoGOST$result[demoGOST$result$source == "GO:MF", ]
#'
#' ## The creation of the network can only be done when Cytoscape 
#' ## is up and running
#' ## A network using GO - Molecular Function enriched terms will be 
#' ## generated and loaded into Cytoscape
#' if (gprofiler2cytoscape:::isCytoscapeRunning()) {
#'     gprofiler2cytoscape:::createCytoscapeNetwork(gostResults=results, 
#'         gostObject=demoGOST, title="Test", collection="Test Collection")
#' }
#' 
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
    
    ## Create the list of genes and terms that will be included in the network
    for (i in seq_len(nrow(gostResults))) {
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
        
        res <- gconvert(query=c(term))
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
    
    ## Create the network using JSON data format and posting it to Cytoscape
    createNetworkFromDataFrames(nodes=nodeInfo, edges=edgeInfo, 
                                title=title, collection=collection)
    
    ## Assign different colors to terms and genes
    column <- 'group'
    control.points <- c("TERM", "GENE")
    setNodeColorMapping(table.column=column, 
        table.column.values=control.points, colors=c('#ffba42', '#99CCFF'), 
        mapping.type="discrete", style.name="default")
    
    ## Override the node labels to use the term descriptions and gene names
    setNodeLabelBypass(node.names=nodeInfo$id, new.labels=nodeInfo$alias)
    
    ## Assign larger node widths to terms and smaller ones to genes
    setNodeWidthMapping(table.column=column, 
        table.column.values=control.points, widths=c(100, 75),
        mapping.type="discrete", style.name="default")
    
    return(TRUE)
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
#' @param removeRoot a \code{logical} that specified if the root terms of 
#' the selected source should be removed (when present). 
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
#'     source="GO:BP", termIDs=NULL, removeRoot=FALSE)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @keywords internal
validateCreateNetworkArguments <- function(gostObject, source, termIDs,
                                                removeRoot) {
    
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
    
    if (!is(removeRoot, "logical")) {
        stop(paste0("The \'removeRoot\' parameter must be the logical ", 
                        "value TRUE or FALSE."))
    }
    
    return(TRUE)   
}


#' @title Remove root term if present in the list of selected terms
#' 
#' @description Remove root term if present in the list of selected terms
#' 
#' @param gostResults a \code{data.frame} containing the terms retained
#' for the creation of the network.
#' 
#' @return a \code{data.frame} of selected terms without the root term.
#' 
#' @examples
#'
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(demoGOST)
#' 
#' ## Only retained the GO - Molecular Function results
#' results <- demoGOST$result[demoGOST$result$source == "WP", ]
#'
#' ## Remove WIKIPATHWAYS root term 
#' gprofiler2cytoscape:::removeRootTerm(gostResult=results)
#' 
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom RCy3 createNetworkFromDataFrames setNodeColorMapping
#' @importFrom RCy3 setNodeLabelBypass setNodeWidthMapping 
#' @encoding UTF-8
#' @keywords internal
removeRootTerm <- function(gostResult) {
    
    source <- c("WP", "KEGG", "REAC", "CORUM", "TF", 
                    "MIRNA", "HPA", "GO:BP", "GO:CC", "GO:MF")
    term_id <- c("WP:000000", "KEGG:00000", "REAC:0000000", 
                    "CORUM:0000000", "TF:M00000", "MIRNA:000000",
                    "HPA:0000000", "GO:0008150", "GO:0005575", "GO:0003674")
    
    for (i in seq_len(10)) {
        test <- which(gostResult$source == source[i] & 
                        gostResult$term_id == term_id[i])
    
        if (length(test) > 0) {
            gostResult <- gostResult[-c(test), ]
        }
    }
    
    rownames(gostResult) <- NULL
    
    return(gostResult)
}
