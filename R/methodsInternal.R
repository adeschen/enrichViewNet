
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
#' ## Only retained the WikiPathways results
#' results <- demoGOST$result[demoGOST$result$source == "WP", ]
#' 
#' \dontrun{
#' 
#' ## The creation of the network can only be done when Cytoscape 
#' ## is up and running
#' ## A network using GO - Molecular Function enriched terms will be 
#' ## generated and loaded into Cytoscape
#' if (gprofiler2cytoscape:::isCytoscapeRunning()) {
#'     gprofiler2cytoscape:::createCytoscapeNetwork(gostResults=results, 
#'         gostObject=demoGOST, title="Test", collection="New Collection")
#' }
#' 
#' }
#' 
#' @author Astrid Deschênes
#' @importFrom RCy3 createNetworkFromDataFrames setNodeColorMapping
#' @importFrom RCy3 setNodeLabelBypass setNodeWidthMapping 
#' @encoding UTF-8
#' @keywords internal
createCytoscapeNetwork <- function(gostResults, gostObject, title, collection) {
    
    ## Extract node and edge information to be used to create network
    if (! "intersection" %in% colnames(gostResults)) {
        info <- extractNodesAndEdgesWhenNoIntersection(gostResults=gostResults,
                                                    gostObject=gostObject)
    } else {
        info <- extractNodesAndEdgesWhenIntersection(gostResults=gostResults,
                                                  gostObject=gostObject)
    }
        
    ## Create the network using JSON data format and posting it to Cytoscape
    createNetworkFromDataFrames(nodes=info$nodes, edges=info$edges, 
                                title=title, collection=collection)
    
    ## Assign different colors to terms and genes
    column <- 'group'
    control.points <- c("TERM", "GENE")
    setNodeColorMapping(table.column=column, 
        table.column.values=control.points, colors=c('#ffba42', '#99CCFF'), 
        mapping.type="discrete", style.name="default")
    
    ## Override the node labels to use the term descriptions and gene names
    setNodeLabelBypass(node.names=info$nodes$id, new.labels=info$nodes$alias)
    
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
            message("Unable to connect to Cytoscape. \n", 
                                "CX JSON file will be created.\n")
            return(FALSE)
        },
        warning=function(cond) {
            message("Unable to connect to Cytoscape. \n", 
                                "CX JSON file will be created.")
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
#' @param fileName a \code{character} string representing the name of the
#' CX JSON file that is created when Cytoscape is not running. The name 
#' must have a '.cx' extension.
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
#'     source="GO:BP", termIDs=NULL, removeRoot=FALSE, fileName="test.cx")
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @importFrom stringr str_ends
#' @keywords internal
validateCreateNetworkArguments <- function(gostObject, source, termIDs,
                                                removeRoot, fileName) {
    
    ## Test that gostObject is a gprofiler2 result 
    if (!("list" %in% class(gostObject) && "result" %in% names(gostObject) &&
            "meta" %in% names(gostObject)))   {
        stop("The gostObject object should be a list with meta ", 
                    "and result as entries corresponding to gprofiler2 ", 
                    "enrichment output.")
    } 
    
    if (source != "TERM_ID") {
        if (sum(gostObject$result$source == source) < 1) {
            stop("There is no enriched term for the selected ", 
                    "source \'", source, "\'.")    
        }
    } else {
        if (is.null(termIDs)) {
            stop("A vector of terms should be given through the ",
                    "\'termIDs\' parameter when source is \'TERM_ID\'.")  
        }
        else {
            if(!all(termIDs %in% gostObject$result$term_id)) {
                stop("Not all listed terms are present in the  ",
                                "enrichment results.")
            }
        }
    }
    
    if (!is(removeRoot, "logical")) {
        stop("The \'removeRoot\' parameter must be the logical ", 
                        "value TRUE or FALSE.")
    }
    
    if (!is(fileName, "character")) {
        stop("The \'fileName\' parameter must a character string.")
    }
    
    if (str_ends(fileName, ".cx", negate = TRUE)) {
        stop("The \'fileName\' parameter must have \'.cx\' extension.")
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
    
    ## Root terms that should be removed
    source <- c("WP", "KEGG", "REAC", "CORUM", "TF", 
                    "MIRNA", "HPA", "GO:BP", "GO:CC", "GO:MF")
    term_id <- c("WP:000000", "KEGG:00000", "REAC:0000000", 
                    "CORUM:0000000", "TF:M00000", "MIRNA:000000",
                    "HPA:0000000", "GO:0008150", "GO:0005575", "GO:0003674")
    
    ## When a root term is present, remove it from the data.frame
    for (i in seq_len(length(source))) {
        test <- which(gostResult$source == source[i] & 
                        gostResult$term_id == term_id[i])
    
        if (length(test) > 0) {
            gostResult <- gostResult[-c(test), ]
        }
    }
    
    ## Reset row names
    rownames(gostResult) <- NULL
    
    return(gostResult)
}

#' @title Create CX JSON text representing the network
#' 
#' @description Create a CX JSON text that represent the network which 
#' includes information about nodes and edges present in the network.
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
#' @return \code{character} string that represent the network in a CX JSON
#' format. 
#' 
#' @examples
#'
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(demoGOST)
#' 
#' ## Only retained the WikiPathways results
#' results <- demoGOST$result[demoGOST$result$source == "WP", ]
#' 
#' jsonFormat <- gprofiler2cytoscape:::createCytoscapeCXJSON(
#'                 gostResults = results, gostObject = demoGOST, 
#'                 title = "WikiPathways")
#' 
#' @author Astrid Deschênes
#' @importFrom jsonlite toJSON
#' @encoding UTF-8
#' @keywords internal
createCytoscapeCXJSON <- function(gostResults, gostObject, title) {
    
    entriesL <- extractNodesAndEdgesInfoForCXJSON(gostResults=gostResults,
        gostObject = gostObject)
    
    networkAttributes <- data.frame(n=c("name"), v=c(title), 
        stringsAsFactors=FALSE)
    
    cyHiddenAttributes <- data.frame(n=c("layoutAlgorithm"),
        y=c("yFiles Circular Layout"), stringsAsFactors=FALSE)
    
    metaData <- createMetaDataSectionCXJSON()
    
    ## Create the network using JSON data format and posting it to Cytoscape
    result <- paste0("[", metaData, ",",
            toJSON(list(networkAttributes=networkAttributes)), ",", 
            toJSON(list(nodes=entriesL$nodes)), ",", 
            toJSON(list(edges=entriesL$edges)), ",",
            toJSON(list(nodeAttributes=entriesL$nodeAttributes)), ",",
            toJSON(list(edgeAttributes=entriesL$edgeAttributes)), ",",
            toJSON(list(cyHiddenAttributes=cyHiddenAttributes)), ",",
            metaData, ",{\"status\":[{\"error\":\"\",\"success\":true}]}]")
    
    return(result)
}


#' @title Create meta data section for the CX JSON file
#' 
#' @description Create meta data section for the CX JSON file that contains
#' the network information
#' 
#' @return a \code{JSON} object that contains the meta data section related 
#' to the network
#' 
#' @examples
#' 
#' ## Create the JSON object that contains the meta data information
#' gprofiler2cytoscape:::createMetaDataSectionCXJSON()
#' 
#' @author Astrid Deschênes
#' @importFrom jsonlite toJSON
#' @encoding UTF-8
#' @keywords internal
createMetaDataSectionCXJSON <- function() {
    
    name <- c("nodes", "edges", "edgeAttributes", "nodeAttributes", 
                "cyHiddenAttributes", "cyNetworkRelations", "cyGroups",
                "networkAttributes", "cyTableColumn", "cySubNetworks")
    
    metaData <- data.frame(name=name, version=rep("1.0", length(name)), 
                            stringsAsFactors=FALSE)
    
    return(toJSON(list(metaData=metaData)))    
}


#' @title Extract information about nodes and edges
#' 
#' @description Extract information about nodes and edges that is necessary
#' to create the CX JSON text representing the network
#' 
#' @param gostResults a \code{data.frame} containing the terms retained
#' for the creation of the network.
#' 
#' @param gostObject a \code{list} created by gprofiler2 that contains
#' the results from an enrichment analysis.
#' 
#' @return a \code{list} containing 4 entries: 
#' \itemize{
#' \item{"nodes"}{a \code{data.frame} containing the information about 
#' the nodes present in the network.}
#' \item{"edges"}{a \code{data.frame} containing the information about 
#' the edges present in the network.}
#' \item{"nodeAttributes"}{a \code{data.frame} containing the attributes 
#' associated to the nodes present in the network.}
#' \item{"edgesAttributes"}{a \code{data.frame} containing the attributes 
#' associated to the edges present in the network}
#' }
#' 
#' @examples
#'
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(demoGOST)
#' 
#' ## Only retained the WikiPathways results
#' results <- demoGOST$result[demoGOST$result$source == "WP", ]
#' 
#' information <- gprofiler2cytoscape:::extractNodesAndEdgesInfoForCXJSON(
#'                 gostResults=results, gostObject=demoGOST)
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @encoding UTF-8
#' @keywords internal
extractNodesAndEdgesInfoForCXJSON <- function(gostResults, gostObject) {
    
    done <- list()   
    nodes <- list()
    edges <- list()
    nodeAttributes <- list()
    edgeAttributes <- list()
    
    id <- 0
    
    ## Create entries for genes and terms that will be included in the network
    for (i in seq_len(nrow(gostResults))) {
        term <- gostResults$term_id[i]
        termName <- gostResults$term_name[i]
        query <- gostResults$query[i]
        
        id <- id + 1
        term_id <- id
        
        ## Create a node entry for the term
        nodes[[length(nodes) + 1]] <- data.frame('@id' = term_id,
            n = c(term), stringsAsFactors = FALSE, check.names = FALSE)
        
        ## Create a node attribute entry for the term alias and the term group
        nodeAttributes[[length(nodeAttributes) + 1]] <- data.frame(
            po = rep(term_id, 2), n = c("alias", "group"),
            v = c(termName, "TERM"), stringsAsFactors = FALSE)
        
        ## Get list of genes associated to the term
        res <- gconvert(query=c(term))
        genes <- gostObject$meta$query_metadata$queries[[query]]
        
        ## Create entries for genes that will be associated with this term
        for (g in genes[genes %in% res$target]) {
            geneName <- res[res$target == g, c("name")]
            ## No need to create node if already created for this gene
            gene_id <- NULL
            if (! g %in% names(done)) {
                id <- id + 1
                gene_id <- id
                
                ## Create a node entry for the gene
                nodes[[length(nodes) + 1]] <- data.frame('@id'=gene_id,
                    n=c(g), stringsAsFactors=FALSE, check.names=FALSE)
                
                ## Create a node attribute entry for the gene alias and group
                nodeAttributes[[length(nodeAttributes) + 1]] <- data.frame(
                    po=rep(gene_id, 2), n=c("alias", "group"),
                    v=c(geneName, "GENE"), stringsAsFactors=FALSE)
                
                done[[g]] <- gene_id
            } else {
                gene_id <- done[[g]] 
            }
            
            id <- id + 1
            
            ## Create an edge connecting term to gene
            edges[[length(edges) + 1]] <- data.frame('@id'=id, s=term_id,
                t=gene_id, i = c("contains"), stringsAsFactors=FALSE, 
                check.names=FALSE)
            
            ## Create edge attributes
            edgeAttributes[[length(edgeAttributes) + 1]] <- data.frame(
                po=rep(id, 3), n=c("name", "source", "target"),
                v=c(paste0(term, " (contains) ", g), term, g),
                stringsAsFactors=FALSE)
        }
    }
    
    ## Transform results into data.frame
    nodeInfo <- do.call(rbind, nodes)
    edgeInfo <- do.call(rbind, edges)
    nodeAttributesInfo <- do.call(rbind, nodeAttributes)
    edgeAttributesInfo <- do.call(rbind, edgeAttributes)
    
    return(list(nodes=nodeInfo, edges=edgeInfo, 
        nodeAttributes=nodeAttributesInfo, edgeAttributes=edgeAttributesInfo))
}



#' @title Extract node and edge information to be used to create Cytoscape
#' network
#' 
#' @description Create a list containing all node and edge information needed 
#' to create the Cytoscape network
#' 
#' @param gostResults a \code{data.frame} containing the terms retained
#' for the creation of the network.
#' 
#' @param gostObject a \code{list} created by gprofiler2 that contains
#' the results from an enrichment analysis.
#' 
#' @return \code{list} containing 2 entries:
#' \itemize{
#' \item{"nodes"}{a \code{data.frame} containing the information about 
#' the nodes present in the network.}
#' \item{"edges"}{a \code{data.frame} containing the information about 
#' the edges present in the network.}
#' }
#' 
#' @examples
#'
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(demoGOST)
#' 
#' ## Only retained the GO Molecular Function results
#' results <- demoGOST$result[demoGOST$result$source == "GO:MF", ]
#' 
#' information <- 
#'     gprofiler2cytoscape:::extractNodesAndEdgesWhenNoIntersection(
#'                 gostResults=results, gostObject=demoGOST)
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @encoding UTF-8
#' @keywords internal
extractNodesAndEdgesWhenNoIntersection <- function(gostResults, gostObject) {

    done <- array()   
    nodes <- list()
    edges <- list()
    
    ## Create the list of genes and terms that will be included in the network
    for (i in seq_len(nrow(gostResults))) {
        term <- gostResults$term_id[i]
        query <- gostResults$query[i]
        termName <- gostResults$term_name[i]
        
        nodes[[length(nodes) + 1]] <- data.frame(id=c(term),
            group=c("TERM"), alias=c(termName), stringsAsFactors=FALSE)
        
        res <- gconvert(query=c(term), 
                            organism=gostObject$meta$query_metadata$organism)
        genes <- gostObject$meta$genes_metadata$query[[query]]$ensgs
        
        for (g in genes) {
            if (g %in% res$target) {
                geneName <- res[res$target == g, c("name")]
                
                if (! g %in% done) {
                    nodes[[length(nodes) + 1]] <- data.frame(id=c(g),
                        group=c("GENE"), alias=c(geneName),
                        stringsAsFactors=FALSE)
                    done <- c(done, g)
                }
                
                edges[[length(edges) + 1]] <- data.frame(source=c(term),
                    target=c(g), interaction=c("contains"), 
                    stringsAsFactors=FALSE)
            }
        }
    }
    
    ## Create data.frame for nodes and edges by merging all extracted info
    nodeInfo <- do.call(rbind, nodes)
    edgeInfo <- do.call(rbind, edges)

    return(list(nodes=nodeInfo, edges=edgeInfo))
}




#' @title Extract node and edge information to be used to create Cytoscape
#' network
#' 
#' @description Create a list containing all node and edge information needed 
#' to create the Cytoscape network
#' 
#' @param gostResults a \code{data.frame} containing the terms retained
#' for the creation of the network. The \code{data.frame} does not contain a   
#' column called "intersection".
#' 
#' @param gostObject a \code{list} created by gprofiler2 that contains
#' the results from an enrichment analysis.
#' 
#' @return \code{list} containing 2 entries:
#' \itemize{
#' \item{"nodes"}{a \code{data.frame} containing the information about 
#' the nodes present in the network.}
#' \item{"edges"}{a \code{data.frame} containing the information about 
#' the edges present in the network.}
#' }
#' 
#' @examples
#'
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(demoGOST)
#' 
#' ## Only retained the GO Molecular Function results
#' results <- demoGOST$result[demoGOST$result$source == "GO:MF", ]
#' 
#' ##information <- 
#' ##      gprofiler2cytoscape:::extractNodesAndEdgesWhenIntersection(
#' ##                gostResults=results, gostObject=demoGOST)
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom stringr str_split
#' @encoding UTF-8
#' @keywords internal
extractNodesAndEdgesWhenIntersection <- function(gostResults, gostObject) {
    
    done <- array()   
    nodes <- list()
    edges <- list()
    
    geneInter <- unique(unlist(str_split(gostResults$intersection, ",")))
    
    res <- gconvert(query=c(geneInter), 
                    organism=gostObject$meta$query_metadata$organism)
    
    ## Create the list of genes and terms that will be included in the network
    for (i in seq_len(nrow(gostResults))) {
        term <- gostResults$term_id[i]
        query <- gostResults$query[i]
        termName <- gostResults$term_name[i]
        
        nodes[[length(nodes) + 1]] <- data.frame(id=c(term),
            group=c("TERM"), alias=c(termName), stringsAsFactors=FALSE)
        
        genes <- res[res$input %in% 
                    unlist(str_split(gostResults$intersection[i], ",")), 
                            c("input", "name") ]
        
        for (j in seq_len(nrow(genes))) {
                g <- genes$input[j]
                geneName <- genes$name[j]
                print(geneName)
                if (! g %in% done) {
                    nodes[[length(nodes) + 1]] <- data.frame(id=c(g),
                                group=c("GENE"), alias=c(geneName),
                                stringsAsFactors=FALSE)
                    done <- c(done, g)
                }
                
                edges[[length(edges) + 1]] <- data.frame(source=c(term),
                                target=c(g), interaction=c("contains"), 
                                stringsAsFactors=FALSE)
        }
    }
    
    ## Create data.frame for nodes and edges by merging all extracted info
    nodeInfo <- do.call(rbind, nodes)
    edgeInfo <- do.call(rbind, edges)
    
    return(list(nodes=nodeInfo, edges=edgeInfo))
}
