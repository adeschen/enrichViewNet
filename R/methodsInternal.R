
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
#' ## The creation of the network can only be done when Cytoscape
#' ## is up and running
#' ## A network using GO - Molecular Function enriched terms will be
#' ## generated and loaded into Cytoscape
#' if (enrichViewNet:::isCytoscapeRunning()) {
#'     enrichViewNet:::createCytoscapeNetwork(gostResults=results,
#'         gostObject=demoGOST, title="Test", collection="New Collection")
#' }
#'
#' @author Astrid Deschênes
#' @importFrom RCy3 createNetworkFromDataFrames setNodeColorMapping
#' @importFrom RCy3 setNodeLabelBypass setNodeWidthMapping
#' @encoding UTF-8
#' @keywords internal
createCytoscapeNetwork <- function(gostResults, gostObject, title, collection) {

    message("Preparing information for Cytoscape.")

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

    message("Done.")

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
#' enrichViewNet:::isCytoscapeRunning()
#'
#' @author Astrid Deschênes
#' @importFrom RCy3 cytoscapePing
#' @encoding UTF-8
#' @keywords internal
isCytoscapeRunning <- function() {

    out <- tryCatch(
        {
            cytoscapePing()
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
#' @param query a \code{character} string that specified the retained query to 
#' generate the network or \code{NULL}. 
#' 
#' @param title a \code{character} string representing the name assigned to 
#' the network.
#' 
#' @param collection a \code{character} string representing the collection 
#' name assigned to the network.
#' 
#' @param fileName a \code{character} string representing the name of the
#' CX JSON file that is created when Cytoscape is not running. The name
#' must have a '.cx' extension.
#' 
#' 
#' @return \code{TRUE} when all arguments are valid
#'
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(demoGOST)
#'
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateNetworkArguments(gostObject=demoGOST,
#'     source="GO:BP", termIDs=NULL, removeRoot=FALSE, query=NULL, 
#'     title="Network graph Test",
#'     collection="test collection", fileName="test.cx")
#'
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @importFrom stringr str_ends
#' @keywords internal
validateCreateNetworkArguments <- function(gostObject, source, termIDs,
                    removeRoot, query, title, collection, fileName) {

    ## Test that gostObject is a gprofiler2 result
    if (!(inherits(gostObject, "list") && "result" %in% names(gostObject) &&
            "meta" %in% names(gostObject)))   {
        stop("The gostObject object should be a list with meta ",
                    "and result as entries corresponding to gprofiler2 ",
                    "enrichment output.")
    }

    if (!is.null(query) & !is.character(query)) {
        stop("The \'query\' parameter must be a character string or \'NULL\'.")
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

    if (!is.logical(removeRoot)) {
        stop("The \'removeRoot\' parameter must be the logical ",
                        "value TRUE or FALSE.")
    }
    
    if (!is.character(title)) {
        stop("The \'title\' parameter must a character string.")
    }
    
    if (!is.character(collection)) {
        stop("The \'collection\' parameter must a character string.")
    }
    
    if (!is.character(fileName)) {
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
#' enrichViewNet:::removeRootTerm(gostResult=results)
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

    term_id <- c("WP:000000", "KEGG:00000", "REAC:0000000",
                    "CORUM:0000000", "TF:M00000", "MIRNA:000000",
                    "HPA:0000000", "GO:0008150", "GO:0005575", "GO:0003674")


    test <- which(gostResult$term_id %in% term_id)
    if (length(test) > 0) {
        gostResult <- gostResult[-c(test), ]
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
#' jsonFormat <- enrichViewNet:::createCytoscapeCXJSON(
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
#' enrichViewNet:::createMetaDataSectionCXJSON()
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
#' information <- enrichViewNet:::extractNodesAndEdgesInfoForCXJSON(
#'                 gostResults=results, gostObject=demoGOST)
#'
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @encoding UTF-8
#' @keywords internal
extractNodesAndEdgesInfoForCXJSON <- function(gostResults, gostObject) {

    listQuery <- unique(gostResults$query)

    listGenes <- unique(gostObject$meta$query_metadata$queries[[listQuery[1]]])
    gostQuery <- gostResults[which(gostResults$query == listQuery[1]),]
    listTerm <- gostQuery[!duplicated(gostQuery$term_id), c("term_id", "term_name")]
    res <- gconvert(query=c(listTerm$term_id))

    query=listQuery[1]

    ## Create a data.frame linking gene and term
    resEdge <- do.call(rbind, lapply(seq_len(nrow(listTerm)),
                FUN=function(x, listTerm, listGenes, resTerm=res, query){
                    tmp <- which(resTerm$input == listTerm$term_id[x] & resTerm$target %in% listGenes)
                    df <- NULL
                    if(length(tmp) > 0){
                        res <- resTerm[tmp,]
                        df <- data.frame(term=rep(listTerm$term_id[x],
                                                    nrow(res)),
                                termName=rep(listTerm$term_name[x],
                                                    nrow(res)),
                                gene=res$target,
                                geneName=res$name,
                                query=query,
                                stringsAsFactors = FALSE)
                    }
                    return(df)
                },
                listTerm=listTerm,
                listGenes=listGenes,
                resTerm=res,
                query=query))

    geneUnique <- resEdge[!duplicated(resEdge$gene),]

    ## Create node entries for the gene
    geneNodes <- data.frame('@id'=seq_len(nrow(geneUnique)),
                        n=geneUnique$gene,
                        stringsAsFactors=FALSE,
                        check.names=FALSE)

    ## Create node attribute entries for the gene alias and group
    geneAttributes  <- data.frame(
        po=rep(seq_len(nrow(geneUnique)), 2),
        n=c(rep("alias", nrow(geneUnique)),
            rep("group", nrow(geneUnique))),
        v=c(geneUnique$geneName,
            rep("GENE", nrow(geneUnique))),
        stringsAsFactors=FALSE)

    rownames(geneNodes) <- geneNodes$n

    # Offset for the id
    termOffSet <- nrow(geneUnique)

    termUnique <- resEdge[!duplicated(resEdge$term),]

    ## Create node entries for the term
    termNodes <- data.frame('@id'=seq_len(nrow(termUnique)) + termOffSet,
                        n=termUnique$term,
                        stringsAsFactors=FALSE,
                        check.names=FALSE)

    ## Create node attribute entries for the term alias and the term group
    termAttributes  <- data.frame(
        po=rep(seq_len(nrow(termUnique)) + termOffSet, 2) ,
        n=c(rep("alias", nrow(termUnique)),
            rep("group", nrow(termUnique))),
        v=c(termUnique$termName, rep("GENE", nrow(termUnique))),
        stringsAsFactors=FALSE)

    rownames(termNodes) <- termNodes$n

    # Offset for the id
    edgeOffSet <- termOffSet + nrow(termUnique)

    ## Create an edge connecting term to gene
    edges <- data.frame('@id'=seq_len(nrow(resEdge)) + edgeOffSet,
                    s=termNodes[resEdge$term, '@id'],
                    t=geneNodes[resEdge$gene, '@id'],
                    i=rep("contains", nrow(resEdge)),
                    stringsAsFactors=FALSE,
                    check.names=FALSE)

    ## Create edge attributes
    edgeAttributes <- data.frame(
        po=rep(edges[,'@id'], 3),
        n=c(rep("name", nrow(resEdge)),
            rep("source", nrow(resEdge)),
            rep("target", nrow(resEdge))),
        v=c(paste0(resEdge$term, " (contains) ", resEdge$gene),
            resEdge$term, resEdge$gene),
        stringsAsFactors=FALSE)

    rownames(geneNodes) <- NULL
    rownames(termNodes) <- NULL

    return(list(nodes=rbind(geneNodes, termNodes), edges=edges,
        nodeAttributes=rbind(geneAttributes, termAttributes),
        edgeAttributes=edgeAttributes))
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
#'     enrichViewNet:::extractNodesAndEdgesWhenNoIntersection(
#'                 gostResults=results, gostObject=demoGOST)
#'
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @encoding UTF-8
#' @keywords internal
extractNodesAndEdgesWhenNoIntersection <- function(gostResults, gostObject) {

    listQuery <- unique(gostResults$query)
    query=listQuery[1]
    listGenes <- unique(gostObject$meta$genes_metadata$query[[query]]$ensgs)
    # genes <- gostObject$meta$genes_metadata$query[[query]]$ensgs
    gostQuery <- gostResults[which(gostResults$query == query),]
    listTerm <- gostQuery[!duplicated(gostQuery$term_id),
                        c("term_id", "term_name")]
    res <- gconvert(query=c(listTerm$term_id),
                       organism=gostObject$meta$query_metadata$organism)
    # res <- gconvert(query=c(listTerm$term_id))



    ## Create a data.frame linking gene and term
    resEdge <- do.call(rbind, lapply(seq_len(nrow(listTerm)),
                                     FUN=function(x, listTerm, listGenes,
                                            resTerm, query){
                                            tmp <- which(resTerm$input ==
                                                             listTerm$term_id[x]
                                                         & resTerm$target %in%
                                                             listGenes)
                                            df <- NULL
                                            if(length(tmp) > 0){
                                                res <- resTerm[tmp,]
                                                df <- data.frame(term=rep(listTerm$term_id[x],
                                                            nrow(res)),
                                                        termName=rep(listTerm$term_name[x],
                                                            nrow(res)),
                                                        gene=res$target,
                                                        geneName=res$name,
                                                        query=query,
                                                        stringsAsFactors = FALSE)
                                            }
                                         return(df)
                                     },
                                     listTerm=listTerm,
                                     listGenes=listGenes,
                                     resTerm=res,
                                     query=query))


    geneUnique <- resEdge[!duplicated(resEdge$gene),]

    ## Create node entries for the gene
    geneNodes <- data.frame(id=geneUnique$gene,
                        group=rep("GENE", nrow(geneUnique)),
                        alias=c(geneUnique$geneName),
                        stringsAsFactors=FALSE)


    termUnique <- resEdge[!duplicated(resEdge$term),]

    ## Create node entries for the term
    termNodes <- data.frame(id=termUnique$term,
                            group=rep("TERM", nrow(termUnique)),
                            alias=termUnique$termName,
                            stringsAsFactors=FALSE)

    ## Create an edge connecting term to gene
    edges <- data.frame(source=resEdge$term,
                    target=resEdge$gene,
                    interaction=rep("contains", nrow(resEdge)),
                    stringsAsFactors=FALSE)

    return(list(nodes=rbind(geneNodes, termNodes), edges=edges))
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
#' data(parentalNapaVsDMSOEnrichment)
#'
#' ## Only retained the GO Molecular Function results
#' results <- parentalNapaVsDMSOEnrichment$result[
#'         parentalNapaVsDMSOEnrichment$result$source == "GO:MF", ]
#'
#' information <-
#'     enrichViewNet:::extractNodesAndEdgesWhenIntersection(
#'         gostResults=results, gostObject=parentalNapaVsDMSOEnrichment)
#'
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom stringr str_split
#' @encoding UTF-8
#' @keywords internal
extractNodesAndEdgesWhenIntersection <- function(gostResults, gostObject) {

    geneInter <- unique(unlist(str_split(gostResults$intersection, ",")))

    res <- gconvert(query=c(geneInter),
                    organism=gostObject$meta$query_metadata$organism)

    resEdge <- do.call(rbind, lapply(seq_len(nrow(gostResults)),
        FUN=function(x, gostResults, res){
            genes <- res[res$input %in%
                        unlist(str_split(gostResults$intersection[x], ",")),
                                c("input", "name") ]
            df <- NULL
            if(gostResults$intersection_size[x] > 0){
                df <- data.frame(term=rep(gostResults$term_id[x],
                                          nrow(genes)),
                                 termName=rep(gostResults$term_name[x],
                                              nrow(genes)),
                                 gene=genes$input,
                                 geneName=genes$name,
                                 query=gostResults$query[x],
                                 stringsAsFactors = FALSE)
            }
            return(df)
    }, gostResults=gostResults,
    res=res))

    geneUnique <- resEdge[!duplicated(resEdge$gene),]

    ## Create node entries for the gene
    geneNodes <- data.frame(id=geneUnique$gene,
                            group=rep("GENE", nrow(geneUnique)),
                            alias=c(geneUnique$geneName),
                            stringsAsFactors=FALSE)


    termUnique <- resEdge[!duplicated(resEdge$term),]

    ## Create node entries for the term
    termNodes <- data.frame(id=termUnique$term,
                            group=rep("TERM", nrow(termUnique)),
                            alias=termUnique$termName,
                            stringsAsFactors=FALSE)

    ## Create an edge connecting term to gene
    edges <- data.frame(source=resEdge$term,
                        target=resEdge$gene,
                        interaction=rep("contains", nrow(resEdge)),
                        stringsAsFactors=FALSE)

    return(list(nodes=rbind(geneNodes, termNodes), edges=edges))
}
