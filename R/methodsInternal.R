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
            "meta" %in% names(gostObject))) {
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


#' @title Filter the results to retain only the selected terms
#'
#' @description Filter the enrichment results to retain only the selected 
#' terms and remove root term if requested.
#'
#' @param gostResults a \code{data.frame} containing the enriched terms that 
#' should be filtered.
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
#' term IDS retained for the creation of the network or \code{NULL}.
#' 
#' @param removeRoot a \code{logical} that specified if the root terms of 
#' the selected source should be removed (when present).
#'
#' @return a \code{data.frame} of filtered terms with or without the root term.
#'
#' @examples
#'
#' ## Loading dataset containing result from an enrichment analysis done with
#' ## gprofiler2
#' data(demoGOST)
#'
#' ## Only retained the GO - Molecular Function results
#' results <- demoGOST$result
#'
#' ## Remove WIKIPATHWAYS root term
#' selected <- enrichViewNet:::filterResults(gostResults=results, source="WP", 
#'     termIDs=NULL, removeRoot=TRUE)
#'
#'
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
filterResults <- function(gostResults, source, termIDs, removeRoot) {
    
    ## Filter results
    if (source == "TERM_ID") {
        gostResults <- gostResults[gostResults$term_id %in% termIDs,]
    } else {
        gostResults <- gostResults[gostResults$source == source,]
    }
    
    ## Remove root term if required
    if (removeRoot) {
        gostResults <- removeRootTerm(gostResults)
    }
    
    return(gostResults)
}


#' @title Remove root term if present in the list of selected terms
#'
#' @description Remove root term if present in the list of selected terms
#'
#' @param gostResult a \code{data.frame} containing the terms retained
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
#' ## Only retained the WikiPathways results
#' results <- demoGOST$result[demoGOST$result$source == "WP", ]
#'
#' ## Remove WIKIPATHWAYS root term
#' enrichViewNet:::removeRootTerm(gostResult=results)
#'
#'
#' @author Astrid Deschênes
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


#' @title Extract node and edge information from the enrichment results
#'
#' @description Create a list containing all node and edge information needed
#' to create the network
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
#' \item{"geneNodes"}{a \code{data.frame} containing the information about
#' the nodes present in the network. The nodes are genes.}
#' \item{"termNodes"}{a \code{data.frame} containing the information about
#' the nodes present in the network. The nodes are terms.}
#' \item{"edges"}{a \code{data.frame} containing the information about
#' the edges present in the network. The edges connect one gene to one term.}
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
#' ## Extract node and edge information
#' information <-
#'     enrichViewNet:::extractNodesAndEdgesInformation(
#'         gostResults=results, gostObject=parentalNapaVsDMSOEnrichment)
#'
#' @author Astrid Deschênes
#' @importFrom stringr str_split
#' @encoding UTF-8
#' @keywords internal
extractNodesAndEdgesInformation <- function(gostResults, gostObject) {

    
    ## Extract node and edge information to be used to create network
    if (! "intersection" %in% colnames(gostResults)) {
        entriesL <- extractInformationWhenNoIntersection(
            gostResults=gostResults, gostObject=gostObject)
    } else {
        entriesL <- extractInformationWhenIntersection(
            gostResults=gostResults)
    }

    return(entriesL)       
}

#' @title Extract node and edge information to be used to create network 
#' when interaction column is present
#'
#' @description Create a list containing all node and edge information needed
#' to create the network
#'
#' @param gostResults a \code{data.frame} containing the terms retained
#' for the creation of the network. The \code{data.frame} must contain a
#' column called "intersection".
#'
#' @return \code{list} containing 2 entries:
#' \itemize{
#' \item{\code{geneNodes}}{  a \code{data.frame} containing the information 
#' about the nodes present in the network. The nodes are genes.}
#' \item{\code{termNodes}}{  a \code{data.frame} containing the information 
#' about the nodes present in the network. The nodes are terms.}
#' \item{\code{edges}}{  a \code{data.frame} containing the information about
#' the edges present in the network. The edges connect one gene to one term.}
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
#' ## Extract node and edge information
#' information <-
#'     enrichViewNet:::extractInformationWhenIntersection(
#'         gostResults=results)
#'
#' @author Astrid Deschênes
#' @importFrom stringr str_split
#' @encoding UTF-8
#' @keywords internal
extractInformationWhenIntersection <- function(gostResults) {
    
    geneInter <- unique(unlist(str_split(gostResults$intersection, ",")))
    
    resEdge <- do.call(rbind, lapply(seq_len(nrow(gostResults)),
        FUN=function(x, gostResults){
            genes <- geneInter[geneInter %in% 
                        unlist(str_split(gostResults$intersection[x], ","))]
            nbEntries <- length(genes)
            df <- NULL
            if (gostResults$intersection_size[x] > 0) {
                df <- data.frame(term=rep(gostResults$term_id[x], nbEntries),
                        termName=rep(gostResults$term_name[x], nbEntries),
                        gene=genes,
                        geneName=genes,
                        query=gostResults$query[x], stringsAsFactors=FALSE)
            }
            return(df)
        }, gostResults=gostResults))
    
    ## Create node entries for the genes
    geneUnique <- resEdge[!duplicated(resEdge$gene),]
    geneNodes <- data.frame(id=geneUnique$gene,
                                group=rep("GENE", nrow(geneUnique)),
                                alias=c(geneUnique$geneName),
                                stringsAsFactors=FALSE)
    
    ## Create node entries for the term
    termUnique <- resEdge[!duplicated(resEdge$term),]
    termNodes <- data.frame(id=termUnique$term,
                                group=rep("TERM", nrow(termUnique)),
                                alias=termUnique$termName,
                                stringsAsFactors=FALSE)
    
    ## Create an edge connecting term (source) to gene (target)
    edges <- data.frame(source=resEdge$term,
                            target=resEdge$gene,
                            interaction=rep("contains", nrow(resEdge)),
                            stringsAsFactors=FALSE)
    
    return(list(geneNodes=geneNodes, termNodes=termNodes, edges=edges))
}

#' @title Extract information about nodes and edges to be used to create 
#' network when interaction column is missing
#'
#' @description Create a list containing all node and edge information needed
#' to create the network
#'
#' @param gostResults a \code{data.frame} containing the terms retained
#' for the creation of the network. The \code{data.frame} does not contain a
#' column called "intersection".
#'
#' @param gostObject a \code{list} created by gprofiler2 that contains
#' the results and meta-data from an enrichment analysis.
#'
#' @return \code{list} containing 2 entries:
#' \itemize{
#' \item{\code{geneNodes}}{  a \code{data.frame} containing the information 
#' about the nodes present in the network. The nodes are genes.}
#' \item{\code{termNodes}}{  a \code{data.frame} containing the information 
#' about the nodes present in the network. The nodes are terms.}
#' \item{\code{edges}}{  a \code{data.frame} containing the information about
#' the edges present in the network. The edges connect one gene to one term.}
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
#' information <- enrichViewNet:::extractInformationWhenNoIntersection(
#'                 gostResults=results, gostObject=demoGOST)
#'
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @encoding UTF-8
#' @keywords internal
extractInformationWhenNoIntersection <- function(gostResults, gostObject) {
    
    listQuery <- unique(gostResults$query)
    
    query <- listQuery[1]
    
    listGenes <- unique(gostObject$meta$query_metadata$queries[[query]])
    
    ## Extract information about the enriched terms from database
    gostQuery <- gostResults[which(gostResults$query == query),]
    listTerm <- gostQuery[!duplicated(gostQuery$term_id), 
                                            c("term_id", "term_name")]
    res <- gconvert(query=c(listTerm$term_id))
    
    ## Create a data.frame linking gene and term
    resEdge <- do.call(rbind, lapply(seq_len(nrow(listTerm)),
        FUN=function(x, listTerm, listGenes, resTerm=res, query){
            tmp <- which(resTerm$input == listTerm$term_id[x] & 
                                resTerm$target %in% listGenes)
            df <- NULL
            if(length(tmp) > 0) {
                res <- resTerm[tmp,]
                df <- data.frame(term=rep(listTerm$term_id[x], nrow(res)),
                            termName=rep(listTerm$term_name[x], nrow(res)),
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
    
    ## Create node entries for the genes
    geneNodes <- data.frame(id=geneUnique$gene, group="GENE",
                                alias=geneUnique$geneName,
                                stringsAsFactors=FALSE,
                                check.names=FALSE)
    
    termUnique <- resEdge[!duplicated(resEdge$term),]
    
    ## Create node entries for the terms
    termNodes <- data.frame(id=termUnique$term, group="TERM",
                                alias=termUnique$termName,
                                stringsAsFactors=FALSE,
                                check.names=FALSE)
    
    ## Create an edge connecting term to gene
    edges <- data.frame(source=resEdge$term, target=resEdge$gene,
                            interaction="contains",
                            stringsAsFactors=FALSE,
                            check.names=FALSE)
    
    rownames(geneNodes) <- NULL
    rownames(termNodes) <- NULL
    rownames(edges) <- NULL
    
    return(list(geneNodes=geneNodes, termNodes=termNodes, edges=edges))
}


#########################################################################
## CX JSON section
#########################################################################


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


#' @title Create network and load it into Cytoscape
#'
#' @description Create network from gprofiler2 results and load it
#' into Cytoscape
#'
#' @param nodeEdgeInfo a TODO
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
#' data(parentalNapaVsDMSOEnrichment)
#'
#' ## Only retained the GO Molecular Function results
#' results <- parentalNapaVsDMSOEnrichment$result[
#'         parentalNapaVsDMSOEnrichment$result$source == "GO:MF", ]
#'
#' ## Extract node and edge information
#' information <- enrichViewNet:::extractInformationWhenIntersection(
#'         gostResults=results)
#'     
#' ## The creation of the network can only be done when Cytoscape
#' ## is up and running
#' ## A network using GO - Molecular Function enriched terms will be
#' ## generated and loaded into Cytoscape
#' if (enrichViewNet:::isCytoscapeRunning()) {
#'     enrichViewNet:::createNetworkForCytoscape(nodeEdgeInfo=information, 
#'         title="Test", collection="New Collection")
#' }
#'
#' @author Astrid Deschênes
#' @importFrom RCy3 createNetworkFromDataFrames setNodeColorMapping
#' @importFrom RCy3 setNodeLabelBypass setNodeWidthMapping
#' @encoding UTF-8
#' @keywords internal
createNetworkForCytoscape <- function(nodeEdgeInfo, title, collection) {
    
    nodes <- rbind(nodeEdgeInfo$termNodes, nodeEdgeInfo$geneNodes)
    
    ## Create the network using JSON data format and posting it to Cytoscape
    createNetworkFromDataFrames(nodes=nodes, 
        edges=nodeEdgeInfo$edges, title=title, collection=collection)
    
    ## Assign different colors to terms and genes
    column <- 'group'
    control.points <- c("TERM", "GENE")
    setNodeColorMapping(table.column=column,
        table.column.values=control.points, colors=c('#ffba42', '#99CCFF'),
        mapping.type="discrete", style.name="default")
    
    ## Override the node labels to use the term descriptions and gene names
    setNodeLabelBypass(node.names=nodes$id, new.labels=nodes$alias)
    
    ## Assign larger node widths to terms and smaller ones to genes
    setNodeWidthMapping(table.column=column,
        table.column.values=control.points, widths=c(100, 75),
        mapping.type="discrete", style.name="default")
    
    return(TRUE)
}


#' @title Transform the node and edge information into a format easy to process 
#' to create a CX JSON text file
#'
#' @description Extract information about nodes and edges that is necessary
#' to create the CX JSON text representing the network
#'
#' @param results a \code{list} containing the information about the 
#' nodes and edges present in the network. TODOS
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
#' data(parentalNapaVsDMSOEnrichment)
#'
#' ## Only retained the GO Molecular Function results
#' results <- parentalNapaVsDMSOEnrichment$result[
#'         parentalNapaVsDMSOEnrichment$result$source == "GO:MF", ]
#'
#' ## Extract node and edge information
#' info <- enrichViewNet:::extractNodesAndEdgesInformation(gostResults=results, 
#'             gostObject=parentalNapaVsDMSOEnrichment)
#'             
#' ## Format node and edge information
#' information <- enrichViewNet:::formatInformationForCXJSON(result=info)
#'
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
formatInformationForCXJSON <- function(results) {

    #####################    
    ## Gene section
    
    nbGenes <- nrow(results$geneNodes)
    results$geneNodes$gene_id <- seq_len(nbGenes)
        
    ## Create node entries for the gene
    geneNodes <- data.frame('@id'=results$geneNodes$gene_id,
        n=results$geneNodes$id, stringsAsFactors=FALSE, check.names=FALSE)
    
    ## Create node attribute entries for the gene alias and group
    geneAttributes  <- data.frame(
        po=rep(seq_len(nbGenes), 2),
        n=c(rep("alias", nbGenes), rep("group", nbGenes)),
        v=c(results$geneNodes$alias, rep("GENE", nbGenes)),
        stringsAsFactors=FALSE)
    rownames(geneNodes) <- geneNodes$n
    
    #####################    
    ## Term section
    
    # Offset for the id
    termOffSet <- nbGenes
    nbTerms <- nrow(results$termNodes)
    results$termNodes$term_id <- seq_len(nbTerms) + termOffSet
    
    ## Create node entries for the term
    termNodes <- data.frame('@id'= results$termNodes$term_id,
        n=results$termNodes$id, stringsAsFactors=FALSE, check.names=FALSE)
    
    ## Create node attribute entries for the term alias and the term group
    termAttributes  <- data.frame(
        po=rep(results$termNodes$term_id, 2) ,
        n=c(rep("alias", nbTerms), rep("group", nbTerms)),
        v=c(results$termNodes$alias, rep("TERM", nbTerms)),
        stringsAsFactors=FALSE)
    rownames(termNodes) <- termNodes$n
    
    #####################    
    ## Edge section
    
    # Offset for the id
    edgeOffSet <- termOffSet + nbTerms
    nbEdges <- nrow(results$edges)
    
    edge_gene <- merge(results$edges, results$geneNodes[, c("id", "gene_id")], 
                        by.x="target", by.y="id", all.x=TRUE)
    edge_term <- merge(results$edges, results$termNodes[, c("id", "term_id")], 
                        by.x="source", by.y="id", all.x=TRUE)
    
    ## Create an edge connecting term to gene
    edges <- data.frame('@id'=seq_len(nbEdges) + edgeOffSet,
                        s=edge_term$term_id,
                        t=edge_gene$gene_id,
                        i=rep("contains", nbEdges),
                        stringsAsFactors=FALSE, check.names=FALSE)
    
    ## Create edge attributes
    edgeAttributes <- data.frame(
        po=rep(edges[,'@id'], 3),
        n=c(rep("name", nbEdges), rep("source", nbEdges),
            rep("target", nbEdges)),
        v=c(paste0(results$edges$source, " (contains) ", results$edges$target),
            results$edges$source, results$edges$target),
        stringsAsFactors=FALSE)
    
    rownames(geneNodes) <- NULL
    rownames(termNodes) <- NULL
    
    return(list(nodes=rbind(geneNodes, termNodes), edges=edges,
                nodeAttributes=rbind(geneAttributes, termAttributes),
                edgeAttributes=edgeAttributes))
}


#' @title Create CX JSON text representing the network
#'
#' @description Create a CX JSON text that represent the network which
#' includes information about nodes and edges present in the network.
#'
#' @param nodeEdgeInfo a \code{list} created by gprofiler2 that contains
#' the results from an enrichment analysis. TODO
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
#' data(parentalNapaVsDMSOEnrichment)
#'
#' ## Only retained the GO Molecular Function results
#' results <- parentalNapaVsDMSOEnrichment$result[
#'         parentalNapaVsDMSOEnrichment$result$source == "GO:MF", ]
#'
#' ## Extract node and edge information
#' information <- enrichViewNet:::extractInformationWhenIntersection(
#'         gostResults=results)
#'
#' jsonFormat <- enrichViewNet:::createCXJSONForCytoscape(
#'                 nodeEdgeInfo=information, title="WikiPathways")
#'
#' @author Astrid Deschênes
#' @importFrom jsonlite toJSON
#' @encoding UTF-8
#' @keywords internal
createCXJSONForCytoscape <- function(nodeEdgeInfo, title) {

    entriesL <- formatInformationForCXJSON(results=nodeEdgeInfo)
    
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
