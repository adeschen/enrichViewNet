
#' @title Using functional enrichment results in  gprofiler2 format to create  
#' an enrichment map
#' 
#' @description User selected enrichment terms are used to create an enrichment 
#' map. The selection of the term can by specifying by the 
#' source of the terms (GO:MF, REAC, TF, etc.) or by listing the selected 
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
#' @param showCategory a positive \code{integer} representing the maximum 
#' number of terms to display.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{NULL}, all terms  
#' will be displayed. Default: \code{30L}.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped. Default: \code{FALSE}.
#' 
#' @param categoryLabel a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1). Default: \code{1}.s
#' 
#' @param categoryNode a positive \code{numeric} representing the amount by 
#' which plotting category nodes should be scaled relative to the default (1).
#' Default: \code{1}.
#'
#' @param line a non-negative \code{numeric} representing the scale of line 
#' width. Default: \code{1}.
#'
#' @param ... additional arguments that will be pass to the 
#' \code{\link[enrichplot]{emapplot}} function. 
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
#'     query=query, source="GO:CC", removeRoot=TRUE)
#' 
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom strex match_arg
#' @encoding UTF-8
#' @export
createEnrichMap <- function(gostObject, query, source=c("TERM_ID", "GO:MF", 
        "GO:CC", "GO:BP", "KEGG", "REAC", "TF", "MIRNA", "HPA", "CORUM", 
        "HP", "WP"), termIDs=NULL, removeRoot=TRUE,  
        showCategory=30L, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, ...) {
    
    ## Validate source is among the possible choices
    source <- match_arg(source, ignore_case=TRUE)
    
    ## Validate parameters
    ## TODO update showCategory
    validateCreateEnrichMapArguments(gostObject=gostObject, query=query, 
        source=source, termIDs=termIDs, removeRoot=removeRoot, 
        showCategory=showCategory, categoryLabel=categoryLabel,
        groupCategory=groupCategory, categoryNode=categoryNode, 
        line=line)
    
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
                backgroundGenes=backgroundGenes, 
                showCategory=showCategory, categoryLabel=categoryLabel,
                groupCategory=groupCategory, categoryNode=categoryNode, 
                significantMethod=significantMethod, line=line,  ...)
    
    return(emap)
}


#' @title Using functional enrichment results in  gprofiler2 format to create  
#' an enrichment map in an igraph format
#' 
#' @description User selected enrichment terms are used to create an enrichment 
#' map. The selection of the term can by specifying by the 
#' source of the terms (GO:MF, REAC, TF, etc.) or by listing the selected 
#' term IDs. The map is only generated when there is at least on 
#' significant term to graph. The output is an \code{igraph}.
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
#' @param showCategory a positive \code{integer} representing the maximum 
#' number of terms to display.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{NULL}, all terms  
#' will be displayed. Default: \code{30L}.
#'  
#' @param similarityCutOff a positive \code{numeric} between 0 and 1 indicating 
#' the minimum level of similarity between two terms to have an edge linking 
#' the terms. Default: \code{0.20}.  
#' 
#' @return a \code{igraph} object which is the enrichment map for enrichment 
#' results. The node have 2 attributes: "name" and "size". The "name" 
#' corresponds to the term description. While the "size" corresponds to the 
#' number of genes found in the specific gene set. The edges have 
#' 3 attributes: "similarity", "width", and "weight". All those 3 attributes 
#' correspond to the Jaccard coefficient.
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
#' mapG <- createEnrichMapAsIgraph(gostObject=parentalNapaVsDMSOEnrichment, 
#'     query=query, source="GO:CC", removeRoot=TRUE, 
#'     showCategory=30L, similarityCutOff=0.20)
#' 
#' ## Required library igraph to show the graph
#' if(requireNamespace("igraph", quietly=TRUE)) {
#'     ## Using library igraph to show the graph
#'     library(igraph)
#'     plot(mapG)
#' }
#' 
#' if (requireNamespace("ggplot2", quietly=TRUE) &&
#'         requireNamespace("igraph", quietly=TRUE) &&
#'         requireNamespace("ggtangle", quietly=TRUE) &&
#'         requireNamespace("ggnetwork", quietly=TRUE)) {
#'     ## Using more complex set of libraries to display personalized graph
#'     library(ggplot2)
#'     library(igraph)
#'     library(ggnetwork)
#'     library(ggtangle)
#'     emapG <- ggplot(mapG, layout=layout_with_fr) +
#'                 geom_edge(color="gray", linewidth=1) + 
#'                 geom_nodes(aes(size=size)) + 
#'                 geom_nodetext(aes(label=name), color="black", size=3) +
#'                 theme_void()
#'     emapG
#' }
#' 
#' @author Astrid Deschênes
#' @importFrom strex match_arg
#' @encoding UTF-8
#' @export
createEnrichMapAsIgraph <- function(gostObject, query, source=c("TERM_ID", 
    "GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "TF", "MIRNA", "HPA", "CORUM", 
    "HP", "WP"), termIDs=NULL, removeRoot=TRUE,  
    showCategory=30L, similarityCutOff=0.20) {
    
    ## Validate source is among the possible choices
    source <- match_arg(source, ignore_case=TRUE)
    
    ## Validate parameters
    validateCreateEnrichMapAsIgraphArg(gostObject=gostObject, query=query, 
        source=source, termIDs=termIDs, removeRoot=removeRoot, 
        showCategory=showCategory, similarityCutOff=similarityCutOff)
    
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
    
    backgroundGenes <- gostObject$meta$query_metadata$queries[[query]]
    
    ## Create basic enrichment map in igraph format
    emapIgraph <- createBasicEmapAsIgraph(gostResults=gostResults, 
                backgroundGenes=backgroundGenes, showCategory=showCategory,
                similarityCutOff=similarityCutOff)
    
    return(emapIgraph)
}


#' @title Using functional enrichment results in gprofiler2 format to create  
#' an enrichment map with multiple groups from different enrichment analyses
#' 
#' @description User selected enrichment terms are used to create an enrichment 
#' map. The selection of the term can by specifying by the 
#' source of the terms (GO:MF, REAC, TF, etc.) or by listing the selected 
#' term IDs. The map is only generated when there is at least on 
#' significant term to graph.
#' 
#' @param gostObjectList a \code{list} of \code{gprofiler2} objects that 
#' contain the results from an enrichment analysis. The list must contain at 
#' least 2 entries. The number of entries must correspond to the number of 
#' entries for the \code{queryList} parameter.
#' 
#' @param queryList a \code{list} of \code{character} strings representing the 
#' names of the queries that are going to be used to generate the graph. 
#' The query names must exist in the associated \code{gostObjectList} objects 
#' and follow the same order. The number of entries must correspond to the 
#' number of entries for the \code{gostObjectList} parameter.
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
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{vector} of terms, 
#' the selected terms will be displayed. Default: \code{30L}.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped. Default: \code{FALSE}.
#' 
#' @param categoryLabel a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1). Default: \code{1}.
#' 
#' @param categoryNode a positive \code{numeric} representing the amount by 
#' which plotting category nodes should be scaled relative to the default (1).
#' Default: \code{1}.
#' 
#' @param line a non-negative \code{numeric} representing the scale of line 
#' width. Default: \code{1}.
#' 
#' @param ... additional arguments that will be pass to the 
#' \code{\link[enrichplot]{emapplot}} function. 
#'
#' @return a \code{ggplot} object which is the enrichment map for enrichment 
#' results.
#' 
#' @examples
#'
#' ## Loading dataset containing results from 2 enrichment analyses done with
#' ## gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' data(rosaNapaVsDMSOEnrichment)
#'
#' ## Extract query information (only one in each dataset)
#' query1 <- unique(parentalNapaVsDMSOEnrichment$result$query)[1]
#' query2 <- unique(rosaNapaVsDMSOEnrichment$result$query)[1]
#' 
#' ## Create graph for KEGG related results from 
#' ## 2 enrichment analyses
#' createEnrichMapMultiBasic(gostObjectList=list(parentalNapaVsDMSOEnrichment, 
#'     rosaNapaVsDMSOEnrichment), 
#'     queryList=list(query1, query2), source="KEGG", removeRoot=TRUE)
#' 
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom strex match_arg
#' @encoding UTF-8
#' @export
createEnrichMapMultiBasic <- function(gostObjectList, queryList, 
    source=c("TERM_ID", "GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "TF", 
    "MIRNA", "HPA", "CORUM", "HP", "WP"), termIDs=NULL, removeRoot=TRUE, 
    showCategory=30L, groupCategory=FALSE, categoryLabel=1, 
    categoryNode=1, line=1, ...) {
    
    ## Validate source is among the possible choices
    source <- match_arg(source, ignore_case=TRUE)
    
    ## Validate parameters
    validateCreateEnrichMapMultiBasicArgs(gostObjectList=gostObjectList, 
        queryList=queryList, source=source, termIDs=termIDs, 
        removeRoot=removeRoot, showCategory=showCategory, 
        categoryLabel=categoryLabel, groupCategory=groupCategory, 
        categoryNode=categoryNode, line=line)
    
    ## Extract results
    gostResultsList <- lapply(gostObjectList, FUN=function(x) {x$result})
    
    ## Retain results associated to query
    gostResultsList <- lapply(seq_len(length(queryList)), 
        FUN=function(i, queryL, gostL) {
            gostL[[i]][gostL[[i]]$query == queryL[[i]], ]}, 
        queryL=queryList, gostL=gostResultsList)
    
    ## Filter results
    if (source == "TERM_ID") {
        gostResultsList <- lapply(gostResultsList,  FUN=function(x, termIDs) {
            x[x$term_id %in% termIDs,]}, termIDs=termIDs)
    } else {
        gostResultsList <- lapply(gostResultsList,  FUN=function(x, source) {
            x[x$source == source,]}, source=source)
    }
    
    ## Remove root term if required
    if (removeRoot) {
        gostResultsList <- lapply(gostResultsList,  FUN=function(x) {
            removeRootTerm(x)})
    }
    
    ## Validate that at least one term is left
    if (sum(unlist(lapply(gostResultsList,  FUN=function(x) {nrow(x)})))  
            == 0) {
        stop("With removal of the root term, there is no ", 
                    "enrichment term left")
    }
    
    ## Create multi categories emap
    emap <- createMultiEmap(gostResultsList=gostResultsList, 
                queryList=queryList, showCategory=showCategory, 
                categoryLabel=categoryLabel, groupCategory=groupCategory, 
                categoryNode=categoryNode, line=line, ...)
    
    return(emap)
}


#' @title Using functional enrichment results in gprofiler2 format to create  
#' an enrichment map with multiple groups from different enrichment analyses as
#' an igraph
#' 
#' @description User selected enrichment terms are used to create an enrichment 
#' map. The selection of the term can by specifying by the 
#' source of the terms (GO:MF, REAC, TF, etc.) or by listing the selected 
#' term IDs. The map is only generated when there is at least on 
#' significant term to graph. The output is an enrichment map in igraph format.
#' 
#' @param gostObjectList a \code{list} of \code{gprofiler2} objects that 
#' contain the results from an enrichment analysis. The list must contain at 
#' least 2 entries. The number of entries must correspond to the number of 
#' entries for the \code{queryList} parameter.
#' 
#' @param queryList a \code{list} of \code{character} strings representing the 
#' names of the queries that are going to be used to generate the graph. 
#' The query names must exist in the associated \code{gostObjectList} objects 
#' and follow the same order. The number of entries must correspond to the 
#' number of entries for the \code{gostObjectList} parameter.
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
#' @param showCategory a positive \code{integer} representing the maximum 
#' number of terms to display.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{NULL}, all terms  
#' will be displayed. Default: \code{30L}.
#'  
#' @param similarityCutOff a positive \code{numeric} between 0 and 1 indicating 
#' the minimum level of similarity between two terms, calculated 
#' using the Jaccard coefficient, to have an edge linking 
#' the terms. Default: \code{0.20}.  
#' 
#' @return a \code{igraph} object which is the enrichment map for enrichment 
#' results. The node have 5 attributes: "name", "size", "pie", "cluster", 
#' and "pieName". The "name" corresponds to the term description. While the 
#' "size" corresponds to the number of unique genes found in the specific 
#' gene set when looking at all the experiments.
#' The edges have 3 attributes: "similarity", "width", and 
#' "weight". All those 3 attributes correspond to the Jaccard coefficient.
#' 
#' @examples
#'
#' ## Loading dataset containing results from 2 enrichment analyses done with
#' ## gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' data(rosaNapaVsDMSOEnrichment)
#'
#' ## Extract query information (only one in each dataset)
#' query1 <- unique(parentalNapaVsDMSOEnrichment$result$query)[1]
#' query2 <- unique(rosaNapaVsDMSOEnrichment$result$query)[1]
#' 
#' ## Create graph for KEGG related results from 
#' ## 2 enrichment analyses
#' emapG <- createEnrichMapMultiBasicAsIgraph(gostObjectList=
#'     list(parentalNapaVsDMSOEnrichment, rosaNapaVsDMSOEnrichment), 
#'     queryList=list(query1, query2), source="KEGG", removeRoot=TRUE, 
#'     similarityCutOff=0.3)
#' 
#' if (requireNamespace("ggplot2", quietly=TRUE) && 
#'         requireNamespace("igraph", quietly=TRUE) && 
#'         requireNamespace("scatterpie", quietly=TRUE) && 
#'         requireNamespace("ggtangle", quietly=TRUE) && 
#'         requireNamespace("ggrepel", quietly=TRUE)) {
#'     ## Create a visual representation of the enrichment map
#'     ## by default
#'     library(igraph)
#'     plot(emapG)
#'     
#'     ## Add see to reproduce the same graph
#'     set.seed(12)
#'     
#'     library(ggplot2)
#'     library(ggtangle)
#'     library(scatterpie)
#'     library(ggrepel)
#'     
#'     emapGraph <- ggplot(emapG, layout=layout_with_fr) + 
#'                 geom_edge(color="gray", linewidth=1)
#'     
#'     pieInfo <- as.data.frame(do.call(rbind, V(emapG)$pie))
#'     colnames(pieInfo) <- V(emapG)$pieName[[1]]
#'     
#'     ## Add information about the groups associated with each node in the 
#'     ## ggplot object so that the node can be colored accordingly
#'     for (i in seq_len(ncol(pieInfo))) {
#'         emapGraph$data[colnames(pieInfo)[i]] <- pieInfo[, i]
#'     }
#'     
#'     ## Using scatterpie, ggtangle and ggrepel to generate the graph
#'     ## geom_scatterpie() allows to have scatter pie plot
#'     ## geom_text_repel() allows to have minimum overlying terms
#'     ## coord_fixed() forces the plot to have a 1:1 aspect ratio
#'     emapGraph + 
#'         geom_scatterpie(aes(x=x, y=y, r=size/50), 
#'             cols=c(colnames(pieInfo)), legend_name = "Group", color=NA) +
#'         geom_scatterpie_legend(radius=emapGraph$data$size/50, n=4, 
#'             x=max(emapGraph$data$x), y=min(emapGraph$data$y),
#'             labeller=function(x) {round(x*50)}, label_position="right") +
#'         geom_text_repel(aes(x=x, y=y, label=label)) +
#'         coord_fixed()
#' }
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom strex match_arg
#' @encoding UTF-8
#' @export
createEnrichMapMultiBasicAsIgraph <- function(gostObjectList, queryList, 
    source=c("TERM_ID", "GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "TF", 
    "MIRNA", "HPA", "CORUM", "HP", "WP"), termIDs=NULL, removeRoot=TRUE, 
    showCategory=30L, similarityCutOff=0.2) {
    
    ## Validate source is among the possible choices
    source <- match_arg(source, ignore_case=TRUE)
    
    ## Validate parameters
    validateCreateEnrichMapMultiBasicAsIgraphArgs(gostObjectList=gostObjectList, 
        queryList=queryList, source=source, termIDs=termIDs, 
        removeRoot=removeRoot, showCategory=showCategory, 
        similarityCutOff=similarityCutOff)
    
    ## Extract results
    gostResultsList <- lapply(gostObjectList, FUN=function(x) {x$result})
    
    ## Retain results associated to query
    gostResultsList <- lapply(seq_len(length(queryList)), 
        FUN=function(i, queryL, gostL) {
            gostL[[i]][gostL[[i]]$query == queryL[[i]], ]}, 
            queryL=queryList, gostL=gostResultsList)
    
    ## Filter results
    if (source == "TERM_ID") {
        gostResultsList <- lapply(gostResultsList,  FUN=function(x, termIDs) {
            x[x$term_id %in% termIDs,]}, termIDs=termIDs)
    } else {
        gostResultsList <- lapply(gostResultsList,  FUN=function(x, source) {
            x[x$source == source,]}, source=source)
    }
    
    ## Remove root term if required
    if (removeRoot) {
        gostResultsList <- lapply(gostResultsList,  FUN=function(x) {
            removeRootTerm(x)})
    }
    
    ## Validate that at least one term is left
    if (sum(unlist(lapply(gostResultsList,  FUN=function(x) {nrow(x)})))  
        == 0) {
        stop("With removal of the root term, there is no ", 
                "enrichment term left")
    }
    
    ## Create multi categories emap
    emap <- createMultiEmapAsIgraph(gostResultsList=gostResultsList, 
        queryList=queryList, showCategory=showCategory, 
        similarityCutOff=similarityCutOff)
    
    return(emap)
}


#' @title Using functional enrichment results in gprofiler2 format to create  
#' an enrichment map with multiple groups from same or different enrichment 
#' analyses
#' 
#' @description User selected enrichment terms are used to create an enrichment 
#' map. The selection of the term can by specifying by the 
#' source of the terms (GO:MF, REAC, TF, etc.) or by listing the selected 
#' term IDs. The map is only generated when there is at least on 
#' significant term to graph.
#' 
#' @param gostObjectList a \code{list} of \code{gprofiler2} objects that 
#' contain the results from an enrichment analysis. The list must contain at 
#' least 2 entries. The number of entries must correspond to the number of 
#' entries for the \code{queryList} parameter.
#' 
#' @param queryInfo a \code{data.frame} contains one row per group being 
#' displayed. The number of rows must correspond to the 
#' number of entries for the \code{gostObjectList} parameter. 
#' The mandatory columns are:
#' \itemize{
#' \item{\code{queryName}: a \code{character} string representing the name 
#' of the query retained for this group). The query names must exist in the 
#' associated \code{gostObjectList} objects and follow the same order. }
#' \item{\code{source}: a \code{character} string representing the selected 
#' source that will be used to generate the network. To hand-pick the terms to 
#' be used, "TERM_ID" should be used and the list of selected term IDs should
#' be passed through the \code{termIDs} parameter. The possible sources are 
#' "GO:BP" for Gene Ontology Biological Process, "GO:CC" for Gene Ontology  
#' Cellular Component, "GO:MF" for Gene Ontology Molecular Function, 
#' "KEGG" for Kegg, "REAC" for Reactome, "TF" for TRANSFAC, "MIRNA" for 
#' miRTarBase, "CORUM" for CORUM database, "HP" for Human phenotype ontology
#' and "WP" for WikiPathways.  Default: "TERM_ID". }
#' \item{\code{removeRoot}: a \code{logical} that specified if the root terms 
#' of the selected source should be removed (when present). }
#' \item{\code{termIDs}: a \code{character} strings that contains the
#' term IDS retained for the creation of the network separated by a comma ',' 
#' when the "TERM_ID" source is selected. Otherwise, it should be a empty 
#' string (""). }
#' \item{\code{groupName}: a \code{character} strings that contains the 
#' name of the group to be shown in the legend. Each group has to have a 
#' unique name. }
#' }
#' 
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{vector} of terms, 
#' the selected terms will be displayed. Default: \code{30L}.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped. Default: \code{FALSE}.
#' 
#' @param categoryLabel a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1). Default: \code{1}.
#' 
#' @param categoryNode a positive \code{numeric} representing the amount by 
#' which plotting category nodes should be scaled relative to the default (1).
#' Default: \code{1}.
#' 
#' @param line a non-negative \code{numeric} representing the scale of line 
#' width. Default: \code{1}.
#' 
#' @param ... additional arguments that will be pass to the 
#' \code{\link[enrichplot]{emapplot}} function. 
#'
#' @return a \code{ggplot} object which is the enrichment map for enrichment 
#' results.
#' 
#' @examples
#'
#' ## Loading dataset containing results from 2 enrichment analyses done with
#' ## gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' data(rosaNapaVsDMSOEnrichment)
#'
#' ## Create list of enrichment objects required to generate the enrichment map
#' gostObjectList=list(parentalNapaVsDMSOEnrichment, 
#'     parentalNapaVsDMSOEnrichment, rosaNapaVsDMSOEnrichment, 
#'     rosaNapaVsDMSOEnrichment)
#'     
#' ## Create data frame containing required information enabling the 
#' ## selection of the retained enriched terms for each enrichment analysis.
#' ## One line per enrichment analyses present in the gostObjectList parameter
#' ## With this data frame, the enrichment results will be split in 4 groups:
#' ## 1) KEGG significant terms from parental napa vs DMSO (no root term)
#' ## 2) REACTOME significant terms from parental napa vs DMSO (no root term)
#' ## 3) KEGG significant terms from rosa napa vs DMSO (no root term)
#' ## 4) REACTOME significant terms from rosa napa vs DMSO (no root term)
#' queryDataFrame <- data.frame(queryName=c("parental_napa_vs_DMSO", 
#'         "parental_napa_vs_DMSO", "rosa_napa_vs_DMSO", "rosa_napa_vs_DMSO"), 
#'     source=c("KEGG", "REAC", "KEGG", "REAC"), 
#'     removeRoot=c(TRUE, TRUE, TRUE, TRUE), termIDs=c("", "", "", ""), 
#'     groupName=c("parental - KEGG", "parental - Reactome", 
#'         "rosa - KEGG", "rosa - Reactome"), stringsAsFactors=FALSE)
#'     
#' ## Create graph for KEGG and REACTOME significant results from 
#' ## 2 enrichment analyses
#' createEnrichMapMultiComplex(gostObjectList=gostObjectList, 
#'     queryInfo=queryDataFrame, line=1.5)
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom strex match_arg
#' @importFrom stringr str_split str_trim
#' @encoding UTF-8
#' @export
createEnrichMapMultiComplex <- function(gostObjectList, queryInfo,  
    showCategory=30L, groupCategory=FALSE, categoryLabel=1, 
    categoryNode=1, line=1, ...) {
    
    ## Validate parameters
    validateCreateEnrichMapMultiComplexArg(gostObjectList=gostObjectList, 
        queryInfo=queryInfo, showCategory=showCategory, 
        categoryLabel=categoryLabel, groupCategory=groupCategory, 
        categoryNode=categoryNode, line=line)
    
    ## Extract results
    gostResultsList <- lapply(gostObjectList, FUN=function(x) {x$result})
    
    ## Retain results associated to query
    gostResultsList <- lapply(seq_len(length(gostObjectList)), 
        FUN=function(i, queryI, gostL) {
            gostL[[i]][gostL[[i]]$query == queryI$queryName[[i]], ]}, 
            queryI=queryInfo, gostL=gostResultsList)
    
    ## Filter results
    gostResultsList <- lapply(seq_len(length(gostObjectList)), 
        FUN=function(i, queryI, gostL) {
            if (queryI$source[i] == "TERM_ID") {
                terms <- unlist(str_split(queryI$termIDs[i], ","))
                terms <- str_trim(terms)
                return(gostL[[i]][gostL[[i]]$term_id %in% terms, ])
            } else {
                return(gostL[[i]][gostL[[i]]$source == queryI$source[i], ])
            }
        }, queryI=queryInfo, gostL=gostResultsList)
    
    
    ## Remove root terms when requested
    gostResultsList <- lapply(seq_len(length(gostObjectList)), 
        FUN=function(i, queryI, gostL) {
            if (queryI$removeRoot[i]) {
                return(removeRootTerm(gostL[[i]]))
            } else {
                return(gostL[[i]])
            }
        }, queryI=queryInfo, gostL=gostResultsList)
    
    ## Validate that at least one term is left
    if (sum(unlist(lapply(gostResultsList,  FUN=function(x) {nrow(x)})))
            == 0) {
        stop("With removal of the root term, there is no ",
                "enrichment term left")
    }
    
    ## Create multi categories emap
    emap <- createMultiEmap(gostResultsList=gostResultsList, 
                queryList=queryInfo$groupName, showCategory=showCategory, 
                categoryLabel=categoryLabel, groupCategory=groupCategory, 
                categoryNode=categoryNode, line=line,  ...)
    
    return(emap)
}


#' @title Using functional enrichment results in gprofiler2 format to create  
#' an enrichment map with multiple groups from same or different enrichment 
#' analyses in an igraph format
#' 
#' @description User selected enrichment terms are used to create an enrichment 
#' map. The selection of the term can by specifying by the 
#' source of the terms (GO:MF, REAC, TF, etc.) or by listing the selected 
#' term IDs. The map is only generated when there is at least on 
#' significant term to graph.
#' 
#' @param gostObjectList a \code{list} of \code{gprofiler2} objects that 
#' contain the results from an enrichment analysis. The list must contain at 
#' least 2 entries. The number of entries must correspond to the number of 
#' entries for the \code{queryList} parameter.
#' 
#' @param queryInfo a \code{data.frame} contains one row per group being 
#' displayed. The number of rows must correspond to the 
#' number of entries for the \code{gostObjectList} parameter. 
#' The mandatory columns are:
#' \itemize{
#' \item{\code{queryName}: a \code{character} string representing the name 
#' of the query retained for this group). The query names must exist in the 
#' associated \code{gostObjectList} objects and follow the same order. }
#' \item{\code{source}: a \code{character} string representing the selected 
#' source that will be used to generate the network. To hand-pick the terms to 
#' be used, "TERM_ID" should be used and the list of selected term IDs should
#' be passed through the \code{termIDs} parameter. The possible sources are 
#' "GO:BP" for Gene Ontology Biological Process, "GO:CC" for Gene Ontology  
#' Cellular Component, "GO:MF" for Gene Ontology Molecular Function, 
#' "KEGG" for Kegg, "REAC" for Reactome, "TF" for TRANSFAC, "MIRNA" for 
#' miRTarBase, "CORUM" for CORUM database, "HP" for Human phenotype ontology
#' and "WP" for WikiPathways.  Default: "TERM_ID". }
#' \item{\code{removeRoot}: a \code{logical} that specified if the root terms 
#' of the selected source should be removed (when present). }
#' \item{\code{termIDs}: a \code{character} strings that contains the
#' term IDS retained for the creation of the network separated by a comma ',' 
#' when the "TERM_ID" source is selected. Otherwise, it should be a empty 
#' string (""). }
#' \item{\code{groupName}: a \code{character} strings that contains the 
#' name of the group to be shown in the legend. Each group has to have a 
#' unique name. }
#' }
#' 
#' @param showCategory a positive \code{integer} or \code{NULL}. 
#' If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{NULL}, all terms will be 
#' displayed. Default: \code{30L}.
#'
#' @param similarityCutOff a positive \code{numeric}, larger than zero and
#' small than 1 that represent the minimum similarity level between two 
#' nodes (terms) to be linked by an edge. Default: \code{0.2}. 
#' 
#' @return  a \code{igraph} object which is the enrichment map for enrichment 
#' results. The node have 5 attributes: "name", "size", "pie", "cluster", 
#' and "pieName". The "name" corresponds to the term description. While the 
#' "size" corresponds to the number of unique genes found in the specific 
#' gene set when looking at all the experiments. 
#' The edges have 3 attributes: "similarity", "width", and 
#' "weight". All those 3 attributes correspond to the Jaccard coefficient.
#' 
#' @examples
#'
#' ## Loading dataset containing results from 2 enrichment analyses done with
#' ## gprofiler2 queries
#' data(parentalNapaVsDMSOEnrichment)
#' data(rosaNapaVsDMSOEnrichment)
#'
#' ## The graph will be split in 4 groups
#' ## Groups 1 and 2 are using the parental Napa vs DMSO dataset
#' ## Group 3 and 4 are using the rosa Napa vs DMSO dataset
#' gostObjectList=list(parentalNapaVsDMSOEnrichment, 
#'     parentalNapaVsDMSOEnrichment, rosaNapaVsDMSOEnrichment, 
#'     rosaNapaVsDMSOEnrichment)
#'     
#' ## Create data frame containing required information enabling the 
#' ## selection of the retained enriched terms for each enrichment analysis.
#' ## One line per enrichment analyses present in the gostObjectList parameter
#' ## With this data frame, the enrichment results will be split in 4 groups:
#' ## 1) KEGG significant terms from parental napa vs DMSO (no root term)
#' ## 2) REACTOME significant terms from parental napa vs DMSO (no root term)
#' ## 3) KEGG significant terms from rosa napa vs DMSO (no root term)
#' ## 4) REACTOME significant terms from rosa napa vs DMSO (no root term)
#' queryDataFrame <- data.frame(queryName=c("parental_napa_vs_DMSO", 
#'         "parental_napa_vs_DMSO", "rosa_napa_vs_DMSO", "rosa_napa_vs_DMSO"), 
#'     source=c("KEGG", "REAC", "KEGG", "REAC"), 
#'     removeRoot=c(TRUE, TRUE, TRUE, TRUE), termIDs=c("", "", "", ""), 
#'     groupName=c("parental - KEGG", "parental - Reactome", 
#'         "rosa - KEGG", "rosa - Reactome"), stringsAsFactors=FALSE)
#'     
#' ## Create graph for KEGG and REACTOME significant results from 
#' ## 2 enrichment analyses in an igraph format
#' emap <- createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjectList, 
#'     queryInfo=queryDataFrame, showCategory=5)
#' 
#' if (requireNamespace("ggplot2", quietly=TRUE) && 
#'         requireNamespace("igraph", quietly=TRUE) && 
#'         requireNamespace("scatterpie", quietly=TRUE) && 
#'         requireNamespace("ggtangle", quietly=TRUE) && 
#'         requireNamespace("ggrepel", quietly=TRUE)) {
#'     ## Create a visual representation of the enrichment map
#'     ## by default
#'     library(igraph)
#'     plot(emap)
#'     
#'     ## Add see to reproduce the same graph
#'     set.seed(12)
#'     
#'     library(ggplot2)
#'     library(ggtangle)
#'     library(scatterpie)
#'     library(ggrepel)
#'     
#'     emapG <- ggplot(emap, layout=layout_with_fr) + 
#'                 geom_edge(color="gray", linewidth=1)
#'     
#'     pieInfo <- as.data.frame(do.call(rbind, V(emap)$pie))
#'     colnames(pieInfo) <- V(emap)$pieName[[1]]
#'     
#'     ## Add information about the groups associated with each node in the 
#'     ## ggplot object so that the node can be colored accordingly
#'     for (i in seq_len(ncol(pieInfo))) {
#'         emapG$data[colnames(pieInfo)[i]] <- pieInfo[, i]
#'     }
#'     
#'     ## Using scatterpie, ggtangle and ggrepel to generate the graph
#'     ## geom_scatterpie() allows to have scatter pie plot
#'     ## geom_text_repel() allows to have minimum overlying terms
#'     ## coord_fixed() forces the plot to have a 1:1 aspect ratio
#'     emapG + 
#'         geom_scatterpie(aes(x=x, y=y, r=size/50), 
#'             cols=c(colnames(pieInfo)), legend_name = "Group", color=NA) +
#'         geom_scatterpie_legend(radius=emapG$data$size/50, n=4, 
#'             x=max(emapG$data$x), y=max(emapG$data$y),
#'             labeller=function(x) {round(x*50)}, label_position="right") +
#'         geom_text_repel(aes(x=x, y = y, label=label)) +
#'         coord_fixed()
#' }
#' 
#' @author Astrid Deschênes
#' @importFrom gprofiler2 gconvert
#' @importFrom strex match_arg
#' @importFrom stringr str_split str_trim
#' @encoding UTF-8
#' @export
createEnrichMapMultiComplexAsIgraph <- function(gostObjectList, queryInfo,  
    showCategory=30L, similarityCutOff=0.2) {
    
    ## Validate parameters
    validateCreateEnrichMapMultiComplexAsIgraphArg(
        gostObjectList=gostObjectList, queryInfo=queryInfo, 
        showCategory=showCategory, similarityCutOff=similarityCutOff)
    
    ## Extract results
    gostResultsList <- lapply(gostObjectList, FUN=function(x) {x$result})
    
    ## Retain results associated to query
    gostResultsList <- lapply(seq_len(length(gostObjectList)), 
            FUN=function(i, queryI, gostL) {
                    gostL[[i]][gostL[[i]]$query == queryI$queryName[[i]], ]}, 
                    queryI=queryInfo, gostL=gostResultsList)
    
    ## Filter results
    gostResultsList <- lapply(seq_len(length(gostObjectList)), 
        FUN=function(i, queryI, gostL) {
            if (queryI$source[i] == "TERM_ID") {
                    terms <- unlist(str_split(queryI$termIDs[i], ","))
                    terms <- str_trim(terms)
                    return(gostL[[i]][gostL[[i]]$term_id %in% terms, ])
            } else {
                    return(gostL[[i]][gostL[[i]]$source == queryI$source[i], ])
            }
        }, queryI=queryInfo, gostL=gostResultsList)
    
    ## Remove root terms when requested
    gostResultsList <- lapply(seq_len(length(gostObjectList)), 
            FUN=function(i, queryI, gostL) {
                if (queryI$removeRoot[i]) {
                    return(removeRootTerm(gostL[[i]]))
                } else {
                    return(gostL[[i]])
                }
        }, queryI=queryInfo, gostL=gostResultsList)
    
    ## Validate that at least one term is left
    if (sum(unlist(lapply(gostResultsList,  FUN=function(x) {nrow(x)})))
        == 0) {
        stop("With removal of the root term, there is no ",
                "enrichment term left")
    }
    
    ## Create multi categories emap in a igraph format
    emap <- createMultiEmapAsIgraph(gostResultsList=gostResultsList, 
        queryList=queryInfo$groupName, showCategory=showCategory, 
        similarityCutOff=similarityCutOff)
    
    return(emap)
}

