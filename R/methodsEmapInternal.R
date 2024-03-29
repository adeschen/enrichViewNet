
#' @title Validate arguments passed to createEnrichMap() function
#' 
#' @description Validate the arguments passed to createEnrichMap() function.
#' First, the object containing the enrichment results must correspond to a 
#' object created by  \code{gprofiler2} software. Second, the selected 
#' source must at least have one enriched term in the results. Then, if the
#' source is 'TERM_ID', the listed terms must be present in the enrichment
#' results.
#' 
#' @param gostObject a \code{list} created by \code{gprofiler2} that contains
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
#' and "WP" for WikiPathways. 
#' 
#' @param termIDs a \code{vector} of \code{character} strings that contains 
#' the term IDs retained for the creation of the network. This parameter is 
#' only used when \code{source} is set to "TERM_ID".
#' 
#' @param removeRoot a \code{logical} that specified if the root terms of 
#' the selected source should be removed (when present). 
#' 
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{vector} of terms, 
#' the selected terms will be displayed.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped.
#' 
#' @param categoryLabel a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1).
#' 
#' @param categoryNode a positive \code{numeric} representing he amount by 
#' which plotting category nodes should be scaled relative to the default (1).
#' 
#' @param force a \code{logical} indicating if the repulsion between 
#' overlapping text labels should be forced.
#'
#' @return \code{TRUE} when all arguments are valid
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(demoGOST)
#' 
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateEnrichMapArguments(gostObject=demoGOST,
#'     query="query_1", source="GO:BP", termIDs=NULL, removeRoot=FALSE, 
#'     showCategory=20, groupCategory=FALSE, 
#'     categoryLabel=1.1, categoryNode=1, force=TRUE)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @importFrom stringr str_ends
#' @keywords internal
validateCreateEnrichMapArguments <- function(gostObject, query, source, 
        termIDs, removeRoot, showCategory, groupCategory, 
        categoryLabel, categoryNode, force) {
    
    ## Test that gostObject is a gprofiler2 result 
    if (!(inherits(gostObject, "list") && "result" %in% names(gostObject) &&
            "meta" %in% names(gostObject)))   {
        stop("The gostObject object should be a list with meta ", 
                "and result as entries corresponding to gprofiler2 ", 
                "enrichment output.")
    } 
    
    if (!is.character(query) || length(query) > 1) {
        stop("The \'query\'must be a character string.")
    }
    
    ## Query must be in gost object
    gostResults <- as.data.frame(gostObject$result)
    if (!(query %in% unique(gostResults$query))) {
        stop("The \'query\' is not present in the results of the gost object.")
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
    
    result <- validateCreateEnrichMapSubSectionArguments(
        showCategory=showCategory, groupCategory=groupCategory, 
        categoryLabel=categoryLabel, categoryNode=categoryNode, force=force)
    
    return(result)     
}


#' @title Validate arguments passed to createEnrichMapMulti() function
#' 
#' @description Validate the arguments passed to createEnrichMapMulti() 
#' function.
#' First, the object containing the enrichment results must correspond to a 
#' object created by  \code{gprofiler2} software. Second, the selected 
#' source must at least have one enriched term in the results. Then, if the
#' source is 'TERM_ID', the listed terms must be present in the enrichment
#' results.
#' 
#' @param gostObjectList a \code{list} of \code{gprofiler2} objects that 
#' contain the results from an enrichment analysis. The list must contain at 
#' least 2 entries. The number of entries must correspond to the number of
#'  entries for the \code{queryList} parameter.
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
#' and "WP" for WikiPathways. 
#' 
#' @param termIDs a \code{vector} of \code{character} strings that contains 
#' the term IDs retained for the creation of the network. This parameter is 
#' only used when \code{source} is set to "TERM_ID".
#' 
#' @param removeRoot a \code{logical} that specified if the root terms of 
#' the selected source should be removed (when present). 
#' 
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{vector} of terms, 
#' the selected terms will be displayed.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped.
#' 
#' @param categoryLabel a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1).
#' 
#' @param categoryNode a positive \code{numeric} representing he amount by 
#' which plotting category nodes should be scaled relative to the default (1).
#' 
#' @param force a \code{logical} indicating if the repulsion between 
#' overlapping text labels should be forced. 
#' 
#' @return \code{TRUE} when all arguments are valid
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' data(rosaNapaVsDMSOEnrichment)
#' 
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateEnrichMapMultiArguments(
#'     gostObjectList=list(parentalNapaVsDMSOEnrichment, 
#'                             rosaNapaVsDMSOEnrichment),
#'     queryList=list("parental_napa_vs_DMSO", "rosa_napa_vs_DMSO"), 
#'     source="GO:BP", termIDs=NULL, removeRoot=FALSE, 
#'     showCategory=20, groupCategory=FALSE, 
#'     categoryLabel=1.1, categoryNode=1, force=FALSE)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @importFrom stringr str_ends
#' @keywords internal
validateCreateEnrichMapMultiArguments <- function(gostObjectList, queryList, 
    source, termIDs, removeRoot, showCategory, groupCategory, 
    categoryLabel, categoryNode, force) {
    
    ## Test that gostObject is a list with minimum of 2 entries
    if (!inherits(gostObjectList, "list") || !(length(gostObjectList) > 1)) {
        stop("The gostObjectList object should be a list of enrichment", 
                " objects. At least 2 enrichment objects are required.")
    }
    
    ## Test that gostObject is a list of gprofiler2 result 
    if (!(all(unlist(lapply(gostObjectList, is.list))) && 
            all(unlist(lapply(gostObjectList, FUN = function(x) {
                "result" %in% names(x) && "meta" %in% names(x)})))))   {
        stop("The gostObjectList should only contain a list of enrichment ", 
                "results. Enrichment results are lists with meta ", 
                "and result as entries corresponding to gprofiler2 ", 
                "enrichment output.")
    } 
    
    ## Test that queryList is a list with 2 entries minimum
    if (!inherits(queryList, "list") || !(length(queryList) > 1) || 
            length(queryList) != length(gostObjectList)) {
        stop("The queryList object should be a list of query names. At least ", 
                "2 query names are required. The number of query names should ", 
                "correspond to the number of enrichment objects.")
    }
    
    ## Test that queryList is a list of character strings
    if (!all(unlist(lapply(queryList, is.character)))) {
        stop("The queryList object should only contain a list of query names", 
                " in character strings.")
    }
    
    ## Validate that all query names are present in the associated results
    if (!all(unlist(lapply(seq_len(length(queryList)), FUN = function(x, 
        queryL, gostL) {res <- as.data.frame(gostL[[x]]$result); 
        queryL[[x]] %in% unique(res$query)}, queryL=queryList, 
        gostL=gostObjectList)))) {
        stop("Each query name present in the \'queryList'\ parameter must be ", 
                "present in the associated enrichment object.")
    }
    
    if (source != "TERM_ID") {
        if (!any(unlist(lapply(gostObjectList, FUN = function(x, source) {
            sum(x$result$source == source) > 1}, source=source)))) {
                stop("There is no enriched term for the selected ", 
                                "source \'", source, "\'.")    
        }
    } else {
        if (is.null(termIDs)) {
            stop("A vector of terms should be given through the ",
                    "\'termIDs\' parameter when source is \'TERM_ID\'.")
        }
    }

    result <- validateCreateEnrichMapSubSectionArguments(
        showCategory=showCategory, groupCategory=groupCategory, 
        categoryLabel=categoryLabel, categoryNode=categoryNode, force=force)
    
    return(result)   
}


#' @title Validate common arguments passed to createEnrichMap() and 
#' createEnrichMapMulti() functions
#' 
#' @description Validate the common arguments passed to createEnrichMap() and 
#' createEnrichMapMulti() functions.
#' 
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{vector} of terms, 
#' the selected terms will be displayed.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped.
#' 
#' @param categoryLabel a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1).
#' 
#' @param categoryNode a positive \code{numeric} representing he amount by 
#' which plotting category nodes should be scaled relative to the default (1).
#' 
#' @param force a \code{logical} indicating if the repulsion between 
#' overlapping text labels should be forced. 
#' 
#' @return \code{TRUE} when all arguments are valid
#' 
#' @examples
#'
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateEnrichMapSubSectionArguments(
#'     showCategory=20, groupCategory=FALSE, categoryLabel=1.1, categoryNode=1, 
#'     force=TRUE)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @importFrom stringr str_ends
#' @keywords internal
validateCreateEnrichMapSubSectionArguments <- function(showCategory, 
    groupCategory, categoryLabel, categoryNode, force) {
    
    if (!is.character(showCategory) && 
        !(is.numeric(showCategory) && (showCategory > 0))) {
        stop("The \'showCategory\' parameter must an positive integer or a ", 
            "vector of character strings representing terms.")
    }

    if (!is.logical(groupCategory)) {
        stop("The \'groupCategory\' parameter must a logical ", 
            "(TRUE or FALSE).")
    }

    if (!is.numeric(categoryLabel) || !(categoryLabel > 0)) {
        stop("The \'categoryLabel\' parameter must be a positive numeric.")
    }

    if (!is.numeric(categoryNode) || !(categoryNode > 0)) {
        stop("The \'categoryNode\' parameter must be a positive numeric.")
    }
    
    if (!is.logical(force)) {
        stop("The \'force\' parameter must a logical (TRUE or FALSE).")
    }
    
    return(TRUE)   
}


#' @title Create a basic enrichment map
#' 
#' @description The function creates a basic enrichment map using functional 
#' enrichment results.
#' 
#' @param gostResults a \code{data.frame} containing the enrichment 
#' results to be plot.
#' 
#' @param backgroundGenes a \code{vector} of \code{character} string 
#' representing the name of the genes present in the request.
#' 
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{vector} of terms, 
#' the selected terms will be displayed.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped.
#' 
#' @param categoryLabel a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1).
#' 
#' @param categoryNode a positive \code{numeric} representing the amount by 
#' which plotting category nodes should be scaled relative to the default (1). 
#' 
#' @param significantMethod a \code{character} string representing the name 
#' of the multiple testing correction method used on the results.
#' 
#' @param force a \code{logical} indicating if the repulsion between 
#' overlapping text labels should be forced. 
#' 
#' @return a \code{ggplot} object representing the enrichment map.
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' 
#' ## Only retain the results section
#' gostResults <- as.data.frame(parentalNapaVsDMSOEnrichment$result)
#' 
#' ## Limit the results to Wikipathways
#' ## and remove the root term
#' gostResults <- gostResults[which(gostResults$source == "WP"),]
#' gostResults <- gostResults[which(gostResults$term_id != "WIKIPATHWAYS"),]
#' 
#' ## Extract meta data information
#' meta <- parentalNapaVsDMSOEnrichment$meta
#' 
#' ## Get all background genes
#' backgroundGenes <- meta$query_metadata$queries[["parental_napa_vs_DMSO"]]
#' 
#' ## Get significant method
#' significantMethod <- meta$query_metadata$significance_threshold_method
#'
#' ## Create basic enrichment map using Wikipathways terms
#' enrichViewNet:::createBasicEmap(gostResults=gostResults, 
#'     backgroundGenes=backgroundGenes, showCategory=30L, 
#'     groupCategory=FALSE, categoryLabel=1, categoryNode=1,
#'     significantMethod=significantMethod)
#'     
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is new
#' @importFrom stringr str_ends str_split str_replace_all
#' @importFrom enrichplot pairwise_termsim emapplot
#' @importClassesFrom DOSE enrichResult
#' @keywords internal
createBasicEmap <- function(gostResults, backgroundGenes, 
        showCategory, groupCategory, categoryLabel, categoryNode, 
        significantMethod, force) {
    
    ## Extract gene list for each term
    geneSets <- lapply(seq_len(nrow(gostResults)), FUN=function(x, gostData) {
        str_split(gostData$intersection[x], pattern=",")[[1]]}, 
        gostData=gostResults
    )
    names(geneSets) <- gostResults$term_id
    
    resultDataFrame <- data.frame(ID=gostResults$term_id, 
        Description=gostResults$term_name,
        GeneRatio=c(paste0(gostResults$intersection_size, "/", 
                                gostResults$query_size)), 
        BgRatio=c(paste0(gostResults$intersection_size, "/", 
                                gostResults$effective_domain_size)), 
        pvalues=gostResults$p_value,
        p.adjust=gostResults$p_value, 
        qvalue=gostResults$p_value, 
        geneID=str_replace_all(gostResults$intersection, ",", "/"),
        Count=c(gostResults$intersection_size), stringsAsFactors=FALSE)
    
    resultDataFrame <- resultDataFrame[order(resultDataFrame$pvalues), ]
    
    rownames(resultDataFrame) <- resultDataFrame$ID
    
    res <- new("enrichResult", result=resultDataFrame, pvalueCutoff=1, 
                pAdjustMethod="UNKNOWN", qvalueCutoff=1, 
                gene=as.character(backgroundGenes), 
                universe=as.character(backgroundGenes), 
                geneSets=geneSets, 
                organism="UNKNOWN", keytype="UNKNOWN", ontology="UNKNOWN", 
                readable=FALSE)
    
    ## Get similarity matrix
    comp <- pairwise_termsim(res)  
    
    graphEmap <- emapplot(x=comp, showCategory=showCategory,
                        cluster.params=list(cluster=groupCategory),
                        cex.params=list(category_node=categoryNode,
                            category_label=categoryLabel), force=force)
    
    return(graphEmap)
}


#' @title Create a basic enrichment map
#' 
#' @description The function creates a basic enrichment map using functional 
#' enrichment results.
#' 
#' @param gostResults a \code{data.frame} containing the enrichment 
#' results to be plot.
#' 
#' @param backgroundGenes a \code{vector} of \code{character} string 
#' representing the name of the genes present in the request.
#' 
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{vector} of terms, 
#' the selected terms will be displayed.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped.
#' 
#' @param categoryLabel a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1).
#' 
#' @param categoryNode a positive \code{numeric} representing the amount by 
#' which plotting category nodes should be scaled relative to the default (1). 
#' 
#' @param significantMethod a \code{character} string representing the name 
#' of the multiple testing correction method used on the results.
#' 
#' @param force a \code{logical} indicating if the repulsion between 
#' overlapping text labels should be forced.
#' 
#' @return a \code{ggplot} object representing the enrichment map.
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' 
#' ## Only retain the results section
#' gostResults <- as.data.frame(parentalNapaVsDMSOEnrichment$result)
#' 
#' ## Limit the results to Wikipathways and KEGG
#' gostResultsWP <- gostResults[which(gostResults$source == "WP"),]
#' gostResultsKEGG <- gostResults[which(gostResults$source == "KEGG"),]
#' 
#' ## Extract meta data information
#' queryList <- list("parental_napa_vs_DMSO - WP", 
#'         "parental_napa_vs_DMSO - KEGG")
#' 
#' ## Create basic enrichment map using Wikipathways terms
#' enrichViewNet:::createMultiEmap(gostResultsList=list(gostResultsWP, 
#'     gostResultsKEGG), queryList=queryList, showCategory=30L, 
#'     groupCategory=FALSE, categoryLabel=1, categoryNode=1, force=TRUE)
#'     
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is new
#' @importFrom stringr str_ends str_split str_replace_all
#' @importFrom enrichplot pairwise_termsim emapplot
#' @importClassesFrom DOSE compareClusterResult
#' @keywords internal
createMultiEmap <- function(gostResultsList, queryList, showCategory, 
    groupCategory, categoryLabel, categoryNode, force) {
    
    resF <- list()
    geneClusters <- list()
    
    for(i in seq_len(length(gostResultsList))) {
        gostResults <- gostResultsList[[i]]
        resultDF <- data.frame(
            Cluster=rep(queryList[[i]], nrow(gostResults)),
            ID=gostResults$term_id,
            Description=gostResults$term_name,
            GeneRatio=c(paste0(gostResults$intersection_size, "/", 
                                gostResults$query_size)), 
            BgRatio=c(paste0(gostResults$intersection_size, "/", 
                                gostResults$effective_domain_size)), 
            pvalues=gostResults$p_value,
            p.adjust=gostResults$p_value, qvalue=gostResults$p_value, 
            geneID=str_replace_all(gostResults$intersection, ",", "/"),
            Count=c(gostResults$intersection_size), stringsAsFactors=FALSE)
        resF[[i]] <- resultDF
        
        geneClusters[[queryList[[i]]]] <- 
            unique(unlist(stringr::str_split(gostResults$intersection, ",")))
    }
    clProfDF <- do.call(rbind, resF)
    
    clProfDF$Cluster = factor(clProfDF$Cluster)
    clProfDF$geneID <- stringr::str_replace_all(string = clProfDF$geneID,
                                                pattern = ",", "/") 
    
    ## Problem when identical description 
    if (any(table(clProfDF$Description) > 1)) {
        test <- which(table(clProfDF$Description) > 1)
        for (i in names(test)) {
            resTest <- clProfDF[clProfDF$Description == i, ]
            if (length(unique(resTest$ID)) != 1) {
                clProfDF$Description[clProfDF$Description == i] <- 
                    paste0(clProfDF$Description[clProfDF$Description == i], 
                            " (", clProfDF$ID[clProfDF$Description == i], ")")   
            }
        }
    }
    
    res <- new("compareClusterResult",
               compareClusterResult = clProfDF,
               geneClusters = geneClusters,
               fun = "createMultiEmap",
               .call = call("createMultiEmap()")
    )
    res@keytype <- "UNKNOWN"
    res@readable <- FALSE
    
    kegg_compar <- pairwise_termsim(res)  
    
    graphEmap <- emapplot(kegg_compar, 
        showCategory=length(table(clProfDF$Description)),
        cluster.params=list(cluster=FALSE),
        cex.params=list(category_node=categoryNode, 
                            category_label=categoryLabel),
        force=force, alpha=0.9, repel=TRUE)
    
    return(graphEmap)
}
