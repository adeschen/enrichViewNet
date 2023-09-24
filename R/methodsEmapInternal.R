
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
#' @param title a \code{character} string representing the name TODO
#' 
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{vector} of terms, 
#' the selected terms will be displayed.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped.
#' 
#' @param cexLabelCategory a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1).
#' 
#' @param cexCategory a positive \code{numeric} representing he amount by which 
#' plotting category nodes should be scaled relative to the default (1).
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
#'     title="GO:BP", showCategory=20, groupCategory=FALSE, 
#'     cexLabelCategory=1.1, cexCategory=1)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @importFrom stringr str_ends
#' @keywords internal
validateCreateEnrichMapArguments <- function(gostObject, query, source, 
        termIDs, removeRoot, title, showCategory, groupCategory, 
        cexLabelCategory, cexCategory) {
    
    ## Test that gostObject is a gprofiler2 result 
    if (!("list" %in% class(gostObject) && "result" %in% names(gostObject) &&
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
    
    if (!is(title, "character")) {
        stop("The \'title\' parameter must a character string.")
    }
    
    if (!is(showCategory, "character") && 
            !(is(showCategory, "numeric") && (showCategory > 0))) {
        stop("The \'showCategory\' parameter must an positive integer or a ", 
                "vector of character strings representing terms.")
    }
    
    if (!is(groupCategory , "logical")) {
        stop("The \'groupCategory\' parameter must a logical ", 
                "(TRUE or FALSE).")
    }
    
    if (!is(cexLabelCategory, "numeric") || !(cexLabelCategory > 0)) {
        stop("The \'cexLabelCategory\' parameter must be a positive numeric.")
    }
    
    if (!is(cexCategory, "numeric") || !(cexCategory > 0)) {
        stop("The \'cexCategory\' parameter must be a positive numeric.")
    }
    
    return(TRUE)   
}

#' @title Create a basic enrichment map
#' 
#' @description TODO
#' 
#' @param gostResults a \code{data.frame} containing the enrichment 
#' results to be plot.
#' 
#' @param backgroundGenes a \code{vector} of \code{character} string 
#' representing the name of the genes present in the request.
#' 
#' @param title a \code{character} string representing TODO
#' 
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{vector} of terms, 
#' the selected terms will be displayed.
#' 
#' @param groupCategory a \code{logical} indicating if the categories should 
#' be grouped.
#' 
#' @param cexLabelCategory a positive \code{numeric} representing the amount by 
#' which plotting category nodes label size should be scaled relative 
#' to the default (1).
#' 
#' @param cexCategory a positive \code{numeric} representing the amount by 
#' which plotting category nodes should be scaled relative to the default (1). 
#' 
#' @param significantMethod a \code{character} string representing the name 
#' of the multiple testing correction method used on the results.
#' 
#' @return TODO
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
#' gostResults <- gostResults[which(gostResults$term_id != "WILIPATHWAYS"),]
#' 
#' ## Extract meta data information
#' meta <- parentalNapaVsDMSOEnrichment$meta
#' 
#' backgroundGenes <- meta$query_metadata$queries[["parental_napa_vs_DMSO"]]
#' 
#' significantMethod <- meta$query_metadata$significance_threshold_method
#'
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @importFrom stringr str_ends
#' @importFrom enrichplot pairwise_termsim emapplot
#' @keywords internal
createBasicEmap <- function(gostResults, backgroundGenes, title, 
                                showCategory, groupCategory, 
                                cexLabelCategory, cexCategory, 
                                significantMethod) {
    # 
    # mapInfo <- data.frame(Cluster=c(rep(title, nrow(gostResults))), 
    #                 ID=c(gostResults$term_id), 
    #                 Description=c(gostResults$term_name), 
    #                 GeneRatio=c(paste0(gostResults$intersection_size, "/", 
    #                                         gostResults$query_size)),
    #                 BgRatio=c(paste0(gostResults$intersection_size, "/", 
    #                                         gostResults$query_size)),
    #                 pvalue=c(gostResults$p_value),
    #                 p.adjust=c(gostResults$p_value),
    #                 qvalue=c(gostResults$p_value),
    #                 geneID=c(gostResults$intersection),
    #                 Count=c(gostResults$intersection_size))
    # 
    
   #  
   #  mapInfo$Cluster <- factor(mapInfo$Cluster, levels=c(title))
   # ## mapInfo$ID <- stringr::str_replace(string = mapInfo$ID , pattern = "KEGG:", "hsa") 
   #  mapInfo$geneID <- stringr::str_replace_all(string = mapInfo$geneID, 
   #                                                  pattern = ",", "/") 
   #  
   #  
   #  geneClusters <- list()
   #  geneClusters[[title]] <- mapInfo$EnsemblID
    # 
    # setClass("compareClusterResult",
    #          representation=representation(
    #              compareClusterResult="data.frame",
    #              geneClusters="list",
    #              fun="character",
    #              gene2Symbol="character",
    #              keytype        = "character",
    #              readable       = "logical",
    #              .call          = "call",
    #              termsim        = "matrix",
    #              method         = "character",
    #              dr             = "list"
    #          )
    # )
    
    # res <- new("compareClusterResult",
    #            compareClusterResult=mapInfo,
    #            geneClusters=geneClusters,
    #            fun="enrichDEMO",
    #            .call = call("compareCluster(geneClusters = x, 
    #                             fun = \"enrichDEMO\", organism = \"hsa\", 
    #                             pvalueCutoff = 0.05)")
    # )
    
    geneSets <- list()
    
    
    for (i in seq_len(nrow(gostResults))) {
        geneSets[[gostResults$term_id[i]]] <- 
            unlist(stringr::str_split(gostResults$intersection[i], 
                                        pattern=","))
    }
    
    resultDataFrame <- data.frame(ID=gostResults$term_id, 
        Description=gostResults$term_name,
        GeneRatio=c(paste0(gostResults$intersection_size, "/", 
                                gostResults$query_size)), 
        BgRatio=c(paste0(gostResults$intersection_size, "/", 
                                gostResults$effective_domain_size)), 
        pvalues=gostResults$p_value,
        p.adjust=gostResults$p_value, 
        qvalue=gostResults$p_value, 
        geneID=stringr::str_replace_all(gostResults$intersection, ",", "/"),
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

    rownames(resultDataFrame) <- resultDataFrame$ID
    
    res <- new("enrichResult", result=resultDataFrame, pvalueCutoff=1, 
             pAdjustMethod=significantMethod, qvalueCutoff=1, 
             gene=as.character(backgroundGenes), 
             universe=as.character(backgroundGenes), 
             geneSets=geneSets, 
             organism="UNKNOWN", keytype="UNKNOWN", ontology="UNKNOWN", 
             readable=FALSE)
    
    ## Get similarity matrix
    comp <- pairwise_termsim(res)  
    
    graphEmap <- emapplot(x=comp, 
                        showCategory=showCategory,
                        cex_label_category=cexLabelCategory,  
                        group_category=groupCategory,
                        cex_category=cexCategory)
    
    return(graphEmap)
}