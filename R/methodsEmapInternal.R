
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
#' @param line a non-negative \code{numeric} representing the scale of line 
#' width. 
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
#'     categoryLabel=1.1, categoryNode=1, line=1.2)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @keywords internal
validateCreateEnrichMapArguments <- function(gostObject, query, source, 
        termIDs, removeRoot, showCategory, groupCategory, 
        categoryLabel, categoryNode, line) {
    
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
        categoryLabel=categoryLabel, categoryNode=categoryNode, line=line)
    
    return(result)     
}

#' @title Validate arguments passed to createEnrichMapAsIgraph() function
#' 
#' @description Validate the arguments passed to createEnrichMapAsIgraph() 
#' function.First, the object containing the enrichment results must 
#' correspond to a 
#' object created by \code{gprofiler2} software. Second, the selected 
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
#' @param similarityCutOff a positive \code{numeric} between 0 and 1 indicating 
#' the minimum level of similarity between two terms to have an edge linking 
#' the terms.
#' 
#' @return \code{TRUE} when all arguments are valid
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(demoGOST)
#' 
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateEnrichMapAsIgraphArg(
#'     gostObject=demoGOST, query="query_1", source="GO:BP", termIDs=NULL, 
#'     removeRoot=FALSE, showCategory=20, similarityCutOff=0.5)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @keywords internal
validateCreateEnrichMapAsIgraphArg <- function(gostObject, query, source, 
        termIDs, removeRoot, showCategory, similarityCutOff) {
    
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
    
    if (!is.logical(removeRoot)) {
        stop("The \'removeRoot\' parameter must a logical (TRUE or FALSE).")
    }
    
    if (!is.null(showCategory) && 
        !(is.numeric(showCategory) && (showCategory > 0))) {
        stop("The \'showCategory\' parameter must an positive integer or", 
                " NULL.")
    }
    
    if (!is.numeric(similarityCutOff) || (similarityCutOff <= 0.0) || 
        (similarityCutOff >= 1.0)) {
        stop("The \'similarityCutOff\' parameter must be a numeric superior",  
             " to zero and inferior to one.")
    }
    
    return(TRUE)     
}


#' @title Validate arguments passed to createEnrichMapMultiBasic() function
#' 
#' @description Validate the arguments passed to createEnrichMapMultiBasic() 
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
#' @param line a non-negative \code{numeric} representing the scale of line 
#' width. 
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
#'     categoryLabel=1.1, categoryNode=1, line=1.2)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @keywords internal
validateCreateEnrichMapMultiArguments <- function(gostObjectList, queryList, 
    source, termIDs, removeRoot, showCategory, groupCategory, 
    categoryLabel, categoryNode, line) {
    
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
                "2 query names are required. The number of query names should", 
                " correspond to the number of enrichment objects.")
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
            stop("Each query name present in the \'queryList'\ parameter ", 
                "must be present in the associated enrichment object.")
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
        categoryLabel=categoryLabel, categoryNode=categoryNode, line=line)
    
    return(result)   
}


#' @title Validate arguments passed to createEnrichMapMultiComplex() 
#' function
#' 
#' @description Validate the arguments passed to 
#' createEnrichMapMultiComplex() function.
#' First, the object containing the enrichment results must correspond to a 
#' object created by  \code{gprofiler2} software. Second, the selected 
#' source must at least have one enriched term in the results. Then, if the
#' source is 'TERM_ID', the listed terms must be present in the enrichment
#' results.
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
#' @param line a non-negative \code{numeric} representing the scale of line 
#' width. 
#' 
#' @return \code{TRUE} when all arguments are valid
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' data(rosaNapaVsDMSOEnrichment)
#' 
#' queryDataFrame <- data.frame(queryName=c("parental_napa_vs_DMSO", 
#'     "rosa_napa_vs_DMSO"), source=c("KEGG", "WP"), removeRoot=c(TRUE, TRUE),
#'     termIDs=c("", ""), groupName=c("parental - KEGG", "rosa - WP"), 
#'     stringsAsFactors=FALSE)
#' 
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateEnrichMapMultiComplexArg(
#'     gostObjectList=list(parentalNapaVsDMSOEnrichment, 
#'                             rosaNapaVsDMSOEnrichment),
#'     queryInfo=queryDataFrame,
#'     showCategory=20, groupCategory=FALSE, 
#'     categoryLabel=1.1, categoryNode=1, line=1.2)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
validateCreateEnrichMapMultiComplexArg <- function(gostObjectList, queryInfo, 
    showCategory, groupCategory, categoryLabel, categoryNode, line) {
    
    validateCreateEnrichMapMultiComplexGostPartOne(
        gostObjectList=gostObjectList, queryInfo=queryInfo)
    
    validateCreateEnrichMapMultiComplexGostPartTwo(
        gostObjectList=gostObjectList, queryInfo=queryInfo)
    
    result <- validateCreateEnrichMapSubSectionArguments(
        showCategory=showCategory, groupCategory=groupCategory, 
        categoryLabel=categoryLabel, categoryNode=categoryNode, line=line)
    
    return(result)   
}

#' @title Validate arguments passed to validateCreateEnrichMapMultiComplex() 
#' function
#' 
#' @description Validate the arguments passed to 
#' validateCreateEnrichMapMultiComplex() function.
#' First, the object containing the enrichment results must correspond to a 
#' object created by  \code{gprofiler2} software. Second, the selected 
#' source must at least have one enriched term in the results. Then, if the
#' source is 'TERM_ID', the listed terms must be present in the enrichment
#' results.
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
#' displayed. 
#' 
#' @param similarityCutOff a positive \code{numeric}, larger than zero and
#' small than 1 that represent the minimum similarity level between two 
#' nodes (terms) to be linked by an edge.
#' 
#' @return \code{TRUE} when all arguments are valid
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' data(rosaNapaVsDMSOEnrichment)
#' 
#' queryDataFrame <- data.frame(queryName=c("parental_napa_vs_DMSO", 
#'     "rosa_napa_vs_DMSO"), source=c("KEGG", "WP"), removeRoot=c(TRUE, TRUE),
#'     termIDs=c("", ""), groupName=c("parental - KEGG", "rosa - WP"), 
#'     stringsAsFactors=FALSE)
#' 
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateEnrichMapMultiComplexAsIgraphArg(
#'     gostObjectList=list(parentalNapaVsDMSOEnrichment, 
#'                             rosaNapaVsDMSOEnrichment),
#'     queryInfo=queryDataFrame, showCategory=20, similarityCutOff=0.2)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
validateCreateEnrichMapMultiComplexAsIgraphArg <- function(gostObjectList, 
        queryInfo, showCategory, similarityCutOff=similarityCutOff) {
    
    validateCreateEnrichMapMultiComplexGostPartOne(
        gostObjectList=gostObjectList, queryInfo=queryInfo)
    
    validateCreateEnrichMapMultiComplexGostPartTwo(
        gostObjectList=gostObjectList, queryInfo=queryInfo)
    
    if (!is.null(showCategory) && 
        !(is.numeric(showCategory) && (showCategory > 0))) {
        stop("The \'showCategory\' parameter must an positive integer or ", 
                "NULL.")
    }
    
    if (!is.numeric(similarityCutOff) || (similarityCutOff <= 0.0) || 
            (similarityCutOff >= 1.0)) {
        stop("The \'similarityCutOff\' parameter must be a numeric superior",  
                " to zero and inferior to one.")
    }
    
    return(TRUE)   
}


#' @title Validate arguments related to GOST passed to 
#' validateCreateEnrichMapMultiComplexArg() function
#' 
#' @description Validate the arguments passed to 
#' validateCreateEnrichMapMultiComplexArg() function.
#' First, the object containing the enrichment results must correspond to a 
#' object created by  \code{gprofiler2} software. Second, the selected 
#' source must at least have one enriched term in the results. Then, if the
#' source is 'TERM_ID', the listed terms must be present in the enrichment
#' results.
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
#' @return \code{TRUE} when all arguments are valid
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' data(rosaNapaVsDMSOEnrichment)
#' 
#' queryDataFrame <- data.frame(queryName=c("parental_napa_vs_DMSO", 
#'     "rosa_napa_vs_DMSO"), source=c("KEGG", "WP"), removeRoot=c(TRUE, TRUE),
#'     termIDs=c("", ""), groupName=c("parental - KEGG", "rosa - WP"), 
#'     stringsAsFactors=FALSE)
#' 
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateEnrichMapMultiComplexGostPartOne(
#'     gostObjectList=list(parentalNapaVsDMSOEnrichment, 
#'                             rosaNapaVsDMSOEnrichment),
#'     queryInfo=queryDataFrame)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @keywords internal
validateCreateEnrichMapMultiComplexGostPartOne <- function(gostObjectList, 
    queryInfo) {
    
    
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
    
    ## Test that queryInfo is a data.frame
    if (!inherits(queryInfo, "data.frame") || !(all(c("queryName", "source", 
            "removeRoot", "termIDs") %in% colnames(queryInfo)))) {
        stop("The queryInfo should a data.frame with those columns: ", 
                "queryName, source, removeRoot and termIDs.")
    }
    
    ## Test that queryInfo has the same number of entries than gostObjectList
    if (nrow(queryInfo) != length(gostObjectList)) {
        stop("The number of rows in queryInfo should ", 
                "correspond to the number of enrichment objects.")
    }
    
    ## Test that source should be character string
    if (!is.character(queryInfo$source)) {
        stop("The \'source'\ column of the \'queryInfo\' data frame ", 
                "should be in a character string format.")
    }
    
    return(TRUE)   
}


#' @title Validate arguments related to GOST passed to 
#' validateCreateEnrichMapMultiComplexArg() function
#' 
#' @description Validate the arguments passed to 
#' validateCreateEnrichMapMultiComplexArg() function.
#' First, the object containing the enrichment results must correspond to a 
#' object created by  \code{gprofiler2} software. Second, the selected 
#' source must at least have one enriched term in the results. Then, if the
#' source is 'TERM_ID', the listed terms must be present in the enrichment
#' results.
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
#' @return \code{TRUE} when all arguments are valid
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' data(rosaNapaVsDMSOEnrichment)
#' 
#' queryDataFrame <- data.frame(queryName=c("parental_napa_vs_DMSO", 
#'     "rosa_napa_vs_DMSO"), source=c("KEGG", "WP"), removeRoot=c(TRUE, TRUE),
#'     termIDs=c("", ""), groupName=c("parental - KEGG", "rosa - WP"), 
#'     stringsAsFactors=FALSE)
#' 
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateEnrichMapMultiComplexGostPartTwo(
#'     gostObjectList=list(parentalNapaVsDMSOEnrichment, 
#'                             rosaNapaVsDMSOEnrichment),
#'     queryInfo=queryDataFrame)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @keywords internal
validateCreateEnrichMapMultiComplexGostPartTwo <- function(gostObjectList, 
                queryInfo) {
    
    ## Test that the source values are valid choices in queryInfo data frame
    sources <- c("TERM_ID", "GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "TF", 
                    "MIRNA", "HPA", "CORUM", "HP", "WP")
    if (!all(unlist(lapply(seq_len(nrow(queryInfo)), FUN=function(x, queryI, 
                sourcesL){ queryI$source[x] %in% sourcesL}, 
                            queryI=queryInfo, sourcesL=sources)))) {
        stop("The values in the \'source\' column of the \'queryInfo\' ", 
                "data frame should be one of those: \"TERM_ID\", \"GO:MF\", ", 
                "\"GO:CC\", \"GO:BP\", \"KEGG\", \"REAC\", \"TF\", ",
                "\"MIRNA\", \"HPA\", \"CORUM\", \"HP\", \"WP\".")
    }
    
    ## Test that queryName should be character string
    if (!is.character(queryInfo$queryName)) {
        stop("The \'queryName'\ column of the \'queryInfo\' data frame ", 
             "should be in a character string format.")
    }
    
    ## Test that source should be character string
    if (!is.character(queryInfo$termIDs)) {
        stop("The \'termIDs'\ column of the \'queryInfo\' data frame ", 
                "should be in a character string format.")
    }
    
    ## Test that removeRoot should be logical
    if (!is.logical(queryInfo$removeRoot)) {
        stop("The \'removeRoot'\ column of the \'queryInfo\' data frame ", 
                "should only contain logical values (TRUE or FALSE).")
    }
    
    ## Test that the query names are present in the associated enrichment
    if (!all(unlist(lapply(seq_len(nrow(queryInfo)), FUN = function(x, queryI, 
                gostL) {res <- as.data.frame(gostL[[x]]$result); 
                            queryI$queryName[x] %in% unique(res$query)}, 
                            queryI=queryInfo, gostL=gostObjectList)))) {
        stop("Each query name present in the \'queryName'\ column of the ", 
                "\'queryInfo\' data frame must be present in the ", 
                "associated enrichment object.")
    }
    
    if (length(which(queryInfo$source == "TERM_ID"))) {
        if(!all(unlist(lapply(which(queryInfo$source == "TERM_ID"), 
                FUN=function(i, queryI) {queryI$termIDs[i] != ""}, 
                                queryI=queryInfo)))) {
            stop("A string of terms should be present in the ",
                     "\'termIDs\' column when source is \'TERM_ID\'.")
        }
    }
    
    if (!is.character(queryInfo$groupName)) {
        stop("The \'groupName'\ column of the \'queryInfo\' data frame ", 
                "should be in a character string format.")
    }
    
    if (any(table(queryInfo$groupName) > 1)) {
        stop("The \'groupName'\ column of the \'queryInfo\' data frame ", 
                "should only contain unique group names.")
    }
    
    return(TRUE)   
}


#' @title Validate common arguments passed to createEnrichMap() and 
#' createEnrichMapMultiBasic() functions
#' 
#' @description Validate the common arguments passed to createEnrichMap() and 
#' createEnrichMapMultiBasic() functions.
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
#' @param line a non-negative \code{numeric} representing the scale of line 
#' width. 
#' 
#' @return \code{TRUE} when all arguments are valid
#' 
#' @examples
#'
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateEnrichMapSubSectionArguments(
#'     showCategory=20, groupCategory=FALSE, categoryLabel=1.1, categoryNode=1, 
#'     line=0.5)
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @keywords internal
validateCreateEnrichMapSubSectionArguments <- function(showCategory, 
    groupCategory, categoryLabel, categoryNode, line) {
    
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

    if (!is.numeric(line) || !(line > 0)) {
        stop("The \'line\' parameter must be a positive numeric.")
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
#' @param line a non-negative \code{numeric} representing the scale of line 
#' width. 
#' 
#' @param ... additional arguments that will be pass to the 
#' \code{\link[enrichplot]{emapplot}} function. 
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
#' gostResults <- gostResults[which(gostResults$term_name != "WIKIPATHWAYS"),]
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
#'     significantMethod=significantMethod, line=1)
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
        significantMethod, line, ...) {
    
    ## Extract gene list for each term
    geneSets <- lapply(seq_len(nrow(gostResults)), FUN=function(x, gostData) {
        str_split(gostData$intersection[x], pattern=",")[[1]]}, 
        gostData=gostResults
    )
    names(geneSets) <- gostResults$term_id
    
    resultDF <- data.frame(ID=gostResults$term_id, 
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
    
    resultDF <- resultDF[order(resultDF$pvalues), ]
    
    ## Problem when identical description 
    if (any(table(resultDF$Description) > 1)) {
        resultDF <- manageNameDuplicationInEmap(clProfDF=resultDF)
    }
    
    rownames(resultDF) <- resultDF$ID
    
    res <- new("enrichResult", result=resultDF, pvalueCutoff=1, 
                pAdjustMethod="UNKNOWN", qvalueCutoff=1, 
                gene=as.character(backgroundGenes), 
                universe=as.character(backgroundGenes), 
                geneSets=geneSets, 
                organism="UNKNOWN", keytype="UNKNOWN", ontology="UNKNOWN", 
                readable=FALSE)
    
    ## Get similarity matrix
    comp <- pairwise_termsim(res)  
    
    graphEmap <- emapplot(x=comp, showCategory=showCategory,
        group=groupCategory, size_category=categoryNode, size_edge=line,  ...)
    
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
#' @param similarityCutOff a positive \code{numeric} between 0 and 1 indicating 
#' the minimum level of similarity between two terms to have an edge linking 
#' the terms.  
#' 
#' @return a \code{igraph} object representing the enrichment map.
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
#' gostResults <- gostResults[which(gostResults$term_name != "WIKIPATHWAYS"),]
#' 
#' ## Extract meta data information
#' meta <- parentalNapaVsDMSOEnrichment$meta
#' 
#' ## Get all background genes
#' backgroundGenes <- meta$query_metadata$queries[["parental_napa_vs_DMSO"]]
#'
#' ## Create basic enrichment map, as an igraph, using Wikipathways terms
#' enrichViewNet:::createBasicEmapAsIgraph(gostResults=gostResults, 
#'     backgroundGenes=backgroundGenes, showCategory=30L, 
#'     similarityCutOff=0.2)
#'     
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is new
#' @importFrom stringr str_ends str_split str_replace_all
#' @importFrom igraph make_empty_graph V add_vertices E add_edges E<- delete_edges
#' @importFrom reshape2 melt
#' @importClassesFrom DOSE enrichResult
#' @keywords internal
createBasicEmapAsIgraph <- function(gostResults, backgroundGenes, 
        showCategory, similarityCutOff) {
    
    ## Extract gene list for each term
    geneSets <- lapply(seq_len(nrow(gostResults)), FUN=function(x, gostData) {
                    str_split(gostData$intersection[x], pattern=",")[[1]]}, 
                    gostData=gostResults)
    names(geneSets) <- gostResults$term_id
    
    resultDF <- data.frame(ID=gostResults$term_id, 
                Description=gostResults$term_name,
                GeneRatio=c(paste0(gostResults$intersection_size, "/", 
                                            gostResults$query_size)), 
                BgRatio=c(paste0(gostResults$intersection_size, "/", 
                                            gostResults$effective_domain_size)), 
                pvalues=gostResults$p_value, 
                p.adjust=gostResults$p_value, 
                geneID=str_replace_all(gostResults$intersection, ",", "/"),
                Count=c(gostResults$intersection_size), stringsAsFactors=FALSE)
    
    resultDF <- resultDF[order(resultDF$pvalues, decreasing=FALSE), ]
    
    ## Select top category to show
    if (is.numeric(showCategory)) {
        if (nrow(resultDF) > showCategory) {
            resultDF <- resultDF[seq_len(showCategory), ]
        }
    }
    
    ## Problem when identical description 
    if (any(table(resultDF$Description) > 1)) {
        resultDF <- manageNameDuplicationInEmap(clProfDF=resultDF)
    }
    
    rownames(resultDF) <- resultDF$ID
    
    if (nrow(resultDF) == 1) {
        ## Create igraph with 1 entry
        g <- make_empty_graph(0, directed=FALSE)
        g <- add_vertices(g, nv=1)
        V(g)$name <- as.character(resultDF$Description)
    } else {
        ## Create igraph with multiple entries
        
        ## Get similarity matrix 
        res <- similarityJaccard(resultDF=resultDF)
        
        ## Generated data frame with similarity information
        simData <- melt(res[as.character(resultDF$Description), 
                    as.character(resultDF$Description)], na.rm=TRUE, 
                    value.name="similarity")
        simData <- simData[which(simData[,1] != simData[,2]), ]
        simData <- simData[which(simData$similarity > 0), ]
        
        ## Each unique name will be a node
        g <- make_empty_graph(n=0, directed=FALSE)
        attrs <- list(name=resultDF$Description, 
                      size=resultDF$Count)
        g <- add_vertices(g, length(attrs$name), attr=attrs)
        
        ## Create links between nodes
        fromN <- as.character(simData[, 1])
        toN <- as.character(simData[, 2])
        edges <- rbind(match(fromN, attrs$name), match(toN, attrs$name))
        
        ## Add attributes (including similarity values) to the edges
        attrsE <- list()
        attrsE[["similarity"]] <- simData$similarity
        attrsE[["width"]] <- simData$similarity
        g <- add_edges(g, edges, attr=attrsE)
        
        ## Remove edges with similarity lower than cut off
        g <- delete_edges(g, E(g)[simData$similarity < similarityCutOff])
    }
    return(g)
}

#' @title Create a basic enrichment map
#' 
#' @description The function creates a basic enrichment map using functional 
#' enrichment results.
#' 
#' @param gostResultsList a \code{list} of \code{data.frame} containing 
#' the enrichment results to be plot with different group identification.
#' 
#' @param queryList a \code{list} of \code{character} string 
#' representing the name of query retained for each enrichment results present 
#' in the \code{gostResultsList} parameter. The query should be present in its 
#' associated enrichment results.
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
#' @param line a non-negative \code{numeric} representing the scale of line 
#' width.
#' 
#' @param ... additional arguments that will be pass to the 
#' \code{\link[enrichplot]{emapplot}} function. 
#'
#' @return a \code{ggplot} object representing the enrichment map with 
#' different colors for each group of enrichment results.
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' 
#' ## Only retain the results section
#' gostResults <- as.data.frame(parentalNapaVsDMSOEnrichment$result)
#' 
#' ## Limit the results subsection of REACTOME and KEGG
#' gostResultsREAC <- gostResults[which(gostResults$source == "REAC"),]
#' gostResultsREAC <- gostResultsREAC[1:13, ]
#' gostResultsKEGG <- gostResults[which(gostResults$source == "KEGG"),]
#' 
#' ## Extract meta data information
#' queryList <- list("parental - REACTOME", "parental - KEGG")
#' 
#' ## Create basic enrichment map using Wikipathways terms
#' enrichViewNet:::createMultiEmap(gostResultsList=list(gostResultsREAC, 
#'     gostResultsKEGG), queryList=queryList, showCategory=30L, 
#'     groupCategory=FALSE, categoryLabel=1, categoryNode=1, line=1.4)
#'     
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is new
#' @importFrom stringr str_split str_replace_all
#' @importFrom enrichplot pairwise_termsim emapplot
#' @importClassesFrom DOSE compareClusterResult
#' @keywords internal
createMultiEmap <- function(gostResultsList, queryList, showCategory, 
    groupCategory, categoryLabel, categoryNode, line,  ...) {
    
    resF <- list()
    geneClusters <- list()
    
    ## Rename query names when duplicated names
    queryList <- manageQueryDuplicationInEmap(queryList=queryList)
    
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
            unique(unlist(str_split(gostResults$intersection, ",")))
    }
    clProfDF <- do.call(rbind, resF)
    
    clProfDF$Cluster  <- factor(clProfDF$Cluster)
    clProfDF$geneID <- str_replace_all(string = clProfDF$geneID,
                                                pattern = ",", "/") 
    
    ## Problem when identical description 
    if (any(table(clProfDF$Description) > 1)) {
        clProfDF <- manageNameDuplicationInEmap(clProfDF=clProfDF)
    }
    
    res <- new("compareClusterResult",
                compareClusterResult = clProfDF,
                geneClusters = geneClusters,
                fun = "createMultiEmap",
                .call = call("createMultiEmap()")
    )
    res@keytype <- "UNKNOWN"
    res@readable <- FALSE
    
    compar <- pairwise_termsim(res)  
    
    graphEmap <- emapplot(compar, showCategory=showCategory, 
        group=groupCategory, size_category=categoryNode, size_edge=line, ...)
    
    return(graphEmap)
}


#' @title Create a complex enrichment map as an igraph
#' 
#' @description The function creates a complex enrichment map, as an igraph, 
#' using functional enrichment results.
#' 
#' @param gostResultsList a \code{list} of \code{data.frame} containing 
#' the enrichment results to be plot with different group identification.
#' 
#' @param queryList a \code{list} of \code{character} string 
#' representing the name of query retained for each enrichment results present 
#' in the \code{gostResultsList} parameter. The query should be present in its 
#' associated enrichment results.
#' 
#' @param showCategory a positive \code{integer} or a \code{vector} of 
#' \code{characters} representing terms.  If a \code{integer}, the first 
#' \code{n} terms will be displayed. If \code{NULL}, 
#' all terms will be displayed.
#' 
#' @param similarityCutOff a positive \code{numeric}, larger than zero and
#' small than 1 that represent the minimum similarity level between two 
#' nodes (terms) to be linked by an edge. 
#' 
#' @return a \code{igraph} object representing the enrichment map with 
#' different colors for each group of enrichment results.
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' 
#' ## Only retain the results section
#' gostResults <- as.data.frame(parentalNapaVsDMSOEnrichment$result)
#' 
#' ## Limit the results subsection of REACTOME and KEGG
#' ## Kegg is replicated to show shared results between queries
#' gostResultsREAC <- gostResults[which(gostResults$source == "REAC"),]
#' gostResultsREAC <- gostResultsREAC[1:13, ]
#' gostResultsKEGG <- gostResults[which(gostResults$source == "KEGG"),]
#' gostResultsKEGG2 <- gostResultsKEGG[1:6,]
#' 
#' ## Extract meta data information
#' queryList <- list("parental - REACTOME", "parental - KEGG - v1", 
#'     "parental - KEGG - v2")
#' 
#' ## Create basic enrichment map using Wikipathways terms
#' igraph <- enrichViewNet:::createMultiEmapAsIgraph(
#'     gostResultsList=list(gostResultsREAC, gostResultsKEGG, 
#'             gostResultsKEGG2), 
#'     queryList=queryList, showCategory=30L, similarityCutOff=0.5)
#'     
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom stringr str_replace_all
#' @importFrom igraph delete_edges V E add_edges V<-
#' @importFrom reshape2 melt
#' @keywords internal
createMultiEmapAsIgraph <- function(gostResultsList, queryList, showCategory, 
                similarityCutOff) {
    
    result <- createCompareResultDataFrame(queryList=queryList, 
                                    gostResultsList=gostResultsList)
    
    ## Limit the number of categories if specified by user
    ## Keep best results with best p-values
    if (!is.null(showCategory)) {
        tmp <- list()
        for (i in unique(as.character(result$Cluster))) {
            resTmp <- result[which(result$Cluster == i), ]
            if (nrow(resTmp) > showCategory) {
                resTmp <- resTmp[order(resTmp$pvalues, decreasing=FALSE), ]
                resTmp <- resTmp[seq_len(showCategory), ]
            }
            tmp[[length(tmp) + 1]] <- resTmp
        }
        result <- do.call(rbind, tmp)
    }
    ## Ensure no result with count to zero
    result <- result[which(as.numeric(result$Count) > 0), ]

    ## Transform description into factors
    result$Description <- as.character(result$Description)

    ## Manage when there are terms present in more than one cluster
    ## Count should include all unique genes 
    resultDF <- result[!duplicated(result$ID), ]
    for (i in seq_len(nrow(resultDF))) {
        tmpRes <- result[which(result$ID == resultDF$ID[i]), ]
        if (nrow(tmpRes) > 1) {
            uniqueGenes <- unique(unlist(strsplit(tmpRes$geneID, "/")))
            resultDF$geneID[i] <- paste0(uniqueGenes, collapse = "/")
            resultDF$Count[i] <- length(uniqueGenes)
        }
    }
    resultDF$Cluster <- NULL
        
    if (nrow(resultDF) == 1) {
        ## Create igraph with 1 entry
        g <- make_empty_graph(0, directed=FALSE)
        g <- add_vertices(g, nv=1)
        V(g)$name <- as.character(resultDF$Description)
    } else {
        ## Create igraph with multiple entries
        
        ## Each unique name will be a node
        g <- make_empty_graph(n=0, directed=FALSE)
        vertex.color <- list()
        vertex.pie <- list()
        vertex.pieName <- list()
        vertex.cluster <- list()
        for (i in seq_len(nrow(resultDF))) {
            e <- result[result$ID== resultDF$ID[i], ]
            vertex.pie[[length(vertex.pie) + 1]] <- 
                    as.integer(unlist(queryList) %in% e$Cluster)
            vertex.pieName[[length(vertex.pieName) + 1]] <-  unlist(queryList)
            vertex.cluster[[length(vertex.cluster) + 1]] <- 
                    as.character(e$Cluster)
        }
        
        attrs <- list(name=resultDF$Description, size=resultDF$Count, 
                        pie=vertex.pie, cluster=vertex.cluster, 
                        pieName=vertex.pieName)
        g <- add_vertices(g, length(attrs$name), attr=attrs)
        
        ## Get similarity matrix 
        res <- similarityJaccard(resultDF=resultDF)
        
        ## Generated data frame with similarity information
        ## Similarity of zero should be removed
        simData <- melt(res[as.character(resultDF$Description), 
                        as.character(resultDF$Description)], na.rm=TRUE, 
                        value.name="similarity")
        simData <- simData[which(simData[,1] != simData[,2]), ]
        simData <- simData[which(simData$similarity > 0), ]
            
        ## Create links between nodes
        fromNode <- as.character(simData[, 1])
        toNode <- as.character(simData[, 2])
        edges <- rbind(match(fromNode, attrs$name), match(toNode, attrs$name))
            
        ## Add attributes (including similarity values) to the edges
        attrsE <- list()
        attrsE[["similarity"]] <- simData$similarity
        attrsE[["width"]] <- simData$similarity
        g <- add_edges(g, edges, attr=attrsE)
            
        ## Remove edges with similarity lower than cut off
        g <- delete_edges(g, E(g)[simData$similarity < similarityCutOff])
    }
    
    return(g)   
}


#' @title Calculate the Jaccard similarity coefficient for the terms present
#' in a data frame
#' 
#' @description The function calculates the Jaccard similarity coefficient for 
#' the enriched terms present i a data frame.  
#' 
#' @param resultDF a \code{data.frame} containing the enrichment 
#' terms that are going to be graphed as an enrichment map. The 
#' \code{data.frame} must contain the \code{geneID} and \code{Description} 
#' columns.
#' 
#' @return a \code{matrix} with the lower triangle section containing the 
#' Jaccard similarity coefficients as \code{numeric}. The upper triangle 
#' section is filled with \code{NA} as well as the diagonal. 
#' 
#' @examples
#'
#' ## Data frame with the enriched terms
#' resultData <- data.frame( 
#'     ID=c("REAC:R-HSA-9031628", "REAC:R-HSA-9648895", "REAC:R-HSA-9614085",
#'             "REAC:R-HSA-166520"),
#'     Description=c("NGF-stimulated transcription", 
#'             "Response of EIF2AK1 (HRI) to heme deficiency",
#'             "FOXO-mediated transcription", "Signaling by NTRKs"),
#'     geneID=c(paste0("ENSG00000120738/ENSG00000122877/ENSG00000125740/",
#'         "ENSG00000135625/ENSG00000170345/ENSG00000171223/ENSG00000173334/", 
#'         "ENSG00000179388/ENSG0000019857"),
#'         paste0("ENSG00000087074/ENSG00000128272/ENSG00000128965/", 
#'         "ENSG00000162772/ENSG00000172216/ENSG00000175197"),
#'         paste0("ENSG00000105327/ENSG00000113916/ENSG00000124762/", 
#'         "ENSG00000133639/ENSG00000136826/ENSG00000153094/ENSG00000175197"),
#'         paste0("ENSG00000105327/ENSG00000113916/ENSG00000124762/", 
#'         "ENSG00000133639/ENSG00000136826/ENSG00000153094/ENSG0000017519/",
#'         "ENSG00000120738/ENSG00000122877/ENSG00000125740/ENSG00000135625/",
#'         "ENSG00000170345/ENSG00000171223/ENSG00000173334/ENSG00000179388/", 
#'         "ENSG00000198576")))
#' 
#' ## Calculate the Jaccard similarity coefficient for all pairwise 
#' ## term combination
#' enrichViewNet:::similarityJaccard(resultDF=resultData)
#'     
#' @author Astrid Deschênes
#' @importFrom stringr str_split
#' @encoding UTF-8
#' @keywords internal
similarityJaccard <- function(resultDF) {
    res <- matrix(data=NA, nrow=nrow(resultDF), ncol=nrow(resultDF))
    for(i in seq_len(nrow(resultDF))) {
        tempI <- unique(unlist(str_split(resultDF$geneID[i], "/")))
        for(j in seq_len(i)) {
            if (i == j) break
            tempJ <- unique(unlist(str_split(resultDF$geneID[j], "/")))
            all <- c(tempI, tempJ)
            if (length(all) < 1) {
                res[i, j] <- 0
            } else {
                res[i, j] <- sum(table(all) > 1)/length(unique(all))
            }
        }
    }
    colnames(res) <- as.character(resultDF$Description)
    rownames(res) <- as.character(resultDF$Description)
    return(res)
}


#' @title Change name description in data frame when more than one term as 
#' the same name
#' 
#' @description The function adds the ID of the term at the name of the name 
#' description when terms with different ID share the same name. For example, 
#' "MAPK pathway" would become "MAPK pathway (KEGG:04010)". 
#' 
#' @param clProfDF a \code{data.frame} containing the enrichment 
#' terms that are going to be graphed as an enrichment map. The 
#' \code{data.frame} must contain at least one case of duplicated name 
#' descriptions for different term IDs.
#' 
#' @return a \code{data.frame} with the name descriptions modified for terms 
#' with duplicated name descriptions but different term IDs.
#' 
#' @examples
#'
#' ## Data frame with duplicated name descriptions for different IDs
#' clustData <- data.frame(Cluster=c("group 1" , "group 1", "group 2"), 
#'     ID=c("WP:WP4925", "WP:WP382", "KEGG:04010"),
#'     Description=c("Unfolded protein response", 
#'                     rep("MAPK signaling pathway", 2)),
#'     GeneRatio=c("4/157", "3/157", "3/157"),
#'     BgRatio=c("4/24022", "3/24022", "3/24022"),
#'     pvalues=c(1.55e-4, 8.13e-8, 4.33e-5),
#'     p.adjust=c(1e-3, 1e-3, 1.4e-3), qvalue=c(1e-3, 1e-3, 1.4e-3), 
#'     geneID=c("ENSG000107968/ENSG000120129/ENSG000123358/ENSG000158050",
#'         "ENSG000107968/ENSG000120129/ENSG000158050",
#'         "ENSG000107968/ENSG000120129/ENSG000158050"),
#'     Count=c(4, 3, 3))
#' 
#' ## Change the name descriptions for the duplicated terms
#' enrichViewNet:::manageNameDuplicationInEmap(clProfDF=clustData)
#'     
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
manageNameDuplicationInEmap <- function(clProfDF) {
    
    ## Problem when identical description 
    test <- which(table(clProfDF$Description) > 1)
    for (i in names(test)) {
        resTest <- clProfDF[clProfDF$Description == i, ]
        if (length(unique(resTest$ID)) != 1) {
            clProfDF$Description[clProfDF$Description == i] <- 
                paste0(clProfDF$Description[clProfDF$Description == i], 
                        " (", clProfDF$ID[clProfDF$Description == i], ")")   
        }
    }
    
    return(clProfDF)
}


#' @title Create a compareClusterResult object from the enrichment results
#' 
#' @description The function creates a compareClusterResult object using the
#' list of queries and list of enrichment results given be the user.
#' 
#' @param gostResultsList a \code{list} of \code{data.frame} containing 
#' the enrichment results to be plot with different group identification.
#' 
#' @param queryList a \code{list} containing the query names.
#' 
#' @return a \code{compareClusterResult} object
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' 
#' ## Only retain the results section
#' gostResults <- as.data.frame(parentalNapaVsDMSOEnrichment$result)
#' 
#' ## Limit the results subsection of REACTOME and KEGG
#' ## Kegg is replicated to show shared results between queries
#' gostResultsREAC <- gostResults[which(gostResults$source == "REAC"),]
#' gostResultsREAC <- gostResultsREAC[1:13, ]
#' gostResultsKEGG <- gostResults[which(gostResults$source == "KEGG"),]
#' gostResultsKEGG2 <- gostResultsKEGG[1:6,]
#' 
#' ## List of query names with duplicated names
#' queryList <- list("parental_vs_DMSO", "rosa_vs_DMSO", "parental_vs_DMSO", 
#'     "rosa_vs_DMSO", "parental_vs_Control", "rosa_vs_DMSO")
#'     
#' ## Extract meta data information
#' queryList <- list("parental - REACTOME", "parental - KEGG - v1", 
#'     "parental - KEGG - v2")
#' 
#' ## Change the query names for the duplicated names
#' enrichViewNet:::createCompareClusterResultObject(queryList=queryList,
#'     gostResultsList=list(gostResultsREAC, gostResultsKEGG, 
#'             gostResultsKEGG2))
#'     
#' @author Astrid Deschênes
#' @importClassesFrom DOSE compareClusterResult
#' @importFrom stringr str_replace_all str_split
#' @encoding UTF-8
#' @keywords internal
createCompareClusterResultObject <- function(queryList, gostResultsList) {
    
    resF <- list()
    geneClusters <- list()
    
    for(i in seq_len(length(gostResultsList))) {
        gostResults <- gostResultsList[[i]]
        resultDF <- data.frame(
            Cluster=rep(queryList[[i]], nrow(gostResults)),
            ID=gostResults$term_id,
            Description=as.character(gostResults$term_name), 
            GeneRatio=as.numeric(gostResults$intersection_size) / 
                as.numeric(gostResults$query_size),
            BgRatio=c(paste0(gostResults$intersection_size, "/", 
                             gostResults$effective_domain_size)), 
            pvalues=as.numeric(gostResults$p_value),
            geneID=str_replace_all(gostResults$intersection, ",", "/"),
            Count=c(as.numeric(gostResults$intersection_size)), 
            stringsAsFactors=FALSE)
        resF[[i]] <- resultDF
        
        geneClusters[[queryList[[i]]]] <- 
            unique(unlist(str_split(gostResults$intersection, ",")))
    }
    clProfDF <- do.call(rbind, resF)
    
    clProfDF$Cluster  <- factor(clProfDF$Cluster)
    clProfDF$geneID <- str_replace_all(string=clProfDF$geneID,
                                       pattern=",", "/") 
    
    ## Problem when identical description (ID added to description)
    if (any(table(clProfDF$Description) > 1)) {
        clProfDF <- manageNameDuplicationInEmap(clProfDF=clProfDF)
    }
    
    res <- new("compareClusterResult",
               compareClusterResult = clProfDF,
               geneClusters = geneClusters,
               fun = "createMultiEmap",
               .call = call("createMultiEmap()")
    )
    res@keytype <- "UNKNOWN"
    res@readable <- FALSE
    
    return(res)
}


#' @title Create a data frame containing the compiled enrichment results from
#' multiple queries
#' 
#' @description The function creates a data frame object using the
#' list of queries and list of enrichment results given be the user.
#' 
#' @param gostResultsList a \code{list} of \code{data.frame} containing 
#' the enrichment results to be plot with different group identification.
#' 
#' @param queryList a \code{list} containing the query names.
#' 
#' @return a \code{data.frame} object
#' 
#' @examples
#'
#' ## Load the result of an enrichment analysis done with gprofiler2
#' data(parentalNapaVsDMSOEnrichment)
#' 
#' ## Only retain the results section
#' gostResults <- as.data.frame(parentalNapaVsDMSOEnrichment$result)
#' 
#' ## Limit the results subsection of REACTOME and KEGG
#' ## Kegg is replicated to show shared results between queries
#' gostResultsREAC <- gostResults[which(gostResults$source == "REAC"),]
#' gostResultsREAC <- gostResultsREAC[1:13, ]
#' gostResultsKEGG <- gostResults[which(gostResults$source == "KEGG"),]
#' gostResultsKEGG2 <- gostResultsKEGG[1:6,]
#' 
#' ## List of query names with duplicated names
#' queryList <- list("parental_vs_DMSO", "rosa_vs_DMSO", "parental_vs_DMSO", 
#'     "rosa_vs_DMSO", "parental_vs_Control", "rosa_vs_DMSO")
#'     
#' ## Extract meta data information
#' queryList <- list("parental - REACTOME", "parental - KEGG - v1", 
#'     "parental - KEGG - v2")
#' 
#' ## Change the query names for the duplicated names
#' enrichViewNet:::createCompareResultDataFrame(queryList=queryList,
#'     gostResultsList=list(gostResultsREAC, gostResultsKEGG, 
#'             gostResultsKEGG2))
#'     
#' @author Astrid Deschênes
#' @importFrom stringr str_replace_all str_split
#' @encoding UTF-8
#' @keywords internal
createCompareResultDataFrame <- function(queryList, gostResultsList) {
    
    resF <- list()
    geneClusters <- list()
    
    for(i in seq_len(length(gostResultsList))) {
        gostResults <- gostResultsList[[i]]
        resultDF <- data.frame(
            Cluster=rep(queryList[[i]], nrow(gostResults)),
            ID=gostResults$term_id,
            Description=as.character(gostResults$term_name), 
            pvalues=as.numeric(gostResults$p_value),
            geneID=str_replace_all(gostResults$intersection, ",", "/"),
            Count=c(as.numeric(gostResults$intersection_size)), 
            stringsAsFactors=FALSE)
        resF[[i]] <- resultDF
        
        geneClusters[[queryList[[i]]]] <- 
            unique(unlist(str_split(gostResults$intersection, ",")))
    }
    allDF <- do.call(rbind, resF)
    
    allDF$Cluster  <- factor(allDF$Cluster)
    allDF$geneID <- str_replace_all(string=allDF$geneID, pattern=",", "/") 
    
    ## Problem when identical description (ID added to description)
    if (any(table(allDF$Description) > 1)) {
        allDF <- manageNameDuplicationInEmap(clProfDF=allDF)
    }
    
    return(allDF)
}

#' @title Change query name when more than one queries have  
#' the same name
#' 
#' @description The function adds an unique ID at the name of the query when 
#' more than one query have the same name.
#' 
#' @param queryList a \code{list} containing the query names.
#' 
#' @return a \code{list} with the query names modified for queries   
#' with duplicated names.
#' 
#' @examples
#'
#' ## List of query names with duplicated names
#' queryList <- list("parental_vs_DMSO", "rosa_vs_DMSO", "parental_vs_DMSO", 
#'     "rosa_vs_DMSO", "parental_vs_Control", "rosa_vs_DMSO")
#' 
#' ## Change the query names for the duplicated names
#' enrichViewNet:::manageQueryDuplicationInEmap(queryList=queryList)
#'     
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
manageQueryDuplicationInEmap <- function(queryList) {
    
    ## Problem when identical query names
    test <- which(table(unlist(queryList)) > 1)
    if (length(test) > 0) {
        for (i in names(test)) {
            pos <- which(queryList == i)
            id <- 1
            for (j in pos) {
                queryList[[j]] <- paste0(queryList[[j]], " (", id, ")")
                id <- id + 1
            }
        }
    }
    
    return(queryList)
}