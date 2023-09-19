
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
#' 
#' ## Check that all arguments are valid
#' enrichViewNet:::validateCreateEnrichMapArguments(gostObject=demoGOST,
#'     source="GO:BP", termIDs=NULL, removeRoot=FALSE, title="GO:BP")
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @importFrom stringr str_ends
#' @keywords internal
validateCreateEnrichMapArguments <- function(gostObject, source, termIDs,
                                           removeRoot, title) {
    
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
    
    if (!is(title, "character")) {
        stop("The \'title\' parameter must a character string.")
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
#' @param title a \code{character} string representing TODO
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
#' 
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @importFrom methods is
#' @importFrom stringr str_ends
#' @importFrom enrichplot pairwise_termsim emapplot
#' @keywords internal
createBasicEmap <- function(gostResults, title) {
    
    ## TODO
    
    
    mapInfo <- data.frame(Cluster=c(rep(title, nrow(gostResults))), 
                    ID=c(gostResults$term_id), 
                    Description=c(gostResults$term_name), 
                    GeneRatio=c(paste0(gostResults$intersection_size, "/", 
                                            gostResults$query_size)),
                    BgRatio=c(paste0(gostResults$intersection_size, "/", 
                                            gostResults$query_size)),
                    pvalue=c(gostResults$p_value),
                    p.adjust=c(gostResults$p_value),
                    qvalue=c(gostResults$p_value),
                    geneID=c(gostResults$intersection),
                    Count=c(gostResults$intersection_size))
    
    
    
    mapInfo$Cluster = factor(mapInfo$Cluster, levels=c(title))
   # mapInfo$ID <- stringr::str_replace(string = mapInfo$ID , pattern = "KEGG:", "hsa") 
    mapInfo$geneID <- stringr::str_replace_all(string = mapInfo$geneID , 
                                                    pattern = ",", "/") 
    
    
    geneClusters <- list()
    geneClusters[[title]] <- mapInfo$EnsemblID
    
    setClass("compareClusterResult",
             representation = representation(
                 compareClusterResult = "data.frame",
                 geneClusters = "list",
                 fun = "character",
                 gene2Symbol    = "character",
                 keytype        = "character",
                 readable       = "logical",
                 .call          = "call",
                 termsim        = "matrix",
                 method         = "character",
                 dr             = "list"
             )
    )
    
    res <- new("compareClusterResult",
               compareClusterResult=mapInfo,
               geneClusters=geneClusters,
               fun="enrichDEMO",
               .call = call("compareCluster(geneClusters = x, 
                                fun = \"enrichDEMO\", organism = \"hsa\", 
                                pvalueCutoff = 0.05)")
    )
    
    res@keytype <- "UNKNOWN"
    res@readable <- FALSE
    
    ## Get similarity matrix
    comp <- pairwise_termsim(res)  
    
    
    graphEmap <- emapplot(x=comp, 
                        showCategory=length(table(mapInfo$Description)),
                        cex_label_category=1,  group_category=F,
                        cex_category=1)
    
    return(graphEmap)
}
