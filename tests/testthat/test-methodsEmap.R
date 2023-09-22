### Unit tests for methodsEmap.R functions

library(enrichViewNet)

data(demoGOST)



### Tests createEnrichMap() results

context("createEnrichMap() results")

test_that("createEnrichMap() must return error when gostObject is a number", {
    
    error_message <- paste0("The gostObject object should be a list ", 
                        "with meta and result as entries corresponding ", 
                        "to gprofiler2 enrichment output.")
    
    expect_error(createEnrichMap(gostObject=33, query="TEST", 
        source="GO:CC", termIDs=NULL, removeRoot=TRUE, title="network", 
        showCategory=30, groupCategory=FALSE, cexLabelCategory=1,
        cexCategory=1), 
        error_message)
})


test_that("createEnrichMap() must return error when gostObject is a string character", {
    
    error_message <- paste0("The gostObject object should be a list ", 
                        "with meta and result as entries corresponding ", 
                        "to gprofiler2 enrichment output.")
    
    expect_error(createEnrichMap(gostObject="TEST", query="TEST", 
        source="GO:CC", termIDs=NULL, removeRoot=TRUE, title="network", 
        showCategory=30, groupCategory=FALSE, cexLabelCategory=1,
        cexCategory=1), 
        error_message)
})


test_that("createEnrichMap() must return error when query is a number", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("The \'query\'must be a character string.")
    
    expect_error(createEnrichMap(gostObject=gostObject, query=33, 
        source="KEGG", termIDs=NULL, removeRoot=TRUE, title="network", 
        showCategory=30, groupCategory=FALSE, cexLabelCategory=1,
        cexCategory=1), error_message)
})


test_that("createEnrichMap() must return error when query is a vector of strings", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("The \'query\'must be a character string.")
    
    expect_error(createEnrichMap(gostObject=gostObject, query=c("1", "2"), 
        source="KEGG", termIDs=NULL, removeRoot=TRUE, title="network", 
        showCategory=30, groupCategory=FALSE, cexLabelCategory=1,
        cexCategory=1), error_message)
})


test_that("createEnrichMap() must return error when query is not in gost", {
  
    error_message <- paste0("The \'query\' is not present in the ", 
                                    "results of the gost object.")
    
    expect_error(createEnrichMap(gostObject=demoGOST, query="CANADA", 
        source="KEGG", termIDs=NULL, removeRoot=TRUE, title="network", 
        showCategory=30, groupCategory=FALSE, cexLabelCategory=1,
        cexCategory=1), error_message)
})


test_that("createEnrichMap() must return error when source is a number", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("Assertion on 'arg' failed: Must be of type ", 
                                "'character', not 'double'.")
    
    expect_error(createEnrichMap(gostObject=gostObject, query="toto", 
        source=333, termIDs=NULL, removeRoot=TRUE, title="network", 
        showCategory=30, groupCategory=FALSE, cexLabelCategory=1,
        cexCategory=1), 
        error_message)
})

test_that("createEnrichMap() must return error when source is a wrong name", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    expect_error(createEnrichMap(gostObject=gostObject, query="toto",  
        source="test", termIDs=NULL, removeRoot=TRUE, title="network", 
        showCategory=30, groupCategory=FALSE, cexLabelCategory=1,
        cexCategory=1))
})


test_that("createEnrichMap() must return error when source is GO", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    expect_error(createEnrichMap(gostObject=gostObject, query="toto", 
        source="GO",
        termIDs=NULL, removeRoot=TRUE, title="gprofiler network", 
        showCategory=30, groupCategory=FALSE, cexLabelCategory=1,
        cexCategory=1))
})


test_that("createEnrichMap() must return error when removeRoot remove last enriched term", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- paste0("With removal of the root term, there is no ", 
                                "enrichment term left")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, title="gprofiler network", 
                    showCategory=30, groupCategory=FALSE, cexLabelCategory=1,
                    cexCategory=1), 
                    error_message)
})


test_that("createEnrichMap() must return error when removeRoot remove last enriched term from term list", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- paste0("With removal of the root term, there is no ", 
                                "enrichment term left")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="TERM_ID",
        termIDs=c("WP:000000"), removeRoot=TRUE, showCategory=30, 
        groupCategory=FALSE, cexLabelCategory=1,
        cexCategory=1), error_message)
})


test_that("createEnrichMap() must return error when showCategory negative value", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'showCategory\' parameter must an positive ", 
            "integer or a vector of character strings representing terms.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=-30, 
                    groupCategory=FALSE, cexLabelCategory=1,
                    cexCategory=1), error_message)
})


test_that("createEnrichMap() must return error when showCategory is boolean", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'showCategory\' parameter must an positive ", 
            "integer or a vector of character strings representing terms.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=TRUE, 
                    groupCategory=FALSE, cexLabelCategory=1,
                    cexCategory=1), error_message)
})


test_that("createEnrichMap() must return error when groupCategory is integer", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'groupCategory\' parameter must a logical ", 
                                "(TRUE or FALSE).")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=30, 
                    groupCategory=22, cexLabelCategory=1,
                    cexCategory=1), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when cexLabelCategory is string", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'cexLabelCategory\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=30, 
                    groupCategory=FALSE, cexLabelCategory="test",
                    cexCategory=1), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when cexLabelCategory is negative", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'cexLabelCategory\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=30, 
                    groupCategory=FALSE, cexLabelCategory=-1.1,
                    cexCategory=1), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when cexCategory is negative", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'cexCategory\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=30, 
                    groupCategory=FALSE, cexLabelCategory=2,
                    cexCategory=-1), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when cexCategory is string", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'cexCategory\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=30, 
                    groupCategory=FALSE, cexLabelCategory=2,
                    cexCategory="te"), error_message, fixed=TRUE)
})