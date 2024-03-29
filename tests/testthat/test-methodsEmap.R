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
        source="GO:CC", termIDs=NULL, removeRoot=TRUE,
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, force=FALSE), error_message)
})


test_that("createEnrichMap() must return error when gostObject is a string character", {
    
    error_message <- paste0("The gostObject object should be a list ", 
                        "with meta and result as entries corresponding ", 
                        "to gprofiler2 enrichment output.")
    
    expect_error(createEnrichMap(gostObject="TEST", query="TEST", 
        source="GO:CC", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, force=FALSE), error_message)
})


test_that("createEnrichMap() must return error when query is a number", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("The \'query\'must be a character string.")
    
    expect_error(createEnrichMap(gostObject=gostObject, query=33, 
        source="KEGG", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, force=FALSE), error_message)
})


test_that("createEnrichMap() must return error when query is a vector of strings", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("The \'query\'must be a character string.")
    
    expect_error(createEnrichMap(gostObject=gostObject, query=c("1", "2"), 
        source="KEGG", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, force=FALSE), error_message)
})


test_that("createEnrichMap() must return error when query is not in gost", {
  
    error_message <- paste0("The \'query\' is not present in the ", 
                                    "results of the gost object.")
    
    expect_error(createEnrichMap(gostObject=demoGOST, query="CANADA", 
        source="KEGG", termIDs=NULL, removeRoot=TRUE,
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, force=FALSE), error_message)
})


test_that("createEnrichMap() must return error when source is a number", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("Assertion on 'arg' failed: Must be of type ", 
                                "'character', not 'double'.")
    
    expect_error(createEnrichMap(gostObject=gostObject, query="toto", 
        source=333, termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, force=FALSE),  error_message)
})

test_that("createEnrichMap() must return error when source is a wrong name", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    expect_error(createEnrichMap(gostObject=gostObject, query="toto",  
        source="test", termIDs=NULL, removeRoot=TRUE, title="network", 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, force=FALSE))
})


test_that("createEnrichMap() must return error when source is GO", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    expect_error(createEnrichMap(gostObject=gostObject, query="toto", 
        source="GO",
        termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, force=FALSE))
})


test_that("createEnrichMap() must return error when removeRoot remove last enriched term", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- paste0("With removal of the root term, there is no ", 
                                "enrichment term left")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, 
                    showCategory=30, groupCategory=FALSE, categoryLabel=1,
                    categoryNode=1, line=1, force=FALSE), error_message)
})


test_that("createEnrichMap() must return error when removeRoot remove last enriched term from term list", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- paste0("With removal of the root term, there is no ", 
                                "enrichment term left")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="TERM_ID",
        termIDs=c("WP:000000"), removeRoot=TRUE, showCategory=30, 
        groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, force=FALSE), error_message)
})


test_that("createEnrichMap() must return error when showCategory negative value", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'showCategory\' parameter must an positive ", 
            "integer or a vector of character strings representing terms.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=-30, 
                    groupCategory=FALSE, categoryLabel=1,
                    categoryNode=1, line=2, force=FALSE), error_message)
})


test_that("createEnrichMap() must return error when showCategory is boolean", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'showCategory\' parameter must an positive ", 
            "integer or a vector of character strings representing terms.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=TRUE, 
                    groupCategory=FALSE, categoryLabel=1,
                    categoryNode=1, line=1, force=TRUE), error_message)
})


test_that("createEnrichMap() must return error when groupCategory is integer", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'groupCategory\' parameter must a logical ", 
                                "(TRUE or FALSE).")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
            source="WP", removeRoot=TRUE, showCategory=30, 
            groupCategory=22, categoryLabel=1,
            categoryNode=1, line=1, force=TRUE), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when categoryLabel is string", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'categoryLabel\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
            source="WP", removeRoot=TRUE, showCategory=30, 
            groupCategory=FALSE, categoryLabel="test",
            categoryNode=1, line=2, force=TRUE), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when cexLabelCategory is negative", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'categoryLabel\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
            source="WP", removeRoot=TRUE, showCategory=30, 
            groupCategory=FALSE, categoryLabel=-1.1,
            categoryNode=1, line=2, force=TRUE), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when categoryNode is negative", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'categoryNode\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=30, 
                    groupCategory=FALSE, categoryLabel=2,
                    categoryNode=-1), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when categoryNode is string", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'categoryNode\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
            source="WP", removeRoot=TRUE, showCategory=30, 
            groupCategory=FALSE, categoryLabel=2,
            categoryNode="te", line=1, force=FALSE), error_message, fixed=TRUE)
})

test_that("createEnrichMap() must return error when line is a string", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'line\' parameter must be a ", 
                                "positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="WP", removeRoot=TRUE,  showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1, line="HI", force=TRUE), 
        error_message, fixed=TRUE)
})

test_that("createEnrichMap() must return error when line is a negative number", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'line\' parameter must be a ", 
                                "positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="WP", removeRoot=TRUE,  showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1, line=-0.3, force=TRUE), 
        error_message, fixed=TRUE)
})

test_that("createEnrichMap() must return error when force is a string", {

    gostTerm <- demoGOST
    
    error_message <- paste0("The \'force\' parameter must a logical ", 
                                "(TRUE or FALSE).")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="WP", removeRoot=TRUE,  showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1, line=1, force="TOTO"), error_message, 
        fixed=TRUE)
})


### Tests createEnrichMapMulti() results

context("createEnrichMapMulti() results")

test_that("createEnrichMapMulti() must return error when gostObjectList is a number", {
    
    error_message <- paste0("The gostObjectList object should be a list ", 
        "of enrichment objects. At least 2 enrichment objects are required.")
    
    expect_error(createEnrichMapMulti(gostObjectList=33, 
        queryList=c("TEST", "Test2"),  source="GO:CC", termIDs=NULL, 
        removeRoot=TRUE, showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1, force=FALSE), error_message)
})


test_that("createEnrichMapMulti() must return error when gostObjectList has only one element", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The gostObjectList object should be a list ", 
        "of enrichment objects. At least 2 enrichment objects are required.")
    
    expect_error(createEnrichMapMulti(gostObjectList=list(gostTerm), 
        queryList=list("TEST"), source="GO:CC", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, force=FALSE), error_message)
})


test_that("createEnrichMapMulti() must return error when queryList is a number", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("he queryList object should be a list of query ", 
        "names. At least 2 query names are required. The number of query ", 
        "names should correspond to the number of enrichment objects.")
    
    expect_error(createEnrichMapMulti(gostObjectList=list(gostTerm, gostTerm), 
        queryList=33, source="GO:CC", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, force=FALSE), error_message)
})


test_that("createEnrichMapMulti() must return error when queryList is longer than gostObjectList", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The queryList object should be a list of query ", 
        "names. At least 2 query names are required. The number of query ", 
        "names should correspond to the number of enrichment objects.")
    
    expect_error(createEnrichMapMulti(gostObjectList=list(gostTerm, gostTerm), 
    queryList=list("TEST", "TEST2", "TEST3"), source="GO:CC", termIDs=NULL, 
    removeRoot=TRUE, showCategory=30, groupCategory=FALSE, categoryLabel=1,
    categoryNode=1, force=FALSE), error_message)
})


test_that("createEnrichMapMulti() must return error when one query in queryList is not in gostObject", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("Each query name present in the ", 
        "\'queryList\' parameter must be present in the associated ", 
        "enrichment object.")
    
    expect_error(createEnrichMapMulti(gostObjectList=list(gostTerm, gostTerm), 
        queryList=list("query_1", "TEST"), source="GO:CC", termIDs=NULL, 
        removeRoot=TRUE, showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, force=FALSE), error_message, fixed=TRUE)
})
