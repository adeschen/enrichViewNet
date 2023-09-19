### Unit tests for methodsEmap.R functions

library(enrichViewNet)

data(demoGOST)



### Tests createEnrichMap() results

context("createEnrichMap() results")

test_that("createEnrichMap() must return error when gostObject is a number", {
    
    error_message <- paste0("The gostObject object should be a list ", 
                        "with meta and result as entries corresponding ", 
                        "to gprofiler2 enrichment output.")
    
    expect_error(createEnrichMap(gostObject=33, source="GO:CC",
        termIDs=NULL, removeRoot=TRUE, title="gprofiler network"), 
        error_message)
})

test_that("createEnrichMap() must return error when gostObject is a string character", {
    
    error_message <- paste0("The gostObject object should be a list ", 
                        "with meta and result as entries corresponding ", 
                        "to gprofiler2 enrichment output.")
    
    expect_error(createEnrichMap(gostObject="TEST", source="GO:CC",
        termIDs=NULL, removeRoot=TRUE, title="gprofiler network"), 
        error_message)
})

test_that("createEnrichMap() must return error when source is a number", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("Assertion on 'arg' failed: Must be of type ", 
                                "'character', not 'double'.")
    
    expect_error(createEnrichMap(gostObject=gostObject, source=333,
        termIDs=NULL, removeRoot=TRUE, title="gprofiler network"), 
        error_message)
})

test_that("createEnrichMap() must return error when source is a wrong name", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    expect_error(createEnrichMap(gostObject=gostObject, source="test",
        termIDs=NULL, removeRoot=TRUE, title="gprofiler network"))
})


test_that("createEnrichMap() must return error when source is GO", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    expect_error(createEnrichMap(gostObject=gostObject, source="GO",
        termIDs=NULL, removeRoot=TRUE, title="gprofiler network"))
})


test_that("createEnrichMap() must return error when removeRoot remove last enriched term", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- paste0("With removal of the root term, there is no ", 
                                "enrichment term left")
    
    expect_error(createEnrichMap(gostObject=gostTerm, source="WP", 
                    removeRoot=TRUE, title="gprofiler network"), 
                    error_message)
})


test_that("createEnrichMap() must return error when removeRoot remove last enriched term from term list", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- paste0("With removal of the root term, there is no ", 
                                "enrichment term left")
    
    expect_error(createEnrichMap(gostObject=gostTerm, source="TERM_ID",
        termIDs=c("WP:000000"), removeRoot=TRUE), error_message)
})


