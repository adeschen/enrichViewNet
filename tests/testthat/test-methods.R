### Unit tests for methods.R functions

library(enrichViewNet)

data("demoGOST")

### Tests createNetwork() results

context("createNetwork() results")

test_that("createNetwork() must return error when gostObject is a number", {

    error_message <- paste0("The gostObject object should be a list ", 
                        "with meta and result as entries corresponding ", 
                        "to gprofiler2 enrichment output.")
    
    expect_error(createNetwork(gostObject=33), error_message)
})

test_that("createNetwork() must return error when gostObject is a string character", {

    error_message <- paste0("The gostObject object should be a list ", 
                        "with meta and result as entries corresponding ", 
                        "to gprofiler2 enrichment output.")
    
    expect_error(createNetwork(gostObject="TEST"), error_message)
})

test_that("createNetwork() must return error when source is a number", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("Assertion on 'arg' failed: Must be of type ", 
                                "'character', not 'double'.")
    
    expect_error(createNetwork(gostObject=gostObject, source=33), 
                    error_message)
})

test_that("createNetwork() must return error when source is a wrong name", {

    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    expect_error(createNetwork(gostObject=gostObject, source="test"))
})


test_that("createNetwork() must return error when source is GO", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    expect_error(createNetwork(gostObject=gostObject, source="GO"))
})


test_that("createNetwork() must return error when query not present in the GOST object", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("With selected query, there is no enrichment term")
    
    expect_error(createNetwork(gostObject=gostTerm, source="WP", 
                               query="Information",
                               removeRoot=TRUE), error_message)
})


test_that("createNetwork() must return error when query is a number", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'query\' parameter must be a character ", 
                                "string or \'NULL\'.")
    
    expect_error(createNetwork(gostObject=gostTerm, source="WP", query=33, 
                               removeRoot=TRUE), error_message)
})


test_that("createNetwork() must return error when removeRoot remove last enriched term", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- paste0("After filtering on the enriched terms, there is ", 
                                    "no enrichment term left")
    
    expect_error(createNetwork(gostObject=gostTerm, source="WP", 
                    removeRoot=TRUE), error_message)
})


test_that("createNetwork() must return error when removeRoot remove last enriched term from term list", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- paste0("After filtering on the enriched terms, there is ", 
                                "no enrichment term left")
    
    expect_error(createNetwork(gostObject=gostTerm, source="TERM_ID",
        termIDs=c("WP:000000"), removeRoot=TRUE), error_message)
})


test_that("createNetwork() must return error when collection is a number", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- "The \'collection\' parameter must a character string."
    
    expect_error(createNetwork(gostObject=gostTerm, source="TERM_ID",
        termIDs=c("WP:000000"), removeRoot=FALSE, collection=33), error_message)
})


test_that("createNetwork() must return error when fileName is a number", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- "The \'fileName\' parameter must a character string."
    
    expect_error(createNetwork(gostObject=gostTerm, source="TERM_ID",
        termIDs=c("WP:000000"), removeRoot=FALSE, fileName=33), error_message)
})


test_that("createNetwork() must return error when fileName is a logical", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- "The \'fileName\' parameter must a character string."
    
    expect_error(createNetwork(gostObject=gostTerm, source="TERM_ID",
        termIDs=c("WP:000000"), removeRoot=FALSE, fileName=FALSE), 
        error_message)
})


test_that("createNetwork() must return error when fileName has .txt extension", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- "The \'fileName\' parameter must have \'.cx\' extension."
    
    expect_error(createNetwork(gostObject=gostTerm, source="TERM_ID",
        termIDs=c("WP:000000"), removeRoot=FALSE, fileName="toto.txt"), 
        error_message)
})


test_that("createNetwork() must return error when no query set and more than one query present", {
    
    gostTerm <- demoGOST
    gostTerm$result$query[1] <- "Test negative"
    
    error_message <- paste0("Multiple queries are present in the results, ",
        "the \'query\' parameter should be used to select one")
    
    expect_error(createNetwork(gostObject=gostTerm, source="CORUM", 
        removeRoot=FALSE, fileName="toto.cx"), error_message)
})


test_that("createNetwork() must return TRUE when all parameters ok and not cytoscape open", {
    
    gostTerm <- demoGOST
    
    tmpDir <- tempdir()
    tmpFile <- paste0(tmpDir, "/toto_0001.cx")
    
    
    if (!isCytoscapeRunning() && !file.exists(tmpFile)) {
        result <- createNetwork(gostObject=gostTerm, source="CORUM", 
                            removeRoot=FALSE, fileName=tmpFile)
    
        expect_true(result)
        expect_true(file.exists(tmpFile))
    }
    
    if (file.exists(tmpFile)) {
        file.remove(tmpFile)
    }
})

