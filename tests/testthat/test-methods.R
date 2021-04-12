### Unit tests for methods.R functions

library(gprofiler2cytoscape)

### Tests prepareInformation() results

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