### Unit tests for methodsInternal.R functions

library(gprofiler2cytoscape)


### Tests isCytoscapeRunning() results

context("isCytoscapeRunning() results")

test_that("isCytoscapeRunning() must return an logical value", {
    
    expect_is(gprofiler2cytoscape:::isCytoscapeRunning(), "logical")
})


### Tests validateCreateNetworkArguments() results

context("validateCreateNetworkArguments() results")

test_that("validateCreateNetworkArguments() must return TRUE when all parameters are good", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                      term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    result <- gprofiler2cytoscape:::validateCreateNetworkArguments(gostObject=gostObj, 
                    source = "GO:BP", termIDs = NULL)
    
    expect_true(result)
})

test_that("validateCreateNetworkArguments() must return error when source is  when all parameters are good", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    source <- "REAC"
    
    error_message <- paste0("There is no enriched term for the selected ", 
                            "source \'", source, "\'.")
    
    expect_error(gprofiler2cytoscape:::validateCreateNetworkArguments(gostObject=gostObj, 
            source=source, termIDs=NULL), error_message)
    
})

test_that("validateCreateNetworkArguments() must return error when source is TERM_ID and termIDs is NULL", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"),2), 
                                      term_id=c("GO:0051171","GO:0010604"))
    gostObj[["meta"]] <- list()
    
    source <- "TERM_ID"
    
    error_message <- paste0("A vector of terms should be given through", 
                    " the \'termIDs\' parameter when source is \'TERM_ID\'.")
    
    expect_error(gprofiler2cytoscape:::validateCreateNetworkArguments(gostObject=gostObj, 
                source=source, termIDs=NULL), error_message)
    
})


