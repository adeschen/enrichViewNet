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
    gostObj[["result"]] <- data.frame(a=c(1,2), b=c(2,4))
    gostObj[["meta"]] <- list()
    
    result <- gprofiler2cytoscape:::validateCreateNetworkArguments(gostObject=gostObj, 
                    source = "REACT", termIDs = NULL)
    
    expect_true(result)
})


 