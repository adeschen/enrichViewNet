### Unit tests for methodsEmapInternal.R functions

library(enrichViewNet)

data(demoGOST)


### Tests validateCreateEnrichMapArguments() results

context("validateCreateEnrichMapArguments() results")


test_that("validateCreateEnrichMapArguments() must return expected result", {
    
    result <- validateCreateEnrichMapArguments(gostObject=demoGOST, 
        query="query_1", source="GO:CC", termIDs=NULL, removeRoot=TRUE, 
        title="network", showCategory=30, groupCategory=FALSE, 
        cexLabelCategory=1, cexCategory=1)
    
    expect_true(result)
})
