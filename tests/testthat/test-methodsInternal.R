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
                    source = "GO:BP", termIDs = NULL, removeRoot=FALSE)
    
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
            source=source, termIDs=NULL, removeRoot=FALSE), error_message)
    
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
                source=source, termIDs=NULL, removeRoot=FALSE), error_message)
    
})

test_that("validateCreateNetworkArguments() must return error when not all term IDS exist", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"),2), 
                                      term_id=c("GO:0051171","GO:0010604"))
    gostObj[["meta"]] <- list()
    
    source <- "TERM_ID"
    
    error_message <- paste0("Not all listed terms are present in the  ",
                                "enrichment results.")
    
    expect_error(gprofiler2cytoscape:::validateCreateNetworkArguments(gostObject=gostObj, 
        source=source, termIDs=c("GO:0051171","GO:0010999"), removeRoot=FALSE), error_message)
})


test_that("validateCreateNetworkArguments() must return error when removeRoot is a string", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"),2), 
                                      term_id=c("GO:0051171","GO:0010604"))
    gostObj[["meta"]] <- list()
    
    source <- "TERM_ID"
    
    error_message <- paste0("The \'removeRoot\' parameter must be the logical ", 
                                "value TRUE or FALSE.")
    
    expect_error(gprofiler2cytoscape:::validateCreateNetworkArguments(gostObject=gostObj, 
        source=source, termIDs=c("GO:0051171"), removeRoot="HI"), error_message)
})


test_that("validateCreateNetworkArguments() must return term IDs are good", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                      term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    result <- gprofiler2cytoscape:::validateCreateNetworkArguments(gostObject=gostObj, 
                 source="TERM_ID", termIDs=c("GO:0051171", "GO:0010604"), removeRoot=FALSE)
    
    expect_true(result)
})


### Tests removeRootTerm() results

context("removeRootTerm() results")

test_that("removeRootTerm() must return same list when not root term present", {
    
    demo <- data.frame(source=c(rep("GO:BP", 3)), 
                        term_id=c("GO:0051171", "GO:0010604", "GO:0014444"))
    
    result <- gprofiler2cytoscape:::removeRootTerm(gostResult=demo)
    
    expect_identical(result, demo)
})


test_that("removeRootTerm() must return KEGG root term", {
    
    demo <- data.frame(source=c(rep("GO:BP", 3), rep("KEGG", 2)), 
                       term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                    "KEGG:00000", "KEGG:00010"))
    
    expected <- data.frame(source=c(rep("GO:BP", 3), "KEGG"), 
                       term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                 "KEGG:00010"))
    
    result <- gprofiler2cytoscape:::removeRootTerm(gostResult=demo)
    
    expect_identical(result, expected)
})


test_that("removeRootTerm() must return MIRNA root term", {
    
    demo <- data.frame(source=c(rep("GO:BP", 3), "MIRNA", rep("KEGG", 2)), 
                       term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                 "MIRNA:000000", "KEGG:00044", "KEGG:00010"))
    
    expected <- data.frame(source=c(rep("GO:BP", 3), rep("KEGG", 2)), 
                           term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                     "KEGG:00044", "KEGG:00010"))
    
    result <- gprofiler2cytoscape:::removeRootTerm(gostResult=demo)
    
    expect_identical(result, expected)
})


test_that("removeRootTerm() must return GO:MF root term", {
    
    demo <- data.frame(source=c(rep("GO:BP", 3), "GO:MF", rep("KEGG", 2)), 
                       term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                 "GO:0003674", "KEGG:00044", "KEGG:00010"))
    
    expected <- data.frame(source=c(rep("GO:BP", 3), rep("KEGG", 2)), 
                           term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                     "KEGG:00044", "KEGG:00010"))
    
    result <- gprofiler2cytoscape:::removeRootTerm(gostResult=demo)
    
    expect_identical(result, expected)
})
