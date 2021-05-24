### Unit tests for methodsInternal.R functions

library(gprofiler2cytoscape)

data(demoGOST)


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


### Tests createMetaDataSectionCXJSON() results

context("createMetaDataSectionCXJSON() results")

test_that("createMetaDataSectionCXJSON() must return expected text", {
    
    expected <- paste0(
        "{\"metaData\":[{\"name\":\"nodes\",\"version\":\"1.0\"},", 
        "{\"name\":\"edges\",\"version\":\"1.0\"},",
        "{\"name\":\"edgeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"nodeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyHiddenAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyNetworkRelations\",\"version\":\"1.0\"},",
        "{\"name\":\"cyGroups\",\"version\":\"1.0\"},", 
        "{\"name\":\"networkAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyTableColumn\",\"version\":\"1.0\"},",
        "{\"name\":\"cySubNetworks\",\"version\":\"1.0\"}]}")
    
    result <- gprofiler2cytoscape:::createMetaDataSectionCXJSON()
    
    expect_identical(as.character(result), expected)
    expect_s3_class(result, "json")
})


### Tests createCytoscapeCXJSON() results

context("createCytoscapeCXJSON() results")

test_that("createCytoscapeCXJSON() must return expected text", {
    
    mirnaDemo <- demoGOST
    mirnaDemo$meta$query_metadata$queries[[1]] <- 
            mirnaDemo$meta$query_metadata$queries[[1]][1:2]
    mirnaData <- demoGOST$result[demoGOST$result$source == "MIRNA", ]
    
    result <- gprofiler2cytoscape:::createCytoscapeCXJSON(
                         gostResults = mirnaData, gostObject = mirnaDemo, 
                         title = "MIRNA")
                         
    expected <- paste0(
        "[{\"metaData\":[{\"name\":\"nodes\",\"version\":\"1.0\"},", 
        "{\"name\":\"edges\",\"version\":\"1.0\"},",
        "{\"name\":\"edgeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"nodeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyHiddenAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyNetworkRelations\",\"version\":\"1.0\"},",
        "{\"name\":\"cyGroups\",\"version\":\"1.0\"},", 
        "{\"name\":\"networkAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyTableColumn\",\"version\":\"1.0\"},",
        "{\"name\":\"cySubNetworks\",\"version\":\"1.0\"}]},",
        "{\"networkAttributes\":[{\"n\":\"name\",\"v\":\"MIRNA\"}]},",
        "{\"nodes\":[{\"@id\":1,\"n\":\"MIRNA:hsa-miR-335-5p\"},",
        "{\"@id\":2,\"n\":\"ENSG00000051108\"},",
        "{\"@id\":4,\"n\":\"MIRNA:hsa-miR-3180-5p\"},",
        "{\"@id\":5,\"n\":\"MIRNA:hsa-miR-759\"},",
        "{\"@id\":6,\"n\":\"ENSG00000051108\"}]},",
        "{\"edges\":[{\"@id\":3,\"s\":1,\"t\":2,\"i\":\"contains\"},",
        "{\"@id\":7,\"s\":5,\"t\":6,\"i\":\"contains\"}]},",
        "{\"nodeAttributes\":[{\"po\":1,\"n\":\"alias\",\"v\":\"hsa-miR-335-5p\"},",
        "{\"po\":1,\"n\":\"group\",\"v\":\"TERM\"},",
        "{\"po\":2,\"n\":\"alias\",\"v\":\"HERPUD1\"},",
        "{\"po\":2,\"n\":\"group\",\"v\":\"GENE\"},",
        "{\"po\":4,\"n\":\"alias\",\"v\":\"hsa-miR-3180-5p\"},",
        "{\"po\":4,\"n\":\"group\",\"v\":\"TERM\"},",
        "{\"po\":5,\"n\":\"alias\",\"v\":\"hsa-miR-759\"},",
        "{\"po\":5,\"n\":\"group\",\"v\":\"TERM\"},",
        "{\"po\":6,\"n\":\"alias\",\"v\":\"HERPUD1\"},",
        "{\"po\":6,\"n\":\"group\",\"v\":\"GENE\"}]},",
        "{\"edgeAttributes\":[{\"po\":3,\"n\":\"name\",\"v\":\"MIRNA:hsa-miR-335-5p (contains) ENSG00000051108\"},",
        "{\"po\":3,\"n\":\"source\",\"v\":\"MIRNA:hsa-miR-335-5p\"},",
        "{\"po\":3,\"n\":\"target\",\"v\":\"ENSG00000051108\"},",
        "{\"po\":7,\"n\":\"name\",\"v\":\"MIRNA:hsa-miR-759 (contains) ENSG00000051108\"},",
        "{\"po\":7,\"n\":\"source\",\"v\":\"MIRNA:hsa-miR-759\"},",
        "{\"po\":7,\"n\":\"target\",\"v\":\"ENSG00000051108\"}]},", 
        "{\"cyHiddenAttributes\":[{\"n\":\"layoutAlgorithm\",\"y\":\"yFiles Circular Layout\"}]},",
        "{\"metaData\":[{\"name\":\"nodes\",\"version\":\"1.0\"},", 
        "{\"name\":\"edges\",\"version\":\"1.0\"},",
        "{\"name\":\"edgeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"nodeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyHiddenAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyNetworkRelations\",\"version\":\"1.0\"},",
        "{\"name\":\"cyGroups\",\"version\":\"1.0\"},", 
        "{\"name\":\"networkAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyTableColumn\",\"version\":\"1.0\"},",
        "{\"name\":\"cySubNetworks\",\"version\":\"1.0\"}]},",
        "{\"status\":[{\"error\":\"\",\"success\":true}]}]")
    
    expect_identical(result, expected)
})





