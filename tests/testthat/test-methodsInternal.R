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
    
    result <- gprofiler2cytoscape:::validateCreateNetworkArguments(
        gostObject=gostObj, source="GO:BP", termIDs=NULL, removeRoot=FALSE,
        fileName="new.cx")
    
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
    
    expect_error(gprofiler2cytoscape:::validateCreateNetworkArguments(
        gostObject=gostObj, source=source, termIDs=NULL, removeRoot=FALSE,
        fileName="hello.cx"), error_message)
    
})

test_that("validateCreateNetworkArguments() must return error when source is TERM_ID and termIDs is NULL", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"),2), 
                                      term_id=c("GO:0051171","GO:0010604"))
    gostObj[["meta"]] <- list()
    
    source <- "TERM_ID"
    
    error_message <- paste0("A vector of terms should be given through", 
                    " the \'termIDs\' parameter when source is \'TERM_ID\'.")
    
    expect_error(gprofiler2cytoscape:::validateCreateNetworkArguments(
        gostObject=gostObj, source=source, termIDs=NULL, removeRoot=FALSE, 
        fileName="toto.cx"), error_message)
    
})

test_that("validateCreateNetworkArguments() must return error when not all term IDS exist", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"),2), 
                                      term_id=c("GO:0051171","GO:0010604"))
    gostObj[["meta"]] <- list()
    
    source <- "TERM_ID"
    
    error_message <- paste0("Not all listed terms are present in the  ",
                                "enrichment results.")
    
    expect_error(gprofiler2cytoscape:::validateCreateNetworkArguments(
        gostObject=gostObj, source=source, 
        termIDs=c("GO:0051171","GO:0010999"), removeRoot=FALSE, 
        fileName="toto.cx"), error_message)
})


test_that("validateCreateNetworkArguments() must return error when removeRoot is a string", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"),2), 
                                      term_id=c("GO:0051171","GO:0010604"))
    gostObj[["meta"]] <- list()
    
    source <- "TERM_ID"
    
    error_message <- paste0("The \'removeRoot\' parameter must be the logical ", 
                                "value TRUE or FALSE.")
    
    expect_error(gprofiler2cytoscape:::validateCreateNetworkArguments(
        gostObject=gostObj, source=source, termIDs=c("GO:0051171"), 
        removeRoot="HI", fileName="test.cx"), error_message)
})


test_that("validateCreateNetworkArguments() must return term IDs are good", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    result <- gprofiler2cytoscape:::validateCreateNetworkArguments(
                gostObject=gostObj, source="TERM_ID", 
                termIDs=c("GO:0051171", "GO:0010604"), fileName="test.cx", 
                removeRoot=FALSE)
    
    expect_true(result)
})


test_that("validateCreateNetworkArguments() must return error when fileName is not a string character", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    error_message <- "The \'fileName\' parameter must a character string."
    
    expect_error(gprofiler2cytoscape:::validateCreateNetworkArguments(
        gostObject = gostObj, source = "GO:BP", termIDs = NULL, 
        removeRoot = FALSE, fileName = 231), error_message)
})


test_that("validateCreateNetworkArguments() must return error when fileName not ending with .cx", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    error_message <- "The \'fileName\' parameter must have \'.cx\' extension."
    
    expect_error(gprofiler2cytoscape:::validateCreateNetworkArguments(
        gostObject=gostObj, source="GO:BP", termIDs=NULL, 
        removeRoot=FALSE, fileName="toto.txt"), error_message)
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
        "{\"@id\":5,\"n\":\"MIRNA:hsa-miR-759\"}]},",
        "{\"edges\":[{\"@id\":3,\"s\":1,\"t\":2,\"i\":\"contains\"},",
        "{\"@id\":6,\"s\":5,\"t\":2,\"i\":\"contains\"}]},",
        "{\"nodeAttributes\":[{\"po\":1,\"n\":\"alias\",\"v\":\"hsa-miR-335-5p\"},",
        "{\"po\":1,\"n\":\"group\",\"v\":\"TERM\"},",
        "{\"po\":2,\"n\":\"alias\",\"v\":\"HERPUD1\"},",
        "{\"po\":2,\"n\":\"group\",\"v\":\"GENE\"},",
        "{\"po\":4,\"n\":\"alias\",\"v\":\"hsa-miR-3180-5p\"},",
        "{\"po\":4,\"n\":\"group\",\"v\":\"TERM\"},",
        "{\"po\":5,\"n\":\"alias\",\"v\":\"hsa-miR-759\"},",
        "{\"po\":5,\"n\":\"group\",\"v\":\"TERM\"}]},",
        "{\"edgeAttributes\":[{\"po\":3,\"n\":\"name\",\"v\":\"MIRNA:hsa-miR-335-5p (contains) ENSG00000051108\"},",
        "{\"po\":3,\"n\":\"source\",\"v\":\"MIRNA:hsa-miR-335-5p\"},",
        "{\"po\":3,\"n\":\"target\",\"v\":\"ENSG00000051108\"},",
        "{\"po\":6,\"n\":\"name\",\"v\":\"MIRNA:hsa-miR-759 (contains) ENSG00000051108\"},",
        "{\"po\":6,\"n\":\"source\",\"v\":\"MIRNA:hsa-miR-759\"},",
        "{\"po\":6,\"n\":\"target\",\"v\":\"ENSG00000051108\"}]},", 
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


### Tests extractNodesAndEdgesInfoForCXJSON() results

context("extractNodesAndEdgesInfoForCXJSON() results")

test_that("extractNodesAndEdgesInfoForCXJSON() must return expected text", {
    
    mirnaDemo <- demoGOST
    
    mirnaDemo$meta$query_metadata$queries[[1]] <- 
        mirnaDemo$meta$query_metadata$queries[[1]][1:4]
    
    mirnaData <- demoGOST$result[demoGOST$result$source == "MIRNA", ]
    
    
    result <- gprofiler2cytoscape:::extractNodesAndEdgesInfoForCXJSON(
        gostResults=mirnaData, gostObject=mirnaDemo)
    
    expected <- list()
    
    expected[["nodes"]] <- data.frame("@id"=c(1, 2, 4, 6, 8, 10),
        "n"=c("MIRNA:hsa-miR-335-5p", "ENSG00000051108", "ENSG00000059728",
            "ENSG00000077616", "MIRNA:hsa-miR-3180-5p", "MIRNA:hsa-miR-759"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["edges"]] <- data.frame("@id"=c(3, 5, 7, 9, 11, 12),
                                "s"=c(1, 1, 1, 8, 10, 10),
                                "t"=c(2, 4, 6, 4, 2, 4),
                                "i"=rep("contains", 6),
                                check.names=FALSE, stringsAsFactors=FALSE)
    
    
    expected[["nodeAttributes"]] <- data.frame(
        "po"=c(1, 1, 2, 2, 4, 4, 6, 6, 8, 8,10, 10),
        "n"=rep(c("alias", "group"), 6),
        "v"=c("hsa-miR-335-5p", "TERM", "HERPUD1", "GENE", "MXD1", "GENE",
            "NAALAD2", "GENE", "hsa-miR-3180-5p", "TERM", 
            "hsa-miR-759", "TERM"),
        check.names=FALSE, stringsAsFactors=FALSE)
   
    expected[["edgeAttributes"]] <- data.frame(
        "po"=c(3, 3, 3, 5, 5, 5, 7, 7, 7, 9, 9, 9, 11, 11, 11, 12, 12, 12),
        "n"=rep(c("name", "source", "target"), 6),
        "v"=c("MIRNA:hsa-miR-335-5p (contains) ENSG00000051108", "MIRNA:hsa-miR-335-5p", "ENSG00000051108", 
            "MIRNA:hsa-miR-335-5p (contains) ENSG00000059728", "MIRNA:hsa-miR-335-5p", "ENSG00000059728",
            "MIRNA:hsa-miR-335-5p (contains) ENSG00000077616", "MIRNA:hsa-miR-335-5p", "ENSG00000077616", 
            "MIRNA:hsa-miR-3180-5p (contains) ENSG00000059728", "MIRNA:hsa-miR-3180-5p", "ENSG00000059728",
            "MIRNA:hsa-miR-759 (contains) ENSG00000051108", "MIRNA:hsa-miR-759", "ENSG00000051108", 
            "MIRNA:hsa-miR-759 (contains) ENSG00000059728", "MIRNA:hsa-miR-759", "ENSG00000059728"),
        check.names=FALSE, stringsAsFactors=FALSE) 
    
    expect_identical(result, expected)
})



### Tests extractNodesAndEdgesWhenNoIntersection() results

context("extractNodesAndEdgesWhenNoIntersection() results")

test_that("extractNodesAndEdgesWhenNoIntersection() must return expected text", {
    
    ccDemo <- demoGOST
    
   # ccDemo$meta$query_metadata$queries[[1]] <- 
    #    ccDemo$meta$query_metadata$queries[[1]][1:2]
    
    ccDemo$meta$genes_metadata$query[[1]]$ensgs <- 
        ccDemo$meta$genes_metadata$query[[1]]$ensgs[1:2]
    
    ccData <- ccDemo$result[ccDemo$result$source == "GO:CC", ]
    
    
    result <- gprofiler2cytoscape:::extractNodesAndEdgesWhenNoIntersection(
        gostResults=ccData, gostObject=ccDemo)
    
    expected <- list()
    
    
    expected[["nodes"]] <- data.frame(
        "id"=c("GO:0005737", "ENSG00000007944", "ENSG00000051108", 
                    "GO:0110165", "GO:0005575", "GO:0005622"),
        "group"=c("TERM", "GENE", "GENE", "TERM", "TERM", "TERM"),
        "alias"=c("cytoplasm", "MYLIP", "HERPUD1", 
                        "cellular anatomical entity", "cellular_component",
                        "intracellular anatomical structure"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["edges"]] <- data.frame(
        "source"=c("GO:0005737", "GO:0005737", "GO:0110165", "GO:0110165", 
                 "GO:0005575", "GO:0005575", "GO:0005622", "GO:0005622"),
         "target"=c("ENSG00000007944", "ENSG00000051108", 
                       "ENSG00000007944", "ENSG00000051108",
                       "ENSG00000007944", "ENSG00000051108",
                       "ENSG00000007944", "ENSG00000051108"),
        "interaction"=rep("contains", 8),
         check.names=FALSE, stringsAsFactors=FALSE)
    
    expect_identical(result, expected)
})




