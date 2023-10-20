### Unit tests for methodsInternal.R functions

library(enrichViewNet)

data(demoGOST)


### Tests isCytoscapeRunning() results

context("isCytoscapeRunning() results")

test_that("isCytoscapeRunning() must return an logical value", {
    
    expect_is(enrichViewNet:::isCytoscapeRunning(), "logical")
})


### Tests validateCreateNetworkArguments() results

context("validateCreateNetworkArguments() results")

test_that("validateCreateNetworkArguments() must return TRUE when all parameters are good", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                      term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    result <- enrichViewNet:::validateCreateNetworkArguments(
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
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
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
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
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
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
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
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
        gostObject=gostObj, source=source, termIDs=c("GO:0051171"), 
        removeRoot="HI", fileName="test.cx"), error_message)
})


test_that("validateCreateNetworkArguments() must return term IDs are good", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    result <- enrichViewNet:::validateCreateNetworkArguments(
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
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
        gostObject = gostObj, source = "GO:BP", termIDs = NULL, 
        removeRoot = FALSE, fileName = 231), error_message)
})


test_that("validateCreateNetworkArguments() must return error when fileName not ending with .cx", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    error_message <- "The \'fileName\' parameter must have \'.cx\' extension."
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
        gostObject=gostObj, source="GO:BP", termIDs=NULL, 
        removeRoot=FALSE, fileName="toto.txt"), error_message)
})


### Tests removeRootTerm() results

context("removeRootTerm() results")

test_that("removeRootTerm() must return same list when not root term present", {
    
    demo <- data.frame(source=c(rep("GO:BP", 3)), 
                        term_id=c("GO:0051171", "GO:0010604", "GO:0014444"))
    
    result <- enrichViewNet:::removeRootTerm(gostResult=demo)
    
    expect_identical(result, demo)
})


test_that("removeRootTerm() must return KEGG root term", {
    
    demo <- data.frame(source=c(rep("GO:BP", 3), rep("KEGG", 2)), 
                       term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                    "KEGG:00000", "KEGG:00010"))
    
    expected <- data.frame(source=c(rep("GO:BP", 3), "KEGG"), 
                       term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                 "KEGG:00010"))
    
    result <- enrichViewNet:::removeRootTerm(gostResult=demo)
    
    expect_identical(result, expected)
})


test_that("removeRootTerm() must return MIRNA root term", {
    
    demo <- data.frame(source=c(rep("GO:BP", 3), "MIRNA", rep("KEGG", 2)), 
                       term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                 "MIRNA:000000", "KEGG:00044", "KEGG:00010"))
    
    expected <- data.frame(source=c(rep("GO:BP", 3), rep("KEGG", 2)), 
                           term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                     "KEGG:00044", "KEGG:00010"))
    
    result <- enrichViewNet:::removeRootTerm(gostResult=demo)
    
    expect_identical(result, expected)
})


test_that("removeRootTerm() must return GO:MF root term", {
    
    demo <- data.frame(source=c(rep("GO:BP", 3), "GO:MF", rep("KEGG", 2)), 
                       term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                 "GO:0003674", "KEGG:00044", "KEGG:00010"))
    
    expected <- data.frame(source=c(rep("GO:BP", 3), rep("KEGG", 2)), 
                           term_id=c("GO:0051171", "GO:0010604", "GO:0014444",
                                     "KEGG:00044", "KEGG:00010"))
    
    result <- enrichViewNet:::removeRootTerm(gostResult=demo)
    
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
    
    result <- enrichViewNet:::createMetaDataSectionCXJSON()
    
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
    
    result <- enrichViewNet:::createCytoscapeCXJSON(
                         gostResults=mirnaData, gostObject=mirnaDemo, 
                         title = "MIRNA")
    
    expected <- paste0(
        "[{\"metaData\":[{\"name\":\"nodes\",\"version\":\"1.0\"},", 
        "{\"name\":\"edges\",\"version\":\"1.0\"},{\"name\":\"edgeAttributes\",",
        "\"version\":\"1.0\"},{\"name\":\"nodeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyHiddenAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyNetworkRelations\",\"version\":\"1.0\"},",
        "{\"name\":\"cyGroups\",\"version\":\"1.0\"},",
        "{\"name\":\"networkAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyTableColumn\",\"version\":\"1.0\"},",
        "{\"name\":\"cySubNetworks\",\"version\":\"1.0\"}]},",
        "{\"networkAttributes\":[{\"n\":\"name\",\"v\":\"MIRNA\"}]},",
        "{\"nodes\":[{\"@id\":1,\"n\":\"ENSG00000051108\"},",
        "{\"@id\":2,\"n\":\"MIRNA:hsa-miR-335-5p\"},",
        "{\"@id\":3,\"n\":\"MIRNA:hsa-miR-759\"}]},",
        "{\"edges\":[{\"@id\":4,\"s\":2,\"t\":1,\"i\":\"contains\"},",
        "{\"@id\":5,\"s\":3,\"t\":1,\"i\":\"contains\"}]},", 
        "{\"nodeAttributes\":[{\"po\":1,\"n\":\"alias\",\"v\":\"HERPUD1\"},",
        "{\"po\":1,\"n\":\"group\",\"v\":\"GENE\"},",
        "{\"po\":2,\"n\":\"alias\",\"v\":\"GENE\"},",
        "{\"po\":3,\"n\":\"group\",\"v\":\"GENE\"},",
        "{\"po\":2,\"n\":\"alias\",\"v\":\"GENE\"},",
        "{\"po\":3,\"n\":\"group\",\"v\":\"GENE\"}]},",
        "{\"edgeAttributes\":[{\"po\":4,\"n\":\"name\",\"v\":\"MIRNA:hsa-miR-335-5p (contains) ENSG00000051108\"},",
        "{\"po\":5,\"n\":\"name\",\"v\":\"MIRNA:hsa-miR-759 (contains) ENSG00000051108\"},",
        "{\"po\":4,\"n\":\"source\",\"v\":\"MIRNA:hsa-miR-335-5p\"},",
        "{\"po\":5,\"n\":\"source\",\"v\":\"MIRNA:hsa-miR-759\"},",
        "{\"po\":4,\"n\":\"target\",\"v\":\"ENSG00000051108\"},",
        "{\"po\":5,\"n\":\"target\",\"v\":\"ENSG00000051108\"}]},",
        "{\"cyHiddenAttributes\":[{\"n\":\"layoutAlgorithm\",\"y\":\"yFiles Circular Layout\"}]},",
        "{\"metaData\":[{\"name\":\"nodes\",\"version\":\"1.0\"},",
        "{\"name\":\"edges\",\"version\":\"1.0\"},{\"name\":\"edgeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"nodeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyHiddenAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyNetworkRelations\",\"version\":\"1.0\"},",
        "{\"name\":\"cyGroups\",\"version\":\"1.0\"},",
        "{\"name\":\"networkAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyTableColumn\",\"version\":\"1.0\"},",
        "{\"name\":\"cySubNetworks\",\"version\":\"1.0\"}]},",
        "{\"status\":[{\"error\":\"\",\"success\":true}]}]")
    
    expect_equal(result, expected)
})


### Tests extractNodesAndEdgesInfoForCXJSON() results

context("extractNodesAndEdgesInfoForCXJSON() results")

test_that("extractNodesAndEdgesInfoForCXJSON() must return expected text", {
    
    mirnaDemo <- demoGOST
    
    mirnaDemo$meta$query_metadata$queries[[1]] <- 
        mirnaDemo$meta$query_metadata$queries[[1]][1:4]
    
    mirnaData <- demoGOST$result[demoGOST$result$source == "MIRNA", ]
    
    
    result <- enrichViewNet:::extractNodesAndEdgesInfoForCXJSON(
        gostResults=mirnaData, gostObject=mirnaDemo)
    
    expected <- list()
    
    expected[["nodes"]] <- data.frame("@id"=c(1, 2, 3, 4, 5, 6),
        "n"=c("ENSG00000059728", "ENSG00000077616", "ENSG00000051108", 
                "MIRNA:hsa-miR-335-5p", 
                "MIRNA:hsa-miR-3180-5p", "MIRNA:hsa-miR-759"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["edges"]] <- data.frame("@id"=c(7, 8, 9, 10, 11, 12),
                                "s"=c(4, 4, 4, 5, 6, 6),
                                "t"=c(1, 2, 3, 1, 3, 1),
                                "i"=rep("contains", 6),
                                check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["nodeAttributes"]] <- data.frame(
        "po"=c(1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6),
        "n"=c("alias", "alias", "alias", "group", "group", "group", 
                  "alias", "alias", "alias", "group", "group", "group"),
        "v"=c("MXD1",  "NAALAD2", "HERPUD1", "GENE", "GENE",  "GENE",  "hsa-miR-335-5p" ,
              "hsa-miR-3180-5p", "hsa-miR-759",  "GENE",  "GENE",  "GENE" ),
        check.names=FALSE, stringsAsFactors=FALSE)
   
    expected[["edgeAttributes"]] <- data.frame(
        "po"=c(7, 8, 9, 10, 11, 12, 7, 8, 9, 10, 11, 12, 7, 8, 9, 10, 11, 12),
        "n"=c(rep(c("name"), 6), rep(c("source"), 6), rep(c("target"), 6)),
        "v"=c("MIRNA:hsa-miR-335-5p (contains) ENSG00000059728", 
            "MIRNA:hsa-miR-335-5p (contains) ENSG00000077616",
            "MIRNA:hsa-miR-335-5p (contains) ENSG00000051108", 
            "MIRNA:hsa-miR-3180-5p (contains) ENSG00000059728",
            "MIRNA:hsa-miR-759 (contains) ENSG00000051108", 
            "MIRNA:hsa-miR-759 (contains) ENSG00000059728",
            "MIRNA:hsa-miR-335-5p", "MIRNA:hsa-miR-335-5p",
            "MIRNA:hsa-miR-335-5p", "MIRNA:hsa-miR-3180-5p",
            "MIRNA:hsa-miR-759", "MIRNA:hsa-miR-759",
            "ENSG00000059728", "ENSG00000077616",
            "ENSG00000051108", "ENSG00000059728",
            "ENSG00000051108","ENSG00000059728"),
        check.names=FALSE, stringsAsFactors=FALSE) 
    
    expect_equal(result$nodes, expected$nodes)
    expect_equal(result$edges, expected$edges)
    expect_equal(result$nodeAttributes, expected$nodeAttributes)
    expect_equal(result$edgeAttributes, expected$edgeAttributes)
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
    
    
    result <- enrichViewNet:::extractNodesAndEdgesWhenNoIntersection(
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




