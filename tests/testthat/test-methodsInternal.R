### Unit tests for methodsInternal.R functions

library(enrichViewNet)

data(demoGOST)
data(parentalNapaVsDMSOEnrichment)


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
        query="test", title="Network", collection="test collection", 
        fileName="new.cx")
    
    expect_true(result)
})

test_that("validateCreateNetworkArguments() must return error when no enriched term and all parameters are good", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    source <- "REAC"
    
    error_message <- paste0("There is no enriched term for the selected ", 
                            "source \'", source, "\'.")
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
        gostObject=gostObj, source=source, termIDs=NULL, removeRoot=FALSE,
        fileName="hello.cx", query=NULL), error_message)
    
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
        fileName="toto.cx", query=NULL), error_message)
    
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
        fileName="toto.cx", query=NULL), error_message)
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
        removeRoot="HI", fileName="test.cx", query=NULL), error_message)
})


test_that("validateCreateNetworkArguments() must return term IDs are good", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    result <- enrichViewNet:::validateCreateNetworkArguments(
            gostObject=gostObj, source="TERM_ID", 
            termIDs=c("GO:0051171", "GO:0010604"), removeRoot=FALSE, query=NULL, 
            collection="Collection", title="Network", fileName="test.cx")
    
    expect_true(result)
})


test_that("validateCreateNetworkArguments() must return error when Collection is not a string character", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    error_message <- "The \'collection\' parameter must a character string."
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
        gostObject = gostObj, source = "GO:BP", termIDs = NULL, 
        removeRoot = FALSE, query=NULL, collection=22,
        title="Test", fileName="file.cx"), error_message)
})


test_that("validateCreateNetworkArguments() must return error when title is not a string character", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                      term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    error_message <- "The \'title\' parameter must a character string."
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
        gostObject = gostObj, source = "GO:BP", termIDs = NULL, 
        removeRoot = FALSE, query=NULL, collection="Test",
        title=111, fileName="file.cx"), error_message)
})


test_that("validateCreateNetworkArguments() must return error when fileName is not a string character", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    error_message <- "The \'fileName\' parameter must a character string."
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
        gostObject = gostObj, source = "GO:BP", termIDs = NULL, 
        removeRoot = FALSE, query=NULL, collection="Collection test",
        title="Test", fileName = 231), error_message)
})


test_that("validateCreateNetworkArguments() must return error when fileName not ending with .cx", {
    
    gostObj <- list()
    gostObj[["result"]] <- data.frame(source=c(rep("GO:BP"), 2), 
                                        term_id=c("GO:0051171", "GO:0010604"))
    gostObj[["meta"]] <- list()
    
    error_message <- "The \'fileName\' parameter must have \'.cx\' extension."
    
    expect_error(enrichViewNet:::validateCreateNetworkArguments(
        gostObject=gostObj, source="GO:BP", termIDs=NULL, 
        removeRoot=FALSE, query=NULL, collection="Collection test", 
        title="Test", fileName="toto.txt"), error_message)
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
    
    set.seed(121)
    
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


### Tests extractNodesAndEdgesWhenNoIntersectionForCXJSON() results

context("extractNodesAndEdgesWhenNoIntersectionForCXJSON() results")

test_that("extractNodesAndEdgesWhenNoIntersectionForCXJSON() must return expected text", {
    
    mirnaDemo <- demoGOST
    
    mirnaDemo$meta$query_metadata$queries[[1]] <- 
        mirnaDemo$meta$query_metadata$queries[[1]][1:4]
    
    mirnaData <- demoGOST$result[demoGOST$result$source == "MIRNA", ]
    
    set.seed(121)
    result <- enrichViewNet:::extractNodesAndEdgesWhenNoIntersectionForCXJSON(
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
        "v"=c("MXD1",  "NAALAD2", "HERPUD1", "GENE", "GENE",  
                "GENE",  "hsa-miR-335-5p" ,
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


### Tests extractNodesAndEdgesInfoWhenIntersection() results

context("extractNodesAndEdgesInfoWhenIntersection() results")

test_that("extractNodesAndEdgesInfoWhenIntersection() must return expected text", {
    
    set.seed(1212)
    
    ccDemo <- parentalNapaVsDMSOEnrichment
    
    ccData <- ccDemo$result[ccDemo$result$term_id %in% 
                                            c("WP:WP1742", "WP:WP410"), ]
    
    result <- enrichViewNet:::extractNodesAndEdgesInfoWhenIntersection(
        gostResults=ccData)
    
    expected <- list()
    
    expected[["geneNodes"]] <- data.frame(
        "id"=c("ENSG00000105327", "ENSG00000124762", "ENSG00000141682",
               "ENSG00000051108", "ENSG00000133639", "ENSG00000141232", 
               "ENSG00000172216"),
        "group"=rep("GENE", 7),
        "alias"=c("ENSG00000105327", "ENSG00000124762", "ENSG00000141682",
                  "ENSG00000051108", "ENSG00000133639", "ENSG00000141232", 
                  "ENSG00000172216"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["termNodes"]] <- data.frame(
        "id"=c("WP:WP1742", "WP:WP410"),
        "group"=c("TERM", "TERM"),
        "alias"=c("TP53 network", "Exercise-induced circadian regulation"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["edges"]] <- data.frame(
        "source"=c("WP:WP1742", "WP:WP1742", "WP:WP1742", "WP:WP410", 
                   "WP:WP410", "WP:WP410", "WP:WP410"),
        "target"=c("ENSG00000105327", "ENSG00000124762", "ENSG00000141682",
                   "ENSG00000051108", "ENSG00000133639", "ENSG00000141232", 
                   "ENSG00000172216"),
        "interaction"=rep("contains", 7),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expect_identical(result, expected)
})


test_that("extractNodesAndEdgesInfoWhenIntersection() must return expected text", {
    
    set.seed(1212)
    
    ccDemo <- parentalNapaVsDMSOEnrichment
    
    ccData <- ccDemo$result[ccDemo$result$term_id %in% 
                                c("WP:WP4879", "WP:WP395"), ]
    
    result <- enrichViewNet:::extractNodesAndEdgesInfoWhenIntersection(
        gostResults=ccData)
    
    expected <- list()
    
    expected[["geneNodes"]] <- data.frame(
        "id"=c("ENSG00000124762", "ENSG00000171223", "ENSG00000172216",
               "ENSG00000245848", "ENSG00000170345", "ENSG00000184557"),
        "group"=rep("GENE", 6),
        "alias"=c("ENSG00000124762", "ENSG00000171223", "ENSG00000172216",
                  "ENSG00000245848", "ENSG00000170345", "ENSG00000184557"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["termNodes"]] <- data.frame(
        "id"=c("WP:WP4879", "WP:WP395"),
        "group"=c("TERM", "TERM"),
        "alias"=c(paste0("Overlap between signal transduction pathways", 
            " contributing to LMNA laminopathies"), "IL-4 signaling pathway"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["edges"]] <- data.frame(
        "source"=c("WP:WP4879", "WP:WP4879", "WP:WP4879", "WP:WP4879", 
                   "WP:WP395", "WP:WP395", "WP:WP395", "WP:WP395"),
        "target"=c("ENSG00000124762", "ENSG00000171223", "ENSG00000172216",
                   "ENSG00000245848", "ENSG00000172216", "ENSG00000245848",
                   "ENSG00000170345", "ENSG00000184557"),
        "interaction"=rep("contains", 8),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expect_identical(result, expected)
})


### Tests extractNodesAndEdgesWhenIntersectionForCXJSON() results

context("extractNodesAndEdgesWhenIntersectionForCXJSON() results")

test_that("extractNodesAndEdgesWhenIntersectionForCXJSON() must return expected text", {
    
    set.seed(1212)
    
    ccDemo <- parentalNapaVsDMSOEnrichment
    
    ccData <- ccDemo$result[ccDemo$result$term_id %in% 
                                c("WP:WP1742", "WP:WP410"), ]
    
    result <- enrichViewNet:::extractNodesAndEdgesWhenIntersectionForCXJSON(
        gostResults=ccData)
    
    expected <- list()
    
    expected[["nodes"]] <- data.frame(
        "@id"=seq_len(9),
        "n"=c("ENSG00000105327", "ENSG00000124762", "ENSG00000141682",
               "ENSG00000051108", "ENSG00000133639", "ENSG00000141232", 
               "ENSG00000172216", "WP:WP1742", "WP:WP410"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["edges"]] <- data.frame(
        "@id"=seq_len(7) + 9,
        "s"=c(8, 8, 8, 9, 9, 9, 9),
        "t"=c(4, 1, 2, 5, 6, 3, 7),
        "i"=rep("contains", 7),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["nodeAttributes"]] <- data.frame(
        "po"=c(seq_len(7), seq_len(7), 8, 9, 8, 9),
        "n"=c(rep("alias", 7), rep("group", 7), "alias", "alias", "group",
                "group"),
        "v"=c("ENSG00000105327", "ENSG00000124762", "ENSG00000141682",
                "ENSG00000051108", "ENSG00000133639", "ENSG00000141232", 
                "ENSG00000172216", rep("GENE", 7), "WP:WP1742", "WP:WP410",
                rep("TERM", 2)),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["edgeAttributes"]] <- data.frame(
        "po"=c(seq_len(7)+9, seq_len(7)+9, seq_len(7)+9),
        "n"=c(rep("name", 7), rep("source", 7), rep("target", 7)),
        "v"=c("WP:WP1742 (contains) ENSG00000105327",
              "WP:WP1742 (contains) ENSG00000124762",
              "WP:WP1742 (contains) ENSG00000141682",
              "WP:WP410 (contains) ENSG00000051108",
              "WP:WP410 (contains) ENSG00000133639",
              "WP:WP410 (contains) ENSG00000141232",
              "WP:WP410 (contains) ENSG00000172216",
              rep("WP:WP1742", 3), rep("WP:WP410", 4),
              "ENSG00000105327", "ENSG00000124762", "ENSG00000141682",
              "ENSG00000051108", "ENSG00000133639", "ENSG00000141232", 
              "ENSG00000172216"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expect_equal(result, expected)
})


### Tests filterResults() results

context("filterResults() results")

test_that("filterResults() must return expected text", {
  
    results <- parentalNapaVsDMSOEnrichment$result
    
    expected <- results[results$source == "WP" & 
                                results$term_id != "WP:000000", ]
    row.names(expected) <- NULL
    
    selected <- enrichViewNet:::filterResults(gostResults=results, source="WP", 
                         termIDs=NULL, removeRoot=TRUE)
    row.names(selected) <- NULL
    
    expect_equal(selected, expected)
})


### Tests extractInformationWhenIntersection() results

context("extractInformationWhenIntersection() results")

test_that("extractInformationWhenIntersection() must return expected text", {
    
    set.seed(1212)
    
    ccDemo <- parentalNapaVsDMSOEnrichment
    
    ccData <- ccDemo$result[ccDemo$result$term_id %in% 
                                    c("WP:WP1742", "WP:WP410"), ]
    
    result <- enrichViewNet:::extractInformationWhenIntersection(
        gostResults=ccData)
    
    expected <- list()
    
    expected[["geneNodes"]] <- data.frame(
        "id"=c("ENSG00000105327", "ENSG00000124762", "ENSG00000141682",
                "ENSG00000051108", "ENSG00000133639", "ENSG00000141232", 
                "ENSG00000172216"),
        "group"=rep("GENE", 7),
        "alias"=c("ENSG00000105327", "ENSG00000124762", "ENSG00000141682",
                    "ENSG00000051108", "ENSG00000133639", "ENSG00000141232", 
                    "ENSG00000172216"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["termNodes"]] <- data.frame(
        "id"=c("WP:WP1742", "WP:WP410"),
        "group"=c("TERM", "TERM"),
        "alias"=c("TP53 network", "Exercise-induced circadian regulation"),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["edges"]] <- data.frame(
        "source"=c("WP:WP1742", "WP:WP1742", "WP:WP1742", "WP:WP410", 
                    "WP:WP410", "WP:WP410", "WP:WP410"),
        "target"=c("ENSG00000105327", "ENSG00000124762", "ENSG00000141682",
                    "ENSG00000051108", "ENSG00000133639", "ENSG00000141232", 
                    "ENSG00000172216"),
        "interaction"=rep("contains", 7),
        check.names=FALSE, stringsAsFactors=FALSE)
    
    expect_identical(result, expected)
})


### Tests extractInformationWhenNoIntersection() results

context("extractInformationWhenNoIntersection() results")

test_that("extractInformationWhenNoIntersection() must return expected text", {
    
    mirnaDemo <- demoGOST
    
    mirnaDemo$meta$query_metadata$queries[[1]] <- 
        mirnaDemo$meta$query_metadata$queries[[1]][1:4]
    
    mirnaData <- demoGOST$result[demoGOST$result$source == "MIRNA", ]
    
    set.seed(121)
    result <- enrichViewNet:::extractInformationWhenNoIntersection(
                    gostResults=mirnaData, gostObject=mirnaDemo)
    
    expected <- list()
    
    expected[["geneNodes"]] <- data.frame("id"=c("ENSG00000059728", 
                    "ENSG00000077616", "ENSG00000051108"),
                    group=rep("GENE", 3), 
                    alias=c("MXD1", "NAALAD2", "HERPUD1"),
                    check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["termNode"]] <- data.frame(id=c("MIRNA:hsa-miR-335-5p", 
                    "MIRNA:hsa-miR-3180-5p", "MIRNA:hsa-miR-759"),
                    target=c("ENSG00000059728", "ENSG00000077616", 
                        "ENSG00000051108", "ENSG00000059728", "ENSG00000051108", 
                        "ENSG00000059728"),
                    interaction=rep("contains", 6),
                    check.names=FALSE, stringsAsFactors=FALSE)
    
    expected[["edges"]] <- data.frame(
        "source"=c("MIRNA:hsa-miR-335-5p", "MIRNA:hsa-miR-335-5p", 
            "MIRNA:hsa-miR-335-5p", "MIRNA:hsa-miR-3180-5p", 
            "MIRNA:hsa-miR-759", "MIRNA:hsa-miR-759"),
        target=c("ENSG00000059728", "ENSG00000077616", "ENSG00000051108", 
                    "ENSG00000059728", "ENSG00000051108", "ENSG00000059728"),
        interaction=rep("contains", 6),
        check.names=FALSE, stringsAsFactors=FALSE)
        
    expect_equal(result$nodes, expected$nodes)
    expect_equal(result$edges, expected$edges)
    expect_equal(result$nodeAttributes, expected$nodeAttributes)
    expect_equal(result$edgeAttributes, expected$edgeAttributes)
})


### Tests createCXJSONForCytoscape() results

context("createCXJSONForCytoscape() results")

test_that("createCXJSONForCytoscape() must return expected text", {
    
    mirnaDemo <- demoGOST
    mirnaDemo$meta$query_metadata$queries[[1]] <- 
        mirnaDemo$meta$query_metadata$queries[[1]][1:2]
    mirnaData <- demoGOST$result[demoGOST$result$source == "MIRNA", ]
    
    set.seed(121)
    
    info <- enrichViewNet:::extractNodesAndEdgesInformation(gostResults=mirnaData, 
                                                gostObject=mirnaDemo)
    result <- enrichViewNet:::createCXJSONForCytoscape(
        nodeEdgeInfo=info, title = "MIRNA")
    
    expected <- paste0(
        "[{\"metaData\":[{\"name\":\"nodes\",\"version\":\"1.0\"},",
        "{\"name\":\"edges\",\"version\":\"1.0\"},{\"name\":\"edgeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"nodeAttributes\",\"version\":\"1.0\"},{\"name\":\"cyHiddenAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyNetworkRelations\",\"version\":\"1.0\"},{\"name\":\"cyGroups\",\"version\":\"1.0\"},",
        "{\"name\":\"networkAttributes\",\"version\":\"1.0\"},{\"name\":\"cyTableColumn\",\"version\":\"1.0\"},",
        "{\"name\":\"cySubNetworks\",\"version\":\"1.0\"}]},{\"networkAttributes\":[{\"n\":\"name\",\"v\":\"MIRNA\"}]},",
        "{\"nodes\":[{\"@id\":1,\"n\":\"ENSG00000051108\"},{\"@id\":2,\"n\":\"MIRNA:hsa-miR-335-5p\"},{\"@id\":3,\"n\":\"MIRNA:hsa-miR-759\"}]},",
        "{\"edges\":[{\"@id\":4,\"s\":2,\"t\":1,\"i\":\"contains\"},{\"@id\":5,\"s\":3,\"t\":1,\"i\":\"contains\"}]},",
        "{\"nodeAttributes\":[{\"po\":1,\"n\":\"alias\",\"v\":\"HERPUD1\"},{\"po\":1,\"n\":\"group\",\"v\":\"GENE\"},",
        "{\"po\":2,\"n\":\"alias\",\"v\":\"hsa-miR-335-5p\"},{\"po\":3,\"n\":\"alias\",\"v\":\"hsa-miR-759\"},",
        "{\"po\":2,\"n\":\"group\",\"v\":\"TERM\"},{\"po\":3,\"n\":\"group\",\"v\":\"TERM\"}]},",
        "{\"edgeAttributes\":[{\"po\":4,\"n\":\"name\",\"v\":\"MIRNA:hsa-miR-335-5p (contains) ENSG00000051108\"},",
        "{\"po\":5,\"n\":\"name\",\"v\":\"MIRNA:hsa-miR-759 (contains) ENSG00000051108\"},",
        "{\"po\":4,\"n\":\"source\",\"v\":\"MIRNA:hsa-miR-335-5p\"},{\"po\":5,\"n\":\"source\",\"v\":\"MIRNA:hsa-miR-759\"},",
        "{\"po\":4,\"n\":\"target\",\"v\":\"ENSG00000051108\"},{\"po\":5,\"n\":\"target\",\"v\":\"ENSG00000051108\"}]},",
        "{\"cyHiddenAttributes\":[{\"n\":\"layoutAlgorithm\",\"y\":\"yFiles Circular Layout\"}]},",
        "{\"metaData\":[{\"name\":\"nodes\",\"version\":\"1.0\"},{\"name\":\"edges\",\"version\":\"1.0\"},",
        "{\"name\":\"edgeAttributes\",\"version\":\"1.0\"},{\"name\":\"nodeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyHiddenAttributes\",\"version\":\"1.0\"},{\"name\":\"cyNetworkRelations\",\"version\":\"1.0\"},",
        "{\"name\":\"cyGroups\",\"version\":\"1.0\"},{\"name\":\"networkAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyTableColumn\",\"version\":\"1.0\"},",
        "{\"name\":\"cySubNetworks\",\"version\":\"1.0\"}]},{\"status\":[{\"error\":\"\",\"success\":true}]}]")
    
    expect_equal(result, expected)
})


test_that("createCXJSONForCytoscape() must return expected text when intersection column present", {
    
    demo <- parentalNapaVsDMSOEnrichment
    demoData <- demo$result[demo$result$source == "REAC", ][c(11, 18, 20),]
    
    set.seed(121)
    
    info <- enrichViewNet:::extractNodesAndEdgesInformation(gostResults=demoData, 
                                                    gostObject=demo)
    
    result <- enrichViewNet:::createCXJSONForCytoscape(
        nodeEdgeInfo = info, title = "DEMO")
    
    expected <- paste0("[{\"metaData\":[{\"name\":\"nodes\",\"version\":\"1.0\"},",
        "{\"name\":\"edges\",\"version\":\"1.0\"},{\"name\":\"edgeAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"nodeAttributes\",\"version\":\"1.0\"},{\"name\":\"cyHiddenAttributes\",\"version\":\"1.0\"},",
        "{\"name\":\"cyNetworkRelations\",\"version\":\"1.0\"},{\"name\":\"cyGroups\",\"version\":\"1.0\"},",
        "{\"name\":\"networkAttributes\",\"version\":\"1.0\"},{\"name\":\"cyTableColumn\",\"version\":\"1.0\"},",
        "{\"name\":\"cySubNetworks\",\"version\":\"1.0\"}]},{\"networkAttributes\":[{\"n\":\"name\",\"v\":\"DEMO\"}]},",
        "{\"nodes\":[{\"@id\":1,\"n\":\"ENSG00000105327\"},{\"@id\":2,\"n\":\"ENSG00000113916\"},",
        "{\"@id\":3,\"n\":\"ENSG00000153094\"},{\"@id\":4,\"n\":\"ENSG00000175197\"},{\"@id\":5,\"n\":\"ENSG00000100292\"},",
        "{\"@id\":6,\"n\":\"ENSG00000124762\"},{\"@id\":7,\"n\":\"ENSG00000170345\"},{\"@id\":8,\"n\":\"ENSG00000171223\"},",
        "{\"@id\":9,\"n\":\"ENSG00000184557\"},{\"@id\":10,\"n\":\"ENSG00000141682\"},{\"@id\":11,\"n\":\"REAC:R-HSA-9614657\"},",
        "{\"@id\":12,\"n\":\"REAC:R-HSA-6785807\"},{\"@id\":13,\"n\":\"REAC:R-HSA-111453\"}]},",
        "{\"edges\":[{\"@id\":14,\"s\":13,\"t\":5,\"i\":\"contains\"},{\"@id\":15,\"s\":13,\"t\":1,\"i\":\"contains\"},",
        "{\"@id\":16,\"s\":13,\"t\":1,\"i\":\"contains\"},{\"@id\":17,\"s\":12,\"t\":2,\"i\":\"contains\"},",
        "{\"@id\":18,\"s\":12,\"t\":2,\"i\":\"contains\"},{\"@id\":19,\"s\":12,\"t\":6,\"i\":\"contains\"},",
        "{\"@id\":20,\"s\":12,\"t\":10,\"i\":\"contains\"},{\"@id\":21,\"s\":12,\"t\":3,\"i\":\"contains\"},",
        "{\"@id\":22,\"s\":12,\"t\":3,\"i\":\"contains\"},{\"@id\":23,\"s\":11,\"t\":7,\"i\":\"contains\"},",
        "{\"@id\":24,\"s\":11,\"t\":8,\"i\":\"contains\"},{\"@id\":25,\"s\":11,\"t\":4,\"i\":\"contains\"},",
        "{\"@id\":26,\"s\":11,\"t\":9,\"i\":\"contains\"}]},",
        "{\"nodeAttributes\":[{\"po\":1,\"n\":\"alias\",\"v\":\"ENSG00000105327\"},{\"po\":2,\"n\":\"alias\",\"v\":\"ENSG00000113916\"},",
        "{\"po\":3,\"n\":\"alias\",\"v\":\"ENSG00000153094\"},{\"po\":4,\"n\":\"alias\",\"v\":\"ENSG00000175197\"},",
        "{\"po\":5,\"n\":\"alias\",\"v\":\"ENSG00000100292\"},{\"po\":6,\"n\":\"alias\",\"v\":\"ENSG00000124762\"},",
        "{\"po\":7,\"n\":\"alias\",\"v\":\"ENSG00000170345\"},{\"po\":8,\"n\":\"alias\",\"v\":\"ENSG00000171223\"},",
        "{\"po\":9,\"n\":\"alias\",\"v\":\"ENSG00000184557\"},{\"po\":10,\"n\":\"alias\",\"v\":\"ENSG00000141682\"},",
        "{\"po\":1,\"n\":\"group\",\"v\":\"GENE\"},{\"po\":2,\"n\":\"group\",\"v\":\"GENE\"},",
        "{\"po\":3,\"n\":\"group\",\"v\":\"GENE\"},{\"po\":4,\"n\":\"group\",\"v\":\"GENE\"},",
        "{\"po\":5,\"n\":\"group\",\"v\":\"GENE\"},{\"po\":6,\"n\":\"group\",\"v\":\"GENE\"},",
        "{\"po\":7,\"n\":\"group\",\"v\":\"GENE\"},{\"po\":8,\"n\":\"group\",\"v\":\"GENE\"},",
        "{\"po\":9,\"n\":\"group\",\"v\":\"GENE\"},{\"po\":10,\"n\":\"group\",\"v\":\"GENE\"},",
        "{\"po\":11,\"n\":\"alias\",\"v\":\"FOXO-mediated transcription of cell death genes\"},", 
        "{\"po\":12,\"n\":\"alias\",\"v\":\"Interleukin-4 and Interleukin-13 signaling\"},",
        "{\"po\":13,\"n\":\"alias\",\"v\":\"BH3-only proteins associate with and inactivate anti-apoptotic BCL-2 members\"},",
        "{\"po\":11,\"n\":\"group\",\"v\":\"TERM\"},",
        "{\"po\":12,\"n\":\"group\",\"v\":\"TERM\"},{\"po\":13,\"n\":\"group\",\"v\":\"TERM\"}]},",
        "{\"edgeAttributes\":[{\"po\":14,\"n\":\"name\",\"v\":\"REAC:R-HSA-9614657 (contains) ENSG00000105327\"},",
        "{\"po\":15,\"n\":\"name\",\"v\":\"REAC:R-HSA-9614657 (contains) ENSG00000113916\"},",
        "{\"po\":16,\"n\":\"name\",\"v\":\"REAC:R-HSA-9614657 (contains) ENSG00000153094\"},",
        "{\"po\":17,\"n\":\"name\",\"v\":\"REAC:R-HSA-9614657 (contains) ENSG00000175197\"},",
        "{\"po\":18,\"n\":\"name\",\"v\":\"REAC:R-HSA-6785807 (contains) ENSG00000113916\"},",
        "{\"po\":19,\"n\":\"name\",\"v\":\"REAC:R-HSA-6785807 (contains) ENSG00000100292\"},",
        "{\"po\":20,\"n\":\"name\",\"v\":\"REAC:R-HSA-6785807 (contains) ENSG00000124762\"},",
        "{\"po\":21,\"n\":\"name\",\"v\":\"REAC:R-HSA-6785807 (contains) ENSG00000170345\"},",
        "{\"po\":22,\"n\":\"name\",\"v\":\"REAC:R-HSA-6785807 (contains) ENSG00000171223\"},",
        "{\"po\":23,\"n\":\"name\",\"v\":\"REAC:R-HSA-6785807 (contains) ENSG00000184557\"},",
        "{\"po\":24,\"n\":\"name\",\"v\":\"REAC:R-HSA-111453 (contains) ENSG00000105327\"},",
        "{\"po\":25,\"n\":\"name\",\"v\":\"REAC:R-HSA-111453 (contains) ENSG00000153094\"},",
        "{\"po\":26,\"n\":\"name\",\"v\":\"REAC:R-HSA-111453 (contains) ENSG00000141682\"},",
        "{\"po\":14,\"n\":\"source\",\"v\":\"REAC:R-HSA-9614657\"},{\"po\":15,\"n\":\"source\",\"v\":\"REAC:R-HSA-9614657\"},",
        "{\"po\":16,\"n\":\"source\",\"v\":\"REAC:R-HSA-9614657\"},{\"po\":17,\"n\":\"source\",\"v\":\"REAC:R-HSA-9614657\"},",
        "{\"po\":18,\"n\":\"source\",\"v\":\"REAC:R-HSA-6785807\"},{\"po\":19,\"n\":\"source\",\"v\":\"REAC:R-HSA-6785807\"},",
        "{\"po\":20,\"n\":\"source\",\"v\":\"REAC:R-HSA-6785807\"},{\"po\":21,\"n\":\"source\",\"v\":\"REAC:R-HSA-6785807\"},",
        "{\"po\":22,\"n\":\"source\",\"v\":\"REAC:R-HSA-6785807\"},{\"po\":23,\"n\":\"source\",\"v\":\"REAC:R-HSA-6785807\"},",
        "{\"po\":24,\"n\":\"source\",\"v\":\"REAC:R-HSA-111453\"},{\"po\":25,\"n\":\"source\",\"v\":\"REAC:R-HSA-111453\"},",
        "{\"po\":26,\"n\":\"source\",\"v\":\"REAC:R-HSA-111453\"},{\"po\":14,\"n\":\"target\",\"v\":\"ENSG00000105327\"},",
        "{\"po\":15,\"n\":\"target\",\"v\":\"ENSG00000113916\"},{\"po\":16,\"n\":\"target\",\"v\":\"ENSG00000153094\"},",
        "{\"po\":17,\"n\":\"target\",\"v\":\"ENSG00000175197\"},{\"po\":18,\"n\":\"target\",\"v\":\"ENSG00000113916\"},",
        "{\"po\":19,\"n\":\"target\",\"v\":\"ENSG00000100292\"},{\"po\":20,\"n\":\"target\",\"v\":\"ENSG00000124762\"},",
        "{\"po\":21,\"n\":\"target\",\"v\":\"ENSG00000170345\"},{\"po\":22,\"n\":\"target\",\"v\":\"ENSG00000171223\"},",
        "{\"po\":23,\"n\":\"target\",\"v\":\"ENSG00000184557\"},{\"po\":24,\"n\":\"target\",\"v\":\"ENSG00000105327\"},",
        "{\"po\":25,\"n\":\"target\",\"v\":\"ENSG00000153094\"},{\"po\":26,\"n\":\"target\",\"v\":\"ENSG00000141682\"}]},",
        "{\"cyHiddenAttributes\":[{\"n\":\"layoutAlgorithm\",\"y\":\"yFiles Circular Layout\"}]},",
                       "{\"metaData\":[{\"name\":\"nodes\",\"version\":\"1.0\"},",
                       "{\"name\":\"edges\",\"version\":\"1.0\"},{\"name\":\"edgeAttributes\",\"version\":\"1.0\"},",
                       "{\"name\":\"nodeAttributes\",\"version\":\"1.0\"},",
                       "{\"name\":\"cyHiddenAttributes\",\"version\":\"1.0\"},{\"name\":\"cyNetworkRelations\",\"version\":\"1.0\"},",
                       "{\"name\":\"cyGroups\",\"version\":\"1.0\"},{\"name\":\"networkAttributes\",\"version\":\"1.0\"},",
                       "{\"name\":\"cyTableColumn\",\"version\":\"1.0\"},{\"name\":\"cySubNetworks\",\"version\":\"1.0\"}]},",
                       "{\"status\":[{\"error\":\"\",\"success\":true}]}]")
    
    expect_equal(result, expected)
})

