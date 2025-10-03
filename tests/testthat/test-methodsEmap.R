### Unit tests for methodsEmap.R functions

library(enrichViewNet)
library(igraph)

data(demoGOST)
data(parentalNapaVsDMSOEnrichment)
data(rosaNapaVsDMSOEnrichment)

### Tests createEnrichMap() results

context("createEnrichMap() results")

test_that("createEnrichMap() must return error when gostObject is a number", {
    
    error_message <- paste0("The gostObject object should be a list ", 
                        "with meta and result as entries corresponding ", 
                        "to gprofiler2 enrichment output.")
    
    expect_error(createEnrichMap(gostObject=33, query="TEST", 
        source="GO:CC", termIDs=NULL, removeRoot=TRUE,
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1), error_message)
})


test_that("createEnrichMap() must return error when gostObject is a string character", {
    
    error_message <- paste0("The gostObject object should be a list ", 
                        "with meta and result as entries corresponding ", 
                        "to gprofiler2 enrichment output.")
    
    expect_error(createEnrichMap(gostObject="TEST", query="TEST", 
        source="GO:CC", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1), error_message)
})


test_that("createEnrichMap() must return error when query is a number", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("The \'query\'must be a character string.")
    
    expect_error(createEnrichMap(gostObject=gostObject, query=33, 
        source="KEGG", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1), error_message)
})


test_that("createEnrichMap() must return error when query is a vector of strings", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("The \'query\'must be a character string.")
    
    expect_error(createEnrichMap(gostObject=gostObject, query=c("1", "2"), 
        source="KEGG", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1), error_message)
})


test_that("createEnrichMap() must return error when query is not in gost", {
  
    error_message <- paste0("The \'query\' is not present in the ", 
                                    "results of the gost object.")
    
    expect_error(createEnrichMap(gostObject=demoGOST, query="CANADA", 
        source="KEGG", termIDs=NULL, removeRoot=TRUE,
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1), error_message)
})


test_that("createEnrichMap() must return error when source is a number", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    error_message <- paste0("Assertion on 'arg' failed: Must be of type ", 
                                "'character', not 'double'.")
    
    expect_error(createEnrichMap(gostObject=gostObject, query="toto", 
        source=333, termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1),  error_message)
})

test_that("createEnrichMap() must return error when source is a wrong name", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    expect_error(createEnrichMap(gostObject=gostObject, query="toto",  
        source="test", termIDs=NULL, removeRoot=TRUE, title="network", 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1))
})


test_that("createEnrichMap() must return error when source is GO", {
    
    gostObject <- list()
    gostObject[["meta"]] <- list()
    gostObject[["result"]] <- list()
    
    expect_error(createEnrichMap(gostObject=gostObject, query="toto", 
        source="GO",
        termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1))
})


test_that("createEnrichMap() must return error when removeRoot remove last enriched term", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- paste0("With removal of the root term, there is no ", 
                                "enrichment term left")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, 
                    showCategory=30, groupCategory=FALSE, categoryLabel=1,
                    categoryNode=1, line=1), error_message)
})


test_that("createEnrichMap() must return error when removeRoot remove last enriched term from term list", {
    
    gostTerm <- demoGOST
    gostTerm$result <- demoGOST$result[54,]
    
    error_message <- paste0("With removal of the root term, there is no ", 
                                "enrichment term left")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="TERM_ID",
        termIDs=c("WP:000000"), removeRoot=TRUE, showCategory=30, 
        groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1), error_message)
})


test_that("createEnrichMap() must return error when showCategory negative value", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'showCategory\' parameter must an positive ", 
            "integer or a vector of character strings representing terms.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=-30, 
                    groupCategory=FALSE, categoryLabel=1,
                    categoryNode=1, line=2), error_message)
})


test_that("createEnrichMap() must return error when showCategory is boolean", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'showCategory\' parameter must an positive ", 
            "integer or a vector of character strings representing terms.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=TRUE, 
                    groupCategory=FALSE, categoryLabel=1,
                    categoryNode=1, line=1), error_message)
})


test_that("createEnrichMap() must return error when groupCategory is integer", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'groupCategory\' parameter must a logical ", 
                                "(TRUE or FALSE).")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
            source="WP", removeRoot=TRUE, showCategory=30, 
            groupCategory=22, categoryLabel=1,
            categoryNode=1, line=1), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when categoryLabel is string", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'categoryLabel\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
            source="WP", removeRoot=TRUE, showCategory=30, 
            groupCategory=FALSE, categoryLabel="test",
            categoryNode=1, line=2), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when cexLabelCategory is negative", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'categoryLabel\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
            source="WP", removeRoot=TRUE, showCategory=30, 
            groupCategory=FALSE, categoryLabel=-1.1,
            categoryNode=1, line=2), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when categoryNode is negative", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'categoryNode\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
                    source="WP", removeRoot=TRUE, showCategory=30, 
                    groupCategory=FALSE, categoryLabel=2,
                    categoryNode=-1), error_message, fixed=TRUE)
})


test_that("createEnrichMap() must return error when categoryNode is string", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'categoryNode\' parameter ", 
                                "must be a positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
            source="WP", removeRoot=TRUE, showCategory=30, 
            groupCategory=FALSE, categoryLabel=2,
            categoryNode="te", line=1), error_message, fixed=TRUE)
})

test_that("createEnrichMap() must return error when line is a string", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'line\' parameter must be a ", 
                                "positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="WP", removeRoot=TRUE,  showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1, line="HI"), 
        error_message, fixed=TRUE)
})

test_that("createEnrichMap() must return error when line is a negative number", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The \'line\' parameter must be a ", 
                                "positive numeric.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="WP", removeRoot=TRUE,  showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1, line=-0.3), 
        error_message, fixed=TRUE)
})

test_that("createEnrichMap() must return error when not term for selected source", {
    
    gostTerm <- demoGOST
    gostTerm$result <- gostTerm$result[gostTerm$result$source != "WP", ]
    
    error_message <- paste0("There is no enriched term for the selected ", 
                                "source \'WP'.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="WP", removeRoot=TRUE,  showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1, line=1), error_message, 
                 fixed=TRUE)
})

test_that("createEnrichMap() must return error when not term id and TERM_ID selected", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("A vector of terms should be given through the ",
                        "\'termIDs\' parameter when source is \'TERM_ID\'.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="TERM_ID", removeRoot=TRUE,  showCategory=30, 
        groupCategory=FALSE, categoryLabel=1, categoryNode=1, line=1), 
        error_message, fixed=TRUE)
})

test_that("createEnrichMap() must return error when not all term ids are present", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("Not all listed terms are present in the  ",
                                    "enrichment results.")
    
    expect_error(createEnrichMap(gostObject=gostTerm, query="query_1", 
        source="TERM_ID", termIDs = c("GO:0051173", "GO:0065004", "GO:1905898"), 
        removeRoot=TRUE,  showCategory=30, 
        groupCategory=FALSE, categoryLabel=1, categoryNode=1, line=1), 
        error_message, fixed=TRUE)
})


### Tests createEnrichMapMultiBasic() results

context("createEnrichMapMultiBasic() results")

test_that("createEnrichMapMultiBasic() must return error when gostObjectList is a number", {
    
    error_message <- paste0("The gostObjectList object should be a list ", 
        "of enrichment objects. At least 2 enrichment objects are required.")
    
    expect_error(createEnrichMapMultiBasic(gostObjectList=33, 
        queryList=c("TEST", "Test2"),  source="GO:CC", termIDs=NULL, 
        removeRoot=TRUE, showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1, line=1), error_message)
})

test_that("createEnrichMapMultiBasic() must return error when gostObjectList has only one element", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The gostObjectList object should be a list ", 
        "of enrichment objects. At least 2 enrichment objects are required.")
    
    expect_error(createEnrichMapMultiBasic(gostObjectList=list(gostTerm), 
        queryList=list("TEST"), source="GO:CC", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1), error_message)
})

test_that("createEnrichMapMultiBasic() must return error when queryList is a number", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The queryList object should be a list of query ", 
        "names. At least 2 query names are required. The number of query ", 
        "names should correspond to the number of enrichment objects.")
    
    expect_error(createEnrichMapMultiBasic(gostObjectList=list(gostTerm, gostTerm), 
        queryList=33, source="GO:CC", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1), error_message)
})

test_that("createEnrichMapMultiBasic() must return error when queryList is longer than gostObjectList", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The queryList object should be a list of query ", 
        "names. At least 2 query names are required. The number of query ", 
        "names should correspond to the number of enrichment objects.")
    
    expect_error(createEnrichMapMultiBasic(gostObjectList=list(gostTerm, gostTerm), 
    queryList=list("TEST", "TEST2", "TEST3"), source="GO:CC", termIDs=NULL, 
    removeRoot=TRUE, showCategory=30, groupCategory=FALSE, categoryLabel=1,
    categoryNode=1), error_message)
})

test_that("createEnrichMapMultiBasic() must return error when one query in queryList is not in gostObject", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("Each query name present in the ", 
        "\'queryList\' parameter must be present in the associated ", 
        "enrichment object.")
    
    expect_error(createEnrichMapMultiBasic(gostObjectList=list(gostTerm, gostTerm), 
        queryList=list("query_1", "TEST"), source="GO:CC", termIDs=NULL, 
        removeRoot=TRUE, showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1), error_message, fixed=TRUE)
})

test_that("createEnrichMapMultiBasic() must return error when one object in gostObjectList in a number", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The gostObjectList should only contain a list ", 
        "of enrichment results. Enrichment results are lists with meta ", 
        "and result as entries corresponding to gprofiler2 ", 
        "enrichment output.")
    
    expect_error(createEnrichMapMultiBasic(gostObjectList=list(gostTerm, 33), 
        queryList=list("query_1", "query_1"), source="GO:CC", termIDs=NULL, 
        removeRoot=TRUE, showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1), error_message, fixed=TRUE)
})

test_that("createEnrichMapMultiBasic() must return error when no enriched term for selected source", {
    
    gostTerm <- demoGOST
    gostTerm$result <- gostTerm$result[which(gostTerm$result$source != "WP"), ]
    
    error_message <- paste0("There is no enriched term for the selected ", 
                                "source \'WP\'.")   
    
    expect_error(createEnrichMapMultiBasic(gostObjectList=list(gostTerm, gostTerm), 
        queryList=list("query_1", "query_1"), source="WP", termIDs=NULL, 
        removeRoot=TRUE, showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1), error_message, fixed=TRUE)
})

test_that("createEnrichMapMultiBasic() must return error when number in queryList", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The queryList object should only contain a list ", 
                                "of query names in character strings.")
    
    expect_error(createEnrichMapMultiBasic(gostObjectList=list(gostTerm, gostTerm), 
        queryList=list("query_1", 33), source="WP", termIDs=NULL, 
        removeRoot=TRUE, showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1), error_message, fixed=TRUE)
})

test_that("createEnrichMapMultiBasic() must return error when number in queryList", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("A vector of terms should be given through the ",
                "\'termIDs\' parameter when source is \'TERM_ID\'.")
    
    expect_error(createEnrichMapMultiBasic(gostObjectList=list(gostTerm, gostTerm), 
        queryList=list("query_1", "query_1"), source="TERM_ID", termIDs=NULL, 
        removeRoot=TRUE, showCategory=30, groupCategory=FALSE, categoryLabel=1,
        categoryNode=1), error_message, fixed=TRUE)
})


### Tests createEnrichMapMultiComplex() results

context("createEnrichMapMultiComplex() results")

test_that("createEnrichMapMultiComplex() must return error when gostObjectList is a number", {
    
    error_message <- paste0("The gostObjectList object should be a list ", 
        "of enrichment objects. At least 2 enrichment objects are required.")
    
    queryDF <- data.frame(queryName=c("parental_napa_vs_DMSO", 
        "rosa_napa_vs_DMSO", "rosa_napa_vs_DMSO"), 
        source=c("GO:CC", "REAC", "GO:CC"), removeRoot=c(TRUE, TRUE, TRUE),
        termIDs=c("", "", ""), stringsAsFactors=FALSE)
    
    expect_error(createEnrichMapMultiComplex(gostObjectList=33, 
        queryInfo=queryDF,  showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1, line=1), error_message)
})

test_that("createEnrichMapMultiComplex() must return error when gostObjectList has only one element", {
    
    gostTerm <- demoGOST
    
    queryDF <- data.frame(queryName=c("query_1"), 
            source=c("GO:CC"), removeRoot=c(TRUE),
            termIDs=c(""), stringsAsFactors=FALSE)
    
    error_message <- paste0("The gostObjectList object should be a list ", 
        "of enrichment objects. At least 2 enrichment objects are required.")
    
    expect_error(createEnrichMapMultiComplex(gostObjectList=list(gostTerm), 
        queryInfo=queryDF, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message)
})

test_that("createEnrichMapMultiComplex() must return error when gostObjectList is a list of numbers", {
    
    error_message <- paste0("The gostObjectList should only contain a list ", 
        "of enrichment results. Enrichment results are lists with meta and", 
        " result as entries corresponding to gprofiler2 enrichment", 
        " output.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(3, 4), 
        queryInfo=33, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message)
})

test_that("createEnrichMapMultiComplex() must return error when queryInfo is a number", {
    
    gostTerm <- demoGOST
    
    error_message <- paste0("The queryInfo should a data.frame with ", 
            "those columns: queryName, source, removeRoot and termIDs.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=33, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message)
})

test_that("createEnrichMapMultiComplex() must return error when queryInfo shorter than gostObjectList", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1"), 
        source=c("KEGG"), removeRoot=c(TRUE), termIDs=c(""), 
        stringsAsFactors=FALSE)
        
    error_message <- paste0("The number of rows in queryInfo should ", 
        "correspond to the number of enrichment objects.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message)
})

test_that("createEnrichMapMultiComplex() must return error when source in queryInfo is not in list of sources", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1", "query_1"), 
        source=c("KEGG", "TEST"), removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The values in the \'source\' column of the \'queryInfo\' ", 
        "data frame should be one of those: \"TERM_ID\", \"GO:MF\", ", 
        "\"GO:CC\", \"GO:BP\", \"KEGG\", \"REAC\", \"TF\", ",
        "\"MIRNA\", \"HPA\", \"CORUM\", \"HP\", \"WP\".")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message)
})

test_that("createEnrichMapMultiComplex() must return error when source in queryInfo is number", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1", "query_1"), 
        source=c(33, 22), removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'source'\ column of the \'queryInfo\' ", 
        "data frame should be in a character string format.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message)
})

test_that("createEnrichMapMultiComplex() must return error when queryName in queryInfo is number", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c(22, 33), 
        source=c("KEGG", "REAC"), removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'queryName'\ column of the \'queryInfo\' ", 
                    "data frame should be in a character string format.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message)
})

test_that("createEnrichMapMultiComplex() must return error when termIDs in queryInfo is number", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1", "query_1"), 
        source=c("KEGG", "REAC"), removeRoot=c(TRUE, TRUE), termIDs=c(33, 22), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'termIDs'\ column of the \'queryInfo\' ", 
    "data frame should be in a character string format.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message)
})

test_that("createEnrichMapMultiComplex() must return error when removeRoot in queryInfo is number", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1", "query_1"), 
        source=c("KEGG", "REAC"), removeRoot=c(33, 22), termIDs=c("", ""), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'removeRoot'\ column of the \'queryInfo\' ", 
        "data frame should only contain logical values (TRUE or FALSE).")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message, fixed=TRUE)
})

test_that("createEnrichMapMultiComplex() must return error when query name in queryInfo not in enrichment object", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1", "query_2"), 
        source=c("KEGG", "REAC"), removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("Each query name present in the \'queryName'\ ", 
        "column of the \'queryInfo\' data frame must be present in the ", 
        "associated enrichment object.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message, fixed=TRUE)
})

test_that("createEnrichMapMultiComplex() must return error when TERM_ID in queryInfo but not term ids", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1", "query_1"), 
        source=c("KEGG", "TERM_ID"), removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("A string of terms should be present in the ",
        "\'termIDs\' column when source is \'TERM_ID\'.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message, fixed=TRUE)
})

test_that("createEnrichMapMultiComplex() must return error when termIDs column in queryInfo contains numbers", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1", "query_1"), 
        source=c("KEGG", "TERM_ID"), removeRoot=c(TRUE, TRUE), termIDs=c(33, 22), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'termIDs\' column of the 'queryInfo' data ", 
        "frame should be in a character string format.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message, 
        fixed=TRUE)
})

test_that("createEnrichMapMultiComplex() must return error when removeRoot column in queryInfo contains numbers", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1", "query_1"), 
        source=c("KEGG", "TERM_ID"), removeRoot=c(33, 22), termIDs=c("", ""), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'removeRoot\' column of the 'queryInfo' ", 
        "data frame should only contain logical values (TRUE or FALSE).")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message, 
        fixed=TRUE)
})

test_that("createEnrichMapMultiComplex() must return error when groupName column in queryInfo contains numbers", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1", "query_1"), 
        source=c("KEGG", "REAC"), removeRoot=c(TRUE, TRUE), 
        termIDs=c("", ""), groupName=c(22, 33),  
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'groupName\' column of the 'queryInfo' ", 
        "data frame should be in a character string format.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message, 
        fixed=TRUE)
})

test_that("createEnrichMapMultiComplex() must return error when groupName column contain identical names", {
    
    gostTerm <- demoGOST
    
    queryDataFrame <- data.frame(queryName=c("query_1", "query_1"), 
        source=c("KEGG", "REAC"), removeRoot=c(TRUE, TRUE), 
        termIDs=c("", ""), groupName=c("REAC", "REAC"),  
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'groupName\' column of the 'queryInfo' ", 
        "data frame should only contain unique group names.")
    
    expect_error(createEnrichMapMultiComplex(
        gostObjectList=list(gostTerm, gostTerm), 
        queryInfo=queryDataFrame, showCategory=30, groupCategory=FALSE, 
        categoryLabel=1, categoryNode=1), error_message, 
        fixed=TRUE)
})


### Tests createEnrichMapMultiComplexAsIgraph() results

context("createEnrichMapMultiComplexAsIgraph() results")

test_that("createEnrichMapMultiComplexAsIgraph() must return error when gostObjectList is numerical list", {
    
    gostObjL <- list(33, 22)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
                        "parental_napa_vs_DMSO"), source=c("KEGG", "REAC"), 
                        removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
                        groupName=c("parental - KEGG", "parental - Reactome"), 
                        stringsAsFactors=FALSE)

    error_message <- paste0("The gostObjectList should only contain a list ", 
        "of enrichment results. Enrichment results are lists with meta and ", 
        "result as entries corresponding to gprofiler2 enrichment output.")
    
    expect_error(createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjL, 
            queryInfo=queryData, showCategory=20, similarityCutOff=0.25), 
            error_message=error_message)
})

test_that("createEnrichMapMultiComplexAsIgraph() must return error when queryInfo missing the groupName field", {
    
    gostObjL <- list(parentalNapaVsDMSOEnrichment, 
                                    parentalNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
                        "parental_napa_vs_DMSO"), source=c("KEGG", "REAC"), 
                        removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
                        test=c("parental - KEGG", "parental - Reactome"), 
                        stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'groupName\' column of the \'queryInfo\' ", 
        "data frame should be in a character string format.")
    
    expect_error(createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjL, 
        queryInfo=queryData, showCategory=20, similarityCutOff=0.25), 
        error_message=error_message)
})

test_that("createEnrichMapMultiComplexAsIgraph() must return error when queryInfo contains more elements than the gostObjectList", {
    
    gostObjL <- list(parentalNapaVsDMSOEnrichment, 
                                parentalNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
        "parental_napa_vs_DMSO", "test"), source=c("KEGG", "REAC", "GO:BP"), 
        removeRoot=c(TRUE, TRUE, TRUE), termIDs=c("", "", ""), 
        groupName=c("parental - KEGG", "parental - Reactome", "toto"), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The number of rows in queryInfo should", 
                        " correspond to the number of enrichment objects.")
    
    expect_error(createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjL, 
        queryInfo=queryData, showCategory=20, similarityCutOff=0.25), 
        error_message=error_message)
})

test_that("createEnrichMapMultiComplexAsIgraph() must return error when showCategory is a string", {
    
    gostObjL <- list(parentalNapaVsDMSOEnrichment, 
                            parentalNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
            "parental_napa_vs_DMSO"), source=c("KEGG", "REAC"), 
            removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
            groupName=c("parental - KEGG", "parental - Reactome"), 
            stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'showCategory\' parameter must an ", 
                                    "positive integer or NULL.")
    
    expect_error(createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjL, 
        queryInfo=queryData, showCategory="TEST", similarityCutOff=0.25), 
        error_message=error_message)
})

test_that("createEnrichMapMultiComplexAsIgraph() must return error when showCategory is negative", {
    
    gostObjL <- list(parentalNapaVsDMSOEnrichment, 
                        parentalNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
                    "parental_napa_vs_DMSO"), source=c("KEGG", "REAC"), 
                    removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
                    groupName=c("parental - KEGG", "parental - Reactome"), 
                    stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'showCategory\' parameter must an ", 
                                    "positive integer or NULL.")
    
    expect_error(createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjL, 
        queryInfo=queryData, showCategory=-2L, similarityCutOff=0.25), 
        error_message=error_message)
})

test_that("createEnrichMapMultiComplexAsIgraph() must return error when similarityCutOff is a string", {
    
    gostObjL <- list(parentalNapaVsDMSOEnrichment, 
                                parentalNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
        "parental_napa_vs_DMSO"), source=c("KEGG", "REAC"), 
        removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
        groupName=c("parental - KEGG", "parental - Reactome"), 
        stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'similarityCutOff\' parameter must be a ", 
        "numeric superior to zero and inferior to one.")
    
    expect_error(createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjL, 
            queryInfo=queryData, showCategory=20L, similarityCutOff="BEST"), 
            error_message=error_message)
})

test_that("createEnrichMapMultiComplexAsIgraph() must return error when not term left", {
    
    gostObjL <- list(parentalNapaVsDMSOEnrichment, 
                        parentalNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
        "parental_napa_vs_DMSO"), source=c("TERM_ID", "TERM_ID"), 
        removeRoot=c(TRUE, TRUE), termIDs=c("WP:000000", "KEGG:00000"), 
        groupName=c("parental - KEGG", "parental - Reactome"), 
            stringsAsFactors=FALSE)
    
    error_message <- paste0("With removal of the root term, there is ", 
                                "no enrichment term left")
    
    expect_error(createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjL, 
        queryInfo=queryData, showCategory=20L, similarityCutOff=0.2), 
        error_message=error_message)
})

test_that("createEnrichMapMultiComplexAsIgraph() must return error when similarityCutOff is negative numeric", {
    
    gostObjL <- list(parentalNapaVsDMSOEnrichment, 
                            parentalNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
            "parental_napa_vs_DMSO"), source=c("KEGG", "REAC"), 
            removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
            groupName=c("parental - KEGG", "parental - Reactome"), 
            stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'similarityCutOff\' parameter must be a ", 
        "numeric superior to zero and inferior to one.")
    
    expect_error(createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjL, 
        queryInfo=queryData, showCategory=20L, similarityCutOff=-0.1), 
        error_message=error_message)
})

test_that("createEnrichMapMultiComplexAsIgraph() must return error when similarityCutOff is superior to 1", {
    
    gostObjL <- list(parentalNapaVsDMSOEnrichment, 
                            parentalNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
            "parental_napa_vs_DMSO"), source=c("KEGG", "REAC"), 
            removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
            groupName=c("parental - KEGG", "parental - Reactome"), 
            stringsAsFactors=FALSE)
    
    error_message <- paste0("The \'similarityCutOff\' parameter must be a ", 
        "numeric superior to zero and inferior to one.")
    
    expect_error(createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjL, 
        queryInfo=queryData, showCategory=20L, similarityCutOff=1.1), 
                 error_message=error_message)
})

test_that("createEnrichMapMultiComplexAsIgraph() must return expected result", {
    
    gostObjectList <- list(parentalNapaVsDMSOEnrichment, 
                           parentalNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
                    "parental_napa_vs_DMSO"), source=c("KEGG", "REAC"), 
                    removeRoot=c(TRUE, TRUE), termIDs=c("", ""), 
                    groupName=c("parental - KEGG", "parental - Reactome"), 
                    stringsAsFactors=FALSE)
    
    result <- createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjectList, 
        queryInfo=queryData, showCategory=20, similarityCutOff=0.25)
    
    exp_size <- c(14, 8, 6, 7, 6, 6, 5, 7, 6, 4, 5, 9, 9, 6, 9, 7, 9, 34, 5, 
                    23, 5, 4, 23, 13, 17, 17, 5, 4, 6, 23, 3)
    exp_name <- c("name",   "size",   "pie",  "cluster", "pieName")
    
    expect_equal(class(result), "igraph")
    expect_equal(length(igraph::V(result)), 31)
    expect_equal(length(igraph::E(result)), 63)
    expect_equal(igraph::V(result)$size, exp_size)
    expect_equal(names(igraph::vertex.attributes(result)), exp_name)
    expect_equal(names(igraph::edge.attributes(result)), c("similarity", 
                                                                    "width"))
})

test_that("createEnrichMapMultiComplexAsIgraph() must return expected result", {
    
    gostObjectList <- list(parentalNapaVsDMSOEnrichment, 
                                rosaNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
                "rosa_napa_vs_DMSO"), source=c("TERM_ID", "TERM_ID"), 
                removeRoot=c(TRUE, TRUE), 
                termIDs=c(c("WP:WP4925,WP:WP3613,WP:WP382,WP:WP395"), 
                        c("WP:WP4925,WP:WP3613,WP:WP382,WP:WP3287")), 
                groupName=c("parental - WP", "rosa - WP"), 
                stringsAsFactors=FALSE)
    
    result <- createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjectList, 
            queryInfo=queryData, showCategory=20, similarityCutOff=0.3)
    
    exp_v_name <- c("Photodynamic therapy-induced unfolded protein response", 
            "Unfolded protein response", "MAPK signaling pathway", 
            "IL-4 signaling pathway", "Overview of nanoparticle effects")
    
    exp_size <- c(7, 6, 10, 4, 2)
    exp_name <- c("name",   "size",   "pie",  "cluster", "pieName")
    
    exp_pie <- list()
    exp_pie[[1]] <- c(1, 1)
    exp_pie[[2]] <- c(1, 1)
    exp_pie[[3]] <- c(1, 1)
    exp_pie[[4]] <- c(1, 0)
    exp_pie[[5]] <- c(0, 1)
    
    expect_equal(class(result), "igraph")
    expect_equal(length(V(result)), 5)
    expect_equal(V(result)$name, exp_v_name)
    expect_equal(V(result)$size, exp_size)
    expect_equal(V(result)$pie, exp_pie)
    expect_equal(length(E(result)), 1)
    expect_equal(E(result)$similarity, c(0.625))
    expect_equal(E(result)$width, c(0.625))
    expect_equal(names(vertex.attributes(result)), exp_name)
    expect_equal(names(edge.attributes(result)), c("similarity", "width"))
})

test_that("createEnrichMapMultiComplexAsIgraph() must return expected result when one node", {
    
    gostObjectList <- list(parentalNapaVsDMSOEnrichment, 
                           rosaNapaVsDMSOEnrichment)
    
    queryData <- data.frame(queryName=c("parental_napa_vs_DMSO", 
            "rosa_napa_vs_DMSO"), source=c("TERM_ID", "TERM_ID"), 
            removeRoot=c(TRUE, TRUE), 
            termIDs=c(c("WP:WP4925"), c("WP:WP4925")), 
            groupName=c("parental - WP", "rosa - WP"), 
            stringsAsFactors=FALSE)
    
    result <- createEnrichMapMultiComplexAsIgraph(gostObjectList=gostObjectList, 
                queryInfo=queryData, showCategory=20, similarityCutOff=0.3)
    
    exp_v_name <- c("Unfolded protein response")
    
    exp_size <- c(6)
    exp_name <- c("name",   "size")
    
    expect_equal(class(result), "igraph")
    expect_equal(length(V(result)), 1)
    expect_equal(V(result)$name, exp_v_name)
    expect_equal(V(result)$size, exp_size)
    expect_equal(length(E(result)),0)
    expect_equal(names(vertex.attributes(result)), exp_name)
})


### Tests createEnrichMapAsIgraph() results

context("createEnrichMapAsIgraph() results")

test_that("createEnrichMapAsIgraph() must return error when gostObjectList is numerical list", {
    
    gostObjL <- list(33, 22)
    
    queryData <- "parental_napa_vs_DMSO"
    
    error_message <- paste0("The gostObject object should be a list with ", 
        "meta and result as entries corresponding to gprofiler2 enrichment ", 
        "output.")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=queryData, showCategory=20, source="GO:CC", 
        similarityCutOff=0.25), error_message=error_message)
})

test_that("createEnrichMapAsIgraph() must return error when similarityCutOff is superior to 1", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "parental_napa_vs_DMSO"
    
    error_message <- paste0("The \'similarityCutOff\' parameter must be a ", 
        "numeric superior to zero and inferior to one.")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=queryData, showCategory=20L, source="GO:CC", removeRoot=TRUE, 
        similarityCutOff=1.1), error_message=error_message)
})

test_that("createEnrichMapAsIgraph() must return error when similarityCutOff is inferior to zero", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "parental_napa_vs_DMSO"
    
    error_message <- paste0("The \'similarityCutOff\' parameter must be a ", 
                        "numeric superior to zero and inferior to one.")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=queryData, showCategory=20L, source="GO:CC", removeRoot=TRUE, 
        similarityCutOff=-0.01), error_message=error_message)
})

test_that("createEnrichMapAsIgraph() must return error when removeRoot is a string", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "parental_napa_vs_DMSO"
    
    error_message <- paste0("The \'removeRoot\' parameter must a ", 
                        "logical (TRUE or FALSE).)")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=queryData, showCategory=20L, source="GO:CC", removeRoot="test", 
        similarityCutOff=0.1), error_message=error_message, fixed=TRUE)
})

test_that("createEnrichMapAsIgraph() must return error when query is a numeric", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    error_message <- paste0("The 'query' must be a character string")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=3, showCategory=20L, source="GO:CC", removeRoot=TRUE, 
        similarityCutOff=0.1), error_message=error_message, fixed=TRUE)
})

test_that("createEnrichMapAsIgraph() must return error when showCategory is a string", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "parental_napa_vs_DMSO"
    
    error_message <- paste0("The 'showCategory' parameter must an ", 
                                "positive integer or NULL.")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=queryData, showCategory="ROMEO", source="GO:CC", removeRoot=TRUE, 
        similarityCutOff=0.1), error_message=error_message, fixed=TRUE)
})

test_that("createEnrichMapAsIgraph() must return error when source has not result", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "parental_napa_vs_DMSO"
    
    error_message <- paste0(" There is no enriched term for the selected ", 
                                "source 'TF'.")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=queryData, showCategory=NULL, source="TF", removeRoot=TRUE, 
        similarityCutOff=0.1), error_message=error_message, fixed=TRUE)
})

test_that("createEnrichMapAsIgraph() must return error when not all listed terms exist", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "parental_napa_vs_DMSO"
    
    errorM <- paste0("Not all listed terms are present in the ", 
                        "enrichment results.")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=queryData, showCategory=NULL, source="TERM_ID",
        termIDs=c("WP:WP3613,REAC:R-HSA-9614085,TEST,GO:0008140"), 
        removeRoot=TRUE, similarityCutOff=0.1), error_message=errorM, 
        fixed=TRUE)
})

test_that("createEnrichMapAsIgraph() must return error when query not in enrichment object", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "ttee"
    
    errorM <- paste0("The 'query' is not present in the results of the", 
                        " gost object.")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=queryData, showCategory=NULL, source="TERM_ID",
        termIDs=c("WP:WP3613,REAC:R-HSA-9614085,GO:0008140"), 
        removeRoot=TRUE, similarityCutOff=0.1), error_message=errorM, 
        fixed=TRUE)
})

test_that("createEnrichMapAsIgraph() must return expected result", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "parental_napa_vs_DMSO"
    
    termIds <- c("WP:WP516", "WP:WP4970", "WP:WP4211", "WP:WP4216")
    
    result <- createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=queryData, showCategory=NULL, source="TERM_ID", 
        termIDs=termIds, similarityCutOff=0.3)
    
    exp_v_name <- c("Transcriptional cascade regulating adipogenesis", 
        "Chromosomal and microsatellite instability in colorectal cancer ", 
        "Galanin receptor pathway", "Hypertrophy model")
    
    exp_size <- c(5, 6, 3, 3)
    exp_name <- c("name",   "size")
    
    exp_pie <- list()
    exp_pie[[1]] <- c(1, 1)
    exp_pie[[2]] <- c(1, 1)
    exp_pie[[3]] <- c(1, 1)
    exp_pie[[4]] <- c(1, 0)
    exp_pie[[5]] <- c(0, 1)
    
    expect_equal(class(result), "igraph")
    expect_equal(length(V(result)), 4)
    expect_equal(V(result)$name, exp_v_name)
    expect_equal(V(result)$size, exp_size)
    expect_equal(length(E(result)), 1)
    expect_equal(E(result)$similarity, c(0.5))
    expect_equal(E(result)$width, c(0.5))
    expect_equal(names(vertex.attributes(result)), exp_name)
    expect_equal(names(edge.attributes(result)), c("similarity", "width"))
})

test_that("createEnrichMapAsIgraph() must return error when query not in enrichment object", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "ttee"
    
    errorM <- paste0("The 'query' is not present in the results of the", 
                        " gost object.")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
            query=queryData, showCategory=NULL, source="TERM_ID",
            termIDs=c("WP:WP3613,REAC:R-HSA-9614085,GO:0008140"), 
            removeRoot=TRUE, similarityCutOff=0.1), error_message=errorM, 
            fixed=TRUE)
})

test_that("createEnrichMapAsIgraph() must return error when query not in enrichment object", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "parental_napa_vs_DMSO"
    
    errorM <- paste0("With removal of the root term, there is no enrichment", 
                            " term left")
    
    expect_error(createEnrichMapAsIgraph(gostObject=gostObjL, 
        query=queryData, showCategory=NULL, source="TERM_ID",
        termIDs=c("WP:000000"), removeRoot=TRUE, similarityCutOff=0.1), 
        error_message=errorM, fixed=TRUE)
})

test_that("createEnrichMapAsIgraph() must return expected result when one entry", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "parental_napa_vs_DMSO"
    
    termIds <- c("WP:WP516")
    
    result <- createEnrichMapAsIgraph(gostObject=gostObjL, 
            query=queryData, showCategory=NULL, source="TERM_ID", 
            termIDs=termIds, similarityCutOff=0.3)
    
    exp_v_name <- c("Hypertrophy model")
    
    exp_size <- c(3)
    exp_name <- c("name",   "size")
    
    expect_equal(class(result), "igraph")
    expect_equal(length(V(result)), 1)
    expect_equal(V(result)$name, exp_v_name)
    expect_equal(V(result)$size, exp_size)
    expect_equal(length(E(result)), 0)
    expect_equal(names(vertex.attributes(result)), exp_name)
})

test_that("createEnrichMapAsIgraph() must return expected result when limit number of nodes", {
    
    gostObjL <- parentalNapaVsDMSOEnrichment
    
    queryData <- "parental_napa_vs_DMSO"
    
    result <- createEnrichMapAsIgraph(gostObject=gostObjL, 
            query=queryData, showCategory=6, source="WP", removeRoot=TRUE,
            termIDs=termIds, similarityCutOff=0.3)
    
    exp_v_name <- c("Photodynamic therapy-induced unfolded protein response",
        "Unfolded protein response" , "VEGFA-VEGFR2 signaling",
        "Transcriptional cascade regulating adipogenesis", 
        "White fat cell differentiation", "Pre-implantation embryo")
    
    exp_size <- c(7, 6, 16, 5, 6, 6)
    exp_name <- c("name",   "size")
    
    expect_equal(class(result), "igraph")
    expect_equal(length(V(result)), 6)
    expect_equal(V(result)$name, exp_v_name)
    expect_equal(V(result)$size, exp_size)
    expect_equal(length(E(result)), 2)
    expect_equal(E(result)$similarity, c(0.62500000000000, 0.8333333333333333))
    expect_equal(names(vertex.attributes(result)), exp_name)
})

