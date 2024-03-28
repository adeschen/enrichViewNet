### Unit tests for methodsEmapInternal.R functions

library(enrichViewNet)

data(parentalNapaVsDMSOEnrichment)
data(rosaNapaVsDMSOEnrichment)


### Tests validateCreateEnrichMapArguments() results

context("validateCreateEnrichMapArguments() results")

test_that("validateCreateEnrichMapArguments() must return expected result", {
    
    result <- validateCreateEnrichMapArguments(
        gostObject=parentalNapaVsDMSOEnrichment, query="parental_napa_vs_DMSO", 
        source="GO:CC", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1, categoryNode=1, 
        force=TRUE)
    
    expect_true(result)
})


### Tests validateCreateEnrichMapMultiArguments() results

context("validateCreateEnrichMapMultiArguments() results")

test_that("validateCreateEnrichMapMultiArguments() must return expected result", {
    
    result <- validateCreateEnrichMapMultiArguments(gostObjectList=list(
        parentalNapaVsDMSOEnrichment, rosaNapaVsDMSOEnrichment), 
        queryList=list("parental_napa_vs_DMSO", "rosa_napa_vs_DMSO"), 
        source="GO:CC", termIDs=NULL, removeRoot=TRUE, showCategory=30, 
        groupCategory=FALSE, categoryLabel=1, categoryNode=1, force=TRUE)
    
    expect_true(result)
})


### Tests validateCreateEnrichMapSubSectionArguments() results

context("validateCreateEnrichMapSubSectionArguments() results")


test_that("validateCreateEnrichMapSubSectionArguments() must return expected result", {
    
    result <- validateCreateEnrichMapSubSectionArguments(
        showCategory=30, groupCategory=FALSE, categoryLabel=1, categoryNode=1,
        force=TRUE)
    
    expect_true(result)
})