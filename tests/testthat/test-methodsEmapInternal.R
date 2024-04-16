### Unit tests for methodsEmapInternal.R functions

library(enrichViewNet)
library(ggplot2)

data(parentalNapaVsDMSOEnrichment)
data(rosaNapaVsDMSOEnrichment)


### Tests validateCreateEnrichMapArguments() results

context("validateCreateEnrichMapArguments() results")

test_that("validateCreateEnrichMapArguments() must return expected result", {
    
    result <- enrichViewNet:::validateCreateEnrichMapArguments(
        gostObject=parentalNapaVsDMSOEnrichment, query="parental_napa_vs_DMSO", 
        source="GO:CC", termIDs=NULL, removeRoot=TRUE, 
        showCategory=30, groupCategory=FALSE, categoryLabel=1, categoryNode=1, 
        line=1, force=TRUE)
    
    expect_true(result)
})


### Tests validateCreateEnrichMapMultiArguments() results

context("validateCreateEnrichMapMultiArguments() results")

test_that("validateCreateEnrichMapMultiArguments() must return expected result", {
    
    result <- enrichViewNet:::validateCreateEnrichMapMultiArguments(
        gostObjectList=list(parentalNapaVsDMSOEnrichment, 
                                        rosaNapaVsDMSOEnrichment), 
        queryList=list("parental_napa_vs_DMSO", "rosa_napa_vs_DMSO"), 
        source="GO:CC", termIDs=NULL, removeRoot=TRUE, showCategory=30, 
        groupCategory=FALSE, categoryLabel=1, categoryNode=1, line=1, 
        force=TRUE)
    
    expect_true(result)
})


### Tests validateCreateEnrichMapMultiComplexArg() results

context("validateCreateEnrichMapMultiComplexArg() results")

test_that("validateCreateEnrichMapMultiComplexArg() must return expected result", {
    
    queryDF <- data.frame(queryName=c("parental_napa_vs_DMSO", 
        "rosa_napa_vs_DMSO", "rosa_napa_vs_DMSO"), 
        source=c("GO:CC", "REAC", "GO:CC"), removeRoot=c(TRUE, TRUE, TRUE),
        termIDs=c("", "", ""), stringsAsFactors=FALSE)
    
    result <- enrichViewNet:::validateCreateEnrichMapMultiComplexArg(
        gostObjectList=list(parentalNapaVsDMSOEnrichment, 
                rosaNapaVsDMSOEnrichment, rosaNapaVsDMSOEnrichment), 
        queryInfo=queryDF, showCategory=30, 
        groupCategory=FALSE, categoryLabel=1, categoryNode=1, line=1, 
        force=TRUE)
    
    expect_true(result)
})


### Tests validateCreateEnrichMapSubSectionArguments() results

context("validateCreateEnrichMapSubSectionArguments() results")

test_that("validateCreateEnrichMapSubSectionArguments() must return expected result", {
    
    result <- enrichViewNet:::validateCreateEnrichMapSubSectionArguments(
        showCategory=30, groupCategory=FALSE, categoryLabel=1, categoryNode=1,
        line=2, force=TRUE)
    
    expect_true(result)
})


### Tests manageNameDuplicationInEmap() results

context("manageNameDuplicationInEmap() results")

test_that("manageNameDuplicationInEmap() must return expected result", {
    
    clustData <- data.frame(Cluster=c("group 1" , "group 1", "group 2", 
                                        "group 2", "group 2", "group 1"), 
        ID=c("WP:WP4925", "WP:WP382", "KEGG:04010",  "KEGG:01010", 
                "KEGG:919191", "KEGG:101010"),
        Description=c("Unfolded protein response", 
            rep("MAPK signaling pathway", 2), "VEGFA-VEGFR2 signaling", 
            rep("Human T-cell pathway", 2)),
        GeneRatio=c("4/157", "3/157", "3/157", rep("2/157", 3)),
        BgRatio=c("4/24022", "3/24022", "3/24022", rep("2/24022", 3)),
        pvalues=c(1.55e-4, 8.13e-8, 4.33e-5, 3.2e-5, 3.1e-5, 3.5e-5),
        p.adjust=c(1e-3, 1e-3, 1.4e-3, 1e-3, 1e-3, 1e-3), 
        qvalue=c(1e-3, 1e-3, 1.4e-3, 1e-3, 1e-3, 1e-3), 
        geneID=c("ENSG000107968/ENSG000120129/ENSG000123358/ENSG000158050",
            "ENSG000107968/ENSG000120129/ENSG000158050",
            "ENSG000107968/ENSG000120129/ENSG000158050", 
            "ENSG000120129/ENSG000158050", "ENSG000120129/ENSG000158050",
            "ENSG000120129/ENSG000158050"),
        Count=c(4, 3, 3, 2, 2, 2))
    
    expected <- c("Unfolded protein response", 
        "MAPK signaling pathway (WP:WP382)", 
        "MAPK signaling pathway (KEGG:04010)", 
        "VEGFA-VEGFR2 signaling",  "Human T-cell pathway (KEGG:919191)", 
        "Human T-cell pathway (KEGG:101010)")
    
    result <- enrichViewNet:::manageNameDuplicationInEmap(clProfDF=clustData)

    expect_equal(result$Description, expected)
})


### Tests manageQueryDuplicationInEmap() results

context("manageQueryDuplicationInEmap() results")

test_that("manageQueryDuplicationInEmap() must return expected result", {
    
    queryList <- list("parental_vs_DMSO", "rosa_vs_DMSO", "parental_vs_DMSO", 
                    "rosa_vs_DMSO", "parental_vs_Control", "rosa_vs_DMSO", 
                    "Mia_vs_DMSO", "parental_vs_Control", "rosa_vs_DMSO")
    
    expected <- list("parental_vs_DMSO (1)", "rosa_vs_DMSO (1)", 
        "parental_vs_DMSO (2)", "rosa_vs_DMSO (2)", "parental_vs_Control (1)", 
        "rosa_vs_DMSO (3)", "Mia_vs_DMSO", "parental_vs_Control (2)", 
        "rosa_vs_DMSO (4)")
    
    result <- enrichViewNet:::manageQueryDuplicationInEmap(queryList=queryList)
    
    expect_equal(result, expected)
})


### Tests createBasicEmap() results

context("createBasicEmap() results")

test_that("createBasicEmap() must return expected result", {
    
    gostResults <- as.data.frame(parentalNapaVsDMSOEnrichment$result)
    gostResults <- gostResults[which(gostResults$source == "KEGG"),]
    gostResults <- gostResults[which(gostResults$term_id != "KEGG:00000"),]
    
    backgroundGenes <- parentalNapaVsDMSOEnrichment$meta$query_metadata$queries[["parental_napa_vs_DMSO"]]
    
    set.seed(111)
    
    graphRes <- enrichViewNet:::createBasicEmap(gostResults=gostResults, 
        backgroundGenes=backgroundGenes, showCategory=30L, 
        groupCategory=FALSE, categoryLabel=1, categoryNode=1,
        significantMethod="FDR", line=1, force=FALSE)
    
    expect_true(is.ggplot(graphRes))
    
    expect_true(all(graphRes$data$name %in% gostResults$term_name))
    
    ##expect_true(all(graphRes$data$size == gostResults$intersection_size))
})


### Tests createMultiEmap() results

context("createMultiEmap() results")

test_that("createMultiEmap() must return expected result when 2 different enrichment analyses", {
    
    gostResults1 <- as.data.frame(parentalNapaVsDMSOEnrichment$result)
    gostResults1 <- gostResults1[which(gostResults1$source == "KEGG"),]
    gostResults1 <- gostResults1[which(gostResults1$term_id != "KEGG:00000"),]
    
    gostResults2 <- as.data.frame(rosaNapaVsDMSOEnrichment$result)
    gostResults2 <- gostResults2[which(gostResults2$source == "KEGG"),]
    gostResults2 <- gostResults2[which(gostResults2$term_id != "KEGG:00000"),]
    
    queryList <- list("parental_napa_vs_DMSO", "rosa_napa_vs_DMSO")
    
    set.seed(111)
    
    graphRes <- enrichViewNet:::createMultiEmap(gostResults=list(gostResults1, 
        gostResults2), queryList=queryList, showCategory=30L, 
        groupCategory=FALSE, categoryLabel=1, categoryNode=1,
        line=1, force=FALSE)
    
    expect_true(is.ggplot(graphRes))
    
    expected_terms <- c("MAPK signaling pathway", 
        "Transcriptional misregulation in cancer", 
        "Parathyroid hormone synthesis, secretion and action",
        "Apoptosis", "Colorectal cancer", "TNF signaling pathway", 
        "p53 signaling pathway", "Human T-cell leukemia virus 1 infection", 
        "Non-alcoholic fatty liver disease", "Amphetamine addiction",
        "Estrogen signaling pathway", "Viral carcinogenesis",
        "IL-17 signaling pathway", "Fluid shear stress and atherosclerosis", 
        "Hepatitis B")
    
    expected_term_size <- c(14, 9, 6, 7, 6, 6, 5, 7, 7, 4, 5, 4, 3, 3, 3)
    
    expect_true(all(graphRes$data$name %in% expected_terms))
    
    ##expect_true(all(graphRes$data$size == expected_term_size))
    
    expect_identical(graphRes$labels$fill, "Cluster")
})

test_that("createMultiEmap() must return expected result when same different enrichment analysis", {
    
    gostResults1 <- as.data.frame(parentalNapaVsDMSOEnrichment$result)
    gostResults1 <- gostResults1[which(gostResults1$source == "KEGG"),]
    gostResults1 <- gostResults1[which(gostResults1$term_id != "KEGG:00000"),]
    
    gostResults2 <- as.data.frame(parentalNapaVsDMSOEnrichment$result)
    gostResults2 <- gostResults2[which(gostResults2$source == "REAC"),]
    gostResults2 <- gostResults2[1:10,]
    
    queryList <- list("parental_napa_vs_DMSO (1)", "parental_napa_vs_DMSO (2)")
    
    set.seed(111)
    
    graphRes <- enrichViewNet:::createMultiEmap(gostResults=list(gostResults1, 
        gostResults2), queryList=queryList, showCategory=30L, 
        groupCategory=FALSE, categoryLabel=1, categoryNode=1,
        line=1, force=FALSE)
    
    expect_true(is.ggplot(graphRes))
    
    expected_terms <- c("MAPK signaling pathway", 
        "Transcriptional misregulation in cancer", 
        "Parathyroid hormone synthesis, secretion and action",
        "Apoptosis", "Colorectal cancer", "TNF signaling pathway", 
        "p53 signaling pathway", "Human T-cell leukemia virus 1 infection", 
        "Non-alcoholic fatty liver disease", "Amphetamine addiction",
        "Estrogen signaling pathway", "NGF-stimulated transcription",
        "Nuclear Events (kinase and transcription factor activation)",
        "Response of EIF2AK1 (HRI) to heme deficiency",                    
        "Signaling by NTRK1 (TRKA)", "FOXO-mediated transcription",                                     
        "Signaling by NTRKs", "Signal Transduction",                                              
        "ATF4 activates genes in response to endoplasmic reticulum  stress", 
        "Generic Transcription Pathway", "PERK regulates gene expression" )
    
    expected_term_size <- c(14, 8, 6, 7, 6, 6, 5, 7, 6, 4, 5, 9, 9, 6, 9, 
                                7, 9, 34, 5, 23, 5)
    
    expect_true(all(graphRes$data$name %in% expected_terms))
    
    ##expect_true(all(graphRes$data$size == expected_term_size))
    
    expect_identical(graphRes$labels$fill, "Cluster")
})
    