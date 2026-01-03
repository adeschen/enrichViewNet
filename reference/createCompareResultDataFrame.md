# Create a data frame containing the compiled enrichment results from multiple queries

The function creates a data frame object using the list of queries and
list of enrichment results given be the user.

## Usage

``` r
createCompareResultDataFrame(queryList, gostResultsList)
```

## Arguments

- queryList:

  a `list` containing the query names.

- gostResultsList:

  a `list` of `data.frame` containing the enrichment results to be plot
  with different group identification.

## Value

a `data.frame` object

## Author

Astrid Deschênes

## Examples

``` r

## Load the result of an enrichment analysis done with gprofiler2
data(parentalNapaVsDMSOEnrichment)

## Only retain the results section
gostResults <- as.data.frame(parentalNapaVsDMSOEnrichment$result)

## Limit the results subsection of REACTOME and KEGG
## Kegg is replicated to show shared results between queries
gostResultsREAC <- gostResults[which(gostResults$source == "REAC"),]
gostResultsREAC <- gostResultsREAC[1:13, ]
gostResultsKEGG <- gostResults[which(gostResults$source == "KEGG"),]
gostResultsKEGG2 <- gostResultsKEGG[1:6,]

## List of query names with duplicated names
queryList <- list("parental_vs_DMSO", "rosa_vs_DMSO", "parental_vs_DMSO", 
    "rosa_vs_DMSO", "parental_vs_Control", "rosa_vs_DMSO")
    
## Extract meta data information
queryList <- list("parental - REACTOME", "parental - KEGG - v1", 
    "parental - KEGG - v2")

## Change the query names for the duplicated names
enrichViewNet:::createCompareResultDataFrame(queryList=queryList,
    gostResultsList=list(gostResultsREAC, gostResultsKEGG, 
            gostResultsKEGG2))
#>                 Cluster                 ID
#> 1   parental - REACTOME REAC:R-HSA-9031628
#> 2   parental - REACTOME  REAC:R-HSA-198725
#> 3   parental - REACTOME REAC:R-HSA-9648895
#> 4   parental - REACTOME  REAC:R-HSA-187037
#> 5   parental - REACTOME REAC:R-HSA-9614085
#> 6   parental - REACTOME  REAC:R-HSA-166520
#> 7   parental - REACTOME  REAC:R-HSA-162582
#> 8   parental - REACTOME  REAC:R-HSA-380994
#> 9   parental - REACTOME  REAC:R-HSA-212436
#> 10  parental - REACTOME  REAC:R-HSA-381042
#> 11  parental - REACTOME REAC:R-HSA-9614657
#> 12  parental - REACTOME   REAC:R-HSA-73857
#> 13  parental - REACTOME REAC:R-HSA-9006934
#> 14 parental - KEGG - v1         KEGG:04010
#> 15 parental - KEGG - v1         KEGG:00000
#> 16 parental - KEGG - v1         KEGG:05202
#> 17 parental - KEGG - v1         KEGG:04928
#> 18 parental - KEGG - v1         KEGG:04210
#> 19 parental - KEGG - v1         KEGG:05210
#> 20 parental - KEGG - v1         KEGG:04668
#> 21 parental - KEGG - v1         KEGG:04115
#> 22 parental - KEGG - v1         KEGG:05166
#> 23 parental - KEGG - v1         KEGG:04932
#> 24 parental - KEGG - v1         KEGG:05031
#> 25 parental - KEGG - v1         KEGG:04915
#> 26 parental - KEGG - v2         KEGG:04010
#> 27 parental - KEGG - v2         KEGG:00000
#> 28 parental - KEGG - v2         KEGG:05202
#> 29 parental - KEGG - v2         KEGG:04928
#> 30 parental - KEGG - v2         KEGG:04210
#> 31 parental - KEGG - v2         KEGG:05210
#>                                                          Description
#> 1                                       NGF-stimulated transcription
#> 2        Nuclear Events (kinase and transcription factor activation)
#> 3                       Response of EIF2AK1 (HRI) to heme deficiency
#> 4                                          Signaling by NTRK1 (TRKA)
#> 5                                        FOXO-mediated transcription
#> 6                                                 Signaling by NTRKs
#> 7                                                Signal Transduction
#> 8  ATF4 activates genes in response to endoplasmic reticulum  stress
#> 9                                      Generic Transcription Pathway
#> 10                                    PERK regulates gene expression
#> 11                   FOXO-mediated transcription of cell death genes
#> 12                                   RNA Polymerase II Transcription
#> 13                            Signaling by Receptor Tyrosine Kinases
#> 14                                            MAPK signaling pathway
#> 15                                                    KEGG root term
#> 16                           Transcriptional misregulation in cancer
#> 17               Parathyroid hormone synthesis, secretion and action
#> 18                                                         Apoptosis
#> 19                                                 Colorectal cancer
#> 20                                             TNF signaling pathway
#> 21                                             p53 signaling pathway
#> 22                           Human T-cell leukemia virus 1 infection
#> 23                                 Non-alcoholic fatty liver disease
#> 24                                             Amphetamine addiction
#> 25                                        Estrogen signaling pathway
#> 26                                            MAPK signaling pathway
#> 27                                                    KEGG root term
#> 28                           Transcriptional misregulation in cancer
#> 29               Parathyroid hormone synthesis, secretion and action
#> 30                                                         Apoptosis
#> 31                                                 Colorectal cancer
#>         pvalues
#> 1  4.150088e-10
#> 2  3.962325e-08
#> 3  1.328872e-07
#> 4  9.278149e-06
#> 5  1.832716e-05
#> 6  3.078300e-05
#> 7  7.557030e-05
#> 8  2.105375e-04
#> 9  2.651522e-04
#> 10 5.500090e-04
#> 11 1.183548e-03
#> 12 1.549768e-03
#> 13 1.713632e-03
#> 14 8.134007e-08
#> 15 2.503078e-06
#> 16 3.367430e-04
#> 17 1.289009e-03
#> 18 1.409759e-03
#> 19 1.588662e-03
#> 20 4.275641e-03
#> 21 6.548809e-03
#> 22 2.187420e-02
#> 23 3.043035e-02
#> 24 3.367925e-02
#> 25 3.560639e-02
#> 26 8.134007e-08
#> 27 2.503078e-06
#> 28 3.367430e-04
#> 29 1.289009e-03
#> 30 1.409759e-03
#> 31 1.588662e-03
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             geneID
#> 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ENSG00000120738/ENSG00000122877/ENSG00000125740/ENSG00000135625/ENSG00000170345/ENSG00000171223/ENSG00000173334/ENSG00000179388/ENSG00000198576
#> 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ENSG00000120738/ENSG00000122877/ENSG00000125740/ENSG00000135625/ENSG00000170345/ENSG00000171223/ENSG00000173334/ENSG00000179388/ENSG00000198576
#> 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ENSG00000087074/ENSG00000128272/ENSG00000128965/ENSG00000162772/ENSG00000172216/ENSG00000175197
#> 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ENSG00000120738/ENSG00000122877/ENSG00000125740/ENSG00000135625/ENSG00000170345/ENSG00000171223/ENSG00000173334/ENSG00000179388/ENSG00000198576
#> 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ENSG00000105327/ENSG00000113916/ENSG00000124762/ENSG00000133639/ENSG00000136826/ENSG00000153094/ENSG00000175197
#> 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ENSG00000120738/ENSG00000122877/ENSG00000125740/ENSG00000135625/ENSG00000170345/ENSG00000171223/ENSG00000173334/ENSG00000179388/ENSG00000198576
#> 7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ENSG00000007944/ENSG00000087074/ENSG00000104140/ENSG00000113070/ENSG00000116741/ENSG00000119630/ENSG00000120129/ENSG00000120738/ENSG00000122877/ENSG00000123358/ENSG00000124216/ENSG00000124762/ENSG00000125740/ENSG00000135625/ENSG00000138166/ENSG00000141582/ENSG00000143333/ENSG00000143878/ENSG00000147437/ENSG00000148926/ENSG00000150991/ENSG00000153094/ENSG00000158050/ENSG00000164056/ENSG00000170345/ENSG00000171223/ENSG00000172602/ENSG00000173334/ENSG00000173530/ENSG00000179388/ENSG00000183691/ENSG00000184545/ENSG00000184557/ENSG00000198576
#> 8                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ENSG00000051108/ENSG00000128272/ENSG00000162772/ENSG00000172216/ENSG00000175197
#> 9                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ENSG00000105327/ENSG00000113916/ENSG00000119508/ENSG00000123358/ENSG00000124762/ENSG00000130766/ENSG00000133639/ENSG00000136826/ENSG00000141582/ENSG00000141682/ENSG00000150991/ENSG00000152433/ENSG00000153094/ENSG00000159388/ENSG00000170345/ENSG00000171223/ENSG00000172216/ENSG00000173530/ENSG00000175197/ENSG00000177873/ENSG00000184557/ENSG00000197566/ENSG00000256223
#> 10                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000051108/ENSG00000128272/ENSG00000162772/ENSG00000172216/ENSG00000175197
#> 11                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000105327/ENSG00000113916/ENSG00000153094/ENSG00000175197
#> 12                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000105327/ENSG00000113916/ENSG00000119508/ENSG00000123358/ENSG00000124762/ENSG00000130766/ENSG00000133639/ENSG00000136826/ENSG00000141582/ENSG00000141682/ENSG00000150991/ENSG00000152433/ENSG00000153094/ENSG00000159388/ENSG00000170345/ENSG00000171223/ENSG00000172216/ENSG00000173530/ENSG00000175197/ENSG00000177873/ENSG00000184557/ENSG00000197566/ENSG00000256223
#> 13                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000113070/ENSG00000119630/ENSG00000120738/ENSG00000122877/ENSG00000125740/ENSG00000135625/ENSG00000150991/ENSG00000164056/ENSG00000170345/ENSG00000171223/ENSG00000173334/ENSG00000179388/ENSG00000198576
#> 14                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000099860/ENSG00000107968/ENSG00000119630/ENSG00000120129/ENSG00000123358/ENSG00000128272/ENSG00000138166/ENSG00000158050/ENSG00000169242/ENSG00000170345/ENSG00000175197/ENSG00000184545/ENSG00000204388/ENSG00000204389
#> 15 ENSG00000007944/ENSG00000051108/ENSG00000087074/ENSG00000099860/ENSG00000100292/ENSG00000102554/ENSG00000105327/ENSG00000107968/ENSG00000108551/ENSG00000113070/ENSG00000113916/ENSG00000116741/ENSG00000119508/ENSG00000119630/ENSG00000120129/ENSG00000120738/ENSG00000122877/ENSG00000123358/ENSG00000124216/ENSG00000124762/ENSG00000125266/ENSG00000125657/ENSG00000125740/ENSG00000128016/ENSG00000128272/ENSG00000128965/ENSG00000130066/ENSG00000130766/ENSG00000131471/ENSG00000133639/ENSG00000136826/ENSG00000138166/ENSG00000139112/ENSG00000141232/ENSG00000141682/ENSG00000142178/ENSG00000143878/ENSG00000144802/ENSG00000147437/ENSG00000148926/ENSG00000150991/ENSG00000152433/ENSG00000153094/ENSG00000157557/ENSG00000158050/ENSG00000159388/ENSG00000163273/ENSG00000169242/ENSG00000170345/ENSG00000171223/ENSG00000172216/ENSG00000172602/ENSG00000173530/ENSG00000175197/ENSG00000176383/ENSG00000177873/ENSG00000179388/ENSG00000183691/ENSG00000184545/ENSG00000184557/ENSG00000197566/ENSG00000198576/ENSG00000204103/ENSG00000204388/ENSG00000204389/ENSG00000245848/ENSG00000256223/ENSG00000257446/ENSG00000275215/ENSG00000278189/ENSG00000278233
#> 16                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000099860/ENSG00000113916/ENSG00000119508/ENSG00000124762/ENSG00000144802/ENSG00000172216/ENSG00000175197/ENSG00000245848
#> 17                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000113070/ENSG00000120738/ENSG00000124762/ENSG00000128272/ENSG00000170345/ENSG00000204103
#> 18                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000099860/ENSG00000105327/ENSG00000128272/ENSG00000141682/ENSG00000153094/ENSG00000170345/ENSG00000175197
#> 19                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000099860/ENSG00000105327/ENSG00000124762/ENSG00000141682/ENSG00000153094/ENSG00000170345
#> 20                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000107968/ENSG00000128272/ENSG00000170345/ENSG00000171223/ENSG00000172216/ENSG00000184557
#> 21                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000099860/ENSG00000105327/ENSG00000124762/ENSG00000130766/ENSG00000141682
#> 22                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000120738/ENSG00000122877/ENSG00000124762/ENSG00000128016/ENSG00000128272/ENSG00000157557/ENSG00000170345
#> 23                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000128272/ENSG00000153094/ENSG00000170345/ENSG00000175197/ENSG00000184557/ENSG00000245848
#> 24                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000125740/ENSG00000128272/ENSG00000170345/ENSG00000198576
#> 25                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000113070/ENSG00000128272/ENSG00000170345/ENSG00000204388/ENSG00000204389
#> 26                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000099860/ENSG00000107968/ENSG00000119630/ENSG00000120129/ENSG00000123358/ENSG00000128272/ENSG00000138166/ENSG00000158050/ENSG00000169242/ENSG00000170345/ENSG00000175197/ENSG00000184545/ENSG00000204388/ENSG00000204389
#> 27 ENSG00000007944/ENSG00000051108/ENSG00000087074/ENSG00000099860/ENSG00000100292/ENSG00000102554/ENSG00000105327/ENSG00000107968/ENSG00000108551/ENSG00000113070/ENSG00000113916/ENSG00000116741/ENSG00000119508/ENSG00000119630/ENSG00000120129/ENSG00000120738/ENSG00000122877/ENSG00000123358/ENSG00000124216/ENSG00000124762/ENSG00000125266/ENSG00000125657/ENSG00000125740/ENSG00000128016/ENSG00000128272/ENSG00000128965/ENSG00000130066/ENSG00000130766/ENSG00000131471/ENSG00000133639/ENSG00000136826/ENSG00000138166/ENSG00000139112/ENSG00000141232/ENSG00000141682/ENSG00000142178/ENSG00000143878/ENSG00000144802/ENSG00000147437/ENSG00000148926/ENSG00000150991/ENSG00000152433/ENSG00000153094/ENSG00000157557/ENSG00000158050/ENSG00000159388/ENSG00000163273/ENSG00000169242/ENSG00000170345/ENSG00000171223/ENSG00000172216/ENSG00000172602/ENSG00000173530/ENSG00000175197/ENSG00000176383/ENSG00000177873/ENSG00000179388/ENSG00000183691/ENSG00000184545/ENSG00000184557/ENSG00000197566/ENSG00000198576/ENSG00000204103/ENSG00000204388/ENSG00000204389/ENSG00000245848/ENSG00000256223/ENSG00000257446/ENSG00000275215/ENSG00000278189/ENSG00000278233
#> 28                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000099860/ENSG00000113916/ENSG00000119508/ENSG00000124762/ENSG00000144802/ENSG00000172216/ENSG00000175197/ENSG00000245848
#> 29                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000113070/ENSG00000120738/ENSG00000124762/ENSG00000128272/ENSG00000170345/ENSG00000204103
#> 30                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000099860/ENSG00000105327/ENSG00000128272/ENSG00000141682/ENSG00000153094/ENSG00000170345/ENSG00000175197
#> 31                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ENSG00000099860/ENSG00000105327/ENSG00000124762/ENSG00000141682/ENSG00000153094/ENSG00000170345
#>    Count
#> 1      9
#> 2      9
#> 3      6
#> 4      9
#> 5      7
#> 6      9
#> 7     34
#> 8      5
#> 9     23
#> 10     5
#> 11     4
#> 12    23
#> 13    13
#> 14    14
#> 15    71
#> 16     8
#> 17     6
#> 18     7
#> 19     6
#> 20     6
#> 21     5
#> 22     7
#> 23     6
#> 24     4
#> 25     5
#> 26    14
#> 27    71
#> 28     8
#> 29     6
#> 30     7
#> 31     6
    
```
