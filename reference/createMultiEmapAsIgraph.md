# Create a complex enrichment map as an igraph

The function creates a complex enrichment map, as an igraph, using
functional enrichment results.

## Usage

``` r
createMultiEmapAsIgraph(
  gostResultsList,
  queryList,
  showCategory,
  similarityCutOff
)
```

## Arguments

- gostResultsList:

  a `list` of `data.frame` containing the enrichment results to be plot
  with different group identification.

- queryList:

  a `list` of `character` string representing the name of query retained
  for each enrichment results present in the `gostResultsList`
  parameter. The query should be present in its associated enrichment
  results.

- showCategory:

  a positive `integer` or a `vector` of `characters` representing terms.
  If a `integer`, the first `n` terms will be displayed. If `NULL`, all
  terms will be displayed.

- similarityCutOff:

  a positive `numeric`, larger than zero and small than 1 that represent
  the minimum similarity level between two nodes (terms) to be linked by
  an edge.

## Value

a `igraph` object representing the enrichment map with different colors
for each group of enrichment results.

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

## Extract meta data information
queryList <- list("parental - REACTOME", "parental - KEGG - v1", 
    "parental - KEGG - v2")

## Create basic enrichment map using Wikipathways terms
igraph <- enrichViewNet:::createMultiEmapAsIgraph(
    gostResultsList=list(gostResultsREAC, gostResultsKEGG, 
            gostResultsKEGG2), 
    queryList=queryList, showCategory=30L, similarityCutOff=0.5)
    
```
