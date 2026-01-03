# Create a basic enrichment map

The function creates a basic enrichment map using functional enrichment
results.

## Usage

``` r
createBasicEmap(
  gostResults,
  backgroundGenes,
  showCategory,
  categoryLabel,
  categoryNode,
  significantMethod,
  line,
  ...
)
```

## Arguments

- gostResults:

  a `data.frame` containing the enrichment results to be plot.

- backgroundGenes:

  a `vector` of `character` string representing the name of the genes
  present in the request.

- showCategory:

  a positive `integer` or a `vector` of `characters` representing terms.
  If a `integer`, the first `n` terms will be displayed. If `vector` of
  terms, the selected terms will be displayed.

- categoryLabel:

  a positive `numeric` representing the amount by which plotting
  category nodes label size should be scaled relative to the default
  (1).

- categoryNode:

  a positive `numeric` representing the amount by which plotting
  category nodes should be scaled relative to the default (1).

- significantMethod:

  a `character` string representing the name of the multiple testing
  correction method used on the results.

- line:

  a non-negative `numeric` representing the scale of line width.

- ...:

  additional arguments that will be pass to the
  [`emapplot`](https://rdrr.io/pkg/enrichplot/man/emapplot.html)
  function.

## Value

a `ggplot` object representing the enrichment map.

## Author

Astrid Deschênes

## Examples

``` r

## Load the result of an enrichment analysis done with gprofiler2
data(parentalNapaVsDMSOEnrichment)

## Only retain the results section
gostResults <- as.data.frame(parentalNapaVsDMSOEnrichment$result)

## Limit the results to Wikipathways
## and remove the root term
gostResults <- gostResults[which(gostResults$source == "WP"),]
gostResults <- gostResults[which(gostResults$term_name != "WIKIPATHWAYS"),]

## Extract meta data information
meta <- parentalNapaVsDMSOEnrichment$meta

## Get all background genes
backgroundGenes <- meta$query_metadata$queries[["parental_napa_vs_DMSO"]]

## Get significant method
significantMethod <- meta$query_metadata$significance_threshold_method

## Create basic enrichment map using Wikipathways terms
enrichViewNet:::createBasicEmap(gostResults=gostResults, 
    backgroundGenes=backgroundGenes, showCategory=30L, 
    categoryLabel=1, categoryNode=1,
    significantMethod=significantMethod, line=1)

    
```
