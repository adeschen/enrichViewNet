# Create a basic enrichment map

The function creates a basic enrichment map using functional enrichment
results.

## Usage

``` r
createBasicEmapAsIgraph(
  gostResults,
  backgroundGenes,
  showCategory,
  similarityCutOff
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

- similarityCutOff:

  a positive `numeric` between 0 and 1 indicating the minimum level of
  similarity between two terms to have an edge linking the terms.

## Value

a `igraph` object representing the enrichment map.

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

## Create basic enrichment map, as an igraph, using Wikipathways terms
enrichViewNet:::createBasicEmapAsIgraph(gostResults=gostResults, 
    backgroundGenes=backgroundGenes, showCategory=30L, 
    similarityCutOff=0.2)
#> IGRAPH 5a08ef0 UNW- 26 53 -- 
#> + attr: name (v/c), size (v/n), similarity (e/n), width (e/n), weight
#> | (e/n)
#> + edges from 5a08ef0 (vertex names):
#> [1] Photodynamic therapy-induced unfolded protein response--Unfolded protein response                                       
#> [2] Photodynamic therapy-induced unfolded protein response--Nonalcoholic fatty liver disease                                
#> [3] Unfolded protein response                             --Chromosomal and microsatellite instability in colorectal cancer 
#> [4] Unfolded protein response                             --Apoptosis modulation and signaling                              
#> + ... omitted several edges
    
```
