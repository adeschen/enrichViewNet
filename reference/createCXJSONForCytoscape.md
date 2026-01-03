# Create CX JSON text representing the network

Create a CX JSON text that represent the network which includes
information about nodes and edges present in the network.

## Usage

``` r
createCXJSONForCytoscape(nodeEdgeInfo, title)
```

## Arguments

- nodeEdgeInfo:

  a `list` created by gprofiler2 that contains the results from an
  enrichment analysis. TODO

- title:

  a `character` string representing the name assigned to the network.

## Value

`character` string that represent the network in a CX JSON format.

## Author

Astrid Deschênes

## Examples

``` r

## Loading dataset containing result from an enrichment analysis done with
## gprofiler2
data(parentalNapaVsDMSOEnrichment)

## Only retained the GO Molecular Function results
results <- parentalNapaVsDMSOEnrichment$result[
        parentalNapaVsDMSOEnrichment$result$source == "GO:MF", ]

## Extract node and edge information
information <- enrichViewNet:::extractInformationWhenIntersection(
        gostResults=results)

jsonFormat <- enrichViewNet:::createCXJSONForCytoscape(
                nodeEdgeInfo=information, title="WikiPathways")
```
