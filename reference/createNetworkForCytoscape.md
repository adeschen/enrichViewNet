# Create network and load it into Cytoscape

Create network from gprofiler2 results and load it into Cytoscape

## Usage

``` r
createNetworkForCytoscape(nodeEdgeInfo, title, collection)
```

## Arguments

- nodeEdgeInfo:

  a TODO

- title:

  a `character` string representing the name assigned to the network.

- collection:

  a `character` string representing the collection name assigned to the
  network.

## Value

`TRUE`

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
    
## The creation of the network can only be done when Cytoscape
## is up and running
## A network using GO - Molecular Function enriched terms will be
## generated and loaded into Cytoscape
if (enrichViewNet:::isCytoscapeRunning()) {
    enrichViewNet:::createNetworkForCytoscape(nodeEdgeInfo=information, 
        title="Test", collection="New Collection")
}
#> Unable to connect to Cytoscape. 
#> CX JSON file will be created.
```
